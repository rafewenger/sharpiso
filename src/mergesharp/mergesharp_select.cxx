/// \file mergesharp_select.cxx
/// Select cubes containing isosurface vertices on sharp edges and corners.

/*
Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iomanip>  
#include <string>  
#include <stdio.h>

#include "ijkcoord.txx"

#include "mergesharp_types.h"
#include "mergesharp_select.h"
#include "mergesharp_debug.h"

using namespace std;
using namespace SHARPISO;
using namespace MERGESHARP;
using namespace IJK;


// *** DEBUG ***
#include "ijkprint.txx"


namespace {

  void reset_covered_isovert_positions
  (const SHARPISO_BOOL_GRID & covered_grid, ISOVERT & isovert);

  bool check_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube0_index);

  void check_and_set_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube_index);

  bool is_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag);

  bool is_facet_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag, int & orth_dir, int & side);

  bool is_facet_on_corner_3x3x3
  (const GRID_CUBE_DATA & c, const SHARPISO_BOOL_GRID & covered_grid,
   const ISOVERT & isovert, int & orth_dir, int & side,
   VERTEX_INDEX & corner_cube_index);

  bool is_dist2_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag);

  bool is_facet_neighbor_covered_B
  (const SHARPISO_BOOL_GRID & covered_grid, const VERTEX_INDEX cube_index0);

  bool is_edge_neighbor_covered_B
  (const SHARPISO_BOOL_GRID & covered_grid, const VERTEX_INDEX cube_index0);

  int get_index_smallest_abs(const COORD_TYPE coord[DIM3]);
};


// **************************************************
// SELECT ROUTINES
// **************************************************

void select_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  isovert.gcube_list[gcube_index].flag = SELECTED_GCUBE;

  selection_data.covered_grid.Set(cube_index, true);
  selection_data.selected_grid.Set(cube_index, true);
  selection_data.mismatch_table.SetMismatchGrid
    (isovert.gcube_list[gcube_index], isovert,
     selection_data.selected_grid, selection_data.mismatch_grid);


  MSDEBUG();
  if (flag_debug) {
    using namespace std;
    cerr << "*** Selecting cube " << cube_index << "  ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << "  Linf_dist: " << isovert.gcube_list[gcube_index].linf_dist;
    cerr << endl;
  }

  bin_grid_insert(scalar_grid, selection_data.BinWidth(), cube_index, 
                  selection_data.bin_grid);

  // mark all the neighbors as covered
  for (int i=0;i < isovert.grid.NumVertexNeighborsC(); i++) {
    VERTEX_INDEX cube_index2 = isovert.grid.VertexNeighborC(cube_index, i);

    selection_data.covered_grid.Set(cube_index2, true);

    NUM_TYPE gcube_index2 = isovert.index_grid.Scalar(cube_index2);
    if(gcube_index2 != ISOVERT::NO_INDEX) {
      isovert.gcube_list[gcube_index2].flag = flag;

      if (isovert.gcube_list[gcube_index2].covered_by ==
          isovert.gcube_list[gcube_index2].cube_index) {
        isovert.gcube_list[gcube_index2].covered_by = cube_index;
      }
    }
  }
}


void check_and_select_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;
  VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;
  VERTEX_INDEX v1, v2;

  if (isovert.gcube_list[gcube_index].linf_dist >= linf_dist_threshold)
    { return; }

  // Check if the sharp vertex is inside a covered cube.
  if (check_covered_point(selection_data.covered_grid, isovert, gcube_index)) {
    isovert.gcube_list[gcube_index].flag = COVERED_POINT;
    return;
  }

  bool triangle_flag =
    creates_triangle_new
    (scalar_grid, isovert, cube_index, isovalue, 
     selection_data.bin_grid, selection_data.BinWidth(), v1, v2);

  if (!triangle_flag) {
    select_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, flag, 
       isovert, selection_data);
  }
  else {

    if (flag_debug) {
      MSDEBUG();
      scalar_grid.PrintIndexAndCoord
        (cerr, "  Cube ", cube_index, " unavailable because of cubes ",
         v1, " and ", v2, "\n");
    }

    isovert.gcube_list[gcube_index].flag = UNAVAILABLE_GCUBE;
	}

  /* DEBUG
  select_cube
    (scalar_grid, isovalue, isovert_param, gcube_index, flag, 
     isovert, selection_data);
  */
}


void check_and_select_cube_no_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;
  VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;
  VERTEX_INDEX v1, v2;

  if (isovert.gcube_list[gcube_index].linf_dist >= linf_dist_threshold)
    { return; }

  // Check if the sharp vertex is inside a covered cube.
  if (isovert.gcube_list[gcube_index].flag_conflict) { return; }

  bool triangle_flag =
    creates_triangle_new
    (scalar_grid, isovert, cube_index, isovalue, 
     selection_data.bin_grid, selection_data.BinWidth(), v1, v2);

  if (!triangle_flag) {
    select_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, flag, 
       isovert, selection_data);
  }
  else {

    if (flag_debug) {
      MSDEBUG();
      scalar_grid.PrintIndexAndCoord
        (cerr, "  Cube ", cube_index, " unavailable because of cubes ",
         v1, " and ", v2, "\n");
    }

    isovert.gcube_list[gcube_index].flag = UNAVAILABLE_GCUBE;
	}
}


/// Unselect cube grid_cube.  Uncover all neighbors.
/// @pre Neigboring cubes are covered only grid_cube.
void unselect_cube(const NUM_TYPE gcube_index, ISOVERT & isovert,
                   SELECTION_DATA & selection_data)
{
  VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE;
  bin_grid_remove(isovert.grid, selection_data.BinWidth(), cube_index, 
                  selection_data.bin_grid);

  if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < isovert.grid.NumVertexNeighborsC(); k++) {

      VERTEX_INDEX cubeB_index = isovert.grid.VertexNeighborC(cube_index, k);
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

      selection_data.covered_grid.Set(cubeB_index, false);
      if (gcubeB_index != ISOVERT::NO_INDEX) {

        if (isovert.gcube_list[gcubeB_index].num_eigenvalues > 1 &&
            !isovert.gcube_list[gcubeB_index].flag_centroid_location) {
          isovert.gcube_list[gcubeB_index].flag = AVAILABLE_GCUBE;
        }
        else {
          isovert.gcube_list[gcubeB_index].flag = SMOOTH_GCUBE;
        }
        isovert.gcube_list[gcubeB_index].covered_by = cubeB_index;
      }
    }
  }
  else {
    // Handle boundary case
  }
}

// Select corner cubes (eigen value 3)
void select_corner_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sortd_ind2gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  PROCEDURE_ERROR error("select_corner_cubes");

  for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    const NUM_TYPE gcube_index = sortd_ind2gcube_list[ind];
    const BOUNDARY_BITS_TYPE boundary_bits =
      isovert.gcube_list[gcube_index].boundary_bits;

    // check boundary
    if (boundary_bits == 0) {
      // select corners first
      if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE &&
          isovert.gcube_list[gcube_index].num_eigenvalues > 2)
        {
          check_and_select_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_CORNER_GCUBE, isovert, selection_data);
        }
    }
  }
}

// Select near corners. neighbor of covered corner
void select_cubes_near_corners
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 vector<NUM_TYPE> sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;


  for (int ind=0; ind < sharp_gcube_list.size(); ind++) {

    const NUM_TYPE gcube_index = sharp_gcube_list[ind];
    const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE) {

      // check boundary
      const BOUNDARY_BITS_TYPE boundary_bits =
        isovert.gcube_list[gcube_index].boundary_bits;

      if (boundary_bits == 0 && 
          isovert.gcube_list[gcube_index].linf_dist <= 0.5 ) {

        bool flag = is_neighbor
          (isovert.gcube_list[gcube_index], isovert,  COVERED_CORNER_GCUBE );

        if (flag) {
          check_and_select_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);
        }

        if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) 
          { isovert.gcube_list[gcube_index].flag_near_corner = true; }
      }
    }
  }
}

void check_and_select_edge_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG gcube_flag,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const BOUNDARY_BITS_TYPE boundary_bits =
    isovert.gcube_list[gcube_index].boundary_bits;

  check_and_set_covered_point
    (selection_data.covered_grid, isovert, gcube_index);

  // check boundary
  if (boundary_bits == 0) {

    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE &&
        isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
        !isovert.gcube_list[gcube_index].flag_conflict)
      {
        check_and_select_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, gcube_flag, 
           isovert, selection_data);
      }
  }
}

/// Select edge cubes (eigenvalue 2)
void select_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    const NUM_TYPE gcube_index = sharp_gcube_list[i];
    check_and_select_edge_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);

    if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {
      const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
      if (selection_data.mismatch_grid.Scalar(cube_index))
        { isovert.gcube_list[gcube_index].flag_ignore_mismatch = true; }
    }
  }
}

void check_and_select_edge_cube_no_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG gcube_flag,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const BOUNDARY_BITS_TYPE boundary_bits =
    isovert.gcube_list[gcube_index].boundary_bits;

  check_and_set_covered_point
    (selection_data.covered_grid, isovert, gcube_index);

  // check boundary
  if (boundary_bits == 0) {

    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE ||
        isovert.gcube_list[gcube_index].flag == COVERED_POINT) {
      if (isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
          !isovert.gcube_list[gcube_index].flag_conflict) {

        check_and_select_cube_no_conflict
          (scalar_grid, isovalue, isovert_param, gcube_index, gcube_flag, 
           isovert, selection_data);
      }
    }
  }
}

/// Select edge cubes (eigenvalue 2) as long as selected point
///   is not in active cube.
void select_edge_cubes_no_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    check_and_select_edge_cube_no_conflict
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);
  }
}

/// Select edge cubes (eigenvalue 2)
void select_edge_cubes_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (selection_data.mismatch_grid.Scalar(cube_index)) { 

      // *** DEBUG ***
      if (flag_debug) {
        MSDEBUG();
        scalar_grid.PrintIndexAndCoord
          (cerr, "Not selecting mismatched cube ", cube_index, "\n");
      }

      continue; 
    }

    check_and_select_edge_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);
  }
}

/// Select edge cubes (eigenvalue 2) with isovert within given distance
///   to cube center.
void select_edge_cubes_within_dist
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 const COORD_TYPE max_distance,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].linf_dist > max_distance)
      { continue; }

    if (selection_data.mismatch_grid.Scalar(cube_index)) { 

      // *** DEBUG ***
      if (flag_debug) {
        MSDEBUG();
        scalar_grid.PrintIndexAndCoord
          (cerr, "Not selecting mismatched cube ", cube_index, "\n");
      }

      continue; 
    }

    check_and_select_edge_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);
  }
}


/// Select edge cubes (eigenvalue 2) near corners
void select_edge_cubes_near_corners
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (selection_data.mismatch_grid.Scalar(cube_index)) { 

      // *** DEBUG ***
      if (flag_debug) {
        MSDEBUG();
        scalar_grid.PrintIndexAndCoord
          (cerr, "Not selecting mismatched cube ", cube_index, "\n");
      }

      continue; 
    }

    //  *** SHOULD REPLACE is_neighbor BY is_edge_neighbor ***
    if (is_dist2_neighbor
        (isovert.gcube_list[gcube_index], isovert,  COVERED_CORNER_GCUBE) ||
        is_neighbor
        (isovert.gcube_list[gcube_index], isovert, COVERED_CORNER_GCUBE)) {

      check_and_select_edge_cube
        (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
         isovert, selection_data);
    }
  }
}

/// Select one edge cube (eigenvalue 2).
/// @param from_list List of gcube indices sorted by increasing
///    distance from isovert_coord to cube centers.
///    Select cube from list from_list.
bool select_one_edge_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & from_list,
 INDEX_DIFF_TYPE & selected_gcube_index,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const int bin_width = isovert_param.bin_width;
  COORD_TYPE linf_dist = isovert_param.linf_dist_thresh_merge_sharp;

  selected_gcube_index = ISOVERT::NO_INDEX;

  for (int ind=0; ind < from_list.size(); ind++) {

    NUM_TYPE gcube_index = from_list[ind];

    GRID_CUBE_DATA c;
    c = isovert.gcube_list[gcube_index];

    if (c.boundary_bits == 0)
      if (c.flag == AVAILABLE_GCUBE &&
          c.linf_dist < linf_dist && c.num_eigenvalues == 2 &&
          !c.flag_conflict)
        {
          check_and_select_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);

          if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {

            MSDEBUG();
            if (flag_debug) {
              cerr << "  Selected new cube: " << c.cube_index << "  ";
              ijkgrid_output_vertex_coord(cerr, scalar_grid, c.cube_index);
              cerr << " Linf dist: " << c.linf_dist << endl;
              cerr << "  isovert_position: ";
              IJK::print_coord3D
                (cerr, isovert.gcube_list[gcube_index].isovert_coord);
              cerr << endl;
            }

            selected_gcube_index = gcube_index;
            return(true); 
          }
        }
  }

  return(false);
}

/// Get corner or edge cubes around given cube.
void get_corner_or_edge_cubes_around_cube
(const SHARPISO_GRID_NEIGHBORS & grid,
 const ISOVERT & isovert,
 const NUM_TYPE gcube_index,
 std::vector<NUM_TYPE> & neighbor_list)
{
  const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  neighbor_list.clear();
  if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < grid.NumVertexNeighborsC(); k++) {

      VERTEX_INDEX adjacent_cube_index = 
        grid.VertexNeighborC(cube_index, k);

      INDEX_DIFF_TYPE adjacent_gcube_index = 
        isovert.GCubeIndex(adjacent_cube_index);

      if (adjacent_gcube_index != ISOVERT::NO_INDEX) {
        if (isovert.gcube_list[adjacent_gcube_index].num_eigenvalues > 1) {
          neighbor_list.push_back(adjacent_gcube_index);
        }
      }
    }
  }
}

/// Reselect two edge cubes around gcube_index.
/// Select one cube from list from_list.
/// @param from_list List of gcube indices sorted by increasing
///    distance from isovert_coord to cube centers.
void reselect_two_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX gcube_index,
 const std::vector<NUM_TYPE> & from_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const int bin_width = isovert_param.bin_width;
  vector<NUM_TYPE> neighbor_list;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);

  get_corner_or_edge_cubes_around_cube
    (isovert.grid, isovert, gcube_index, neighbor_list);
  sort(neighbor_list.begin(), neighbor_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) {
    scalar_grid.PrintIndexAndCoord
      (cerr, "  From list: ", isovert.CubeIndex(gcube_index), ": ");
    for (int i = 0; i < from_list.size(); i++) {
      scalar_grid.PrintIndexAndCoord
        (cerr, " ", isovert.CubeIndex(from_list[i]), "");
    }
    cerr << endl;
    scalar_grid.PrintIndexAndCoord
      (cerr, "  Neighbors to cube: ", isovert.CubeIndex(gcube_index), ": ");
    for (int i = 0; i < neighbor_list.size(); i++) {
      scalar_grid.PrintIndexAndCoord
        (cerr, " ", isovert.CubeIndex(neighbor_list[i]), "");
    }
    cerr << endl;
  }

  unselect_cube(gcube_index, isovert, selection_data);

  INDEX_DIFF_TYPE gcubeB_index, gcubeC_index;
  if (select_one_edge_cube(scalar_grid, isovalue, isovert_param, 
                           from_list, gcubeB_index, isovert, selection_data)) {

    bool flag = select_one_edge_cube
      (scalar_grid, isovalue, isovert_param, neighbor_list, 
       gcubeC_index, isovert, selection_data);
  }
  else {

    // Reselect cube gcube_index
    select_cube(scalar_grid, isovalue, isovert_param, gcube_index, 
                COVERED_A_GCUBE, isovert, selection_data);
  }
}

/// Reselect edge cubes around gcube_index.
/// @param nearby_sharp_list List of nearby sharp cubes.
void reselect_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX gcube_index,
 const std::vector<VERTEX_INDEX> & nearby_selected_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
  GRID_COORD_TYPE cube_coord[DIM3], cubeB_coord[DIM3], cubeC_coord[DIM3];
  GRID_COORD_TYPE rmin[DIM3], rmax[DIM3];
  std::vector<NUM_TYPE> from_list;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);

  scalar_grid.ComputeCoord(cube_index, cube_coord);

  int max_overlap_dim = -1;
  VERTEX_INDEX max_overlap_cube = 0;
  for (NUM_TYPE i = 0; i < nearby_selected_list.size(); i++) {
    VERTEX_INDEX cubeB_index = nearby_selected_list[i];

    if (cube_index != cubeB_index) {

      scalar_grid.ComputeCoord(cubeB_index, cubeB_coord);

      int overlap_dim;
      find_3x3x3_overlap(cube_coord, cubeB_coord, rmin, rmax, overlap_dim);
      if (overlap_dim > max_overlap_dim) { 
        max_overlap_dim = overlap_dim; 
        max_overlap_cube = cubeB_index;
      }
    }
  }

  if (max_overlap_dim >= 3) { 
    // Some cube is covered by both gcube_index and another selected cube.
    return; 
  }
  else {

    INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(max_overlap_cube);
    if (gcubeB_index == ISOVERT::NO_INDEX) { return; }

    const COORD_TYPE * isovert_coordA =
      isovert.gcube_list[gcube_index].isovert_coord;
    const COORD_TYPE * directionA =
      isovert.gcube_list[gcube_index].direction;
    const COORD_TYPE * isovert_coordB =
      isovert.gcube_list[gcubeB_index].isovert_coord;
    const COORD_TYPE * directionB =
      isovert.gcube_list[gcubeB_index].direction;

    // Check that line containing sharp edges in each cube pass
    //   through/near other cube.

    if (isovert.gcube_list[gcube_index].num_eigenvalues == 2) {
      COORD_TYPE distB_to_lineA;
      IJK::compute_distance_to_line_3D
        (isovert_coordB, isovert_coordA, directionA, distB_to_lineA);
      if (distB_to_lineA > isovert_param.max_dist_to_sharp_edge)
        { return; };
    }


    if (isovert.gcube_list[gcubeB_index].num_eigenvalues == 2) {
      COORD_TYPE distA_to_lineB;
      IJK::compute_distance_to_line_3D
        (isovert_coordA, isovert_coordB, directionB, distA_to_lineB);
      if (distA_to_lineB > isovert_param.max_dist_to_sharp_edge)
        { return; };
    }

    MSDEBUG();
    if (flag_debug) {
      cerr << "Reselecting around cube " << cube_index << ".  ";
      IJK::print_coord3D(cerr, cube_coord);
      cerr << "  linf_dist: " << isovert.gcube_list[gcube_index].linf_dist
           << endl;
      cerr << "  Overlaps cube: " << max_overlap_cube << "  ";
      ijkgrid_output_vertex_coord(cerr, isovert.grid, max_overlap_cube);
      cerr << "  overlap_dim: " << max_overlap_dim << endl;

      cerr << "  Cube isovertex: ";
      IJK::print_coord3D(cerr, isovert_coordA);
      cerr << "  edge dir: ";
      IJK::print_coord3D(cerr, directionA);
      cerr << endl;
      cerr << "  Overlap cube isovertex: ";
      IJK::print_coord3D(cerr, isovert_coordB);
      cerr << "  edge dir: ";
      IJK::print_coord3D(cerr, directionB);
      cerr << endl;
    }

    if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

      scalar_grid.ComputeCoord(max_overlap_cube, cubeB_coord);

      for (NUM_TYPE k = 0; k < isovert.grid.NumVertexNeighborsC(); k++) {

        VERTEX_INDEX cubeC_index = isovert.grid.VertexNeighborC(cube_index, k);
        INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(cubeC_index);

        if (gcubeC_index != ISOVERT::NO_INDEX) {

          scalar_grid.ComputeCoord(cubeC_index, cubeC_coord);

          int overlap_dim;
          find_3x3x3_overlap(cubeB_coord, cubeC_coord, rmin, rmax, overlap_dim);
          if (overlap_dim >= 3) 
            { from_list.push_back(gcubeC_index); }
        }
      }

      if (from_list.size() > 0) {
        sort(from_list.begin(), from_list.end(), gcube_compare);

        reselect_two_edge_cubes
          (scalar_grid, isovalue, isovert_param, gcube_index, from_list, 
           isovert, selection_data);
      }
    }
    else {
      // Handle boundary case
    }
  }

}

/// Reselect edge cubes around gcube_index.
/// @param nearby_selected_list List of nearby sharp cubes.
void reselect_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  const int bin_width = isovert_param.bin_width;
  std::vector<VERTEX_INDEX> nearby_selected_list;

  for (NUM_TYPE i = 0; i < sharp_gcube_list.size(); i++) {

    NUM_TYPE gcube_index = sharp_gcube_list[i];
    if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE &&
        isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
        !isovert.gcube_list[gcube_index].flag_near_corner) {

      VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;

      get_selected(scalar_grid, cube_index, selection_data.bin_grid, 
                   selection_data.BinWidth(), nearby_selected_list);

      reselect_edge_cubes
        (scalar_grid, isovalue, isovert_param, gcube_index, 
         nearby_selected_list, isovert, selection_data);
    }
  }
}

// Set gcube_list[i].covered_by to gcube_list[i].cube_index for each i.
void initialize_covered_by(ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) 
    { isovert.gcube_list[i].covered_by = isovert.gcube_list[i].cube_index; }
}

// *** NOT TESTED ***
/// Version of select sharp for hermite data.
void MERGESHARP::select_sharp_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const int dimension = scalar_grid.Dimension();
  const int bin_width = isovert_param.bin_width;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
  std::vector<NUM_TYPE> sharp_gcube_list;
  SELECTION_DATA selection_data(scalar_grid, isovert_param);

  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  initialize_covered_by(isovert);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting corner cubes." << endl; 
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting near corner cubes." << endl; 
  }

  // select cubes near corners
  select_cubes_near_corners
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

  // select edge cubes.
  select_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // reselect edge cubes
  reselect_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, isovert,
     selection_data);
}


void MERGESHARP::select_sharp_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const int dimension = scalar_grid.Dimension();
  const int bin_width = isovert_param.bin_width;
  const bool flag_recompute = 
    isovert_param.flag_recompute_isovert &&
    isovert_param.flag_recompute_changing_gradS_offset;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
  std::vector<NUM_TYPE> sharp_gcube_list;
  SELECTION_DATA selection_data(scalar_grid, isovert_param);


  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  initialize_covered_by(isovert);

  MSDEBUG();
  flag_debug = false;
  if (flag_debug) {
    cerr << endl << "*** Selecting corner cubes." << endl;
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  if (flag_recompute) {

    recompute_covered_point_positions
      (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
       isovert_param, isovert);

    reset_covered_isovert_positions(selection_data.covered_grid, isovert);
  }

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting near corner cubes." << endl; 
  }

  // select cubes near corners
  select_cubes_near_corners
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  if (flag_recompute) {

    recompute_covered_point_positions
      (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
       isovert_param, isovert);

    reset_covered_isovert_positions(selection_data.covered_grid, isovert);
  }

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

  // select edge cubes.
  select_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  if (flag_recompute) {

    MSDEBUG();
    if (flag_debug) 
      { cerr << endl << "*** Recomputing covered point positions." << endl; }

    recompute_covered_point_positions
      (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
       isovert_param, isovert);

    reset_covered_isovert_positions(selection_data.covered_grid, isovert);
  }

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl 
           << "*** Selecting edge cubes with isovert location outside cube." 
           << endl; }

  // select edge cubes which do not conflict with active cubes.
  select_edge_cubes_no_conflict
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Reselecting edge cubes." << endl; }

  // reselect edge cubes
  reselect_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug) {
    cerr << endl << "*** Selecting edge cubes (again.)" << endl; 
  }

  // Retry selecting edge cubes.
  select_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  if (flag_recompute) {

    MSDEBUG();
    if (flag_debug) 
      { cerr << endl << "*** Recomputing covered point positions." << endl; }

    recompute_covered_point_positions
      (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
       isovert_param, isovert);
  }

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
  set_cover_type(isovert);

  MSDEBUG();
  flag_debug = false;
}


// **************************************************
// MISMATCH GRID ROUTINES
// **************************************************


/// Set mismatch flags for cubes near corners.
void set_mismatch_near_corner
(const ISOVERT & isovert, const GRID_CUBE_DATA & gcubeA, 
 const SELECTION_DATA & selection_data, SHARPISO_BOOL_GRID & mismatch_grid)
{
  int orth_dir, side;
  GRID_COORD_TYPE distance2boundary;
  PROCEDURE_ERROR error("set_mismatch_near_corner");

  isovert.grid.ComputeCubeDistanceToGridBoundary
    (gcubeA.cube_coord, distance2boundary);

  if (distance2boundary < 2) { return; }

  VERTEX_INDEX corner_cube_index;
  if (is_facet_on_corner_3x3x3
      (gcubeA, selection_data.covered_grid, isovert, 
       orth_dir, side, corner_cube_index)) {

    const VERTEX_INDEX cubeA_index = gcubeA.cube_index;
    const INDEX_DIFF_TYPE corner_gcube_index = 
      isovert.GCubeIndex(corner_cube_index, error);

    if (corner_gcube_index == ISOVERT::NO_INDEX) { throw error; }

    for (int i = 1; i < DIM3; i++) {

      const int d1 = (orth_dir+i)%DIM3;
      const int d2 = (orth_dir+2*i)%DIM3;

      for (int j1 = 0; j1 < 2; j1++) {
        for (int j2 = -1; j2 < 2; j2++) {

          const VERTEX_INDEX cubeB_index =
            cubeA_index + 2*isovert.grid.AxisIncrement(d1)*(2*j1-1) +
            isovert.grid.AxisIncrement(d2)*j2;
          const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

          if (gcubeB_index == ISOVERT::NO_INDEX) { continue; }
          if (isovert.NumEigenvalues(gcubeB_index) == 1)
            { continue; }

          COORD_TYPE linf_distance;
          const GRID_COORD_TYPE * cubeA_coord = gcubeA.cube_coord;
          const GRID_COORD_TYPE * cubeB_coord =
            isovert.gcube_list[gcubeB_index].cube_coord;

          IJK::compute_Linf_distance
            (DIM3, cubeA_coord, cubeB_coord, linf_distance);

          if (linf_distance > 2) { continue; }

          MSDEBUG();
          if (flag_debug) {
            isovert.grid.PrintIndexAndCoord
              (cerr, "  Cube ", cubeA_index,
               " near corner cube ", corner_cube_index, "\n.");
            isovert.grid.PrintIndexAndCoord
              (cerr, "      does not match ", cubeB_index, "\n");
          }

          mismatch_grid.Set(cubeA_index, true);
        }
      }
    }
  }
}


/// Set mismatch flags for cubes near corners.
void set_mismatch_near_corner
(const ISOVERT & isovert,  const vector<NUM_TYPE> & sharp_gcube_list,
 const SELECTION_DATA & selection_data, SHARPISO_BOOL_GRID & mismatch_grid)
{
  for (int i = 0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];

    if (isovert.NumEigenvalues(gcube_index) == 2) {
      set_mismatch_near_corner
        (isovert, isovert.gcube_list[gcube_index], selection_data,
         mismatch_grid);
    }
  }
}

/// Set mismatch flags for cubes near corners.
/// Version using selection_data.mismatch_grid
void set_mismatch_near_corner
(const ISOVERT & isovert,  const vector<NUM_TYPE> & sharp_gcube_list,
 SELECTION_DATA & selection_data)
{
  set_mismatch_near_corner
    (isovert, sharp_gcube_list, selection_data, selection_data.mismatch_grid);
}

/// Compute mismatch grid.
void compute_mismatch_grid
(const ISOVERT & isovert,  const vector<NUM_TYPE> & sharp_gcube_list,
 const SELECTION_DATA & selection_data, SHARPISO_BOOL_GRID & mismatch_grid)
{
  mismatch_grid.SetAll(false);

  for (int i = 0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];

    if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {
      selection_data.mismatch_table.SetMismatchGrid
        (isovert.gcube_list[gcube_index], isovert,
         selection_data.selected_grid, mismatch_grid);
    }

    if (isovert.NumEigenvalues(gcube_index) == 2) {
      set_mismatch_near_corner
        (isovert, isovert.gcube_list[gcube_index], 
         selection_data, mismatch_grid);
    }
  }

}

/// Compute mismatch grid.
/// Version using selection_data.mismatch_grid
void compute_mismatch_grid
(const ISOVERT & isovert,  const vector<NUM_TYPE> & sharp_gcube_list,
 SELECTION_DATA & selection_data)
{
  compute_mismatch_grid
    (isovert, sharp_gcube_list, selection_data, selection_data.mismatch_grid);
}


// **************************************************
// SELECT MOD K ROUTINES
// **************************************************

/// Return true if k coordinates are congruent to zero.
template <const int k, const int modulus>
bool are_k_coord_cong2zero(const GRID_COORD_TYPE coord[DIM3])
{
  int num_zero = 0;
  if ((coord[0]%modulus) == 0) { num_zero++; }
  if ((coord[1]%modulus) == 0) { num_zero++; }
  if ((coord[2]%modulus) == 0) { num_zero++; }
  if (num_zero == k) { return(true); }
  return(false);
}

/// Return true if k coordinates are congruent to zero.
/// Specialization for k = 3, modulus = 3.
template <>
bool are_k_coord_cong2zero<3,3>(const GRID_COORD_TYPE coord[DIM3])
{
  const int modulus = 3;

  if ((coord[0]%modulus) != 0) { return(false); }
  if ((coord[1]%modulus) != 0) { return(false); }
  if ((coord[2]%modulus) != 0) { return(false); }
  return(true);
}

/// Return true if k coordinates are congruent to zero.
/// Specialization for k = 3, modulus = 6.
template <>
bool are_k_coord_cong2zero<3,6>(const GRID_COORD_TYPE coord[DIM3])
{
  const int modulus = 6;

  if ((coord[0]%modulus) != 0) { return(false); }
  if ((coord[1]%modulus) != 0) { return(false); }
  if ((coord[2]%modulus) != 0) { return(false); }
  return(true);
}

/// Return true if k coordinates are congruent to zero.
/// Specialization for k = 0, modulus = 3.
template <>
bool are_k_coord_cong2zero<0,3>(const GRID_COORD_TYPE coord[DIM3])
{
  const int modulus = 3;

  if ((coord[0]%modulus) == 0) { return(false); }
  if ((coord[1]%modulus) == 0) { return(false); }
  if ((coord[2]%modulus) == 0) { return(false); }
  return(true);
}

/// Compute (coord[] mod modulus)
template <const int modulus>
inline void compute_coord_mod
(const GRID_COORD_TYPE coord[DIM3], GRID_COORD_TYPE coord_mod[DIM3])
{
  coord_mod[0] = coord[0]%modulus;
  coord_mod[1] = coord[1]%modulus;
  coord_mod[2] = coord[2]%modulus;
}

template <const int modulus>
bool is_fadj2post(const GRID_COORD_TYPE coord[DIM3])
{
  int num_diff = 0;
  for (int d = 0; d < DIM3; d++) {
    int r = coord[d]%modulus;
    
    if (r == 1 || r == modulus-1)
      { num_diff++; }
    else if (r != 0)
      { return(false); }
  }

  if (num_diff == 1) { return(true); }

  return(false);
}

template <const int modulus>
bool is_eadj2post(const GRID_COORD_TYPE coord[DIM3])
{
  int num_diff = 0;
  for (int d = 0; d < DIM3; d++) {
    int r = coord[d]%modulus;
    
    if (r == 1 || r == modulus-1)
      { num_diff++; }
    else if (r != 0)
      { return(false); }
  }

  if (num_diff == 2) { return(true); }

  return(false);
}

/// Select edge cubes (eigenvalue 2) with k coord cong to 0 mod modulus
template <const int k, const int modulus>
void select_cubes_with_k_coord_cong2zero
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with " << k
           << " coord cong to 0 mod " << modulus << "." << endl; }

  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    const GRID_COORD_TYPE * cube_coord = 
      isovert.gcube_list[gcube_index].cube_coord;

    if (are_k_coord_cong2zero<k,modulus>(cube_coord)) {

      check_and_select_edge_cube
        (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
         isovert, selection_data);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}


// **************************************************
// SELECT MOD 3 ROUTINES
// **************************************************

template <int d>
void remove_mod3_remainder_d
(const GRID_COORD_TYPE coord0[DIM3], GRID_COORD_TYPE coord1[DIM3],
 int & num_equal, int & mod3_coord, int & not_mod3_coord)
{
  coord1[d] = (coord0[d]-(int(coord0[d])%3));
  if (coord1[d] == coord0[d]) {
    num_equal++;
    mod3_coord = d;
  }
  else {
    not_mod3_coord = d;
  }
}

void remove_mod3_remainder
(const GRID_COORD_TYPE coord0[DIM3], GRID_COORD_TYPE coord1[DIM3],
 int & num_equal, int & mod3_coord, int & not_mod3_coord)
{
  num_equal = 0;
  remove_mod3_remainder_d<0>
    (coord0, coord1, num_equal, mod3_coord, not_mod3_coord);
  remove_mod3_remainder_d<1>
    (coord0, coord1, num_equal, mod3_coord, not_mod3_coord);
  remove_mod3_remainder_d<2>
    (coord0, coord1, num_equal, mod3_coord, not_mod3_coord);
}

/// Return number of zero coordinates.
inline int count_num_zero(const GRID_COORD_TYPE coord[DIM3])
{
  int num_zero = 0;
  if (coord[0] == 0) { num_zero++; }
  if (coord[1] == 0) { num_zero++; }
  if (coord[2] == 0) { num_zero++; }
  return(num_zero);
}

/// Snap coordinates to +/- 1.
void snap_coord_to_plus_minus_one
(const COORD_TYPE coord[DIM3], GRID_COORD_TYPE pm_one_coord[DIM3])
{
  for (int d = 0; d < DIM3; d++) {
    if (coord[d] >= 0) { pm_one_coord[d] = 1; }
    else {pm_one_coord[d] = -1; }
  }
}

/// Return true if edge dir points to an adjacent 0 mod3 cube.
bool does_edge_dir_point_to_adjacent_zero_mod3
(const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE edge_dir[DIM3])
{
  GRID_COORD_TYPE pm_one_edge_dir[DIM3];
  GRID_COORD_TYPE temp[DIM3];

  snap_coord_to_plus_minus_one(edge_dir, pm_one_edge_dir);

  for (int d = 0; d < DIM3; d++) 
    { temp[d] = (cube_coord[d]+pm_one_edge_dir[d])%3; }

  if (temp[0] == 0 && temp[1] == 0 && temp[2] == 0)
    { return(true); }

  for (int d = 0; d < DIM3; d++) 
    { temp[d] = (cube_coord[d]-pm_one_edge_dir[d])%3; }

  if (temp[0] == 0 && temp[1] == 0 && temp[2] == 0)
    { return(true); }

  return(false);
}


// Return true if nearby cube whose coordinates are 0 mod3 is selected.
bool is_nearby_mod3_selected
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 VERTEX_INDEX cube0_index)
{
  const INDEX_DIFF_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);

  if (gcube0_index == ISOVERT::NO_INDEX) { return(false); }

  const GRID_COORD_TYPE * cube0_coord = 
    isovert.gcube_list[gcube0_index].cube_coord;
  GRID_COORD_TYPE cubeL_coord[DIM3];   // Lower cube coord
  GRID_COORD_TYPE cubeH_coord[DIM3];   // Higher cube coord
  GRID_COORD_TYPE x[DIM3], diff[DIM3];
  bool flag_zero_mod3;                 // true, if some coordinate is 0 mod 3.

  flag_zero_mod3 = false;
  for (int d = 0; d < DIM3; d++) {
    cubeL_coord[d] = (cube0_coord[d]-(int(cube0_coord[d])%3));
    if (cube0_coord[d] == cubeL_coord[d]) {
      flag_zero_mod3 = true;
      cubeH_coord[d] = cubeL_coord[d];
    }
    else if (cube0_coord[d] + 5 >= grid.AxisSize(d)) {
      cubeH_coord[d] = cubeL_coord[d];
    }
    else {
      cubeH_coord[d] = cubeL_coord[d]+3;
    }
  }

  for (x[0] = cubeL_coord[0]; x[0] <= cubeH_coord[0]; x[0]+=3) {
    diff[0] = abs(x[0] - cube0_coord[0]);

    for (x[1] = cubeL_coord[1]; x[1] <= cubeH_coord[1]; x[1]+=3) {
      diff[1] = abs(x[1] - cube0_coord[1]);

      for (x[2] = cubeL_coord[2]; x[2] <= cubeH_coord[2]; x[2]+=3) {
        diff[2] = abs(x[2] - cube0_coord[2]);

        if (flag_zero_mod3) {
          if (diff[0]+diff[1]+diff[2] == 4) { continue; }
        }
        else {
          if (diff[0]+diff[1]+diff[2] == 6) { continue; }
        }

        VERTEX_INDEX cube1_index = grid.ComputeVertexIndex(x);

        INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
        if (gcube1_index == ISOVERT::NO_INDEX) { continue; }

        if (isovert.gcube_list[gcube1_index].flag == SELECTED_GCUBE)
          { return(true); }
      }
    }
  }

  return(false);
}

/// Return true if cube in 3x3x3 interior is selected.
bool is_cube_in_3x3x3_interior_selected
(const ISOVERT & isovert, const GRID_COORD_TYPE coord0[DIM3],
 VERTEX_INDEX & cube_index)
{
  GRID_COORD_TYPE coord1[DIM3];
  GRID_COORD_TYPE x[DIM3];

  for (int d = 0; d < DIM3; d++) {
    coord1[d] = coord0[d]+3;
    if (coord1[d] >= isovert.grid.AxisSize(d)) { return(false); }
  }

  for (x[0] = coord0[0]+1; x[0] < coord1[0]; x[0]++) {
    for (x[1] = coord0[1]+1; x[1] < coord1[1]; x[1]++) {
      for (x[2] = coord0[2]+1; x[2] < coord1[2]; x[2]++) {

        cube_index = isovert.grid.ComputeVertexIndex(x);

        INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);
        if (gcube_index == ISOVERT::NO_INDEX) { continue; }

        if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) 
          { return(true); }
      }
    }
  }

  return(false);
}

/// Return true if cube in 3x3x3 interior is selected.
bool is_cube_in_3x3x3_interior_selected
(const ISOVERT & isovert, const GRID_COORD_TYPE coord0[DIM3])
{
  VERTEX_INDEX cube_index;

  return(is_cube_in_3x3x3_interior_selected(isovert, coord0, cube_index));
}

/// Return true if cube in 3x3 region is sharp. 
///   (Note: 3x3 region, not 3x3x3 region.)
bool is_cube_in_3x3_region_sharp
(const ISOVERT & isovert, const GRID_COORD_TYPE coord0[DIM3], 
 const int orth_dir, VERTEX_INDEX & cube_index)
{
  const int d1 = (orth_dir+1)%DIM3;
  const int d2 = (orth_dir+2)%DIM3;
  GRID_COORD_TYPE coord1[DIM3];
  GRID_COORD_TYPE x[DIM3];

  IJK::copy_coord_3D(coord0, coord1);
  coord1[d1] += 3;
  coord1[d2] += 3;
  if (coord1[d1] >= isovert.grid.AxisSize(d1)) { return(false); }
  if (coord1[d2] >= isovert.grid.AxisSize(d2)) { return(false); }

  x[orth_dir] = coord0[orth_dir];
  for (x[d1] = coord0[d1]; x[d1] <= coord1[d1]; x[d1]++) {
    for (x[d2] = coord0[d2]; x[d2] <= coord1[d2]; x[d2]++) {

        cube_index = isovert.grid.ComputeVertexIndex(x);

        INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);
        if (gcube_index == ISOVERT::NO_INDEX) { continue; }

        if (isovert.NumEigenvalues(gcube_index) > 1 &&
            !isovert.gcube_list[gcube_index].flag_conflict)
          { return(true); }
    }
  }

  return(false);
}

/// Return true if cube in 3x3 region is sharp. 
///   (Note: 3x3 region, not 3x3x3 region.)
bool is_cube_in_3x3_region_sharp
(const ISOVERT & isovert, const GRID_COORD_TYPE coord0[DIM3], 
 const int orth_dir)
{
  VERTEX_INDEX cube_index;
  return(is_cube_in_3x3_region_sharp(isovert, coord0, orth_dir, cube_index));
}

/// Return true if nearby cube whose coordinates are not 0 mod3 is selected.
/// @pre Exactly one coordinate of cube0_index equals 0 mod 3. 
bool is_nearby_not_mod3_selected
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 VERTEX_INDEX cube0_index)
{
  const INDEX_DIFF_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);
  PROCEDURE_ERROR error("is_nearby_not_mod3_selected");

  if (gcube0_index == ISOVERT::NO_INDEX) { return(false); }

  const GRID_COORD_TYPE * cube0_coord = 
    isovert.gcube_list[gcube0_index].cube_coord;
  GRID_COORD_TYPE region_coord[DIM3];
  int num_zero;
  int mod3_coord, not_mod3_coord;

  remove_mod3_remainder
    (cube0_coord, region_coord, num_zero, mod3_coord, not_mod3_coord);

  if (num_zero == 1) {
    if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
      { return(true); }

    if (region_coord[mod3_coord] >= 3) {
      region_coord[mod3_coord] -= 3;
      if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
        { return(true); }
    }
  }
  else if (num_zero == 2) {
    int d1 = (not_mod3_coord+1)%DIM3;
    int d2 = (not_mod3_coord+2)%DIM3;

    if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
      { return(true); }

    if (region_coord[d1] >= 3) {
      region_coord[d1] -= 3;
      if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
        { return(true); }

      if (region_coord[d2] >= 3) {
        region_coord[d2] -= 3;
        if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
          { return(true); }

        region_coord[d1] += 3;
        if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
          { return(true); }
      }
    }
    else {
      if (region_coord[d2] >= 3) {
        region_coord[d2] -= 3;
        if (is_cube_in_3x3x3_interior_selected(isovert, region_coord))
          { return(true); }
      }
    }
  }

  return(false);
}

/// Select edge cubes (eigenvalue 2) whose coordinates are all NOT
///   congruent to 0 mod 3 and which are near a selected mod 3 cube.
void select_cubes_not_cong_zero_mod3_A
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Selecting cubes not cong to 0 mod 3 near a selected mod 3 cube." << endl; 
  }

  // Select from cubes whose coordinates are not congruent to 0 mod 3
  //   and who are near a selected mod 3 cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    const GRID_COORD_TYPE * cube_coord = 
      isovert.gcube_list[gcube_index].cube_coord;

    if (are_k_coord_cong2zero<0,3>(cube_coord)) {

      if (is_nearby_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
           isovert, selection_data);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) whose coordinates are all NOT
///   congruent to 0 mod 3 and which are facet adjacent to a covered cube.
void select_cubes_not_cong_zero_mod3_fadj2covered
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 const bool flag_require_cube_contains_isovert,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Selecting cubes not cong to 0 mod 3 which are facet adjacent to a covered cube." << endl; 
  }

  // Select from cubes whose coordinates are not congruent to 0 mod 3
  //   and who are facet adjacent to a covered cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    // On first pass, always require that cube contains isovert.
    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      if (is_facet_neighbor_covered_B
          (selection_data.covered_grid, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE,
           isovert, selection_data);
      }
    }
  }

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Removing requirement that cube contains isovert."
         << endl;
  }

  if (!flag_require_cube_contains_isovert) {
    // Redo selection, but without requiring that cube contains isovert.

    // Select from cubes whose coordinates are not congruent to 0 mod 3
    //   and who are facet adjacent to a covered cube.
    for (int i=0; i < sharp_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sharp_gcube_list[i];
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

      compute_coord_mod<3>
        (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

      if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
          cube_coord_mod3[2] != 0) {

        if (is_facet_neighbor_covered_B
            (selection_data.covered_grid, cube_index)) {
          check_and_select_edge_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);
        }
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) with one coord congruent to 0 mod 3
///   and which are near a selected mod 3 cube.
void select_cubes_with_one_coord_cong_zero_mod3_A
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 const bool flag_require_cube_contains_isovert,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 1 coord cong to 0 mod 3 near a selected mod 3 cube." << endl; }

  // Select from cubes with one coord congruent to 0 mod 3 near a selected mod 3 cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 1) {

      if (is_nearby_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
           isovert, selection_data);
      }
    }
  }

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Removing requirement that cube contains isovert."
         << endl;
  }

  if (!flag_require_cube_contains_isovert) {
    // Redo selection, but without requiring that cube contains isovert.

    // Select from cubes with one coord congruent to 0 mod 3 near a selected mod 3 cube.
    for (int i=0; i < sharp_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sharp_gcube_list[i];
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

      compute_coord_mod<3>
        (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

      int num_zero = count_num_zero(cube_coord_mod3);
      if (num_zero == 1) {

        if (is_nearby_mod3_selected(scalar_grid, isovert, cube_index)) {
          check_and_select_edge_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);
        }
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) with one coord congruent to 0 mod 3
///   and which are near a selected NOT mod 3 cube.
void select_cubes_with_one_coord_cong_zero_mod3_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 1 coord cong to 0 mod 3 near a selected NOT mod 3 cube." << endl; }

  // Select from cubes with one coord congruent to 0 mod 3 near a selected NOT mod 3 cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 1) {
      if (is_nearby_not_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
           isovert, selection_data);
      }
    }
  }

}


/// Select edge cubes (eigenvalue 2) with one coordinate congruent to 0 mod 3
///   and which are facet adjacent to a covered cube.
void select_cubes_with_one_coord_cong_zero_mod3_fadj2covered
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 const bool flag_require_cube_contains_isovert,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Selecting cubes with 1 coord cong to 0 mod 3 which are facet adjacent to a covered cube." << endl; 
  }

  // Select from cubes with one coordinate congruent to 0 mod 3
  //   and who are facet adjacent to a covered cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    // On first pass, always require that cube contains isovert.
    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 1) {
      if (is_facet_neighbor_covered_B
          (selection_data.covered_grid, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, 
           COVERED_A_GCUBE, isovert, selection_data);
      }
    }
  }

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Removing requirement that cube contains isovert."
         << endl;
  }

  if (!flag_require_cube_contains_isovert) {
    // Redo selection, but without requiring that cube contains isovert.

    // Select from cubes with one coordinate congruent to 0 mod 3
    //   and who are facet adjacent to a covered cube.
    for (int i=0; i < sharp_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sharp_gcube_list[i];
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

      compute_coord_mod<3>
        (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

      int num_zero = count_num_zero(cube_coord_mod3);
      if (num_zero == 1) {

        if (is_facet_neighbor_covered_B
            (selection_data.covered_grid, cube_index)) {
          check_and_select_edge_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);
        }
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) with two coord congruent to 0 mod 3
///   and which are near a selected NOT mod 3 cube.
void select_cubes_with_two_coord_cong_zero_mod3_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 const bool flag_require_cube_contains_isovert,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 2 coord cong to 0 mod 3 near a selected NOT mod 3 cube." << endl; }

  // Select from cubes with two coord congruent to 0 mod 3 near a selected NOT mod 3 cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 2) {
      if (is_nearby_not_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
           isovert, selection_data);
      }
    }
  }

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Removing requirement that cube contains isovert."
         << endl;
  }

  if (!flag_require_cube_contains_isovert) {
    // Redo selection, but without requiring that cube contains isovert.

    // Select from cubes with two coord congruent to 0 mod 3 near a selected NOT mod 3 cube.
    for (int i=0; i < sharp_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sharp_gcube_list[i];
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

      compute_coord_mod<3>
        (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

      int num_zero = count_num_zero(cube_coord_mod3);
      if (num_zero == 2) {
        if (is_nearby_not_mod3_selected(scalar_grid, isovert, cube_index)) {
          check_and_select_edge_cube
            (scalar_grid, isovalue, isovert_param, gcube_index, 
             COVERED_A_GCUBE, isovert, selection_data);
        }
      }
    }
  }
}

/// Select edge cubes (eigenvalue 2) with two coordinates congruent to 0 mod 3
///   where the sharp edge does not "twist" in the 3x3x3 regions.
void select_cubes_with_two_coord_cong_zero_mod3_C
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE region_coord[DIM3];
  int num_equal;
  int mod3_coord, not_mod3_coord;

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 2 coord cong to 0 mod 3, no edge twist." << endl; }

  // Select from cubes with two coordinates congruent to 0 mod 3.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    if (isovert.NumEigenvalues(gcube_index) != 2)
      { continue; }

    const GRID_COORD_TYPE * cube_coord = 
      isovert.gcube_list[gcube_index].cube_coord;

    remove_mod3_remainder
      (cube_coord, region_coord, num_equal, mod3_coord, not_mod3_coord);

    if (num_equal == 2) {

      if (cube_coord[not_mod3_coord] - region_coord[not_mod3_coord] == 1) {
        region_coord[not_mod3_coord] += 3;
      }

      VERTEX_INDEX cube2_index;

      if (region_coord[not_mod3_coord] < 
          isovert.grid.AxisSize(not_mod3_coord)) {

        if (is_cube_in_3x3_region_sharp
            (isovert, region_coord, not_mod3_coord, cube2_index)) {

          MSDEBUG();
          if (flag_debug) {
            isovert.grid.PrintIndexAndCoord
              (cerr, "  Cube ", cube_index, " conflicts with ", 
               cube2_index, "\n");
          }

          continue; 
        }

        const int d1 = (not_mod3_coord+1)%DIM3;
        const int d2 = (not_mod3_coord+2)%DIM3;
        if (region_coord[d1] >= 3) {
          region_coord[d1] -= 3;
          if (is_cube_in_3x3_region_sharp
              (isovert, region_coord, not_mod3_coord, cube2_index)) {

            MSDEBUG();
            if (flag_debug) {
              isovert.grid.PrintIndexAndCoord
                (cerr, "  Cube ", cube_index, " conflicts with ", 
                 cube2_index, "\n");
            }

            continue; 
          }

          if (region_coord[d2] >= 3) {
            region_coord[d2] -= 3;
            if (is_cube_in_3x3_region_sharp
                (isovert, region_coord, not_mod3_coord, cube2_index)) {

              MSDEBUG();
              if (flag_debug) {
                isovert.grid.PrintIndexAndCoord
                  (cerr, "  Cube ", cube_index, " conflicts with ", 
                   cube2_index, "\n");
              }

              continue; 
            }

            region_coord[d1] += 3;
            if (is_cube_in_3x3_region_sharp
                (isovert, region_coord, not_mod3_coord)) {

              if (flag_debug) {
                isovert.grid.PrintIndexAndCoord
                  (cerr, "  Cube ", cube_index, " conflicts with ", 
                   cube2_index, "\n");
              }

              continue; 
            }
          }
        }
        else {
          if (region_coord[d2] >= 3) {
            region_coord[d2] -= 3;
            if (is_cube_in_3x3_region_sharp
                (isovert, region_coord, not_mod3_coord)) {

              MSDEBUG();
              if (flag_debug) {
                isovert.grid.PrintIndexAndCoord
                  (cerr, "  Cube ", cube_index, " conflicts with ", 
                   cube2_index, "\n");
              }

              continue; 
            }
          }
        }
      }

      check_and_select_edge_cube
        (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
         isovert, selection_data);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}


/// Select edge cubes (eigenvalue 2) using mod 3 algorithm
void select_edge_cubes_mod3	
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert,
 SELECTION_DATA & selection_data)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];
  bool require_cube_contains_isovert;

  // Select edge cubes (eigenvalue 2) whose coord are all congruent to 0 mod 3.
  select_cubes_with_k_coord_cong2zero<3,3>
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select edge cubes whose coord are all NOT congruent to 0 mod 3
  //   and which are near a selected 0 mod 3 cube.
  select_cubes_not_cong_zero_mod3_A
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select edge cubes with one coord congruent to 0 mod 3
  //   and which are near a selected 0 mod 3 cube.
  require_cube_contains_isovert = true;
  select_cubes_with_one_coord_cong_zero_mod3_A
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     require_cube_contains_isovert, isovert, selection_data);

  // Select edge cubes with one coord congruent to 0 mod 3
  //   and which are near a selected NOT 0 mod 3 cube.
  select_cubes_with_one_coord_cong_zero_mod3_B
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select edge cubes with one coord congruent to 0 mod 3
  //   and which are near a selected 0 mod 3 cube.
  require_cube_contains_isovert = false;
  select_cubes_with_one_coord_cong_zero_mod3_A
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     require_cube_contains_isovert, isovert, selection_data);

  // Select edge cubes with two coord congruent to 0 mod 3
  //   and which are near a selected NOT 0 mod 3 cube.
  require_cube_contains_isovert = false;
  select_cubes_with_two_coord_cong_zero_mod3_B
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     require_cube_contains_isovert, isovert, selection_data);

  // Select edge cubes with two coord congruent to 0 mod 3
  //   where sharp edge does not "twist".
  select_cubes_with_two_coord_cong_zero_mod3_C
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select edge cubes with two coord congruent to 0 mod 3
  select_cubes_with_k_coord_cong2zero<2,3>
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select from cubes whose coordinates are not congruent to 0 mod 3
  //   and who are facet adjacent to a covered cube.
  require_cube_contains_isovert = false;
  select_cubes_not_cong_zero_mod3_fadj2covered
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     require_cube_contains_isovert, isovert, selection_data);

  // Select from cubes with one coord congruent to 0 mod 3
  //   and who are facet adjacent to a covered cube.
  require_cube_contains_isovert = false;
  select_cubes_with_one_coord_cong_zero_mod3_fadj2covered
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     require_cube_contains_isovert, isovert, selection_data);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with coord not cong to 0 mod 3." << endl; }

  // Select from cubes whose coordinates are not congruent to 0 mod 3
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      check_and_select_edge_cube
        (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
         isovert, selection_data);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 1 coord cong to 0 mod 3 adjacent to a covered cube." << endl; }

  // Select from cubes with one coordinate congruent to 0 mod 3
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_congruent(0);
    for (int d=0; d < 3; d++) {
      if (cube_coord_mod3[d] == 0) { num_congruent++; }
    }

    if (num_congruent == 1) {

      if (is_facet_neighbor_covered_B
          (selection_data.covered_grid, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
           isovert, selection_data);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select edge cubes with one coord congruent to 0 mod 3
  select_cubes_with_k_coord_cong2zero<1,3>
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  // Select edge cubes with two coord congruent to 0 mod 3
  select_cubes_with_k_coord_cong2zero<2,3>
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes containing isovert coord." << endl; }

  // Select from any cube which contains its isovert_coord.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    check_and_select_edge_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting from any cube." << endl; }

  // Select from cubes whose isovert coordinates are outside their cubes.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];

    compute_coord_mod<3>
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    check_and_select_edge_cube
      (scalar_grid, isovalue, isovert_param, gcube_index, COVERED_A_GCUBE, 
       isovert, selection_data);
  }
}

void MERGESHARP::select_sharp_isovert_mod3
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const int dimension = scalar_grid.Dimension();
  const int bin_width = isovert_param.bin_width;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
  std::vector<NUM_TYPE> sharp_gcube_list;
  SELECTION_DATA selection_data(scalar_grid, isovert_param);


  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  initialize_covered_by(isovert);

  MSDEBUG();
  flag_debug = false;
  if (flag_debug) {
    cerr << endl << "*** Selecting corner cubes." << endl;
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) {
    cerr << endl << "*** Selecting edge cubes (mod3)." 
         << endl;
  }

  // select edge cubes.
  select_edge_cubes_mod3
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Recomputing covered point positions." << endl; }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  // Select edge cubes (again)
  select_edge_cubes_mod3
    (scalar_grid, gradient_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  MSDEBUG();
  flag_debug = false;
}


// **************************************************
// SELECT MOD 6 ROUTINES
// **************************************************

/// Select edge cubes (eigenvalue 2) using mod 6 algorithm
void select_edge_cubes_mod6
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert,
 SELECTION_DATA_MOD6 & selection_data)
{
  // Select (0,0,0) mod 6
  select_edge_cubes_B
    (scalar_grid, gradient_grid, isovalue, isovert_param,
     selection_data.gcube_lists_mod6[0][0][0], isovert, selection_data);

  // Select (k,0,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.k_0_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_B
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_k_0_0(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_B
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_0(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,kc) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_kc_mod6_list.NumEntries(); 
       i++) {
    select_edge_cubes_B
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_kc(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) using mod 6 algorithm
/// Only select cubes whose isovert are within given distance to cube center.
void select_edge_cubes_within_dist_mod6
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const COORD_TYPE max_distance,
 ISOVERT & isovert,
 SELECTION_DATA_MOD6 & selection_data)
{
  // Select (0,0,0) mod 6
  select_edge_cubes_within_dist
    (scalar_grid, gradient_grid, isovalue, isovert_param,
     selection_data.gcube_lists_mod6[0][0][0], max_distance,
     isovert, selection_data);

  // Select (k,0,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.k_0_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_within_dist
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_k_0_0(i), max_distance,
       isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_within_dist
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_0(i), max_distance,
       isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,kc) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_kc_mod6_list.NumEntries(); 
       i++) {
    select_edge_cubes_within_dist
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_kc(i), max_distance,
       isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

// Select cubes near corner cubes.
void select_cubes_near_corners_mod6
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert,
 SELECTION_DATA_MOD6 & selection_data)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;


  // Select (0,0,0) mod 6
  select_edge_cubes_near_corners
    (scalar_grid, gradient_grid, isovalue, isovert_param,
     selection_data.gcube_lists_mod6[0][0][0], isovert, selection_data);

  // Select (k,0,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.k_0_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_near_corners
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_k_0_0(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,0) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_0_mod6_list.NumEntries(); i++) {
    select_edge_cubes_near_corners
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_0(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select (ka,kb,kc) mod 6
  for (NUM_TYPE i = 0; i < selection_data.ka_kb_kc_mod6_list.NumEntries(); 
       i++) {
    select_edge_cubes_near_corners
      (scalar_grid, gradient_grid, isovalue, isovert_param,
       selection_data.GCubeLists_ka_kb_kc(i), isovert, selection_data);
  }

  // Recomputing covered point positions
  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, isovalue, 
     isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);
}

void MERGESHARP::select_sharp_isovert_mod6
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const int dimension = scalar_grid.Dimension();
  const int bin_width = isovert_param.bin_width;
  const COORD_TYPE half_cube_distance = 0.51;
  COORD_TYPE max_distance;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
  std::vector<NUM_TYPE> sharp_gcube_list;
  SELECTION_DATA_MOD6 selection_data(scalar_grid, isovert_param);

  // Set selection_data.mismatch_table.
  selection_data.mismatch_table.AddSignedPermutations(3, 2, 2, 2);
  selection_data.mismatch_table.AddSignedPermutations(3, 2, 1, 2);
  selection_data.mismatch_table.AddSignedPermutations(3, 2, 0, 2);

  // *** TEST ***
  selection_data.mismatch_table.AddSignedPermutations(5, 3, 2, 4);

  selection_data.mismatch_table.ComputeAllBlockingLocations();

  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  initialize_covered_by(isovert);

  selection_data.AddEdgeCubesToGCubeLists(isovert.gcube_list);

  MSDEBUG();
  flag_debug = false;
  if (flag_debug) 
    { cerr << endl << "*** Selecting corner cubes." << endl; }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, 
     isovalue, isovert_param, isovert);
  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  set_mismatch_near_corner(isovert, sharp_gcube_list, selection_data);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Selecting cubes near corners (mod6)." << endl; }

  select_cubes_near_corners_mod6
    (scalar_grid, gradient_grid,isovalue, isovert_param, 
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "*** Recomputing mismatch grid. ***" << endl; }

  // Recompute mismatch grid so that selected near corner cubes can block
  //   effects of corner cubes.
  compute_mismatch_grid(isovert, sharp_gcube_list, selection_data);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Selecting edge cubes (mod6)."  << endl; }

  // select edge cubes.
  select_edge_cubes_within_dist_mod6
    (scalar_grid, gradient_grid, isovalue, isovert_param, half_cube_distance,
     isovert, selection_data);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Recomputing covered point positions." << endl; }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, 
     isovalue, isovert_param, isovert);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "*** Recomputing mismatch grid. ***" << endl; }

  // Recompute mismatch grid so that selected cubes can block effects
  //   of other cubes.
  compute_mismatch_grid(isovert, sharp_gcube_list, selection_data);

  // Select edge cubes (again). distance 0.6.
  max_distance = 0.6;
  select_edge_cubes_within_dist_mod6
    (scalar_grid, gradient_grid, isovalue, isovert_param, max_distance,
     isovert, selection_data);

  // Select edge cubes (again). distance 0.7.
  max_distance = 0.7;
  select_edge_cubes_within_dist_mod6
    (scalar_grid, gradient_grid, isovalue, isovert_param, max_distance,
     isovert, selection_data);

  // Select edge cubes (again). distance 0.8.
  max_distance = 0.8;
  select_edge_cubes_within_dist_mod6
    (scalar_grid, gradient_grid, isovalue, isovert_param, max_distance,
     isovert, selection_data);

  // Recompute mismatch grid one last time.
  compute_mismatch_grid(isovert, sharp_gcube_list, selection_data);

  // Select edge cubes (again).  (No distance restrictions.)
  select_edge_cubes_mod6
    (scalar_grid, gradient_grid, isovalue, isovert_param,
     isovert, selection_data);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, selection_data.covered_grid, 
     isovalue, isovert_param, isovert);

  reset_covered_isovert_positions(selection_data.covered_grid, isovert);

  // Select edge cubes (again)
  MSDEBUG();
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  // Select edge cubes.  Ignore mismatch restrictions.
  select_edge_cubes
    (scalar_grid, isovalue, isovert_param, sharp_gcube_list, 
     isovert, selection_data);

  MSDEBUG();
  flag_debug = false;
}


// **************************************************
// RESET ISOVERT POSITIONS
// **************************************************

namespace {

  // Change isovert positions which lie in covered cubes.
  void reset_covered_isovert_positions
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert)
  {
    std::vector<NUM_TYPE> gcube_sharp_list;

    get_corner_or_edge_cubes(isovert.gcube_list, gcube_sharp_list);

    for (NUM_TYPE i = 0; i < gcube_sharp_list.size(); i++) {

      NUM_TYPE gcubeA_index = gcube_sharp_list[i];

      if (isovert.gcube_list[gcubeA_index].flag_conflict) {

        VERTEX_INDEX conflicting_cube = 
          isovert.gcube_list[gcubeA_index].cube_containing_isovert;
        INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(conflicting_cube);

        if (gcubeB_index != ISOVERT::NO_INDEX) {

          if (isovert.gcube_list[gcubeB_index].flag != SELECTED_GCUBE) {

            VERTEX_INDEX cubeB_index = isovert.CubeIndex(gcubeB_index);
            VERTEX_INDEX cubeC_index =
              isovert.gcube_list[gcubeB_index].cube_containing_isovert;

            if (cubeC_index != cubeB_index) {


              if (covered_grid.Scalar(cubeC_index)) {

                MSDEBUG();
                if (flag_debug) {
                  using namespace std;
                  VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
                  cerr << "Copying position from " << cubeA_index << " ";
                  ijkgrid_output_vertex_coord(cerr, isovert.grid, cubeA_index);
                  cerr << " to " << cubeB_index << " ";
                  ijkgrid_output_vertex_coord(cerr, isovert.grid, cubeB_index);
                  cerr << endl;
                  cerr << "  Old position: ";
                  IJK::print_coord3D(cerr, isovert.gcube_list[gcubeB_index].isovert_coord);
                  cerr << ".  New position: ";
                  IJK::print_coord3D(cerr, isovert.gcube_list[gcubeA_index].isovert_coord);
                  cerr << endl;
                }

                copy_isovert_position(isovert.grid, gcubeA_index, gcubeB_index, isovert);
              }
            }
          }
        }
        else {
          // Should never happen, but...
          IJK::PROCEDURE_ERROR error("reset_covered_isovert_positions");
          error.AddMessage("Programming error.  Illegal conflicting cube index: ",
                           conflicting_cube, " for cube ",
                           isovert.CubeIndex(gcubeA_index), ".");
          throw error;

        }
      }
    }

  }

}


// **************************************************
// Local routines
// **************************************************

namespace {

  /// Return true if sharp vertex is in a cube which is covered.
  bool check_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube0_index)
  {
    const VERTEX_INDEX cube1_index =
      isovert.gcube_list[gcube0_index].cube_containing_isovert;
  
    if (covered_grid.Scalar(cube1_index)) { return true; }

    // *** IS THIS NECESSARY? ***
    // *** SHOULD USE IsCoveredOrSelected() ***
    // if the sharp vertex is not in the same cube. 
    if (isovert.GCubeIndex(cube1_index) != ISOVERT::NO_INDEX) {
      if (isovert.gcube_list[gcube0_index].cube_index != cube1_index){
        if (isovert.isFlag(cube1_index, COVERED_A_GCUBE))
          { return true; }
      }
    }

    return false;
  }

  /// If sharp vertex is in covered cube, set cube to COVERED_POINT.
  void check_and_set_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube_index)
  {
    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE) {

      if (isovert.gcube_list[gcube_index].cube_index !=
          isovert.gcube_list[gcube_index].cube_containing_isovert) {

        if (check_covered_point
            (covered_grid, isovert, gcube_index)) {

          isovert.gcube_list[gcube_index].flag = COVERED_POINT;
        }
      }
    }
  }

  // Return true if a neighbor of c has given flag.
  bool is_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag)
  {
    VERTEX_INDEX neighbor_cube_index;

    // Check facet neighbors
    for (NUM_TYPE j = 0; j < isovert.grid.NumCubeNeighborsF(); j++) {

      // *** SHOULD CHECK boundary_bits
      neighbor_cube_index = isovert.grid.CubeNeighborF(c.cube_index, j);
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.index_grid.Scalar(neighbor_cube_index);

      if(neighbor_gcube_index == ISOVERT::NO_INDEX) { continue; }

      if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
        { return true; }
    }

    // Check edge neighbors
    for (NUM_TYPE j = 0; j < isovert.grid.NumCubeNeighborsE(); j++) {

      neighbor_cube_index = isovert.grid.CubeNeighborE(c.cube_index, j);
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.index_grid.Scalar(neighbor_cube_index);

      if(neighbor_gcube_index == ISOVERT::NO_INDEX) { continue; }

      if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
        { return true; }
    }

    // Check vertex neighbors
    for (NUM_TYPE j = 0; j < isovert.grid.NumCubeNeighborsV(); j++) {
      neighbor_cube_index = isovert.grid.CubeNeighborV(c.cube_index, j);

      // *** SHOULD CHECK boundary_bits

      // neighbor_gcube_index is an entry into the gcube list 
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.index_grid.Scalar(neighbor_cube_index);

      if (neighbor_gcube_index == ISOVERT::NO_INDEX)
        { continue; }
      if (isovert.gcube_list[neighbor_gcube_index].flag == flag)
        { return true; }
    }

    return false;
  }

  // Return true if cube at (c0,c1,c2) +/- (1,0,0) or (0,1,0) or (0,0,1)
  // @param[out] orth_dir Direction orthogonal to facet separating neighbor.
  // @param[out] side Side of facet separating neighbor.
  bool is_facet_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag, int & orth_dir, int & side)
  {
    const VERTEX_INDEX cubeA_index = c.cube_index;

    if (c.boundary_bits == 0) {

      for (orth_dir = 0; orth_dir < DIM3; orth_dir++) {
        for (side = 0; side < 2; side++) {

          VERTEX_INDEX cubeB_index = 
            isovert.grid.AdjacentVertex(cubeA_index, orth_dir, side);

          INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

          if (gcubeB_index != ISOVERT::NO_INDEX) {
            if (isovert.gcube_list[gcubeB_index].flag == flag)
              { return true; }
          }
        }
      }
    }

    return(false);
  }

  // Return true if cube facet lies on a 3x3x3 region around a corner cube.
  bool is_facet_on_corner_3x3x3
  (const GRID_CUBE_DATA & c, const SHARPISO_BOOL_GRID & covered_grid,
   const ISOVERT & isovert, int & orth_dir, int & side,
   VERTEX_INDEX & corner_cube_index)
  {
    const VERTEX_INDEX cubeA_index = c.cube_index;
    GRID_COORD_TYPE distance2boundary;

    corner_cube_index = 0;

    isovert.grid.ComputeCubeDistanceToGridBoundary
      (c.cube_coord, distance2boundary);

    if (distance2boundary >= 2) {

      for (orth_dir = 0; orth_dir < DIM3; orth_dir++) {
        for (side = 0; side < 2; side++) {

          VERTEX_INDEX cubeB_index = 
            isovert.grid.AdjacentVertex(cubeA_index, orth_dir, side);

          if (!covered_grid.Scalar(cubeB_index)) { continue; }

          cubeB_index = isovert.grid.AdjacentVertex
            (cubeB_index, orth_dir, side);

          int d1 = (orth_dir+1)%DIM3;
          int d2 = (orth_dir+2)%DIM3;

          for (int k1 = -1; k1 < 2; k1++) {
            for (int k2 = -1; k2 < 2; k2++) {

              const VERTEX_INDEX cubeC_index = 
                cubeB_index + isovert.grid.AxisIncrement(d1)*k1
                + isovert.grid.AxisIncrement(d2)*k2;

              const INDEX_DIFF_TYPE gcubeC_index = 
                isovert.GCubeIndex(cubeC_index);

              if (gcubeC_index == ISOVERT::NO_INDEX) { continue; }

              if (isovert.gcube_list[gcubeC_index].flag == SELECTED_GCUBE &&
                  isovert.NumEigenvalues(gcubeC_index) == 3) {
                corner_cube_index = cubeC_index;
                return(true);
              }
            }
          }
        }
      }
    }

    return(false);
  }

  // Return true if cube at (c0,c1,c2) +/- (a,0,0) or (0,a,0) or (0,0,a)
  //   for a = 1 or 2 has given flag.
  bool is_dist2_neighbor
  (const GRID_CUBE_DATA & c, const ISOVERT & isovert, 
   const GRID_CUBE_FLAG flag)
  {
    const GRID_COORD_TYPE max_dist = 2;
    const VERTEX_INDEX cubeA_index = c.cube_index;
    GRID_COORD_TYPE distance2boundary;

    isovert.grid.ComputeCubeDistanceToGridBoundary
      (c.cube_coord, distance2boundary);

    if (distance2boundary >= max_dist) {
      for (int d = 0; d < DIM3; d++) {
        for (int j = 0; j < 2; j++) {

          VERTEX_INDEX cubeB_index = cubeA_index;
          for (int i = 0; i < max_dist; i++) {

            cubeB_index = isovert.grid.AdjacentVertex(cubeB_index, d, j);

            INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

            if (gcubeB_index != ISOVERT::NO_INDEX) {
              if (isovert.gcube_list[gcubeB_index].flag == flag)
                { return true; }
            }
          }
        }
      }
    }

    return(false);
  }

  /// Return true if some facet neighbor is covered by a sharp cube.
  /// Version based on covered_grid.
  /// @pre cube0 is an internal cube.
  bool is_facet_neighbor_covered_B
  (const SHARPISO_BOOL_GRID & covered_grid, const VERTEX_INDEX cube_index0)
  {
    for (int d = 0; d < DIM3; d++) {
      for (int iside = 0; iside < 2; iside++) {
        VERTEX_INDEX cube_index1 = 
          covered_grid.AdjacentVertex(cube_index0, d, iside);

        if (covered_grid.Scalar(cube_index1)) { return(true); }
      }
    }
    return(false);
  }

  // *** UNUSED ***
  /// Return true if some edge neighbor is covered by a sharp cube.
  /// Version based on covered_grid.
  /// @pre cube0 is an internal cube.
  bool is_edge_neighbor_covered_B
  (const SHARPISO_BOOL_GRID & covered_grid, const VERTEX_INDEX cube_index0)
  {
    for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {
      int d1 = (edge_dir+1)%DIM3;
      int d2 = (edge_dir+2)%DIM3;

      for (int d1 = 0; d1 < DIM3; d1++) {
        for (int j1 = 0; j1 < 2; j1++) {
          for (int d2 = 0; d2 < DIM3; d2++) {
            for (int j2 = 0; j2 < 2; j2++) {

              VERTEX_INDEX cube_index1 = 
                cube_index0 + (2*j1-1)*covered_grid.AxisIncrement(d1)
                + (2*j2-1)*covered_grid.AxisIncrement(d2);

              if (covered_grid.Scalar(cube_index1)) { 
                return(true); 
              }
            }
          }
        }
      }
    }
    return(false);
  }

  // Return index of coordinate with smallest absolute value.
  int get_index_smallest_abs(const COORD_TYPE coord[DIM3])
  {
    int index = 0;
    for (int d = 1; d < DIM3; d++) {
      if (abs(coord[d]) < abs(coord[index]))
        { index = d; }
    }

    return(index);
  }

}


// **************************************************
// MISMATCH ROUTINES
// **************************************************

namespace {

  // Return true if 0 <= c_d+k < axis_d_size.
  inline bool is_cube_coord_d_in_grid
  (const AXIS_SIZE_TYPE axis_d_size, const GRID_COORD_TYPE c_d, const int k)
  {    
    if (c_d + k < 0) { return(false); }
    if (c_d + k + 1 >= axis_d_size) { return(false); }

    return(true);
  }

  // Check if (cube_coord[0]+k0, cube_coord[1]+k1, cube_coord[2]+k2)
  //   is in grid.
  bool is_cube_coord_in_grid
  (const SHARPISO_GRID & grid, const GRID_COORD_TYPE cube_coord[DIM3],
   const int k0, const int k1, const int k2)
  {
    if (!is_cube_coord_d_in_grid(grid.AxisSize(0), cube_coord[0], k0))
      { return(false); }
    if (!is_cube_coord_d_in_grid(grid.AxisSize(1), cube_coord[1], k1))
      { return(false); }
    if (!is_cube_coord_d_in_grid(grid.AxisSize(2), cube_coord[2], k2))
      { return(false); }

    return(true);
  }

  // Check if (cube_coord[d]+offset_coord[d]) is in grid.
  bool is_cube_coord_in_grid
  (const SHARPISO_GRID & grid, const GRID_COORD_TYPE cube_coord[DIM3],
   const GRID_COORD_DIFF_TYPE offset_coord[DIM3])
  {
    return(is_cube_coord_in_grid
           (grid, cube_coord, offset_coord[0], offset_coord[1], 
            offset_coord[2]));
  }

  // Returns absolute value and sign.
  // Note: Returns sign of +1 if x = 0.
  template <typename T, typename Tabs, typename Tsign>
  void get_abs_and_sign(const T x, Tabs & x_abs, Tsign & x_sign)
  {
    if (x >= 0) { 
      x_abs = x;
      x_sign = 1;
    }
    else {
      x_abs = -x;
      x_sign = -1;
    }
  }

}


// **************************************************
// Class MISMATCH_TABLE member functions
// **************************************************


MISMATCH_TABLE::MISMATCH_TABLE(const SHARPISO_GRID & sharpiso_grid)
{
  grid.SetSize(sharpiso_grid);
}

NUM_TYPE MISMATCH_TABLE::AddEntry
(const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
 const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance)
{
  const NUM_TYPE ientry = entry.size();
  PROCEDURE_ERROR error("MISMATCH_TABLE::AddEntry");

  entry.resize(ientry+1);

  entry[ientry].offset_coord[0] = c0;
  entry[ientry].offset_coord[1] = c1;
  entry[ientry].offset_coord[2] = c2;

  for (int d = 0; d < DIM3; d++) {

    if (entry[ientry].offset_coord[d] > MaxOffsetCoord() ||
        -entry[ientry].offset_coord[d] > MaxOffsetCoord()) {

      error.AddMessage("Programming error.  Offset coordinate too large.");
      error.AddMessage("  Offset coordinate = ", 
                       entry[ientry].offset_coord[d], ".");
      error.AddMessage("  Maximum allowable value = ", 
                       MaxOffsetCoord(), ".");
      throw error;
    }
  }

  entry[ientry].offset_increment =
    grid.AxisIncrement(0)*c0 + grid.AxisIncrement(1)*c1 + 
    grid.AxisIncrement(2)*c2;

  entry[ientry].max_distance = max_distance;

  return(ientry);
}

void MISMATCH_TABLE::AddBlockingLocation
(const NUM_TYPE ientry, const GRID_COORD_DIFF_TYPE c0, 
 const GRID_COORD_DIFF_TYPE c1, const GRID_COORD_DIFF_TYPE c2)
{
  INDEX_DIFF_TYPE blocking_location_increment =
    grid.AxisIncrement(0)*c0 + grid.AxisIncrement(1)*c1 +
    grid.AxisIncrement(2)*c2;

  entry[ientry].blocking_location_increment.push_back
    (blocking_location_increment);
}

void MISMATCH_TABLE::ComputeBlockingLocations(const NUM_TYPE ientry)
{
  GRID_COORD_TYPE abs_c[DIM3];
  int sign_c[DIM3];

  for (int d = 0; d < DIM3; d++) 
    { get_abs_and_sign(entry[ientry].offset_coord[d], abs_c[d], sign_c[d]); }

  for (GRID_COORD_TYPE i0 = 0; i0 <= abs_c[0]; i0++) {

    if (abs_c[0] >= 3 &&
        (i0 == 0 || i0 == abs_c[0])) { continue; }

    for (GRID_COORD_TYPE i1 = 0; i1 <= abs_c[1]; i1++) {

      if (abs_c[1] >= 3 &&
          (i1 == 0 || i1 == abs_c[1])) { continue; }

      for (GRID_COORD_TYPE i2 = 0; i2 <= abs_c[2]; i2++) {

        if (abs_c[2] >= 3 &&
            (i2 == 0 || i2 == abs_c[2])) { continue; }

        if (i0 < 2 && i1 < 2 && i2 < 2) {
          // Cube is adjacent to 0,0,0
          continue;
        }

        if (i0 > abs_c[0]-2 && i1 > abs_c[1]-2 && i2 > abs_c[2]-2) {
          // Cube is adjacent to (abs_c[0],abs_c[1],abs_c[2])
          continue;
        }

        AddBlockingLocation(ientry, i0*sign_c[0], i1*sign_c[1], i2*sign_c[2]);
      }
    }
  }

}

void MISMATCH_TABLE::ComputeAllBlockingLocations()
{
  for (NUM_TYPE ientry = 0; ientry < entry.size(); ientry++)
    { ComputeBlockingLocations(ientry); }
}

// Add all permutations of c0, c1 and c2.
void MISMATCH_TABLE::AddPermutations
(const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
 const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance)
{
  AddEntry(c0, c1, c2, max_distance);
  if (c1 != c2) {
    AddEntry(c0, c2, c1, max_distance);
  }

  if (c0 != c1) {
    AddEntry(c1, c0, c2, max_distance);

    if (c0 != c2) {
      AddEntry(c1, c2, c0, max_distance);
    }
  }

  if (c0 != c2 && c1 != c2) {
    AddEntry(c2, c0, c1, max_distance);

    if (c0 != c1) {
      AddEntry(c2, c1, c0, max_distance);
    }
  }
}


// Add all permutations of +/- c0, +/- c1 and +/- c2.
void MISMATCH_TABLE::AddSignedPermutations
(const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
 const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance)
{
  AddPermutations(c0, c1, c2, max_distance);
  if (c0 != 0) { AddPermutations(-c0, c1, c2, max_distance); }
  if (c1 != 0) { AddPermutations(c0, -c1, c2, max_distance); }
  if (c2 != 0) { AddPermutations(c0, c1, -c2, max_distance); }
  if (c0 != 0 && c1 != 0) { AddPermutations(-c0, -c1, c2, max_distance); }
  if (c0 != 0 && c2 != 0) { AddPermutations(-c0, c1, -c2, max_distance); }
  if (c1 != 0 && c2 != 0) { AddPermutations(c0, -c1, -c2, max_distance); }
  if (c0 != 0 && c1 != 0 && c2 != 0) 
    { AddPermutations(-c0, -c1, -c2, max_distance); }
}


bool MISMATCH_TABLE::IsBlocked
(const NUM_TYPE ientry, const VERTEX_INDEX cubeA_index,
 const SHARPISO_BOOL_GRID & selected_grid) const
{
  for (NUM_TYPE j = 0; j < entry[ientry].NumBlockingLocations(); j++) {

    VERTEX_INDEX cubeB_index = 
      cubeA_index + entry[ientry].blocking_location_increment[j];

    if (selected_grid.Scalar(cubeB_index)) { return(true); }
  }

  return(false);
}


void MISMATCH_TABLE::SetMismatchGrid
(const GRID_CUBE_DATA & gcubeA,
 const ISOVERT & isovert,
 const SHARPISO_BOOL_GRID & selected_grid,
 SHARPISO_BOOL_GRID & mismatch_grid) const
{
  const VERTEX_INDEX cubeA_index = gcubeA.cube_index;
  GRID_COORD_TYPE distance2boundary;

  mismatch_grid.ComputeCubeDistanceToGridBoundary
    (gcubeA.cube_coord, distance2boundary);

  if (distance2boundary >= MaxOffsetCoord()) {

    for (NUM_TYPE ientry = 0; ientry < entry.size(); ientry++) {
      VERTEX_INDEX cubeB_index = cubeA_index+entry[ientry].offset_increment;

      if (!IsBlocked(ientry, cubeA_index, selected_grid)) {

        const COORD_TYPE max_distance = entry[ientry].max_distance;

        if ((gcubeA.num_eigenvalues == 3) ||
            does_sharp_edge_point_to_cube
            (grid, isovert, cubeA_index, cubeB_index, max_distance)) {

          // *** DEBUG ***
          MSDEBUG();
          if (flag_debug) {
            if (!mismatch_grid.Scalar(cubeB_index)) {
              isovert.grid.PrintIndexAndCoord
                (cerr, "  Cube ", cubeB_index,
                 " does not match cube ", cubeA_index, "\n");
            }
          }

          mismatch_grid.Set(cubeB_index, true); 
        }
      }
    }
  }
  else {
    for (NUM_TYPE ientry = 0; ientry < entry.size(); ientry++) {

      if (is_cube_coord_in_grid
          (grid, gcubeA.cube_coord, entry[ientry].offset_coord)) {
        VERTEX_INDEX cubeB_index = cubeA_index+entry[ientry].offset_increment;

        if (!IsBlocked(ientry, cubeA_index, selected_grid)) {

          const COORD_TYPE max_distance = entry[ientry].max_distance;

          if (does_sharp_edge_point_to_cube
              (grid, isovert, cubeA_index, cubeB_index, max_distance)) {
            mismatch_grid.Set(cubeB_index, true); 
          }
        }
      }
    }
  }
}


// **************************************************
// Class MOD6_LIST member functions
// **************************************************

void MOD6_LIST::AddEntry(const int ka, const int kb, const int kc)
{
  NUM_TYPE ientry = entry.size();

  entry.resize(ientry+1);

  entry[ientry].coord[0] = ka;
  entry[ientry].coord[1] = kb;
  entry[ientry].coord[2] = kc;
}

void MOD6_LIST::AddPermutations(const int ka, const int kb, const int kc)
{
  AddEntry(ka, kb, kc);
  if (kb != kc) {
    AddEntry(ka, kc, kb);
  }

  if (ka != kb) {
    AddEntry(kb, ka, kc);

    if (ka != kc) {
      AddEntry(kb, kc, ka);
    }
  }

  if (ka != kc && kb != kc) {
    AddEntry(kc, ka, kb);

    if (ka != kb) {
      AddEntry(kc, kb, ka);
    }
  }
}

namespace {

  inline bool not_0_3(const int k)
  {
    if (k == 0 || k == 3) { return(false); }
    else { return(true); }
  }

  inline bool not_0_3(const int ka, const int kb)
  {
    if (not_0_3(ka) && not_0_3(kb)) { return(true); }
    else { return(false); }
  }

  inline bool not_0_3(const int ka, const int kb, const int kc)
  {
    if (not_0_3(ka) && not_0_3(kb) && not_0_3(kc)) { return(true); }
    else { return(false); }
  }

}

void MOD6_LIST::AddSignedPermutations
(const int ka, const int kb, const int kc)
{
  AddPermutations(ka, kb, kc);
  if (not_0_3(ka)) { AddPermutations(6-ka, kb, kc); }
  if (not_0_3(kb) && ka != kb) { AddPermutations(ka, 6-kb, kc); }
  if (not_0_3(kc) && kc != ka && kc != kb) 
    { AddPermutations(ka, kb, 6-kc); }
  if (not_0_3(ka,kb)) { AddPermutations(6-ka, 6-kb, kc); }
  if (not_0_3(ka,kc) && kc != kb) { AddPermutations(6-ka, kb, 6-kc); }
  if (not_0_3(kb,kc) && kb != ka && kc != ka) 
    { AddPermutations(ka, 6-kb, 6-kc); }
  if (not_0_3(ka,kb,kc)) { AddPermutations(6-ka, 6-kb, 6-kc); }
}


// **************************************************
// Class SELECTION_DATA member functions
// **************************************************

// Constructor
SELECTION_DATA::SELECTION_DATA
(const SHARPISO_GRID & grid, const SHARP_ISOVERT_PARAM & isovert_param):
  mismatch_table(grid)
{
  // Initialize bin grid.
  bin_width = isovert_param.bin_width;
  init_bin_grid(grid, bin_width, bin_grid);

  // Initialize covered grid.
  covered_grid.SetSize(grid);
  covered_grid.SetAll(false);

  // Initialize selected grid.
  selected_grid.SetSize(grid);
  selected_grid.SetAll(false);

  // Initialize mismatch grid.
  mismatch_grid.SetSize(grid);
  mismatch_grid.SetAll(false);
}

// Constructor
SELECTION_DATA_MOD6::SELECTION_DATA_MOD6
(const SHARPISO_GRID & grid, const SHARP_ISOVERT_PARAM & isovert_param):
  SELECTION_DATA(grid, isovert_param)
{
  SetMod6Lists();
}

// Add gcube_index to gcube_list_mod6[][][]/
void SELECTION_DATA_MOD6::AddCube
(const GRID_COORD_TYPE coord_mod6[DIM3], const VERTEX_INDEX gcube_index)
{
  const GRID_COORD_TYPE x0 = coord_mod6[0];
  const GRID_COORD_TYPE x1 = coord_mod6[1];
  const GRID_COORD_TYPE x2 = coord_mod6[2];

  gcube_lists_mod6[x0][x1][x2].push_back(gcube_index);
}

// Set cube_list_mod6[][][]
void SELECTION_DATA_MOD6::AddEdgeCubesToGCubeLists
(const std::vector<GRID_CUBE_DATA> & gcube_list)
{
  GRID_COORD_TYPE coord_mod6[DIM3];

  for (NUM_TYPE i=0; i<gcube_list.size(); i++) {

    if (gcube_list[i].num_eigenvalues == 2 &&
        !gcube_list[i].flag_centroid_location) {

      compute_coord_mod<6>(gcube_list[i].cube_coord, coord_mod6);
      AddCube(coord_mod6, i);
    }
  }
}

// Set up mod 6 lists.
void SELECTION_DATA_MOD6::SetMod6Lists()
{
  // Set up k_0_0_mod6_list.
  k_0_0_mod6_list.AddSignedPermutations(1, 0, 0);
  k_0_0_mod6_list.AddSignedPermutations(2, 0, 0);
  k_0_0_mod6_list.AddSignedPermutations(3, 0, 0);

  // Set up ka_kb_0_mod6_list
  ka_kb_0_mod6_list.AddSignedPermutations(1, 1, 0);
  ka_kb_0_mod6_list.AddSignedPermutations(1, 2, 0);
  ka_kb_0_mod6_list.AddSignedPermutations(1, 3, 0);
  ka_kb_0_mod6_list.AddSignedPermutations(2, 2, 0);
  ka_kb_0_mod6_list.AddSignedPermutations(2, 3, 0);
  ka_kb_0_mod6_list.AddSignedPermutations(3, 3, 0);

  // Set up ka_kb_kc_mod6_list
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 1, 1);
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 1, 2);
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 1, 3);
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 2, 2);
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 2, 3);
  ka_kb_kc_mod6_list.AddSignedPermutations(1, 3, 3);
  ka_kb_kc_mod6_list.AddSignedPermutations(2, 2, 2);
  ka_kb_kc_mod6_list.AddSignedPermutations(2, 2, 3);
  ka_kb_kc_mod6_list.AddSignedPermutations(2, 3, 3);
  ka_kb_kc_mod6_list.AddSignedPermutations(3, 3, 3);
}

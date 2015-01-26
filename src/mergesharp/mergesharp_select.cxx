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
  (const SHARPISO_GRID_NEIGHBORS & grid, 
   const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert);

  bool check_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube0_index);

  void check_and_set_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube_index);

  bool is_neighbor
  (const GRID_CUBE_DATA & c,
   const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   SHARPISO_GRID_NEIGHBORS & grid,
   ISOVERT & isovert, GRID_CUBE_FLAG flag);

  bool is_facet_neighbor_covered_B
  (const SHARPISO_BOOL_GRID & covered_grid, const VERTEX_INDEX cube_index0);

  int get_index_smallest_abs(const COORD_TYPE coord[DIM3]);

};


// **************************************************
// SELECT ROUTINES
// **************************************************

void select_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert
 )
{
  const int bin_width = isovert_param.bin_width;
  const VERTEX_INDEX cube_index = 
    isovert.gcube_list[gcube_index].cube_index;

  isovert.gcube_list[gcube_index].flag = SELECTED_GCUBE;

  covered_grid.Set(cube_index, true);

  MSDEBUG();
  if (flag_debug) {
    using namespace std;
    cerr << "*** Selecting cube " << cube_index << "  ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << "  Linf_dist: " << isovert.gcube_list[gcube_index].linf_dist;
    cerr << endl;
  }

  bin_grid_insert(scalar_grid, bin_width, cube_index, bin_grid);

  // mark all the neighbors as covered
  for (int i=0;i < gridn.NumVertexNeighborsC(); i++) {
    VERTEX_INDEX cube_index2 = gridn.VertexNeighborC(cube_index, i);

    covered_grid.Set(cube_index2, true);

    NUM_TYPE gcube_index2 = isovert.sharp_ind_grid.Scalar(cube_index2);
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
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;
  VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;
  VERTEX_INDEX v1, v2;

  if (isovert.gcube_list[gcube_index].linf_dist >= linf_dist_threshold)
    { return; }

  // Check if the sharp vertex is inside a covered cube.
  if (check_covered_point(covered_grid, isovert, gcube_index)) {
    isovert.gcube_list[gcube_index].flag = COVERED_POINT;
    return;
  }

  bool triangle_flag =
    creates_triangle_new(scalar_grid, isovert, cube_index,
                         isovalue, bin_grid, bin_width, v1, v2);

  if (!triangle_flag) {
    select_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
       isovert_param, gcube_index, flag, isovert);
  }
  else
	{
      isovert.gcube_list[gcube_index].flag = UNAVAILABLE_GCUBE;
	}

}


/// Unselect cube grid_cube.  Uncover all neighbors.
/// @pre Neigboring cubes are covered only grid_cube.
void unselect_cube
(const SHARPISO_GRID_NEIGHBORS & grid,
 const AXIS_SIZE_TYPE bin_width,
 const NUM_TYPE gcube_index,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<int> & bin_grid,
 ISOVERT & isovert)
{
  VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE;
  bin_grid_remove(grid, bin_width, cube_index, bin_grid);

  if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < grid.NumVertexNeighborsC(); k++) {

      VERTEX_INDEX cubeB_index = grid.VertexNeighborC(cube_index, k);
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

      covered_grid.Set(cubeB_index, false);
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
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sortd_ind2gcube_list,
 ISOVERT &isovert)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;
  const int bin_width = isovert_param.bin_width;

  for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    const NUM_TYPE gcube_index = sortd_ind2gcube_list[ind];
    const BOUNDARY_BITS_TYPE boundary_bits =
      isovert.gcube_list[gcube_index].boundary_bits;

    // check boundary
    if (boundary_bits == 0)

      // select corners first
      if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE &&
          isovert.gcube_list[gcube_index].num_eigenvalues > 2)
        {
          check_and_select_cube
            (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_CORNER_GCUBE, isovert);
        }
  }

}

// Select near corners. neighbor of covered corner
void select_cubes_near_corners
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID &covered_grid,
 BIN_GRID<VERTEX_INDEX> &bin_grid,
 SHARPISO_GRID_NEIGHBORS &gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 vector<NUM_TYPE> sharp_gcube_list,
 ISOVERT & isovert)
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
          (isovert.gcube_list[gcube_index], scalar_grid, gridn,  
           isovert,  COVERED_CORNER_GCUBE );

        if (flag) {
          check_and_select_cube
            (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_A_GCUBE, isovert );
        }

        if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) 
          { isovert.gcube_list[gcube_index].flag_near_corner = true; }
      }
    }
  }
}

void check_and_select_edge_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const NUM_TYPE gcube_index,
 const GRID_CUBE_FLAG flag,
 ISOVERT & isovert)
{
  const BOUNDARY_BITS_TYPE boundary_bits =
    isovert.gcube_list[gcube_index].boundary_bits;

  check_and_set_covered_point(covered_grid, isovert, gcube_index);

  // check boundary
  if (boundary_bits == 0) {

    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE &&
        isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
        !isovert.gcube_list[gcube_index].flag_conflict)
      {
        check_and_select_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
  }
}

/// Select edge cubes (eigenvalue 2)
void select_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    check_and_select_edge_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
       isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
  }
}

/// Select one edge cube (eigenvalue 2).
/// @param from_list List of gcube indices sorted by increasing
///    distance from isovert_coord to cube centers.
///    Select cube from list from_list.
bool select_one_edge_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & from_list,
 INDEX_DIFF_TYPE & selected_gcube_index,
 ISOVERT & isovert)
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
            (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_A_GCUBE, isovert);

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
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX gcube_index,
 const std::vector<NUM_TYPE> & from_list,
 ISOVERT & isovert)
{
  const int bin_width = isovert_param.bin_width;
  vector<NUM_TYPE> neighbor_list;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);

  get_corner_or_edge_cubes_around_cube
    (gridn, isovert, gcube_index, neighbor_list);
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

  unselect_cube
    (gridn, bin_width, gcube_index, covered_grid, bin_grid, isovert);

  INDEX_DIFF_TYPE gcubeB_index, gcubeC_index;
  if (select_one_edge_cube(scalar_grid, covered_grid, bin_grid, gridn, isovalue,
                           isovert_param, from_list, gcubeB_index, isovert)) {

    bool flag = select_one_edge_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
       isovert_param, neighbor_list, gcubeC_index, isovert);
  }
  else {

    // Reselect cube gcube_index
    select_cube(scalar_grid, covered_grid, bin_grid, gridn, isovalue,
                isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
  }
}

/// Reselect edge cubes around gcube_index.
/// @param nearby_sharp_list List of nearby sharp cubes.
void reselect_edge_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX gcube_index,
 const std::vector<VERTEX_INDEX> & nearby_selected_list,
 ISOVERT & isovert)
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
      ijkgrid_output_vertex_coord(cerr, gridn, max_overlap_cube);
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

      for (NUM_TYPE k = 0; k < gridn.NumVertexNeighborsC(); k++) {

        VERTEX_INDEX cubeC_index = gridn.VertexNeighborC(cube_index, k);
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
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, from_list, isovert);
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
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  const int bin_width = isovert_param.bin_width;
  std::vector<VERTEX_INDEX> nearby_selected_list;

  for (NUM_TYPE i = 0; i < sharp_gcube_list.size(); i++) {

    NUM_TYPE gcube_index = sharp_gcube_list[i];
    if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE &&
        isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
        !isovert.gcube_list[gcube_index].flag_near_corner) {

      VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;

      get_selected(scalar_grid, cube_index, bin_grid, bin_width, 
                   nearby_selected_list);

      reselect_edge_cubes
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
         isovert_param, gcube_index, nearby_selected_list, isovert);
    }
  }
}

// Set gcube_list[i].covered_by to gcube_list[i].cube_index for each i.
void initialize_covered_by(ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) 
    { isovert.gcube_list[i].covered_by = isovert.gcube_list[i].cube_index; }
}

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

  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  BIN_GRID<VERTEX_INDEX> bin_grid;
  init_bin_grid(scalar_grid, bin_width, bin_grid);

  initialize_covered_by(isovert);

  SHARPISO_GRID_NEIGHBORS gridn;
  gridn.SetSize(scalar_grid);

  //Keeps track of all cubes which are covered.
  SHARPISO_BOOL_GRID covered_grid;
  covered_grid.SetSize(scalar_grid);
  covered_grid.SetAll(false);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting corner cubes." << endl; 
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting near corner cubes." << endl; 
  }

  // select cubes near corners
  select_cubes_near_corners
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  MSDEBUG();
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

  // select edge cubes.
  select_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param,
     sharp_gcube_list, isovert);

  // reselect edge cubes
  reselect_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);
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
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
  std::vector<NUM_TYPE> sharp_gcube_list;

  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  BIN_GRID<VERTEX_INDEX> bin_grid;
  init_bin_grid(scalar_grid, bin_width, bin_grid);

  initialize_covered_by(isovert);

  SHARPISO_GRID_NEIGHBORS gridn;
  gridn.SetSize(scalar_grid);

  //Keeps track of all cubes which are covered.
  SHARPISO_BOOL_GRID covered_grid;
  covered_grid.SetSize(scalar_grid);
  covered_grid.SetAll(false);

  MSDEBUG();
  flag_debug = false;
  if (flag_debug) {
    cerr << endl << "*** Selecting corner cubes." << endl;
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { 
    cerr << endl << "*** Selecting near corner cubes." << endl; 
  }

  // select cubes near corners
  select_cubes_near_corners
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

  // select edge cubes.
  select_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param,
     sharp_gcube_list, isovert);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Recomputing covered point positions." << endl; }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Reselecting edge cubes." << endl; }

  // reselect edge cubes
  reselect_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  MSDEBUG();
  if (flag_debug) {
    cerr << endl << "*** Selecting edge cubes (again.)" << endl; 
  }

  // Retry selecting edge cubes.
  select_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param,
     sharp_gcube_list, isovert);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Recomputing covered point positions." << endl; }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  flag_debug = false;
}


// **************************************************
// SELECT MOD3 ROUTINES
// **************************************************

/// Compute (coord[] mod 3)
inline void compute_coord_mod3
(const GRID_COORD_TYPE coord[DIM3], GRID_COORD_TYPE coord_mod3[DIM3])
{
  coord_mod3[0] = coord[0]%3;
  coord_mod3[1] = coord[1]%3;
  coord_mod3[2] = coord[2]%3;
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
  GRID_COORD_TYPE cubeL_coord[DIM3];   // Lower cube coord
  GRID_COORD_TYPE cubeH_coord[DIM3];   // Higher cube coord
  GRID_COORD_TYPE x[DIM3], diff[DIM3];
  bool flag_zero_mod3;                 // true, if some coordinate is 0 mod 3.
  int mod3_coord;

  flag_zero_mod3 = false;
  for (int d = 0; d < DIM3; d++) {
    cubeL_coord[d] = (cube0_coord[d]-(int(cube0_coord[d])%3));
    if (cube0_coord[d] == cubeL_coord[d]) {

      if (flag_zero_mod3) {
        error.AddMessage
          ("Programming error.  More than one coord is 0 mod 3.");
        throw error;
      }

      flag_zero_mod3 = true;
      mod3_coord = 0;
      cubeH_coord[d] = cubeL_coord[d];
    }
    else if (cube0_coord[d] + 5 >= grid.AxisSize(d)) {
      cubeH_coord[d] = cubeL_coord[d];
    }
    else {
      cubeH_coord[d] = cubeL_coord[d]+3;
    }
  }

  if (!flag_zero_mod3) {

    MSDEBUG();
    grid.PrintIndexAndCoord(cerr, "Error. Cube: ", cube0_index, "\n");
    cerr << "  Coord: ";
    print_coord3D(cerr, cube0_coord);
    cerr << "\n";

    error.AddMessage("Programming error.  No coord is 0 mod 3.");
    throw error;
  }

  // Return false if cube0 is near boundary.
  if (cube0_coord[mod3_coord] <= 1 || 
      cube0_coord[mod3_coord]+4 >= grid.AxisSize(mod3_coord))
    { return(false); }

  int d1 = (mod3_coord+1)%DIM3;
  int d2 = (mod3_coord+2)%DIM3;

  for (x[mod3_coord] = cube0_coord[mod3_coord]-2;
       x[mod3_coord] <= cube0_coord[mod3_coord]+2;
       x[mod3_coord]++) {

    if (x[mod3_coord] == cube0_coord[mod3_coord]) { continue; }

    for (x[d1] = cubeL_coord[d1]+1; x[d1] < cubeH_coord[d1]; x[d1]++) {
      for (x[d2] = cubeL_coord[2]+1; x[d2] < cubeH_coord[d2]; x[d2]++) {

        VERTEX_INDEX cube1_index = grid.ComputeVertexIndex(x);

        INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
        if (gcube1_index == ISOVERT::NO_INDEX) { continue; }

        if (isovert.gcube_list[gcube1_index].flag == SELECTED_GCUBE) {

          MSDEBUG();
          if (flag_debug) {
            grid.PrintIndexAndCoord
              (cerr, "  Cube ", cube0_index, " is near not mod 3 ");
            grid.PrintIndexAndCoord(cerr, "cube ", cube1_index, "\n");
          }
                                  
          return(true);
        }
      }
    }
  }

  return(false);
}

/// Select edge cubes (eigenvalue 2) whose coordinates are all 
///   congruent to 0 mod 3.
void select_cubes_cong_zero_mod3
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes cong to 0 mod 3." << endl; }

  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] == 0 && cube_coord_mod3[1] == 0 &&
        cube_coord_mod3[2] == 0) {

      check_and_select_edge_cube
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
         isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) whose coordinates are all NOT
///   congruent to 0 mod 3 and which are near a selected mod 3 cube.
void select_cubes_not_cong_zero_mod3_A
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

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

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      if (is_nearby_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) whose coordinates are all NOT
///   congruent to 0 mod 3 and whose isovert edge dir does NOT point
///   to an adjacent 0 mod 3 cube.
void select_cubes_not_cong_zero_mod3_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug) {
    cerr << endl 
         << "--- Selecting cubes not cong to 0 mod 3 whose isovert edge dir does not point to an adjacent mod 3 cube." << endl; 
  }

  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      if (!does_edge_dir_point_to_adjacent_zero_mod3
          (cube_coord_mod3, isovert.gcube_list[gcube_index].direction)) {

        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) with one coord congruent to 0 mod 3
///   and which are near a selected mod 3 cube.
void select_cubes_with_one_coord_cong_zero_mod3_A
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
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

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 1) {

      if (is_nearby_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) with one coord congruent to 0 mod 3
///   and which are near a selected NOT mod 3 cube.
void select_cubes_with_one_coord_cong_zero_mod3_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
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

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_zero = count_num_zero(cube_coord_mod3);
    if (num_zero == 1) {
      if (is_nearby_not_mod3_selected(scalar_grid, isovert, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

}

/// Select edge cubes (eigenvalue 2) with two coordinates congruent to 0 mod 3
///   where the two coordinates are closest to the isovert edge direction.
void select_cubes_with_two_coord_cong_zero_mod3_A
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 2 coord closest to isovert edge direction cong to 0 mod 3." << endl; }

  // Select from cubes with two coordinates congruent to 0 mod 3.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (isovert.gcube_list[gcube_index].num_eigenvalues != 2)
      { continue; }

    int d0 = 
      get_index_smallest_abs(isovert.gcube_list[gcube_index].direction);
    int d1 = (d0+1)%DIM3;
    int d2 = (d0+2)%DIM3;

    if (cube_coord_mod3[d1] == 0 && cube_coord_mod3[d2] == 0) {
      check_and_select_edge_cube
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
         isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);
}

/// Select edge cubes (eigenvalue 2) using mod 3 algorithm
void select_edge_cubes_mod3	
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<VERTEX_INDEX> & bin_grid,
 SHARPISO_GRID_NEIGHBORS & gridn,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const vector<NUM_TYPE> & sharp_gcube_list,
 ISOVERT & isovert)
{
  GRID_COORD_TYPE cube_coord_mod3[DIM3];

  // Select edge cubes (eigenvalue 2) whose coord are all congruent to 0 mod 3.
  select_cubes_cong_zero_mod3
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, isovalue,
     isovert_param, sharp_gcube_list, isovert);

  // Select edge cubes whose coord are all NOT congruent to 0 mod 3
  //   and which are near a selected 0 mod 3 cube.
  select_cubes_not_cong_zero_mod3_A
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, isovalue,
     isovert_param, sharp_gcube_list, isovert);

  // Select edge cubes with one coord congruent to 0 mod 3
  //   and which are near a selected 0 mod 3 cube.
  select_cubes_with_one_coord_cong_zero_mod3_A
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, isovalue,
     isovert_param, sharp_gcube_list, isovert);

  // Select edge cubes with one coord congruent to 0 mod 3
  //   and which are near a selected NOT 0 mod 3 cube.
  select_cubes_with_one_coord_cong_zero_mod3_B
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, isovalue,
     isovert_param, sharp_gcube_list, isovert);

  // Select edge cubes with two coord congruent to 0 mod 3
  //   where the two coord are closest to the isovert edge direction.
  select_cubes_with_two_coord_cong_zero_mod3_A
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, isovalue,
     isovert_param, sharp_gcube_list, isovert);


  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with coord not cong to 0 mod 3 adjacent to a covered cube." 
           << endl; }

  // Select (again) from cubes whose coordinates are not congruent to 0 mod 3
  //   and who are facet adjacent to a covered cube.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      if (is_facet_neighbor_covered_B(covered_grid, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with coord not cong to 0 mod 3." << endl; }

  // Select from cubes whose coordinates are not congruent to 0 mod 3
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    if (cube_coord_mod3[0] != 0 && cube_coord_mod3[1] != 0 &&
        cube_coord_mod3[2] != 0) {

      check_and_select_edge_cube
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
         isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 1 coord cong to 0 mod 3 adjacent to a covered cube." << endl; }

  // Select from cubes with one coordinate congruent to 0 mod 3
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_congruent(0);
    for (int d=0; d < 3; d++) {
      if (cube_coord_mod3[d] == 0) { num_congruent++; }
    }

    if (num_congruent == 1) {

      if (is_facet_neighbor_covered_B(covered_grid, cube_index)) {
        check_and_select_edge_cube
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
      }
    }
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 1 coord cong to 0 mod 3." << endl; }

  // Select from cubes with one coordinate congruent to 0 mod 3
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_congruent(0);
    for (int d=0; d < 3; d++) {
      if (cube_coord_mod3[d] == 0) { num_congruent++; }
    }

    if (num_congruent == 1) {

      check_and_select_edge_cube
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
         isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
    }
  }

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes with 2 coord cong to 0 mod 3." << endl; }

  // Select from cubes with two coordinates congruent to 0 mod 3.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    int num_congruent(0);
    for (int d=0; d < 3; d++) {
      if (cube_coord_mod3[d] == 0) { num_congruent++; }
    }

    if (num_congruent == 2) {
      check_and_select_edge_cube
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
         isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
    }
  }


  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting cubes containing isovert coord." << endl; }

  // Select from any cube which contains its isovert_coord.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

    if (isovert.gcube_list[gcube_index].cube_containing_isovert != cube_index)
      { continue; }

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    check_and_select_edge_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
       isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
  }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);
  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  if (flag_debug)
    { cerr << endl << "--- Selecting from any cube." << endl; }

  // Select from cubes whose isovert coordinates are outside their cubes.
  for (int i=0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];

    compute_coord_mod3
      (isovert.gcube_list[gcube_index].cube_coord, cube_coord_mod3);

    check_and_select_edge_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
       isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
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

  // get corner or edge cubes
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  BIN_GRID<VERTEX_INDEX> bin_grid;
  init_bin_grid(scalar_grid, bin_width, bin_grid);

  initialize_covered_by(isovert);

  SHARPISO_GRID_NEIGHBORS gridn;
  gridn.SetSize(scalar_grid);

  //Keeps track of all cubes which are covered.
  SHARPISO_BOOL_GRID covered_grid;
  covered_grid.SetSize(scalar_grid);
  covered_grid.SetAll(false);

  MSDEBUG();
  flag_debug = false;
  if (flag_debug) {
    cerr << endl << "*** Selecting corner cubes." << endl;
  }

  // select corner cubes
  select_corner_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  MSDEBUG();
  if (flag_debug) {
    cerr << endl << "*** Selecting edge cubes (mod3)." 
         << endl;
  }

  // select edge cubes.
  select_edge_cubes_mod3
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, 
     isovalue, isovert_param, sharp_gcube_list, isovert);

  MSDEBUG();
  if (flag_debug) 
    { cerr << endl << "*** Recomputing covered point positions." << endl; }

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  // Select edge cubes (again)
  select_edge_cubes_mod3
    (scalar_grid, gradient_grid, covered_grid, bin_grid, gridn, 
     isovalue, isovert_param, sharp_gcube_list, isovert);

  recompute_covered_point_positions
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param,
     isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  MSDEBUG();
  flag_debug = false;
}


// **************************************************
// RESET ISOVERT POSITIONS
// **************************************************

namespace {

  // Change isovert positions which lie in covered cubes.
  void reset_covered_isovert_positions
  (const SHARPISO_GRID_NEIGHBORS & grid, 
   const SHARPISO_BOOL_GRID & covered_grid,
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
                  ijkgrid_output_vertex_coord(cerr, grid, cubeA_index);
                  cerr << " to " << cubeB_index << " ";
                  ijkgrid_output_vertex_coord(cerr, grid, cubeB_index);
                  cerr << endl;
                  cerr << "  Old position: ";
                  IJK::print_coord3D(cerr, isovert.gcube_list[gcubeB_index].isovert_coord);
                  cerr << ".  New position: ";
                  IJK::print_coord3D(cerr, isovert.gcube_list[gcubeA_index].isovert_coord);
                  cerr << endl;
                }

                copy_isovert_position(grid, gcubeA_index, gcubeB_index, isovert);
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
  (const GRID_CUBE_DATA & c,
   const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   SHARPISO_GRID_NEIGHBORS & gridn,
   ISOVERT &isovert,
   GRID_CUBE_FLAG flag)
  {
    VERTEX_INDEX neighbor_cube_index;

    // Check facet neighbors
    for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) {

      neighbor_cube_index = gridn.CubeNeighborF(c.cube_index, j);
      //neighbor_gcube_index is an entry into the gcube list 
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

      if(neighbor_gcube_index == ISOVERT::NO_INDEX)
        {
          continue;
        }
      if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
        {
          return true;
        }
    }

    // Check edge neighbors
    for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) {
      neighbor_cube_index = gridn.CubeNeighborE(c.cube_index, j);
      //neighbor_gcube_index is an entry into the gcube list 
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

      if(neighbor_gcube_index == ISOVERT::NO_INDEX)
        {
          continue;
        }
      if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
        {
          return true;
        }
    }

    // Check vertex neighbors
    for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) {
      neighbor_cube_index = gridn.CubeNeighborV(c.cube_index, j);
      //neighbor_gcube_index is an entry into the gcube list 
      INDEX_DIFF_TYPE neighbor_gcube_index
        = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

      if(neighbor_gcube_index == ISOVERT::NO_INDEX)
        {
          continue;
        }
      if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
        {
          return true;
        }
    }

    return false;
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

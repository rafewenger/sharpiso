/// \file mergesharp_check_map.cxx
/// Check if mapping one isosurface vertex to another distorts or reverses
///   any isosurface triangles.


/*
Copyright (C) 2015 Arindam Bhattacharya and Rephael Wenger

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
#include <math.h>

#include "ijkcoord.txx"
#include "ijkgrid_macros.h"

#include "sharpiso_array.txx"
#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"
#include "mergesharp_check_map.h"

// *** DEBUG ***
#include "mergesharp_debug.h"
#include "ijkprint.txx"


using namespace MERGESHARP;


// **************************************************
//  Forward declarations
// **************************************************

namespace {

  bool is_strictly_between
  (const COORD_TYPE x, const COORD_TYPE a0, const COORD_TYPE a1);

  bool is_covered_A_or_corner(const GRID_CUBE_FLAG flag);

  bool is_in_3x3x3_region(const GRID_CUBE_FLAG flag);

  bool are_separated
  (const SHARPISO_GRID & grid,
   const GRID_CUBE_DATA & gcubeA, const COORD_TYPE coordA[DIM3],
   const GRID_CUBE_DATA & gcubeB, const COORD_TYPE coordB[DIM3],
   const int orth_dir);

  void subtract_project_normalize
  (const COORD_TYPE pA[DIM3], const COORD_TYPE pB[DIM3], 
   const int orth_dir, const COORD_TYPE max_small_magnitude,
   COORD_TYPE w[DIM3], COORD_TYPE & magnitude, bool & flag_zero);

  void determine_quad_maps
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_vertex_maps_to_cube1,
   bool & flag_vertex_maps_to_cube2,
   bool & flag_quad_maps_to_012X,
   bool & flag_quad_maps_to_021X_or_02X1,
   bool & flag_edge_maps_to_02);

  void determine_quad_edge_map
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_vertex_maps_to_cube0,
   bool & flag_vertex_maps_to_cube1,
   bool & flag_edge_maps_to_both_cubes);

  void determine_quad_map_to_both_cubes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_maps_to_cube1,
   bool & flag_maps_to_both_cubes);

  void temp_map_determine_quad_map_to_both_cubes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const NUM_TYPE from_gcube_index,
   const NUM_TYPE to_gcube_index,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_maps_to_cube1,
   bool & flag_maps_to_both_cubes);

  bool temp_map_check_single_edge_manifold
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube0_index,
   const VERTEX_INDEX to_cube0_index,
   const VERTEX_INDEX from_cube1_index,
   const VERTEX_INDEX cubeA_index,
   const VERTEX_INDEX cubeB_index,
   const ISOVERT & isovert,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);

  bool temp_map_check_edge_XA_manifoldII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube0_index,
   const VERTEX_INDEX to_cube0_index,
   const VERTEX_INDEX from_cube1_index,
   const VERTEX_INDEX to_cube1_index,
   const VERTEX_INDEX from_cube2_index,
   const VERTEX_INDEX to_cube2_index,
   const VERTEX_INDEX cubeA_index,
   const MERGESHARP::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);

	bool are_connected_by_iso_quad
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX & cube0_index,
   const VERTEX_INDEX & cube1_index,
   const SCALAR_TYPE isovalue,
   const MERGESHARP::ISOVERT & isovert, 
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);

}


// **************************************************
//  Check distortion
// **************************************************

/// Return true if mapping of from_cube to to_cube does not reverse any
///   triangles incident on vertex in from_cube 
///   or create any degenerate triangles.
/// @param flag_strict If true, use strict tests.
bool MERGESHARP::check_distortion
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const bool flag_strict,
 const MERGE_PARAM & param)
{
  bool flag;

  if (flag_strict) {
    flag = 
      check_distortion_strict
      (scalar_grid, isovalue, isovert, gcube_map, from_cube, to_cube, param);
  }
  else {
    flag =
      check_distortion_loose
      (scalar_grid, isovalue, isovert, gcube_map, from_cube, to_cube, param);
  }

  return(flag);
}


/// Return true if mapping of from_cube to to_cube does not reverse any
///   triangles incident on vertex in from_cube 
///   or create any degenerate triangles.
bool MERGESHARP::check_distortion_strict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const MERGE_PARAM & param)
{
  IJK::PROCEDURE_ERROR error("check_distortion_strict");
  const INDEX_DIFF_TYPE from_gcube_index = 
    isovert.GCubeIndex(from_cube, error);
  BOUNDARY_BITS_TYPE boundary_bits;

  if (from_gcube_index == ISOVERT::NO_INDEX) { throw error; }

  boundary_bits = isovert.gcube_list[from_gcube_index].boundary_bits;
  if (boundary_bits == 0) {

    for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

      const int dir1 = (edge_dir+1)%DIM3;
      const int dir2 = (edge_dir+2)%DIM3;

      for (int j1 = 0; j1 < 2; j1++) {
        for (int j2 = 0; j2 < 2; j2++) {

          VERTEX_INDEX iend0 = from_cube + j1*scalar_grid.AxisIncrement(dir1)
            + j2*scalar_grid.AxisIncrement(dir2);

          if (!check_quad_distortion_strict
              (scalar_grid, isovalue, isovert, gcube_map, from_cube, to_cube,
               edge_dir, iend0, j1, j2, param)) 
            { return(false); }
        }
      }
    }
  }
  else {
    // HANDLE BOUNDARY CASE
  }

  return(true);
}

/// Return true if mapping of from_cube to to_cube does not reverse any
///   triangles incident on vertex in from_cube 
///   or create any degenerate triangles.
bool MERGESHARP::check_distortion_loose
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const MERGE_PARAM & param)
{
  IJK::PROCEDURE_ERROR error("check_distortion_loose");
  const INDEX_DIFF_TYPE from_gcube_index = 
    isovert.GCubeIndex(from_cube, error);
  BOUNDARY_BITS_TYPE boundary_bits;

  if (from_gcube_index == ISOVERT::NO_INDEX) { throw error; }

  boundary_bits = isovert.gcube_list[from_gcube_index].boundary_bits;
  if (boundary_bits == 0) {

    for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

      const int dir1 = (edge_dir+1)%DIM3;
      const int dir2 = (edge_dir+2)%DIM3;

      for (int j1 = 0; j1 < 2; j1++) {
        for (int j2 = 0; j2 < 2; j2++) {

          VERTEX_INDEX iend0 = from_cube + j1*scalar_grid.AxisIncrement(dir1)
            + j2*scalar_grid.AxisIncrement(dir2);

          if (!check_quad_distortion_loose
              (scalar_grid, isovalue, isovert, gcube_map, from_cube, to_cube,
               edge_dir, iend0, j1, j2, param)) 
            { return(false); }
        }
      }
    }
  }
  else {
    // HANDLE BOUNDARY CASE
  }

  return(true);
}


/// Return true if simultaneous mapping of from_cube0 to to_cube0
///   and from_cube1 to to_cube1 does not reverse any triangles 
///   incident on vertices in from_cube0 or from_cube1 or create 
///    any degenerate triangles.
bool MERGESHARP::check_distortionII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube0, const VERTEX_INDEX to_cube0,
 const VERTEX_INDEX from_cube1, const VERTEX_INDEX to_cube1,
 const bool flag_strict,
 const MERGE_PARAM & param)
{
  IJK::PROCEDURE_ERROR error("check_distortionII");
  const INDEX_DIFF_TYPE from_gcube0_index = 
    isovert.GCubeIndex(from_cube0, error);
  const INDEX_DIFF_TYPE from_gcube1_index = 
    isovert.GCubeIndex(from_cube1, error);
  const INDEX_DIFF_TYPE to_gcube0_index = 
    isovert.GCubeIndex(to_cube0, error);
  const INDEX_DIFF_TYPE to_gcube1_index = 
    isovert.GCubeIndex(to_cube1, error);
  NUM_TYPE store_map[2];
  bool flag;

  if (from_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
  if (from_gcube1_index == ISOVERT::NO_INDEX) { throw error; }
  if (to_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
  if (to_gcube1_index == ISOVERT::NO_INDEX) { throw error; }

  store_map[0] = gcube_map[from_gcube0_index];
  store_map[1] = gcube_map[from_gcube1_index];

  // temporarily set gcube_map[from_gcube1_index] to to_cube1
  gcube_map[from_gcube1_index] = to_gcube1_index;

  flag = check_distortion
    (scalar_grid, isovalue, isovert, gcube_map, from_cube0, to_cube0, 
     flag_strict, param);

  // restore gcube_map[from_gcube1_index]
  gcube_map[from_gcube1_index] = store_map[1];

  if (!flag) { return(false); }

  // temporarily set gcube_map[from_gcube0_index] to to_cube0
  gcube_map[from_gcube0_index] = to_gcube0_index;

  flag = check_distortion
    (scalar_grid, isovalue, isovert, gcube_map, from_cube1, to_cube1, 
     flag_strict, param);

  // restore gcube_map[from_gcube0_index]
  gcube_map[from_gcube0_index] = store_map[0];

  if (!flag) { return(false); }

  return(true);
}


/// Return true if simultaneous mapping of three cubes to to_cube
///   does not reverse any triangles incident on vertices on cubes
///   or create any degenerate triangles.
bool MERGESHARP::check_distortionIII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX cube_index[3], const VERTEX_INDEX to_cube,
 const bool flag_strict,
 const MERGE_PARAM & param)
{
  INDEX_DIFF_TYPE gcube_index[3], to_gcube;
  IJK::PROCEDURE_ERROR error("check_distortionIII");
  NUM_TYPE store_map[3];
  bool flag;

  to_gcube = isovert.GCubeIndex(to_cube, error);
  if (to_gcube == ISOVERT::NO_INDEX) { throw error; }
  for (int i = 0; i < 3; i++) {
    gcube_index[i] = isovert.GCubeIndex(cube_index[i], error);
    if (gcube_index[i] == ISOVERT::NO_INDEX) { throw error; }

    store_map[i] = gcube_map[gcube_index[i]];
  }

  for (int i0 = 0; i0 < 3; i0++) {

    int i1 = (i0+1)%3;
    int i2 = (i0+2)%3;

    // temporarily set gcube_map[] for gcube_index[i1] and gcube_index[i2]
    //   to to_cube
    gcube_map[gcube_index[i1]] = to_gcube;
    gcube_map[gcube_index[i2]] = to_gcube;

    flag = check_distortion
      (scalar_grid, isovalue, isovert, gcube_map, cube_index[i0], 
       to_cube, flag_strict, param);

    // restore gcube_map[]
    gcube_map[gcube_index[i1]] = store_map[i1];
    gcube_map[gcube_index[i2]] = store_map[i2];

    if (!flag) { return(false); }
  }

  return(true);
}


/// Return true if mapping of from_cube to to_cube does not reverse/distort
///   triangles FAB or FBC on quad (FABC) dual to (iend0, iend1).
bool MERGESHARP::check_quad_distortion_strict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const int edge_dir,
 const VERTEX_INDEX iend0,
 const int j1,
 const int j2,
 const MERGE_PARAM & param)
{
  const VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

  if (!is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {
    // Edge is not bipolar
    return(true);
  }

  const int dir1 = (edge_dir+1)%DIM3;
  const int dir2 = (edge_dir+2)%DIM3;

  const VERTEX_INDEX icubeA = scalar_grid.AdjacentVertex(from_cube, dir1, j1);
  const VERTEX_INDEX icubeB = scalar_grid.AdjacentVertex(icubeA, dir2, j2);
  const VERTEX_INDEX icubeC = scalar_grid.AdjacentVertex(from_cube, dir2, j2);

  if (!check_tri_distortion_mapA
      (scalar_grid, isovert, gcube_map, from_cube, icubeA, icubeB, 
       to_cube, dir1, dir2, j1, j2, param))
    { return(false); }

  if (!check_tri_distortion_mapA
      (scalar_grid, isovert, gcube_map, from_cube, icubeC, icubeB, 
       to_cube, dir2, dir1, j2, j1, param))
    { return(false); }

  return(true);
}

/// Return true if mapping of from_cube to to_cube does not reverse/distort
///   triangles FAB or FBC on quad (FABC) dual to (iend0, iend1)
///   or does not reverse/distort triangle FAC.
bool MERGESHARP::check_quad_distortion_loose
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const int edge_dir,
 const VERTEX_INDEX iend0,
 const int j1,
 const int j2,
 const MERGE_PARAM & param)
{
  const VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);
  IJK::PROCEDURE_ERROR error("check_quad_distortion_loose");

  if (!is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {
    // Edge is not bipolar
    return(true);
  }

  const int dir1 = (edge_dir+1)%DIM3;
  const int dir2 = (edge_dir+2)%DIM3;

  const VERTEX_INDEX icubeA = scalar_grid.AdjacentVertex(from_cube, dir1, j1);
  const VERTEX_INDEX icubeB = scalar_grid.AdjacentVertex(icubeA, dir2, j2);
  const VERTEX_INDEX icubeC = scalar_grid.AdjacentVertex(from_cube, dir2, j2);

  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA, error);
  const INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(icubeC, error);
  const INDEX_DIFF_TYPE to_gcube_index = isovert.GCubeIndex(to_cube, error);
  if ((gcubeA_index == ISOVERT::NO_INDEX) || (gcubeC_index == ISOVERT::NO_INDEX) ||
      (to_gcube_index == ISOVERT::NO_INDEX))
    { throw error; }

  if (gcube_map[gcubeA_index] != to_gcube_index && 
      gcube_map[gcubeC_index] != to_gcube_index) {
    if (check_tri_distortion_mapB
        (scalar_grid, isovert, gcube_map, icubeA, from_cube, icubeC,
         to_cube, dir2, dir1, j2, j1, param)) { 
      // Does not reverse FAC.
      return(true);
    }
  }

  if (!check_tri_distortion_mapA
      (scalar_grid, isovert, gcube_map, from_cube, icubeA, icubeB, 
       to_cube, dir1, dir2, j1, j2, param)) {
    // Reverses FAB
    return(false); 
  }

  if (!check_tri_distortion_mapA
      (scalar_grid, isovert, gcube_map, from_cube, icubeC, icubeB, 
       to_cube, dir2, dir1, j2, j1, param)) {
    // Reverses FBC
    return(false); 
  }

  return(true);
}


// Return true if mapping of icubeA to to_cube does not reverse/distort
//   triangle with vertices in (icubeA,icubeB,icubeC).
// @param Cube icubeB shares a facet with icubeA.
// @param Cube icubeC shares an edge with icubeA.
// @param relposAB Relative positions of cubes icubeA and icubeB
//         If relposAB = 0, then icubeB precedes icubeA in direction dirAB.
//         If relposAB = 1, then icubeB follows icubeA in direction dirAB.
// @param relposAB Relative positions of cubes icubeA and icubeB.
//         If relposAB = 0, then icubeC precedes icubeB in direction dirAB.
//         If relposAB = 1, then icubeC follows icubeB in direction dirAB.
bool MERGESHARP::check_tri_distortion_mapA
(const SHARPISO_GRID & grid, const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB, 
 const VERTEX_INDEX icubeC,
 const VERTEX_INDEX to_cube,
 const int dirAB, const int dirBC, const int relposAB, const int relposBC,
 const MERGE_PARAM & param)
{
  const COORD_TYPE * spacing = grid.SpacingPtrConst();
  IJK::PROCEDURE_ERROR error("check_tri_distortion_mapA");
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA, error);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB, error);
  const INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(icubeC, error);
  const INDEX_DIFF_TYPE to_gcube_index = isovert.GCubeIndex(to_cube, error);
  COORD_TYPE unscaled_to_coord[DIM3];
  COORD_TYPE unscaled_coordA[DIM3], unscaled_coordB[DIM3], 
    unscaled_coordC[DIM3];
  bool flag_small_magnitude;

  // Orientation (-1,1) of icubeA, icubeB and icubeC around shared edge.
  const int cube_orient = (2*relposAB-1)*(2*relposBC-1);

  if ((gcubeA_index == ISOVERT::NO_INDEX) ||
      (gcubeB_index == ISOVERT::NO_INDEX) || 
      (gcubeC_index == ISOVERT::NO_INDEX) || 
      (to_gcube_index == ISOVERT::NO_INDEX))
    { throw error; }

  if ((gcube_map[gcubeB_index] == to_gcube_index) || 
      (gcube_map[gcubeC_index] == to_gcube_index) ||
      (gcube_map[gcubeB_index] == gcube_map[gcubeC_index])) 
    { return(true); }

  GRID_CUBE_FLAG gcubeB_flag = isovert.gcube_list[gcube_map[gcubeB_index]].flag;
  GRID_CUBE_FLAG gcubeC_flag = isovert.gcube_list[gcube_map[gcubeC_index]].flag;

  if (is_covered_A_or_corner(gcubeB_flag)) { return(true); }
  if (is_covered_A_or_corner(gcubeC_flag)) { return(true); }

  if (gcubeB_flag == UNAVAILABLE_GCUBE || gcubeC_flag == UNAVAILABLE_GCUBE)
    { return(true); }

  // Reverse scaling on coordA[], coordB[], coordC[]
  const COORD_TYPE * to_coord = isovert.IsoVertCoord(to_gcube_index);
  const COORD_TYPE * coordA = isovert.IsoVertCoord(gcubeA_index);
  const COORD_TYPE * coordB = isovert.IsoVertCoord(gcube_map[gcubeB_index]);
  const COORD_TYPE * coordC = isovert.IsoVertCoord(gcube_map[gcubeC_index]);
  IJK::reverse_coord_scaling (DIM3, spacing, to_coord, unscaled_to_coord);
  IJK::reverse_coord_scaling (DIM3, spacing, coordA, unscaled_coordA);
  IJK::reverse_coord_scaling(DIM3, spacing, coordB, unscaled_coordB);
  IJK::reverse_coord_scaling(DIM3, spacing, coordC, unscaled_coordC);

  // Triangle angle test.
  const COORD_TYPE min_dist = 0.01;
  const COORD_TYPE cos_min_triangle_angle = param.CosMinTriangleAngle();
  const COORD_TYPE cos_min_sharp_cube_triangle_angle = 
    param.CosMinSharpCubeTriangleAngle();
  COORD_TYPE cos_angle_ABC, cos_angle_ACB;

  compute_cos_triangle_angles
    (grid, unscaled_to_coord, unscaled_coordB, unscaled_coordC,
     min_dist, cos_angle_ABC, cos_angle_ACB, flag_small_magnitude);

  // *** DEBUG ***
  if (flag_debug) {
    MSDEBUG();
    grid.PrintIndexAndCoord
      (cerr, "To cube: ", to_cube, " cube B: ", icubeB, 
       " cube C: ", icubeC,"\n");
    cerr << "  To coord: ";
    IJK::print_coord3D(cerr, unscaled_to_coord);
    cerr << "  Coord B: ";
    IJK::print_coord3D(cerr, unscaled_coordB);
    cerr << "  Coord C: ";
    IJK::print_coord3D(cerr, unscaled_coordC);
    cerr << endl;

    cerr << "  angle_ABC: " << acos(cos_angle_ABC)*180.0/M_PI
         << "  angle_ACB: " << acos(cos_angle_ACB)*180.0/M_PI
         << "  flag_small_magnitude: " << int(flag_small_magnitude) << endl;
  }

  if (!flag_small_magnitude) {

    if (cos_angle_ABC > cos_min_triangle_angle ||
        cos_angle_ACB > cos_min_triangle_angle) {

      // *** DEBUG ***
      if (flag_debug) {
        MSDEBUG();
        cerr << "xxx Failed angle test." 
             << "  angle_ABC: " << acos(cos_angle_ABC) * 180.0/M_PI
             << "  angle_ACB: " << acos(cos_angle_ACB) * 180.0/M_PI
             << endl;
      }

      return(false);
    }
  }

  // Normal angle test
  const COORD_TYPE cos_min_normal_angle = param.CosMinNormalAngle();

  COORD_TYPE cos_normal_angle;
  compute_cos_dihedral_angle
    (grid, unscaled_to_coord, unscaled_coordA, unscaled_coordB, unscaled_coordC,
     min_dist, cos_normal_angle, flag_small_magnitude);

  if (!flag_small_magnitude) {
    if (cos_normal_angle >= cos_min_normal_angle) {

      if (flag_debug) {
        MSDEBUG();
        // Numerical error could cause cos_normal_angle to be out of bounds.
        if (cos_normal_angle > 1) { cos_normal_angle = 1; }
        if (cos_normal_angle < -1) { cos_normal_angle = -1; }
        cerr << "    Passed normal angle test.  cos_normal_angle: "
             << cos_normal_angle
             << ".  Normal angle: "
             << acos(cos_normal_angle)*180.0/M_PI << endl;
      }

      return(true);
    }
  }

  // Cube orientation test
  const COORD_TYPE min_proj_angle = 1;
  const COORD_TYPE cos_min_proj_angle = cos(min_proj_angle*M_PI/180.0);

  // edge_dir = Direction of edge shared by gcubeA, gcubeB and gcubeC.
  const int edge_dir = (2*(dirAB+dirBC))%DIM3;

  int orient;
  bool flag_small_BC;
  compute_projected_tri_orientation
    (grid, unscaled_to_coord, unscaled_coordB, unscaled_coordC,
     cos_min_proj_angle, min_dist, edge_dir, dirAB, dirBC, 
     orient, flag_small_BC);


  // *** DEBUG ***
  if (flag_debug) {

    MSDEBUG();
    if (orient == cube_orient) {

      cerr << "xxx Matches cube orientations. Fails normal test." << endl;
      grid.PrintIndexAndCoord
        (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
      grid.PrintIndexAndCoord
        (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
      if (!flag_small_magnitude) {
        if (cos_normal_angle < -1) { cos_normal_angle = -1; }
        if (cos_normal_angle > 1) { cos_normal_angle = 1; }
        cerr << "    normal angle: " << acos(cos_normal_angle)*180.0/M_PI
             << endl;
      }
      else {
        cerr << "    magnitude < " << min_dist << endl;
      }
    }
  }


  if (orient == cube_orient) { return(true); }

  // Cube separation test.
  if (flag_small_BC) {

    if (are_separated
        (grid, isovert.gcube_list[gcubeA_index], unscaled_coordA, 
         isovert.gcube_list[gcubeB_index], unscaled_coordB, dirAB)) {
      if (are_separated
          (grid, isovert.gcube_list[gcubeA_index], unscaled_coordA, 
           isovert.gcube_list[gcubeC_index], unscaled_coordC, dirAB)) {

        // Since direction from coordB[] and coordC[] to coordA[]
        //   matches direction from gcubeB and gcubeC to gcubeA,
        //   new triangle orientation matches original.

        if (flag_debug) {
          MSDEBUG();
          cerr << "zzz Passes cube separation test.  Fails normal & cube orient tests." << endl;
          grid.PrintIndexAndCoord
            (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
          grid.PrintIndexAndCoord
            (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
          cerr << "    normal angle: " << acos(cos_normal_angle*M_PI/180.0)
               << endl;
        }

        return(true);
      }
    }

    return(false);
  }
  else { 

    if (flag_debug) {
      MSDEBUG();
      cerr << "    Triangle: ";
      IJK::print_coord3D(cerr, unscaled_coordA);
      IJK::print_coord3D(cerr, unscaled_coordB);
      IJK::print_coord3D(cerr, unscaled_coordC);
      cerr << endl;
      grid.PrintIndexAndCoord
        (cerr, "    Cubes ", icubeA, " ", icubeB, " ", icubeC, "\n");
      cerr << "    orth_dir: " << edge_dir 
           << "  orient: " << orient << endl;
    }

    return(false); 
  }
}


// Return true if mapping of icubeB to to_cube does not reverse/distort
//   triangle with vertices in (icubeA,icubeB,icubeC).
// @param Cube icubeB shares a facet with icubeA.
// @param Cube icubeC shares an edge with icubeA.
// @param relposAB Relative positions of cubes icubeA and icubeB
//         If relposAB = 0, then icubeB precedes icubeA in direction dirAB.
//         If relposAB = 1, then icubeB follows icubeA in direction dirAB.
// @param relposAB Relative positions of cubes icubeA and icubeB.
//         If relposAB = 0, then icubeC precedes icubeB in direction dirAB.
//         If relposAB = 1, then icubeC follows icubeB in direction dirAB.
bool MERGESHARP::check_tri_distortion_mapB
(const SHARPISO_GRID & grid, const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB, 
 const VERTEX_INDEX icubeC,
 const VERTEX_INDEX to_cube,
 const int dirAB, const int dirBC, const int relposAB, const int relposBC,
 const MERGE_PARAM & param)
{
  const COORD_TYPE * spacing = grid.SpacingPtrConst();
  IJK::PROCEDURE_ERROR error("check_tri_distortion_mapB");
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA, error);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB, error);
  const INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(icubeC, error);
  const INDEX_DIFF_TYPE to_gcube_index = isovert.GCubeIndex(to_cube, error);
  COORD_TYPE unscaled_to_coord[DIM3];
  COORD_TYPE unscaled_coordA[DIM3], unscaled_coordB[DIM3], 
    unscaled_coordC[DIM3];
  bool flag_small_magnitude;

  // Orientation (-1,1) of icubeA, icubeB and icubeC around shared edge.
  const int cube_orient = (2*relposAB-1)*(2*relposBC-1);

  if ((gcubeA_index == ISOVERT::NO_INDEX) ||
      (gcubeB_index == ISOVERT::NO_INDEX) || 
      (gcubeC_index == ISOVERT::NO_INDEX) || 
      (to_gcube_index == ISOVERT::NO_INDEX))
    { throw error; }

  if ((gcube_map[gcubeA_index] == to_gcube_index) || 
      (gcube_map[gcubeC_index] == to_gcube_index) ||
      (gcube_map[gcubeA_index] == gcube_map[gcubeC_index])) 
    { return(true); }

  GRID_CUBE_FLAG gcubeA_flag = isovert.gcube_list[gcube_map[gcubeA_index]].flag;
  GRID_CUBE_FLAG gcubeC_flag = isovert.gcube_list[gcube_map[gcubeC_index]].flag;

  if (is_covered_A_or_corner(gcubeA_flag)) { return(true); }
  if (is_covered_A_or_corner(gcubeC_flag)) { return(true); }

  if (gcubeA_flag == UNAVAILABLE_GCUBE || gcubeC_flag == UNAVAILABLE_GCUBE)
    { return(true); }

  // Reverse scaling on coordA[], coordB[], coordC[]
  const COORD_TYPE * to_coord = isovert.IsoVertCoord(to_gcube_index);
  const COORD_TYPE * coordA = isovert.IsoVertCoord(gcube_map[gcubeA_index]);
  const COORD_TYPE * coordB = isovert.IsoVertCoord(gcubeB_index);
  const COORD_TYPE * coordC = isovert.IsoVertCoord(gcube_map[gcubeC_index]);
  IJK::reverse_coord_scaling (DIM3, spacing, to_coord, unscaled_to_coord);
  IJK::reverse_coord_scaling (DIM3, spacing, coordA, unscaled_coordA);
  IJK::reverse_coord_scaling(DIM3, spacing, coordB, unscaled_coordB);
  IJK::reverse_coord_scaling(DIM3, spacing, coordC, unscaled_coordC);

  // Triangle angle test.
  const COORD_TYPE min_dist = 0.01;
  const COORD_TYPE cos_min_triangle_angle = param.CosMinTriangleAngle();
  const COORD_TYPE cos_min_sharp_cube_triangle_angle = 
    param.CosMinSharpCubeTriangleAngle();
  COORD_TYPE cos_angle_BAC, cos_angle_BCA;

  compute_cos_triangle_angles
    (grid, unscaled_to_coord, unscaled_coordA, unscaled_coordC,
     min_dist, cos_angle_BAC, cos_angle_BCA, flag_small_magnitude);

  // *** DEBUG ***
  if (flag_debug) {
    MSDEBUG();
    grid.PrintIndexAndCoord
      (cerr, "Cube A: ", icubeA, " cube B: ", icubeB, " cube C: ", icubeC,"\n");
    cerr << "  Coord A: ";
    IJK::print_coord3D(cerr, unscaled_coordA);
    cerr << "  Coord B: ";
    IJK::print_coord3D(cerr, unscaled_coordB);
    cerr << "  Coord C: ";
    IJK::print_coord3D(cerr, unscaled_coordC);
    cerr << endl;

    cerr << "  angle_BAC: " << acos(cos_angle_BAC)*180.0/M_PI
         << "  angle_BCA: " << acos(cos_angle_BCA)*180.0/M_PI
         << "  flag_small_magnitude: " << int(flag_small_magnitude) << endl;
  }

  if (!flag_small_magnitude) {

    if (cos_angle_BAC > cos_min_triangle_angle ||
        cos_angle_BCA > cos_min_triangle_angle) {

      // *** DEBUG ***
      if (flag_debug) {
        MSDEBUG();
        cerr << "xxx Failed angle test." 
             << "  angle_BAC: " << acos(cos_angle_BAC) * 180.0/M_PI
             << "  angle_BCA: " << acos(cos_angle_BCA) * 180.0/M_PI
             << endl;
      }
      
      return(false);
    }
  }


  // Normal angle test
  const COORD_TYPE cos_min_normal_angle = param.CosMinNormalAngle();

  COORD_TYPE cos_normal_angle;
  compute_cos_dihedral_angle
    (grid, unscaled_to_coord, unscaled_coordB, unscaled_coordA, unscaled_coordC,
     min_dist, cos_normal_angle, flag_small_magnitude);

  if (!flag_small_magnitude) {
    if (cos_normal_angle >= cos_min_normal_angle)
      { return(true); }
  }

  // Cube orientation test
  const COORD_TYPE min_proj_angle = 1;
  const COORD_TYPE cos_min_proj_angle = cos(min_proj_angle*M_PI/180.0);

  // edge_dir = Direction of edge shared by gcubeA, gcubeB and gcubeC.
  const int edge_dir = (2*(dirAB+dirBC))%DIM3;

  int orient;
  bool flag_small_BC;
  compute_projected_tri_orientation
    (grid, unscaled_to_coord, unscaled_coordB, unscaled_coordC,
     cos_min_proj_angle, min_dist, edge_dir, dirAB, dirBC, 
     orient, flag_small_BC);


  // *** DEBUG ***
  if (flag_debug) {

    MSDEBUG();
    if (orient == cube_orient) {

      cerr << "xxx Matches cube orientations. Fails normal test." << endl;
      grid.PrintIndexAndCoord
        (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
      grid.PrintIndexAndCoord
        (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
      if (!flag_small_magnitude) {
        if (cos_normal_angle < -1) { cos_normal_angle = -1; }
        if (cos_normal_angle > 1) { cos_normal_angle = 1; }
        cerr << "    normal angle: " << acos(cos_normal_angle)*180.0/M_PI
             << endl;
      }
      else {
        cerr << "    magnitude < " << min_dist << endl;
      }
    }
  }

  if (orient == cube_orient) { return(true); }
  else { return(false); }
}


// Compute projected triangle orientation.
// @param orth_dir Direction of orthogonal projection onto plane.
// @param orient[out] Orientation, 1,-1 or 0.  
//       0 if cos projected angle ABC or ACB is below cos_min_proj_angle
//       or if distance between projected points is less than 
//       or equal to min_dist.
// @param flag_small_BC[out] If true, distance between projected 
//       coordB and coordC is less than or equal to min_dist.
// @pre (orth_dir, dirAB, dirBC) are distinct.
void MERGESHARP::compute_projected_tri_orientation
(const SHARPISO_GRID & grid, 
 const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3],
 const COORD_TYPE coordC[DIM3],
 const COORD_TYPE cos_min_proj_angle, const COORD_TYPE min_dist,
 const int orth_dir, const int dirAB, const int dirBC,
 int & orient, bool & flag_small_BC)
{
  COORD_TYPE wBA[DIM3], wBC[DIM3], wAC[DIM3];
  COORD_TYPE cos_angle;
  COORD_TYPE magnitude;
  bool flag_zero;
  IJK::PROCEDURE_ERROR error("compute_projected_tri_orientation");

 
  orient = 0;
  flag_small_BC = false;

  // edge_dir = Direction of edge shared by gcubeA, gcubeB and gcubeC.
  const int edge_dir = (2*(dirAB+dirBC))%DIM3;

  subtract_project_normalize
    (coordA, coordB, orth_dir, min_dist, wBA, magnitude, flag_zero); 
  if (flag_zero) { return; }

  subtract_project_normalize
    (coordC, coordA, orth_dir, min_dist, wAC, magnitude, flag_zero); 
  if (flag_zero) { return; }

  subtract_project_normalize
    (coordC, coordB, orth_dir, min_dist, wBC, magnitude, flag_zero); 
  if (flag_zero) { 
    flag_small_BC = true;
    return;
  }

  IJK::compute_inner_product_3D(wBA, wBC, cos_angle);
  if (cos_angle > cos_min_proj_angle) { return; }

  IJK::compute_inner_product_3D(wBC, wAC, cos_angle);
  if (cos_angle > cos_min_proj_angle) { return; }

  // Now check orientation.
  double D;
  IJK::determinant_2x2
    (wBC[dirAB], wBC[dirBC], wBA[dirAB], wBA[dirBC], D);

  if (D > 0) { orient = 1; }
  else if (D < 0) { orient = -1; }
}


// Compute cosine of the normal angle between (A0,B,C) and (A1,B,C)
void MERGESHARP::compute_cos_dihedral_angle
(const SHARPISO_GRID & grid, 
 const COORD_TYPE coordA0[DIM3], const COORD_TYPE coordA1[DIM3],
 const COORD_TYPE coordB[DIM3], const COORD_TYPE coordC[DIM3],
 const COORD_TYPE min_dist,
 COORD_TYPE & cos_dihedral, bool & flag_small_magnitude)
{
  COORD_TYPE wA0B[DIM3], wA1B[DIM3], wBC[DIM3];
  COORD_TYPE u0[DIM3], u1[DIM3];
  COORD_TYPE magnitude;
  bool flag_zero;

  cos_dihedral = -1;

  IJK::compute_unit_vector_3D
    (coordB, coordC, min_dist, wBC, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_unit_vector_3D
    (coordA0, coordB, min_dist, wA0B, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_unit_vector_3D
    (coordA1, coordB, min_dist, wA1B, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_cross_product_3D(wBC, wA0B, u0);
  IJK::compute_cross_product_3D(wBC, wA1B, u1);

  IJK::normalize_vector
    (DIM3, u0, min_dist, u0, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::normalize_vector
    (DIM3, u1, min_dist, u1, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_inner_product_3D(u0, u1, cos_dihedral);

  // *** DEBUGXXX ***
  if (flag_debug) {
    MSDEBUG();
    IJK::print_coord3D(cerr, "    A0: ", coordA0, "  A1: ", coordA1, "\n");
    IJK::print_coord3D(cerr, "    wBC: ", wBC, "\n");
    IJK::print_coord3D(cerr, "    wA0B: ", wA0B, "  wA1B:", wA1B, "\n");
    IJK::print_coord3D(cerr, "    u0: ", u0, "  u1: ", u1, "\n");
  }
}


// Compute cos triangle angles.
void MERGESHARP::compute_cos_triangle_angles
(const SHARPISO_GRID & grid, 
 const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3],
 const COORD_TYPE coordC[DIM3],
 const COORD_TYPE min_dist,
 COORD_TYPE & cos_angle_ABC, COORD_TYPE & cos_angle_ACB,
 bool & flag_small_magnitude)
{
  COORD_TYPE wBA[DIM3], wBC[DIM3], wAC[DIM3];
  COORD_TYPE magnitude;

  cos_angle_ABC = 1;
  cos_angle_ACB = 1;

  IJK::compute_unit_vector_3D
    (coordB, coordA, min_dist, wBA, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_unit_vector_3D
    (coordB, coordC, min_dist, wBC, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_unit_vector_3D
    (coordA, coordC, min_dist, wAC, magnitude, flag_small_magnitude);
  if (flag_small_magnitude) { return; }

  IJK::compute_inner_product_3D(wBA, wBC, cos_angle_ABC);
  IJK::compute_inner_product_3D(wAC, wBC, cos_angle_ACB);
}


// **************************************************
//  Check for violation of manifold conditions
// **************************************************

/// Check for manifold violations of 4 or more triangles on edge (cubeA, cubeB)
///   caused by mapping from from_cube to cubeA.
/// Return false if violation found.
bool check_single_edge_manifold
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX from_cube_index,
 const VERTEX_INDEX cubeA_index,
 const VERTEX_INDEX cubeB_index,
 const MERGESHARP::ISOVERT & isovert, 
 const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{
  if (flag_debug) {
    MSDEBUG();
    cerr << "In " << __func__ << endl;
    scalar_grid.PrintIndexAndCoord
      (cerr, "    Checking if mapping ", from_cube_index, " to ",
       cubeA_index, "\n");
    scalar_grid.PrintIndexAndCoord
      (cerr, "      creates non-manifold edge ", cubeA_index, " ", 
       cubeB_index, "\n");
  }

  if (cubeA_index == cubeB_index) { return(true); }

  bool flag_vertex_maps_to_cubeA, flag_vertex_maps_to_cubeB;
  bool flag_quad_maps_to_FABX;
  bool flag_quad_maps_to_FBAX_or_FBXA;
  bool flag_edge_maps_to_FB;
  determine_quad_maps
    (scalar_grid, isovalue, from_cube_index, cubeA_index, cubeB_index,
     isovert, gcube_map, flag_vertex_maps_to_cubeA, flag_vertex_maps_to_cubeB,
     flag_quad_maps_to_FABX, flag_quad_maps_to_FBAX_or_FBXA, 
     flag_edge_maps_to_FB);

  // *** DEBUGXXX ***
  if (flag_debug) {
    MSDEBUG();
    cerr << "            flag_quad_maps_to_FBAX_or_FBXA: "
         << int(flag_quad_maps_to_FBAX_or_FBXA) << endl;
    cerr << "            flag_quad_maps_to_FABX: "
         << int(flag_quad_maps_to_FABX) << endl;
    cerr << "            flag_edge_maps_to_FB: "
         << int(flag_edge_maps_to_FB) << endl;
  }

  if (!flag_vertex_maps_to_cubeA) {
    // Vertex in from_cube not connected to vertex in cubeA.
    return(false);
  }
  else if (flag_vertex_maps_to_cubeB) {

    if (flag_quad_maps_to_FBAX_or_FBXA) {
      return(true);
    } 
    else if (flag_quad_maps_to_FABX) {
      // Quad maps to quad with diagonal (from_cube, cubeB).
      if (flag_edge_maps_to_FB) {
        // (from_cube, cubeB) is an edge of some other quad.
        // Mapping from_cube to cubeA would create edge in 4 polygons.


        MSDEBUG();
        if (flag_debug) {
          scalar_grid.PrintIndexAndCoord
            (cerr, "  Cubes ", from_cube_index, "");
          scalar_grid.PrintIndexAndCoord
            (cerr, " and ", cubeB_index, " are diagonals of isoquad\n");
          scalar_grid.PrintIndexAndCoord
            (cerr, "    containing ", cubeA_index, "\n");
          scalar_grid.PrintIndexAndCoord
            (cerr, "  Edge ", from_cube_index, "");
          scalar_grid.PrintIndexAndCoord
            (cerr, " ", cubeB_index, " already exists.\n");
          scalar_grid.PrintIndexAndCoord
            (cerr, "  Mapping of ", from_cube_index, "");
          scalar_grid.PrintIndexAndCoord
            (cerr, " to ", cubeA_index, " fails.\n");
        }

        return(false);
      }
    }
    else if (are_connected_by_iso_quad
             (scalar_grid, cubeA_index, cubeB_index,
              isovalue, isovert, gcube_map, flag_extended)) {

      MSDEBUG();
      if (flag_debug) {
        scalar_grid.PrintIndexAndCoord
          (cerr, "  Cubes ", cubeA_index, "");
        scalar_grid.PrintIndexAndCoord
          (cerr, " and ", cubeB_index, " are connected by isoquad\n");
        scalar_grid.PrintIndexAndCoord
          (cerr, "    No isoquad connects ", cubeA_index,
           " ", cubeB_index, " and ", from_cube_index, "\n");
        scalar_grid.PrintIndexAndCoord
          (cerr, "  Mapping of ", from_cube_index, "");
        scalar_grid.PrintIndexAndCoord
          (cerr, " to ", cubeA_index, " fails.\n");
      }

      return(false); 
    }
  }

  return(true);
}


/// Check for manifold violations of 4 or more triangles on an edge
///   caused by mapping from from_cube to to_cube.
bool MERGESHARP::check_edge_manifold
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const MERGESHARP::ISOVERT & isovert, 
 const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{
  static CUBE_CONNECTED_ARRAY connected_sharp;

  find_connected_sharp
    (scalar_grid, isovalue, from_cube, isovert, gcube_map, connected_sharp);

  // *** DEBUGXXX ***
  /*
  if (from_cube == 1153131 || from_cube == 1153132) {
    MSDEBUG();
    scalar_grid.PrintIndexAndCoord
      (cerr, "$$$ NOT Mapping ", from_cube, " to " , to_cube, "\n");
    return(false);
  }
  */

  for (NUM_TYPE i = 0; i < connected_sharp.NumElements(); i++) {

    if (!check_single_edge_manifold
        (scalar_grid, isovalue, from_cube, to_cube, connected_sharp[i],
         isovert, gcube_map, flag_extended)) 
      { return(false); }
  }

  return(true);
}

// Check for manifold violations caused by edges between to_cube0 and cubeA
//   or to_cube1 and cubeA caused by mapping from_cube0 to to_cube0 
//   and from_cube1 to to_cube1.
bool check_edge_XA_manifoldII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX from_cube0_index,
 const VERTEX_INDEX to_cube0_index,
 const VERTEX_INDEX from_cube1_index,
 const VERTEX_INDEX to_cube1_index,
 const VERTEX_INDEX cubeA_index,
 const MERGESHARP::ISOVERT & isovert, 
 std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{

  // *** DEBUGXXX ***
  if (flag_debug) {
    MSDEBUG();
    cerr << "  In " << __func__ << endl;
  }

  // *** DEBUGXXX ***
  /*
  if (from_cube0_index == 1755446 || from_cube0_index == 755447 ||
      from_cube0_index == 777946 || from_cube0_index == 777947 ||
      from_cube0_index == 755296) {
    MSDEBUG();
    scalar_grid.PrintIndexAndCoord
      (cerr, "$$$ NOT Mapping ", from_cube0_index, " to " , 
       to_cube0_index, "\n");
    return(false);
  }

  if (from_cube1_index == 1755446 || from_cube1_index == 755447 ||
      from_cube1_index == 777946 || from_cube1_index == 777947 ||
      from_cube1_index == 755296) {
    MSDEBUG();
    scalar_grid.PrintIndexAndCoord
      (cerr, "$$$ NOT Mapping ", from_cube1_index, " to " , to_cube1_index, "\n");
    return(false);
  }
  */


  if (check_single_edge_manifold
      (scalar_grid, isovalue, from_cube0_index, to_cube0_index, cubeA_index,
       isovert, gcube_map, flag_extended)) {
    
    if (temp_map_check_single_edge_manifold
        (scalar_grid, isovalue, from_cube0_index, to_cube0_index, 
         from_cube1_index, to_cube1_index, cubeA_index, isovert, 
         gcube_map, flag_extended)) {

      return(true);
    }
  }

  if (check_single_edge_manifold
      (scalar_grid, isovalue, from_cube1_index, to_cube1_index, cubeA_index,
       isovert, gcube_map, flag_extended)) {

    if (temp_map_check_single_edge_manifold
        (scalar_grid, isovalue, from_cube1_index, to_cube1_index, 
         from_cube0_index, to_cube0_index, cubeA_index, isovert, 
         gcube_map, flag_extended)) {

      return(true);
    }
  }

  return(false);
}


// Check for manifold violations caused by edges between sharp cubes.
//   caused by mapping from_cube0 to to_cube0 and from_cube1 to to_cube1.
bool MERGESHARP::check_edge_manifoldII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX from_cube0_index,
 const VERTEX_INDEX to_cube0_index,
 const VERTEX_INDEX from_cube1_index,
 const VERTEX_INDEX to_cube1_index,
 const MERGESHARP::ISOVERT & isovert, 
 std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{
  IJK::PROCEDURE_ERROR error("check_edge_between_sharp_cubesII");
  static CUBE_CONNECTED_ARRAY connected_sharp;
  INDEX_DIFF_TYPE from_gcube0_index = 
    isovert.GCubeIndex(from_cube0_index, error);
  INDEX_DIFF_TYPE from_gcube1_index = 
    isovert.GCubeIndex(from_cube1_index, error);
  INDEX_DIFF_TYPE to_gcube0_index = 
    isovert.GCubeIndex(to_cube0_index, error);
  INDEX_DIFF_TYPE to_gcube1_index = 
    isovert.GCubeIndex(to_cube1_index, error);


  if (from_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
  if (from_gcube1_index == ISOVERT::NO_INDEX) { throw error; }
  if (to_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
  if (to_gcube1_index == ISOVERT::NO_INDEX) { throw error; }

  // *** DEBUGXXX ***
  if (flag_debug) {
    MSDEBUG();
    cerr << endl << "In " << __func__ << endl;
  }

  if (gcube_map[from_gcube0_index] == to_gcube0_index) {
    // gcube0 already mapped to to_gcube0.
    return(check_edge_manifold
           (scalar_grid, isovalue, from_cube1_index, to_cube1_index,
            isovert, gcube_map, flag_extended));
  }

  if (gcube_map[from_gcube1_index] == to_gcube1_index) {
    // gcube1 already mapped to to_gcube1.
    return(check_edge_manifold
           (scalar_grid, isovalue, from_cube0_index, to_cube0_index,
            isovert, gcube_map, flag_extended));
  }

  find_connected_sharp
    (scalar_grid, isovalue, from_cube0_index, isovert, 
     gcube_map, connected_sharp);
  add2list_connected_sharp
    (scalar_grid, isovalue, from_cube1_index, isovert, 
     gcube_map, connected_sharp);

  // *** DEBUGXXX ***
  if (flag_debug) {
    MSDEBUG();
    cerr << "  In " << __func__ << endl;
    cerr << "    connected_sharp: ";
    for (int i = 0; i < connected_sharp.NumElements(); i++)
      { cerr << " " << connected_sharp[i]; }
    cerr << endl;
  }

  for (NUM_TYPE i = 0; i < connected_sharp.NumElements(); i++) {

    const VERTEX_INDEX cubeA_index = connected_sharp[i];

    if (to_cube0_index == to_cube1_index) {
      if (to_cube0_index == cubeA_index) { 
        // Case I: to_cube0 == to_cube1 == cubeA_index
        continue; 
      }
    }

    if (!check_edge_XA_manifoldII
        (scalar_grid, isovalue, from_cube0_index, to_cube0_index, 
         from_cube1_index, to_cube1_index, cubeA_index,
         isovert, gcube_map, flag_extended)) {

      MSDEBUG();
      if (flag_debug) {
        isovert.grid.PrintIndexAndCoord
          (cerr, "xxx  Both maps fail. From: ", from_cube0_index, " to ", 
           to_cube0_index, "\n");
        isovert.grid.PrintIndexAndCoord
          (cerr, "                  and ", from_cube1_index, " to ", 
           to_cube1_index, "\n");
      }

      return(false);
    }

  }

  return(true);
}


/* OBSOLETE
namespace {

  // check for manifold violations on a single edge with endpoints
  //   in (to_cube,cubeA_index) caused by mapping three cubes to to_cube.
  bool check_single_edge_manifoldIII_old
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube_index[3],
   const VERTEX_INDEX to_cube,
   const VERTEX_INDEX cubeA_index,
   const MERGESHARP::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    INDEX_DIFF_TYPE gcube_index[3], to_gcube;
    NUM_TYPE store_map[3];
    bool flag_maps_to_cubeA[3];
    bool flag_maps_to_cubeT[3];
    bool flag_edge_maps_to_both_cubes[3];

    IJK::PROCEDURE_ERROR error("check_single_edge_manifoldIII");

    to_gcube = isovert.GCubeIndex(to_cube, error);
    for (int i = 0; i < 3; i++) {
      gcube_index[i] = isovert.GCubeIndex(cube_index[i], error); 
      if (gcube_index[i] == ISOVERT::NO_INDEX) { throw error; }
      store_map[i] = gcube_map[gcube_index[i]];
    }

    for (int i0 = 0; i0 < 3; i0++) {

      determine_quad_edge_map
        (scalar_grid, isovalue, cube_index[i0], cubeA_index, to_cube,
         isovert, gcube_map, flag_maps_to_cubeA[0], flag_maps_to_cubeT[0],
         flag_edge_maps_to_both_cubes[0]);

      if (!flag_maps_to_cubeA[0] || flag_edge_maps_to_both_cubes[0]) {

        gcube_map[gcube_index[i0]] = to_gcube;

        for (int k = 1; k <= 2; k++) {

          int i1 = (i0+k)%3;
          int i2 = (i0+2*k)%3;

          determine_quad_edge_map
            (scalar_grid, isovalue, cube_index[i1], cubeA_index, to_cube,
             isovert, gcube_map, flag_maps_to_cubeA[1], flag_maps_to_cubeT[1],
             flag_edge_maps_to_both_cubes[1]);

          if (!flag_maps_to_cubeA[1] || flag_edge_maps_to_both_cubes[1]) {
            gcube_map[gcube_index[i1]] = to_gcube;

            determine_quad_edge_map
              (scalar_grid, isovalue, cube_index[i2], cubeA_index, to_cube,
               isovert, gcube_map, flag_maps_to_cubeA[2], flag_maps_to_cubeT[2],
               flag_edge_maps_to_both_cubes[2]);

            if (!flag_maps_to_cubeA[2] || flag_edge_maps_to_both_cubes[2]) {

              // restore gcube_map
              gcube_map[gcube_index[i0]] = store_map[i0];
              gcube_map[gcube_index[i1]] = store_map[i1];
              return(true);
            }
          }

          // restore gcube_map[i1]
          gcube_map[gcube_index[i1]] = store_map[i1];        
        }
      }

      // restore gcube_map
      gcube_map[gcube_index[i0]] = store_map[i0];
    }

    return(false);
  }

}
*/

// check for manifold violations on a single edge with endpoints
//   in (to_cube,cubeA_index) caused by mapping three cubes to to_cube.
bool check_single_edge_manifoldIII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index[3],
 const VERTEX_INDEX to_cube,
 const VERTEX_INDEX cubeA_index,
 const MERGESHARP::ISOVERT & isovert, 
 std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{
  if (cubeA_index == to_cube) { return(true); }

  for (int i0 = 0; i0 < 3; i0++) {

    if (check_single_edge_manifold
        (scalar_grid, isovalue, cube_index[i0], to_cube, cubeA_index,
         isovert, gcube_map, flag_extended)) {

      NUM_TYPE i1 = (i0+1)%3;
      NUM_TYPE i2 = (i0+2)%3;

      if (temp_map_check_edge_XA_manifoldII
          (scalar_grid, isovalue, cube_index[i0], to_cube, 
           cube_index[i1], to_cube, cube_index[i2], to_cube, cubeA_index,
           isovert, gcube_map, flag_extended))
        { return(true); }
    }
  }

  return(false);
}


// check for manifold violations caused by edges between sharp cubes.
//   caused by mapping three cubes to to_cube.
bool MERGESHARP::check_edge_manifoldIII
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index[3],
 const VERTEX_INDEX to_cube,
 const MERGESHARP::ISOVERT & isovert, 
 std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
 const bool flag_extended)
{
  INDEX_DIFF_TYPE gcube_index[3], to_gcube_index;
  static CUBE_CONNECTED_ARRAY connected_sharp;
  IJK::PROCEDURE_ERROR error("check_edge_between_sharp_cubesIII");

  to_gcube_index = isovert.GCubeIndex(to_cube, error);
  for (int i = 0; i < 3; i++) {
    gcube_index[i] = isovert.GCubeIndex(cube_index[i], error); 
    if (gcube_index[i] == ISOVERT::NO_INDEX) { throw error; }
  }

  find_connected_sharp
    (scalar_grid, isovalue, cube_index[0], isovert, gcube_map, 
     connected_sharp);
  add2list_connected_sharp
    (scalar_grid, isovalue, cube_index[1], isovert, 
     gcube_map, connected_sharp);
  add2list_connected_sharp
    (scalar_grid, isovalue, cube_index[2], isovert, 
     gcube_map, connected_sharp);

  for (NUM_TYPE i = 0; i < connected_sharp.NumElements(); i++) {

    const VERTEX_INDEX cubeA_index = connected_sharp[i];

    if (cubeA_index == to_cube) { continue; }

    /* OBSOLETE
    if (!are_connected_by_iso_quad
        (scalar_grid, cubeA_index, to_cube, isovalue, isovert, gcube_map, 
         flag_extended)) { continue; }
    */

    if (!check_single_edge_manifoldIII
        (scalar_grid, isovalue, cube_index, to_cube, cubeA_index,
         isovert, gcube_map, flag_extended))
      { return(false); }
  }

  return(true);
}


// **************************************************
//  Misc routines
// **************************************************

/// Get cubes in clockwise order around edge.
/// @pre edge is in the interior of the grid.
void MERGESHARP::get_cubes_around_edge
(const SHARPISO_GRID & grid, const VERTEX_INDEX iend0,
 const int edge_dir, VERTEX_INDEX cube_index[NUM_VERT_PER_QUAD])
{
  const int d1 = (edge_dir+1)%DIM3;
  const int d2 = (edge_dir+2)%DIM3;

  cube_index[0] =
    iend0 - grid.AxisIncrement(d1) - grid.AxisIncrement(d2);
  cube_index[1] = iend0 - grid.AxisIncrement(d2);
  cube_index[2] = iend0;
  cube_index[3] = iend0 - grid.AxisIncrement(d1);
}


// **************************************************
//  Local routines
// **************************************************

namespace {

  bool is_strictly_between
  (const COORD_TYPE x, const COORD_TYPE a0, const COORD_TYPE a1)
  {
    if ((a0 < x) && (x < a1)) { return(true); }
    if ((a0 > x) && (x > a1)) { return(true); }
    return(false);
  }

  bool is_between
  (const COORD_TYPE x, const COORD_TYPE a0, const COORD_TYPE a1)
  {
    if ((a0 <= x) && (x <= a1)) { return(true); }
    if ((a0 >= x) && (x >= a1)) { return(true); }
    return(false);
  }

  bool is_covered_A_or_corner(const GRID_CUBE_FLAG flag)
  {
    if (flag == COVERED_A_GCUBE || flag == COVERED_CORNER_GCUBE)
      { return(true); }
    else
      { return(false); }
  }

  // If a cube is an 3x3x3 region around a selected cube,
  //   then the cube is COVERED_A_GCUBE or COVERED_CORNER_GCUBE
  //   or SELECTED_GCUBE.
  bool is_in_3x3x3_region(const GRID_CUBE_FLAG flag)
  {
    if (flag == COVERED_A_GCUBE || flag == COVERED_CORNER_GCUBE ||
        flag == SELECTED_GCUBE)
      { return(true); }
    else
      { return(false); }
  }

  bool are_separated
  (const SHARPISO_GRID & grid,
   const GRID_CUBE_DATA & gcubeA, const COORD_TYPE coordA[DIM3],
   const GRID_CUBE_DATA & gcubeB, const COORD_TYPE coordB[DIM3],
   const int orth_dir) 
  {
    if (gcubeA.cube_coord[orth_dir] < gcubeB.cube_coord[orth_dir]) {
      COORD_TYPE x = gcubeB.cube_coord[orth_dir];
      if (coordA[orth_dir] < x && x < coordB[orth_dir])
        { return(true); }
    }
    else if (gcubeA.cube_coord[orth_dir] > gcubeB.cube_coord[orth_dir]) {
      COORD_TYPE x = gcubeB.cube_coord[orth_dir]+1;
      if (coordA[orth_dir] > x && x > coordB[orth_dir])
        { return(true); }
    }

    return(false);
  }

  /// Compute (coordA-coordB), project to plane orthogonal to orth_dir,
  ///   and normalize vector.
  /// @param[out] w[] Resulting vector.
  void subtract_project_normalize
  (const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3], 
   const int orth_dir, const COORD_TYPE max_small_magnitude,
   COORD_TYPE w[DIM3], COORD_TYPE & magnitude, bool & flag_zero)
  {
    IJK::subtract_coord_3D(coordA, coordB, w);

    // Project onto plane orthogonal to orth_dir.
    w[orth_dir] = 0;

    IJK::normalize_vector
      (DIM3, w, max_small_magnitude, w, magnitude, flag_zero);
  }

	/// Compute the overlap region between two cube indices
  /// @param dist2boundary Distance from cube to region boundary
  //         region size = (2*dist2boundary+1)
	bool find_overlap(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2,
    const GRID_COORD_TYPE dist2boundary,
		GRID_COORD_TYPE rmin[],
		GRID_COORD_TYPE rmax[])
	{
		COORD_TYPE coord1[DIM3], coord2[DIM3];
		scalar_grid.ComputeCoord(cube_index1, coord1);
		scalar_grid.ComputeCoord(cube_index2, coord2);
		for (int d=0;d<DIM3;d++){
			rmin[d] = std::max(coord1[d]-dist2boundary, coord2[d]-dist2boundary);
      rmin[d] = std::max(rmin[d], 0);
			rmax[d] = std::min(coord1[d]+dist2boundary+1, coord2[d]+dist2boundary+1);
      rmax[d] = std::min(rmax[d], scalar_grid.AxisSize(d)-1);
		}

		// track the dimension of the tracked,
		// if the tracked regions has at least 2 dimension
		int dim_of_overlap=0;
		for (int d=0;d<DIM3;d++)
		{
			if(rmin[d] > rmax[d])
			{ return false; }
			if(rmin[d] < rmax[d])
				dim_of_overlap++;
		}
		if (dim_of_overlap>=2)
			return true;
		else
			return false;
	}

	bool find_overlap(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2,
		GRID_COORD_TYPE rmin[],
		GRID_COORD_TYPE rmax[])
  {
    return(find_overlap
           (scalar_grid, cube_index1, cube_index2, 1, rmin, rmax));
  }

  /// Determine if isosurface quad dual to grid edge (iend0, edge_dir)
  ///   maps to cube0 or cube1.
  /// @param iend0 Lower endpoint of edge.
  /// @param edge_dir Edge direction.
  /// @pre Edge (iend0, edge_dir) is an internal edge.
  /// @pre Edge (iend0, edge_dir) is bipolar.
  void determine_quad_map_to_cubes
  (const VERTEX_INDEX iend0, const int edge_dir,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_maps_to_cube0,
   bool & flag_maps_to_cube1)
  {
    const NUM_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);
    const NUM_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);

    // Initialize
    flag_maps_to_cube0 = false;
    flag_maps_to_cube1 = false;

    int d1 = (edge_dir+1)%DIM3;
    int d2 = (edge_dir+2)%DIM3;

    // *** DEBUGXXX ***
    MSDEBUG();
    VERTEX_INDEX iend1 = isovert.grid.NextVertex(iend0, edge_dir);
    if (flag_debug) {
      cerr << "      In " << __func__ << endl;
      isovert.grid.PrintIndexAndCoord
        (cerr, "        Checking quad around edge ", iend0, " ", iend1, "\n");
    }

    VERTEX_INDEX iv0 = 
      iend0 - isovert.grid.AxisIncrement(d1) - isovert.grid.AxisIncrement(d2);
    for (int j = 0; j < NUM_VERT_PER_QUAD; j++) {
      VERTEX_INDEX cube2_index = isovert.grid.FacetVertex(iv0, edge_dir, j);
      NUM_TYPE gcube2_index = isovert.GCubeIndex(cube2_index);

      NUM_TYPE to_gcube = gcube_map[gcube2_index];
      if (to_gcube == gcube0_index) { 
        flag_maps_to_cube0 = true; 

        // *** DEBUGXXX ***
        if (flag_debug) {
          MSDEBUG();
          isovert.grid.PrintIndexAndCoord
            (cerr, "          Cube ", cube2_index, " maps to ", 
             cube0_index, "\n");
        }

      }
      else if (to_gcube == gcube1_index) { 
        flag_maps_to_cube1 = true; 
      }
    }

  }

  /// Determine if some isosurface quad dual to an edge of cube0
  ///   maps to both cube1 and cube2.
  void determine_quad_map_to_both_cubes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_maps_to_cube1,
   bool & flag_maps_to_both_cubes)
  {
    const INDEX_DIFF_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);

    flag_maps_to_cube1 = false;
    flag_maps_to_both_cubes = false;

    if (gcube0_index == ISOVERT::NO_INDEX) { return; }

    if (flag_debug) {
      MSDEBUG();
      cerr << "      In " << __func__ << endl;
      scalar_grid.PrintIndexAndCoord
        (cerr, "        cube0: ", cube0_index, "\n");
    }

    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[gcube0_index].boundary_bits;

    if (boundary_bits == 0) {
      for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

        for (NUM_TYPE j = 0; j < NUM_CUBE_FACET_VERTICES3D; j++) {

          VERTEX_INDEX iend0 = 
            scalar_grid.FacetVertex(cube0_index, edge_dir, j);
          VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

          if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

            bool flagB_maps_to_cube1;
            bool flagB_maps_to_cube2;
            determine_quad_map_to_cubes
              (iend0, edge_dir, cube1_index, cube2_index,
               isovert, gcube_map, flagB_maps_to_cube1, flagB_maps_to_cube2);

            if (flagB_maps_to_cube1) {
              flag_maps_to_cube1 = true;

              if (flagB_maps_to_cube2) {
                flag_maps_to_both_cubes = true;
                return;
              }
            }
          }
        }
      }
    }
    else {
      // Handle boundary case.
    }
  }

  /// Determine if an edge of the isosurface quad dual 
  ///   to grid edge (iend0, edge_dir) maps to both cube0 and cube1.
  /// @param iend0 Lower endpoint of edge.
  /// @param edge_dir Edge direction.
  /// @pre Edge (iend0, edge_dir) is an internal edge.
  /// @pre Edge (iend0, edge_dir) is bipolar.
  void determine_single_quad_edge_map
  (const VERTEX_INDEX iend0, const int edge_dir,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_vertex_maps_to_cube0,
   bool & flag_vertex_maps_to_cube1,
   bool & flag_edge_maps_to_both_cubes)
  {
    IJK::PROCEDURE_ERROR error("determine_single_quad_edge_map");
    const NUM_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);
    const NUM_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
    VERTEX_INDEX cube_index[NUM_VERT_PER_QUAD];
    INDEX_DIFF_TYPE gcube_index[NUM_VERT_PER_QUAD];
    NUM_TYPE to_gcube_index[NUM_VERT_PER_QUAD];

    // Initialize
    flag_vertex_maps_to_cube0 = false;
    flag_vertex_maps_to_cube1 = false;
    flag_edge_maps_to_both_cubes = false;

    // *** DEBUGXXX ***
    MSDEBUG();    
    VERTEX_INDEX iend1 = isovert.grid.NextVertex(iend0, edge_dir);
    if (flag_debug) {
      cerr << "      In " << __func__ << endl;
      isovert.grid.PrintIndexAndCoord
        (cerr, "        Checking quad around edge ", iend0, " ", iend1, "\n");
    }

    get_cubes_around_edge(isovert.grid, iend0, edge_dir, cube_index);

    for (int j = 0; j < NUM_VERT_PER_QUAD; j++) {
      gcube_index[j] = isovert.GCubeIndex(cube_index[j], error);
      if (gcube_index[j] == ISOVERT::NO_INDEX) { throw error; }

      to_gcube_index[j] = gcube_map[gcube_index[j]];
    }

    for (int j = 0; j < NUM_VERT_PER_QUAD; j++) {

      if (to_gcube_index[j] == gcube0_index)
        { flag_vertex_maps_to_cube0 = true; }

      if (to_gcube_index[j] == gcube1_index)
        { flag_vertex_maps_to_cube1 = true; }
    }

    for (int j0 = 0; j0 < NUM_VERT_PER_QUAD; j0++) {

      int j1 = (j0+1)%NUM_VERT_PER_QUAD;

      if (to_gcube_index[j0] == gcube0_index && 
          to_gcube_index[j1] == gcube1_index) {

        // *** DEBUGXXX **
        if (flag_debug) {
          MSDEBUG();

          isovert.grid.PrintIndexAndCoord
            (cerr, "        Cube ", cube_index[j0], " maps to ", 
             cube0_index, "\n");
          isovert.grid.PrintIndexAndCoord
            (cerr, "        Cube ", cube_index[j1], " maps to ", 
             cube1_index, "\n");
        }

        flag_edge_maps_to_both_cubes = true;
      }
      else if (to_gcube_index[j0] == gcube1_index && 
               to_gcube_index[j1] == gcube0_index) {

        // *** DEBUGXXX **
        if (flag_debug) {
          MSDEBUG();

          isovert.grid.PrintIndexAndCoord
            (cerr, "        Cube ", cube_index[j0], " maps to ", 
             cube1_index, "\n");
          isovert.grid.PrintIndexAndCoord
            (cerr, "        Cube ", cube_index[j1], " maps to ", 
             cube0_index, "\n");
        }


        flag_edge_maps_to_both_cubes = true;
      }
    }

  }

  /// Determine if the endpoints of some edge of an isosurface quad dual 
  ///   to an edge of cube0 map to both cube1 and cube2.
  void determine_quad_edge_map
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_vertex_maps_to_cube1,
   bool & flag_vertex_maps_to_cube2,
   bool & flag_edge_maps_to_both_cubes)
  {
    const INDEX_DIFF_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);

    flag_vertex_maps_to_cube1 = false;
    flag_vertex_maps_to_cube2 = false;
    flag_edge_maps_to_both_cubes = false;

    if (gcube0_index == ISOVERT::NO_INDEX) { return; }

    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[gcube0_index].boundary_bits;

    if (boundary_bits == 0) {
      for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

        for (NUM_TYPE j = 0; j < NUM_CUBE_FACET_VERTICES3D; j++) {

          VERTEX_INDEX iend0 = 
            scalar_grid.FacetVertex(cube0_index, edge_dir, j);
          VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

          if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

            bool flagB_maps_to_cube1;
            bool flagB_maps_to_cube2;
            bool flagB_edge_maps_to_both_cubes;
            determine_single_quad_edge_map
              (iend0, edge_dir, cube1_index, cube2_index,
               isovert, gcube_map, flagB_maps_to_cube1, flagB_maps_to_cube2,
               flagB_edge_maps_to_both_cubes);

              // *** DEBUGXXX ***
              MSDEBUG();
              if (flag_debug) {
                isovert.grid.PrintIndexAndCoord
                  (cerr, "          Flag maps to ", cube1_index, " : ");
                cerr << int(flagB_maps_to_cube1) << endl;
                isovert.grid.PrintIndexAndCoord
                  (cerr, "          Flag maps to ", cube2_index, " : ");
                cerr << int(flagB_maps_to_cube2) << endl;
              }

            if (flagB_maps_to_cube1)
              { flag_vertex_maps_to_cube1 = true; }
            if (flagB_maps_to_cube2)
              { flag_vertex_maps_to_cube2 = true; }
            if (flagB_edge_maps_to_both_cubes) {
              flag_edge_maps_to_both_cubes = true; 

              // *** DEBUG ***
              MSDEBUG();
              if (flag_debug) {
                cerr << "    In " << __func__
                     << " maps to cube1 " << int(flag_vertex_maps_to_cube1)
                     << " maps to both cubes " 
                     << int(flag_edge_maps_to_both_cubes) << endl;
              }

              return;
            }
          }
        }
      }
    }
    else {
      // Handle boundary case.
    }
  }

  /// Temporarily map from_gcube to to_gcube and then determine if some 
  ///   isosurface quad dual to an edge of cube0 maps to both cube1 and cube2.
  void temp_map_determine_quad_map_to_both_cubes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const NUM_TYPE from_gcube_index,
   const NUM_TYPE to_gcube_index,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_maps_to_cube1,
   bool & flag_maps_to_both_cubes)
  {
    const NUM_TYPE store_from_map = gcube_map[from_gcube_index];

    // *** DEBUGXXX ***
    if (flag_debug) {
      MSDEBUG();
      isovert.grid.PrintIndexAndCoord
        (cerr, "    Temporarily mapping ",
         isovert.CubeIndex(from_gcube_index), " to ",
         isovert.CubeIndex(to_gcube_index), "\n");
    }

    // temporarily set gcube_map[from_gcube_index] to to_gcube_index.
    gcube_map[from_gcube_index] = to_gcube_index;
    determine_quad_map_to_both_cubes
      (scalar_grid, isovalue, cube0_index, cube1_index, cube2_index,
       isovert, gcube_map, flag_maps_to_cube1, flag_maps_to_both_cubes);

    // restore gcube_map[from_gcube_index]
    gcube_map[from_gcube_index] = store_from_map;
  }

  /// Temporarily map from_cube0 to to_cube0 and 
  ///   then check for manifold violations of 4 or more triangles 
  ///   on edge (cubeA, cubeB) caused by mapping from from_cube1 to cubeA.
  bool temp_map_check_single_edge_manifold
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube0_index,
   const VERTEX_INDEX to_cube0_index,
   const VERTEX_INDEX from_cube1_index,
   const VERTEX_INDEX cubeA_index,
   const VERTEX_INDEX cubeB_index,
   const ISOVERT & isovert,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended)
  {
    IJK::PROCEDURE_ERROR error("temp_map_check_edge_manifold");
    INDEX_DIFF_TYPE from_gcube0_index = 
      isovert.GCubeIndex(from_cube0_index, error);
    INDEX_DIFF_TYPE to_gcube0_index = 
      isovert.GCubeIndex(to_cube0_index, error);

    if (from_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
    if (to_gcube0_index == ISOVERT::NO_INDEX) { throw error; }

    const NUM_TYPE store_from0_map = gcube_map[from_gcube0_index];

    // *** DEBUGXXX ***
    if (flag_debug) {
      MSDEBUG();
      isovert.grid.PrintIndexAndCoord
        (cerr, "      Temporarily mapping ",
         from_cube0_index, " to ", to_cube0_index, "\n");
    }

    // temporarily set gcube_map[from_gcube_index] to to_gcube_index.
    gcube_map[from_gcube0_index] = to_gcube0_index;

    bool flag = 
      check_single_edge_manifold
      (scalar_grid, isovalue, from_cube1_index, cubeA_index, cubeB_index,
       isovert, gcube_map, flag_extended);

    // restore gcube_map[from_gcube0_index]
    gcube_map[from_gcube0_index] = store_from0_map;

    return(flag);
  }

  /// Temporarily map from_cube0 to to_cube0 and 
  ///   then check for manifold violations caused by edges 
  ///   between to_cube1 and cubeA or to_cube2 and cubeA caused 
  ///   by mapping from_cube1 to to_cube1 and from_cube2 to to_cube2.
  bool temp_map_check_edge_XA_manifoldII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube0_index,
   const VERTEX_INDEX to_cube0_index,
   const VERTEX_INDEX from_cube1_index,
   const VERTEX_INDEX to_cube1_index,
   const VERTEX_INDEX from_cube2_index,
   const VERTEX_INDEX to_cube2_index,
   const VERTEX_INDEX cubeA_index,
   const MERGESHARP::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended)
  {
    IJK::PROCEDURE_ERROR error("temp_map_check_edge_manifold");
    INDEX_DIFF_TYPE from_gcube0_index = 
      isovert.GCubeIndex(from_cube0_index, error);
    INDEX_DIFF_TYPE to_gcube0_index = 
      isovert.GCubeIndex(to_cube0_index, error);

    if (from_gcube0_index == ISOVERT::NO_INDEX) { throw error; }
    if (to_gcube0_index == ISOVERT::NO_INDEX) { throw error; }

    const NUM_TYPE store_from0_map = gcube_map[from_gcube0_index];

    // *** DEBUGXXX ***
    if (flag_debug) {
      MSDEBUG();
      isovert.grid.PrintIndexAndCoord
        (cerr, "      Temporarily mapping ",
         from_cube0_index, " to ", to_cube0_index, "\n");
    }

    // temporarily set gcube_map[from_gcube_index] to to_gcube_index.
    gcube_map[from_gcube0_index] = to_gcube0_index;

    bool flag = 
      check_edge_XA_manifoldII
      (scalar_grid, isovalue, from_cube1_index, to_cube1_index, 
       from_cube2_index, to_cube2_index, cubeA_index,
       isovert, gcube_map, flag_extended);

    // restore gcube_map[from_gcube0_index]
    gcube_map[from_gcube0_index] = store_from0_map;

    return(flag);
  }

	// Return true if two cubes are connected by an isosurface quadrilaterals.
	bool are_connected_by_iso_quad
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX & cube0_index,
   const VERTEX_INDEX & cube1_index,
   const SCALAR_TYPE isovalue,
   const MERGESHARP::ISOVERT & isovert, 
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended)
  {
    int dist2boundary;      // Distance from cube to region boundary
		GRID_COORD_TYPE rmin[DIM3], rmax[DIM3];

    // *** DEBUGXXX ***
    MSDEBUG();
    if (flag_debug) {
      cerr << "      In " << __func__ << endl;
    }

    if (flag_extended) 
      { dist2boundary = 2; }
    else
      { dist2boundary = 1; }

		// find the overlap region
		bool is_overlap =
			find_overlap(scalar_grid, cube0_index, cube1_index, 
                   dist2boundary, rmin, rmax);

		if (!is_overlap) { return false; }

		VERTEX_INDEX vbase = scalar_grid.ComputeVertexIndex(rmin);

    for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

      int d1 = (edge_dir+1)%DIM3;
      int d2 = (edge_dir+2)%DIM3;

      VERTEX_INDEX v0 = vbase;
      for (GRID_COORD_TYPE x0 = rmin[edge_dir]; x0 < rmax[edge_dir]; x0++) {

        VERTEX_INDEX v1 = v0;
        for (GRID_COORD_TYPE x1 = rmin[d1]; x1 <= rmax[d1]; x1++) {

          if (x1 == 0 || x1 >= scalar_grid.AxisSize(d1)) { continue; }

          VERTEX_INDEX v2 = v1;
          for (GRID_COORD_TYPE x2 = rmin[d2]; x2 <= rmax[d2]; x2++) {

            if (x2 == 0 || x2 >= scalar_grid.AxisSize(d2)) { continue; }

            const VERTEX_INDEX iend0 = v2;
            const VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

            if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

              bool flag_maps_to_cube0;
              bool flag_maps_to_cube1;

              determine_quad_map_to_cubes
                (iend0, edge_dir, cube0_index, cube1_index,
                 isovert, gcube_map, flag_maps_to_cube0, flag_maps_to_cube1);

              // *** DEBUGXXX ***
              MSDEBUG();
              if (flag_debug) {
                isovert.grid.PrintIndexAndCoord
                  (cerr, "          Flag maps to ", cube0_index, " : ");
                cerr << int(flag_maps_to_cube0) << endl;
                isovert.grid.PrintIndexAndCoord
                  (cerr, "          Flag maps to ", cube1_index, " : ");
                cerr << int(flag_maps_to_cube1) << endl;
              }

              if (flag_maps_to_cube0 && flag_maps_to_cube1) {
                return(true); 
              }
            }

            v2 = scalar_grid.NextVertex(v2, d2);
          }

          v1 = scalar_grid.NextVertex(v1, d1);
        }

        v0 = scalar_grid.NextVertex(v0, edge_dir);
      }
    }

    return(false);
  }


  /// Determine if isosurface quad dual to grid edge (iend0, edge_dir)
  ///   maps to cube0, cube1, or cube2
  /// @param iend0 Lower endpoint of edge.
  /// @param edge_dir Edge direction.
  /// @pre Edge (iend0, edge_dir) is an internal edge.
  /// @pre Edge (iend0, edge_dir) is bipolar.
  void determine_quad_map_to_cubes
  (const VERTEX_INDEX iend0, const int edge_dir,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool flag_maps_to_cube[3])
  {
    IJK::PROCEDURE_ERROR error("determine_quad_map_to_cubes");
    INDEX_DIFF_TYPE gcube_index[3];

    // Initialize
    for (int i = 0; i < 3; i++) 
      { flag_maps_to_cube[i] = false; }

    gcube_index[0] = isovert.GCubeIndex(cube0_index, error);
    gcube_index[1] = isovert.GCubeIndex(cube1_index, error);
    gcube_index[2] = isovert.GCubeIndex(cube2_index, error);

    for (int i = 0; i < 3; i++) {
      if (gcube_index[i] == ISOVERT::NO_INDEX) { throw error; }
    }

    int d1 = (edge_dir+1)%DIM3;
    int d2 = (edge_dir+2)%DIM3;

    // *** DEBUGXXX ***
    MSDEBUG();
    VERTEX_INDEX iend1 = isovert.grid.NextVertex(iend0, edge_dir);
    if (flag_debug) {
      cerr << "      In " << __func__ << endl;
      isovert.grid.PrintIndexAndCoord
        (cerr, "        Checking quad around edge ", iend0, " ", iend1, "\n");
    }

    VERTEX_INDEX iv0 = 
      iend0 - isovert.grid.AxisIncrement(d1) - isovert.grid.AxisIncrement(d2);
    for (int j = 0; j < NUM_VERT_PER_QUAD; j++) {
      VERTEX_INDEX cubeA_index = isovert.grid.FacetVertex(iv0, edge_dir, j);
      NUM_TYPE gcubeA_index = isovert.GCubeIndex(cubeA_index);

      NUM_TYPE to_gcube = gcube_map[gcubeA_index];
      for (int i = 0; i < 3; i++) {
        if (to_gcube == gcube_index[i]) { 
          flag_maps_to_cube[i] = true;
        }
      }
    }
  }

  /// Determine if isosurface quad dual to grid edge (iend0, edge_dir)
  ///   maps to cube1, or cube2
  /// Set flag_edge_maps_to_cube2 to true if some edge dual to a facet
  ///   of cube0 maps to cube2.
  /// @param iend0 Lower endpoint of edge.
  /// @param edge_dir Edge direction.
  /// @pre Edge (iend0, edge_dir) is an internal edge.
  /// @pre Edge (iend0, edge_dir) is bipolar.
  void determine_single_quad_map
  (const VERTEX_INDEX iend0, const int edge_dir,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_quad_vertex_maps_to_cube1,
   bool & flag_quad_vertex_maps_to_cube2,
   bool & flag_edge_endpoint_maps_to_cube2)
  {
    IJK::PROCEDURE_ERROR error("determine_quad_map_to_cubes");
    VERTEX_INDEX quad_cube_index[NUM_VERT_PER_QUAD];
    NUM_TYPE quad_gcube_index[NUM_VERT_PER_QUAD];
    INDEX_DIFF_TYPE gcube1_index, gcube2_index;

    // Initialize
    flag_quad_vertex_maps_to_cube1 = false;
    flag_quad_vertex_maps_to_cube2 = false;
    flag_edge_endpoint_maps_to_cube2 = false;

    gcube1_index = isovert.GCubeIndex(cube1_index, error);
    gcube2_index = isovert.GCubeIndex(cube2_index, error);

    if (gcube1_index == ISOVERT::NO_INDEX) { throw error; }
    if (gcube2_index == ISOVERT::NO_INDEX) { throw error; }

    // *** DEBUGXXX ***
    MSDEBUG();
    VERTEX_INDEX iend1 = isovert.grid.NextVertex(iend0, edge_dir);
    if (flag_debug) {
      cerr << "      In " << __func__ << endl;
      isovert.grid.PrintIndexAndCoord
        (cerr, "        Checking quad around edge ", iend0, " ", iend1, "\n");
    }

    get_cubes_around_edge(isovert.grid, iend0, edge_dir, quad_cube_index);

    for (int i = 0; i < NUM_VERT_PER_QUAD; i++) {
      quad_gcube_index[i] = isovert.GCubeIndex(quad_cube_index[i]);
    }

    for (int j0 = 0; j0 < NUM_VERT_PER_QUAD; j0++) {

      const int j1 = (j0+1)%NUM_VERT_PER_QUAD;

      const NUM_TYPE gcubeA_index = quad_gcube_index[j0];
      const NUM_TYPE gcubeB_index = quad_gcube_index[j1];
      const NUM_TYPE to_gcubeA = gcube_map[gcubeA_index];
      const NUM_TYPE to_gcubeB = gcube_map[gcubeB_index];

      if (gcube1_index == to_gcubeA) { flag_quad_vertex_maps_to_cube1 = true; }
      if (gcube2_index == to_gcubeA) { flag_quad_vertex_maps_to_cube2 = true; }

      if (cube0_index == quad_cube_index[j0] &&
          gcube2_index == to_gcubeB)
        { flag_edge_endpoint_maps_to_cube2 = true; }

      if (cube0_index == quad_cube_index[j1] &&
          gcube2_index == to_gcubeA)
        { flag_edge_endpoint_maps_to_cube2 = true; }
    }

  }


  /// Determine if some isosurface quad dual to an edge of cube0
  ///   maps to cube1 and cube2.
  /// @param flag_quad_maps_to_012X True if vertices of some quad
  ///   map to cube0, cube1, cube2 and no edge of that quad has endpoints
  ///   mapping to cube0, cube2.
  /// @param flag_quad_maps_to_021X_or_02X1 True if vertices of some quad
  ///   map to cube0, cube1, cube2 and endpoints of some edge of that quad
  ///   map to cube0, cube2.
  void determine_quad_maps
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube0_index,
   const VERTEX_INDEX cube1_index,
   const VERTEX_INDEX cube2_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   bool & flag_vertex_maps_to_cube1,
   bool & flag_vertex_maps_to_cube2,
   bool & flag_quad_maps_to_012X,
   bool & flag_quad_maps_to_021X_or_02X1,
   bool & flag_edge_maps_to_02)
  {
    const INDEX_DIFF_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);

    flag_vertex_maps_to_cube1 = false;
    flag_vertex_maps_to_cube2 = false;
    flag_quad_maps_to_012X = false;
    flag_quad_maps_to_021X_or_02X1 = false;
    flag_edge_maps_to_02 = false;

    if (gcube0_index == ISOVERT::NO_INDEX) { return; }


    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[gcube0_index].boundary_bits;

    if (boundary_bits == 0) {
      for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

        for (NUM_TYPE j = 0; j < NUM_CUBE_FACET_VERTICES3D; j++) {

          VERTEX_INDEX iend0 = 
            scalar_grid.FacetVertex(cube0_index, edge_dir, j);
          VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

          if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

            bool flagB_quad_vertex_maps_to_cube1;
            bool flagB_quad_vertex_maps_to_cube2;
            bool flagB_edge_endpoint_maps_to_cube2;
            determine_single_quad_map
              (iend0, edge_dir, cube0_index, cube1_index, cube2_index,
               isovert, gcube_map, flagB_quad_vertex_maps_to_cube1,
               flagB_quad_vertex_maps_to_cube2,
               flagB_edge_endpoint_maps_to_cube2);

            if (flagB_quad_vertex_maps_to_cube1) 
              { flag_vertex_maps_to_cube1 = true; }

            if (flagB_quad_vertex_maps_to_cube2) 
              { flag_vertex_maps_to_cube2 = true; }
              
            if (flagB_edge_endpoint_maps_to_cube2)
              { flag_edge_maps_to_02 = true; }

            if (flagB_quad_vertex_maps_to_cube1 &&
                flagB_quad_vertex_maps_to_cube2) {
              if (flagB_edge_endpoint_maps_to_cube2) {
                flag_quad_maps_to_021X_or_02X1 = true;
              }
              else {
                flag_quad_maps_to_012X = true;
              }
            }

            // *** DEBUGXXX ***
            MSDEBUG();
            if (flag_debug) {
              isovert.grid.PrintIndexAndCoord
                (cerr, "          Flag maps to ", cube1_index, " : ");
              cerr << int(flagB_quad_vertex_maps_to_cube1) << endl;
              isovert.grid.PrintIndexAndCoord
                (cerr, "          Flag maps to ", cube2_index, " : ");
              cerr << int(flagB_quad_vertex_maps_to_cube2) << endl;
              cerr << "          flagB_edge_endpoint_maps_to_cube2: " 
                   << int(flagB_edge_endpoint_maps_to_cube2) << endl;
              cerr << "          flag_quad_maps_to_012X: " 
                   << int(flag_quad_maps_to_012X) << endl;
              cerr << "          flag_quad_maps_to_021X_or_02X1: " 
                   << int(flag_quad_maps_to_021X_or_02X1) << endl;
            }

          }
        }
      }
    }
    else {
      // Handle boundary case.
    }

  }

}

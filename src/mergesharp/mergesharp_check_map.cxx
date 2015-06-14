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

  bool are_separated
  (const SHARPISO_GRID & grid,
   const GRID_CUBE_DATA & gcubeA, const COORD_TYPE coordA[DIM3],
   const GRID_CUBE_DATA & gcubeB, const COORD_TYPE coordB[DIM3],
   const int orth_dir);

  void subtract_project_normalize
  (const COORD_TYPE pA[DIM3], const COORD_TYPE pB[DIM3], 
   const int orth_dir, const COORD_TYPE max_small_magnitude,
   COORD_TYPE w[DIM3], COORD_TYPE & magnitude, bool & flag_zero);
}


// **************************************************
//  Check distortion
// **************************************************


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

  flag = check_distortion_strict
    (scalar_grid, isovalue, isovert, gcube_map, from_cube0, to_cube0, param);

  // restore gcube_map[from_gcube1_index]
  gcube_map[from_gcube1_index] = store_map[1];

  if (!flag) { return(false); }

  // temporarily set gcube_map[from_gcube0_index] to to_cube0
  gcube_map[from_gcube0_index] = to_gcube0_index;

  flag = check_distortion_strict
    (scalar_grid, isovalue, isovert, gcube_map, from_cube1, to_cube1, param);

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

    flag = check_distortion_strict
      (scalar_grid, isovalue, isovert, gcube_map, cube_index[i0], 
       to_cube, param);

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

  if (gcubeB_flag == COVERED_A_GCUBE || gcubeB_flag == COVERED_CORNER_GCUBE) {
    return(true);
  }
  if (gcubeC_flag == COVERED_A_GCUBE || gcubeC_flag == COVERED_CORNER_GCUBE) {
    return(true);
  }

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
  const COORD_TYPE min_triangle_angle = param.min_triangle_angle;
  const COORD_TYPE cos_min_triangle_angle =
    cos(min_triangle_angle*M_PI/180.0);
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
        cerr << "--- Failed angle test." 
             << "  angle_ABC: " << acos(cos_angle_ABC) * 180.0/M_PI
             << "  angle_ACB: " << acos(cos_angle_ACB) * 180.0/M_PI
             << endl;
      }

      return(false);
    }
  }

  // Dihedral angle test
  const COORD_TYPE min_dihedral_angle = param.min_dihedral_angle;
  const COORD_TYPE cos_min_dihedral_angle = 
    cos(min_dihedral_angle*M_PI/180.0);

  COORD_TYPE cos_dihedral_angle;
  compute_cos_dihedral_angle
    (grid, unscaled_to_coord, unscaled_coordA, unscaled_coordB, unscaled_coordC,
     min_dist, cos_dihedral_angle, flag_small_magnitude);

  if (!flag_small_magnitude) {
    if (cos_dihedral_angle >= cos_min_dihedral_angle)
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

      cerr << "xxx Matches cube orientations. Fails dihedral test." << endl;
      grid.PrintIndexAndCoord
        (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
      grid.PrintIndexAndCoord
        (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
      if (!flag_small_magnitude) {
        if (cos_dihedral_angle < -1) { cos_dihedral_angle = -1; }
        if (cos_dihedral_angle > 1) { cos_dihedral_angle = 1; }
        cerr << "    dihedral angle: " << acos(cos_dihedral_angle)*180.0/M_PI
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
          cerr << "zzz Passes cube separation test.  Fails dihedral & cube orient tests." << endl;
          grid.PrintIndexAndCoord
            (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
          grid.PrintIndexAndCoord
            (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
          cerr << "    dihedral angle: " << acos(cos_dihedral_angle*M_PI/180.0)
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
  if (gcubeA_flag == COVERED_A_GCUBE || gcubeA_flag == COVERED_CORNER_GCUBE) {
    return(true);
  }
  if (gcubeC_flag == COVERED_A_GCUBE || gcubeC_flag == COVERED_CORNER_GCUBE) {
    return(true);
  }

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
  const COORD_TYPE min_triangle_angle = param.min_triangle_angle;
  const COORD_TYPE cos_min_triangle_angle =
    cos(min_triangle_angle*M_PI/180.0);
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

      if (flag_debug) {
        MSDEBUG();
        cerr << "--- Failed angle test." 
             << "  angle_BAC: " << acos(cos_angle_BAC) * 180.0/M_PI
             << "  angle_BCA: " << acos(cos_angle_BCA) * 180.0/M_PI
             << endl;

        return(false);
      }
    }
  }


  // Dihedral angle test
  const COORD_TYPE min_dihedral_angle = param.min_dihedral_angle;
  const COORD_TYPE cos_min_dihedral_angle = 
    cos(min_dihedral_angle*M_PI/180.0);

  COORD_TYPE cos_dihedral_angle;
  compute_cos_dihedral_angle
    (grid, unscaled_to_coord, unscaled_coordB, unscaled_coordA, unscaled_coordC,
     min_dist, cos_dihedral_angle, flag_small_magnitude);

  if (!flag_small_magnitude) {
    if (cos_dihedral_angle >= cos_min_dihedral_angle)
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

      cerr << "xxx Matches cube orientations. Fails dihedral test." << endl;
      grid.PrintIndexAndCoord
        (cerr, "    From: ", icubeA, " to ", to_cube, "\n");
      grid.PrintIndexAndCoord
        (cerr, "    cubeB: ", icubeB, " cubeC: ", icubeC, "\n");
      if (!flag_small_magnitude) {
        if (cos_dihedral_angle < -1) { cos_dihedral_angle = -1; }
        if (cos_dihedral_angle > 1) { cos_dihedral_angle = 1; }
        cerr << "    dihedral angle: " << acos(cos_dihedral_angle)*180.0/M_PI
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


// Compute cosine of the dihedral angle between (A0,B,C) and (A1,B,C)
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

}

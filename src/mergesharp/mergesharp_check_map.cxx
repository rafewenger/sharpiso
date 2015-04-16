/// \file mergesharp_check_map.cxx
/// Check if mapping one isosurface vertex to another distorts or reverses
///   any isosurface triangles.


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

  bool is_between
  (const COORD_TYPE x, const COORD_TYPE a0, const COORD_TYPE a1);

  bool are_separated
  (const SHARPISO_GRID & grid,
   const GRID_CUBE_DATA & gcubeA, const COORD_TYPE coordA[DIM3],
   const GRID_CUBE_DATA & gcubeB, const COORD_TYPE coordB[DIM3],
   const int orth_dir);
}


// **************************************************
//  Check distortion
// **************************************************


/// Return true if mapping of from_cube to to_cube does not reverse any
///   triangles incident on vertex in from_cube 
///   or create any degenerate triangles.
bool MERGESHARP::check_distortion
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube)
{
  IJK::PROCEDURE_ERROR error("check_distortion");
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

          if (!check_quad_distortion
              (scalar_grid, isovalue, isovert, gcube_map, from_cube, to_cube,
              edge_dir, iend0, j1, j2)) 
            { return(false); }
        }
      }
    }
  }
  else {
    // HANDLE BOUNDARY CASE
  }

  // *** DEBUG ***
  /*
  MSDEBUG();
  scalar_grid.PrintIndexAndCoord
    (cerr, "  Passed distortion test!.  From cube: ", from_cube,
     " to cube ", to_cube, "\n");
  */

  return(true);
}



/// Return true if mapping of from_cube to to_cube does not reverse/distort
///   any triangles on quad dual to (iend0,iend1).
bool MERGESHARP::check_quad_distortion
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube,
 const VERTEX_INDEX to_cube,
 const int edge_dir,
 const VERTEX_INDEX iend0,
 const int j1,
 const int j2)
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

  if (!check_tri_distortion
      (scalar_grid, isovert, gcube_map, from_cube, to_cube,
       icubeA, icubeB, dir1, dir2))
    { return(false); }

  if (!check_tri_distortion
      (scalar_grid, isovert, gcube_map, from_cube, to_cube,
       icubeC, icubeB, dir2, dir1))
    { return(false); }

  // *** DEBUG ***
  /*
  MSDEBUG();
  scalar_grid.PrintIndexAndCoord
    (cerr, "  NO quad distortion.  From cube: ", from_cube,
     " to cube ", to_cube, ". ");
  scalar_grid.PrintIndexAndCoord
    (cerr, "  Edge: [", iend0, ",", iend1, "]\n");
  */

  return(true);
}

/// Return true if mapping of from_cube to to_cube does not reverse/distort
///   triangle with vertices in (from_cube,ivA,ivB)
/// @param Cube icubeA shares a facet with from_cube.
/// @param Cube icubeB shares an edge with from_cube.
/// @param dirFA Direction (0,1,2) of from_cube to icubeA.
/// @param dirAB Direction (0,1,2) of icubeA to icubeB.
bool MERGESHARP::check_tri_distortion
(const SHARPISO_GRID & grid, const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map,
 const VERTEX_INDEX from_cube, const VERTEX_INDEX to_cube,
 const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB,
 const int dirFA, const int dirAB)
{
  const COORD_TYPE * spacing = grid.SpacingPtrConst();
  IJK::PROCEDURE_ERROR error("check_tri_distortion");
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA, error);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB, error);
  const INDEX_DIFF_TYPE to_gcube_index = isovert.GCubeIndex(to_cube, error);
  const INDEX_DIFF_TYPE from_gcube_index = isovert.GCubeIndex(from_cube, error);

  if ((gcubeA_index == ISOVERT::NO_INDEX) || 
      (gcubeB_index == ISOVERT::NO_INDEX) || 
      (to_gcube_index == ISOVERT::NO_INDEX) || 
      (from_gcube_index == ISOVERT::NO_INDEX))
    { throw error; }

  if ((gcube_map[gcubeA_index] == to_gcube_index) || 
      (gcube_map[gcubeB_index] == to_gcube_index) ||
      (gcube_map[gcubeA_index] == gcube_map[gcubeB_index])) 
    { return(true); }

  bool return_flag =
    check_tri_distortion
    (grid, isovert, from_cube, icubeA, icubeB,
     isovert.IsoVertCoord(to_gcube_index),
     isovert.IsoVertCoord(gcube_map[gcubeA_index]),
     isovert.IsoVertCoord(gcube_map[gcubeB_index]),
     dirFA, dirAB);

  return(return_flag);
}


// Return true if mapping to given coordinates does not
//   reverse/distort triangle with vertices in (icubeA, icubeB, icubeC)
// @param icubeA Cube icubeA shares an edge with cube icubeC.
// @param icubeB Cube icubeB shares facets with cubes icubeA and icubeC.
// @param icubeC Cube icubeC shares an edge with cube icubeA.
// @param dirAB Direction (0,1,2) of icubeA to icubeB.
// @param dirBC Direction (0,1,2) of icubeB to icubeC.
bool MERGESHARP::check_tri_distortion
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB, 
 const VERTEX_INDEX icubeC, 
 const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3],
 const COORD_TYPE coordC[DIM3],
 const int dirAB, const int dirBC)
{
  const COORD_TYPE * spacing = grid.SpacingPtrConst();
  IJK::PROCEDURE_ERROR error("check_tri_distortion");
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA, error);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB, error);
  const INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(icubeC, error);


  if ((gcubeA_index == ISOVERT::NO_INDEX) || 
      (gcubeB_index == ISOVERT::NO_INDEX) || 
      (gcubeC_index == ISOVERT::NO_INDEX))
    { throw error; }

  COORD_TYPE unscaled_coordA[DIM3], unscaled_coordB[DIM3], 
    unscaled_coordC[DIM3];

  // Reverse scaling on coordA[], coordB[], coordC[]
  IJK::reverse_coord_scaling(DIM3, spacing, coordA, unscaled_coordA);
  IJK::reverse_coord_scaling(DIM3, spacing, coordB, unscaled_coordB);
  IJK::reverse_coord_scaling(DIM3, spacing, coordC, unscaled_coordC);

  // *** DEBUG ***
  /*
  MSDEBUG();
  grid.PrintIndexAndCoord
    (cerr, "  Triangle cubes: ", icubeA, " ", icubeB, " ", icubeC, "\n");
  cerr << "  Triangle vertices: ";
  IJK::print_coord3D(cerr, coordA);
  IJK::print_coord3D(cerr, coordB);
  IJK::print_coord3D(cerr, coordC);
  cerr << endl;
  */

  if (are_separated(grid, isovert.gcube_list[gcubeA_index], unscaled_coordA,
                    isovert.gcube_list[gcubeB_index], unscaled_coordB, dirAB)) {

    if (are_separated(grid, isovert.gcube_list[gcubeB_index], unscaled_coordB, 
                      isovert.gcube_list[gcubeC_index], unscaled_coordC, 
                      dirBC)) {

      // *** DEBUG ***
      /*
      cerr << "    Passed first two separation tests." << endl;
      */

      if (are_separated(grid, isovert.gcube_list[gcubeA_index], unscaled_coordA,
                        isovert.gcube_list[gcubeC_index], unscaled_coordC, 
                        dirAB)) {

        // *** DEBUG ***
        /*
        cerr << "    Passed separation test in direction " << dirAB << endl;
        cerr << "    coord[" << dirBC << "] : "
             << coordB[dirBC] << " " << coordA[dirBC]
             << " " << coordC[dirBC] << endl;
        */

        if (is_between(unscaled_coordA[dirBC], unscaled_coordB[dirBC], 
                       unscaled_coordC[dirBC]))
          { return(true); }      
      }

      if (are_separated(grid, isovert.gcube_list[gcubeA_index], unscaled_coordA,
                        isovert.gcube_list[gcubeC_index], unscaled_coordC, 
                        dirBC)) {

        // *** DEBUG ***
        /*
        cerr << "    Passed separation test in direction " << dirBC << endl;
        cerr << "    coord[" << dirAB << "] : "
             << coordA[dirAB] << " " << coordC[dirAB]
             << " " << coordB[dirAB] << endl;
        */

        if (is_between(coordC[dirAB], coordA[dirAB], coordB[dirAB]))
          { return(true); }
      }
    }
  }

  return(false);
}


// **************************************************
//  Local routines
// **************************************************

namespace {

  bool is_between
  (const COORD_TYPE x, const COORD_TYPE a0, const COORD_TYPE a1)
  {
    if ((a0 < x) && (x < a1)) { return(true); }
    if ((a0 > x) && (x > a1)) { return(true); }
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

}

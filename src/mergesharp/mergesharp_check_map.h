/// \file mergesharp_check_map.h
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


#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"


namespace MERGESHARP {

  // **************************************************
  //  Check distortion
  // **************************************************

  /// Return true if mapping of from_cube to to_cube does not reverse any
  ///   triangles incident on vertex in from_cube 
  ///   or create any degenerate triangles.
  /// @param flag_strict If true, use strict tests.
  bool check_distortion
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const bool flag_strict,
   const MERGE_PARAM & param);

  /// Return true if mapping of from_cube to to_cube does not reverse any
  ///   triangles incident on vertex in from_cube 
  ///   or create any degenerate triangles.
  /// Use check_quad_distortion_strict.
  bool check_distortion_strict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const MERGE_PARAM & param);

  /// Return true if mapping of from_cube to to_cube does not reverse any
  ///   triangles incident on vertex in from_cube 
  ///   or create any degenerate triangles.
  /// Use check_quad_distortion_loose.
  bool check_distortion_loose
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const MERGE_PARAM & param);

  /// Return true if simultaneous mapping of from_cube0 to to_cube0
  ///   and from_cube1 to to_cube1 does not reverse any triangles 
  ///   incident on vertices in from_cube0 or from_cube1 or create 
  ///    any degenerate triangles.
  bool check_distortionII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube0, const VERTEX_INDEX to_cube0,
   const VERTEX_INDEX from_cube1, const VERTEX_INDEX to_cube1,
   const bool flag_strict,
   const MERGE_PARAM & param);

  // Return true if simultaneous mapping of three cubes to to_cube
  //   does not reverse any triangles incident on vertices on cubes
  //   or create any degenerate triangles.
  bool check_distortionIII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX cube_index[3], const VERTEX_INDEX to_cube,
   const bool flag_strict,
   const MERGE_PARAM & param);

  /// Return true if mapping of from_cube to to_cube does not reverse/distort
  ///   any triangles on quad dual to (iend0,iend1).
  bool check_quad_distortion_strict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const int edge_direction,
   const VERTEX_INDEX iend0,
   const int j1,
   const int j2,
   const MERGE_PARAM & param);

  /// Return true if mapping of from_cube to to_cube does not reverse/distort
  ///   triangles FAB or FBC on quad (FABC) dual to (iend0, iend1).
  bool check_quad_distortion_loose
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const int edge_direction,
   const VERTEX_INDEX iend0,
   const int j1,
   const int j2,
   const MERGE_PARAM & param);


  /// Return true if mapping of icubA to to_cube does not reverse/distort
  ///   triangles FAB or FBC on quad (FABC) dual to (iend0, iend1)
  ///   or does not reverse/distort triangle FAC.
  /// @param Cube icubeB shares a facet with icubeA.
  /// @param Cube icubeC shares an edge with icubeA.
  /// @param relposAB Relative positions of cubes icubeA and icubeB
  ///         If relposAB = 0, then icubeB precedes icubeA in direction dirAB.
  ///         If relposAB = 1, then icubeB follows icubeA in direction dirAB.
  /// @param relposAB Relative positions of cubes icubeA and icubeB.
  ///         If relposAB = 0, then icubeC precedes icubeB in direction dirAB.
  ///         If relposAB = 1, then icubeC follows icubeB in direction dirAB.
  bool check_tri_distortion_mapA
  (const SHARPISO_GRID & grid, const ISOVERT & isovert, 
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB, 
   const VERTEX_INDEX icubeC,
   const VERTEX_INDEX to_cube,
   const int dirAB, const int dirBC, const int relposAB, const int relposBC,
   const MERGE_PARAM & param);

  /// Return true if mapping of icubeB to to_cube does not reverse/distort
  ///   triangle with vertices in (icubeA,icubeB,icubeC).
  /// @param Cube icubeB shares a facet with icubeA.
  /// @param Cube icubeC shares an edge with icubeA.
  /// @param relposAB Relative positions of cubes icubeA and icubeB
  ///         If relposAB = 0, then icubeB precedes icubeA in direction dirAB.
  ///         If relposAB = 1, then icubeB follows icubeA in direction dirAB.
  /// @param relposAB Relative positions of cubes icubeA and icubeB.
  ///         If relposAB = 0, then icubeC precedes icubeB in direction dirAB.
  ///         If relposAB = 1, then icubeC follows icubeB in direction dirAB.
  bool check_tri_distortion_mapB
  (const SHARPISO_GRID & grid, const ISOVERT & isovert, 
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB, 
   const VERTEX_INDEX icubeC,
   const VERTEX_INDEX to_cube,
   const int dirAB, const int dirBC, const int relposAB, const int relposBC,
   const MERGE_PARAM & param);


  /// Compute cos triangle angles.
  void compute_cos_triangle_angles
  (const SHARPISO_GRID & grid, 
   const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3],
   const COORD_TYPE coordC[DIM3],
   const COORD_TYPE min_dist,
   COORD_TYPE & cos_angle_ABC, COORD_TYPE & cos_angle_ACB,
   bool & flag_small_magnitude);

  /// Compute projected triangle orientation.
  /// @param orth_dir Direction of orthogonal projection onto plane.
  /// @param orient[out] Orientation, 1,-1 or 0.  
  ///       0 if cos projected angle ABC or ACB is below cos_min_proj_angle
  ///       or if distance between projected points is less than 
  ///       or equal to min_dist.
  /// @param flag_small_BC[out] If true, distance between projected 
  ///       coordB and coordC is less than or equal to min_dist.
  /// @pre (orth_dir, dirAB, dirBC) are distinct.
  void compute_projected_tri_orientation
  (const SHARPISO_GRID & grid, 
   const COORD_TYPE coordA[DIM3], const COORD_TYPE coordB[DIM3],
   const COORD_TYPE coordC[DIM3],
   const COORD_TYPE cos_min_proj_angle, const COORD_TYPE min_dist,
   const int orth_dir, const int dirAB, const int dirBC,
   int & orient, bool & flag_small_BC);

  /// Compute cosine of the dihedral angle between (A0,B,C) and (A1,B,C).
  void compute_cos_dihedral_angle
  (const SHARPISO_GRID & grid, 
   const COORD_TYPE coordA0[DIM3], const COORD_TYPE coordA1[DIM3],
   const COORD_TYPE coordB[DIM3], const COORD_TYPE coordC[DIM3],
   const COORD_TYPE min_dist,
   COORD_TYPE & cos_dihedral, bool & flag_small_magnitude);


  // **************************************************
  //  Check for violation of manifold conditions
  // **************************************************

  /// Check for manifold violations caused by edges between sharp cubes.
  bool check_edge_manifold
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const MERGESHARP::ISOVERT & isovert, 
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);

  /// Check for manifold violations caused by edges between sharp cubes.
  ///   caused by mapping from_cube0 to to_cube0 and from_cube1 to to_cube1.
  bool check_edge_manifoldII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX from_cube0,
   const VERTEX_INDEX to_cube0,
   const VERTEX_INDEX from_cube1,
   const VERTEX_INDEX to_cube1,
   const MERGESHARP::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);

  /// Check for manifold violations caused by edges between sharp cubes.
  ///   caused by mapping three cubes to to_cube.
  bool check_edge_manifoldIII
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube_index[3],
   const VERTEX_INDEX to_cube,
   const MERGESHARP::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
   const bool flag_extended);


  // **************************************************
  //  Misc routines
  // **************************************************

  /// Get cubes in clockwise order around edge.
  /// @pre edge is in the interior of the grid.
  void get_cubes_around_edge
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iend0,
   const int edge_dir, VERTEX_INDEX cube_index[NUM_VERT_PER_QUAD]);

}

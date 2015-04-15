/// \file mergesharp_check_map.h
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


#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"


namespace MERGESHARP {

  // **************************************************
  //  Check distortion
  // **************************************************

  /// Return true if mapping of from_cube to to_cube does not reverse any
  ///   triangles incident on vertex in from_cube 
  ///   or create any degenerate triangles.
  bool check_distortion
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube);

  /// Return true if mapping of from_cube to to_cube does not reverse/distort
  ///   any triangles on quad dual to (iend0,iend1).
  bool check_quad_distortion
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube,
   const VERTEX_INDEX to_cube,
   const int edge_direction,
   const VERTEX_INDEX iend0,
   const int j1,
   const int j2);

  /// Return true if mapping of from_cube to to_cube does not reverse/distort
  ///   triangle with vertices in (from_cube,ivA,ivB)
  /// @param Cube icubeA shares a facet with from_cube.
  /// @param Cube icubeB shares an edge with from_cube.
  /// @param dirFA Direction (0,1,2) of from_cube to icubeA.
  /// @param dirAB Direction (0,1,2) of icubeA to icubeB.
  bool check_tri_distortion
  (const SHARPISO_GRID & grid, const ISOVERT & isovert, 
   const std::vector<VERTEX_INDEX> & gcube_map,
   const VERTEX_INDEX from_cube, const VERTEX_INDEX to_cube,
   const VERTEX_INDEX icubeA, const VERTEX_INDEX icubeB,
   const int dirFA, const int dirAB);
}

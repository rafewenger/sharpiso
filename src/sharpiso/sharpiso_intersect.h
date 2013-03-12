/// \file sharpiso_intersect.h
/// Use sharp formulat to compute intersection of isosurface and grid edge.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2012-2013 Rephael Wenger
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public License
 (LGPL) as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _SHARPISO_INTERSECT_
#define _SHARPISO_INTERSECT_

#include "sharpiso_types.h"
#include "sharpiso_grids.h"

namespace SHARPISO {

  /// Compute intersections of isosurface and all grid edges.
  /// @param[out] edgeI_coord[] Coordinates of isosurface-edge intersections
  /// @param[out] edgeI_normal_coord[] Coordinates of normal vectors
  ///             at each location in edgeI_coord[].
  void compute_all_edgeI
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_magnitude,
   std::vector<COORD_TYPE> & edgeI_coord,
   std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord);

  /// Compute intersections of isosurface and all grid edges 
  ///   using linear interpolation.
  /// Note: This is NOT the recommended method for computing intersections
  ///   of isosurface and grid edges when gradient data is available.
  ///   This routine is provided for testing/comparison purposes.
  void compute_all_edgeI_linear_interpolate
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_magnitude,
   std::vector<COORD_TYPE> & edgeI_coord,
   std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord);

  /// Compute intersection of isosurface and grid edge using sharp formula.
  void compute_isosurface_grid_edge_intersection
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
   COORD_TYPE p[DIM3]);

  /// Compute intersection of isosurface and grid edge and 
  ///    normal at the intersection point using sharp formula.
  void compute_isosurface_grid_edge_intersection
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
   const GRADIENT_COORD_TYPE max_small_magnitude,
   COORD_TYPE p[DIM3],
   GRADIENT_COORD_TYPE normal[DIM3]);

  /// Compute intersection of edge and plane determined by gradient g, scalar s.
  /// Intersection point is v0+t*dir, 0 <= t <= 1.
  /// If plane does not intersect edge in a single point, then t < 0 or t > 1.
  void compute_edge_intersection
  (const SCALAR_TYPE s, const GRADIENT_COORD_TYPE g,
   const SCALAR_TYPE isovalue, SCALAR_TYPE & t);
};

#endif


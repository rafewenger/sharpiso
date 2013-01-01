/// \file isodual3D_extract.h
/// Subroutines for extracting dual isosurface mesh

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

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

#ifndef _ISODUAL3D_EXTRACT_
#define _ISODUAL3D_EXTRACT_

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"
#include "isodual3D_isovert.h"


namespace ISODUAL3D {

  /// Extract dual isosurface polytopes.
  /// Returns list of isosurface polytope vertices.
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polytope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  void extract_dual_isopoly
    (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
     ISODUAL_INFO & isodual_info);

  /// Extract dual isosurface polytopes.
  /// Returns list of isosurface polytope vertices.
  /// Return locations of isosurface vertices on each facet.
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_cube[] = cubes containing isosurface polytope vertices.
  ///   iso_cube[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param facet_vertex = Location of iso vertex on facet.
  void extract_dual_isopoly
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_cube,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   ISODUAL_INFO & isodual_info);

  /// Extract dual isosurface polytopes from list of edges.
  /// Returns list of isosurface polytope vertices.
  /// Return locations of isosurface vertices on each facet.
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param edge_list = List of edges. Polytopes are dual to edges.
  /// @pre Each edge is an internal grid edge.
  /// @param iso_cube[] = cubes containing isosurface polytope vertices.
  ///   iso_cube[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param facet_vertex = Location of iso vertex on facet.
  void extract_dual_isopoly_from_list
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<EDGE_INDEX> & edge_list,
   std::vector<ISO_VERTEX_INDEX> & iso_cube,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   ISODUAL_INFO & isodual_info);


  // **************************************************
  // MAP TO ISOPOLY VERTICES
  // **************************************************

  /// Map cube indices in iso_poly_vert to isosurface vertices
  ///   in isovert.gcube_list.
  void map_isopoly_vert
  (const ISOVERT & isovert, std::vector<ISO_VERTEX_INDEX> & iso_poly_vert);

}

#endif

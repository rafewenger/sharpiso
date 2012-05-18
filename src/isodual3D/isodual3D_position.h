/// \file isodual3D_position.h
/// Position dual isosurface vertices

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


#ifndef _ISODUAL3D_POSITION_
#define _ISODUAL3D_POSITION_

#include <vector>

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"

#include "ijkdualtable.h"

namespace ISODUAL3D {


  // **************************************************
  // Position using only scalar data
  // **************************************************

  /// Position dual isosurface vertices in cube centers
  void position_dual_isovertices_cube_center
    (const ISODUAL_GRID & grid,
    const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in cube centers
  void position_dual_isovertices_cube_center
    (const ISODUAL_GRID & grid,
    const std::vector<ISO_VERTEX_INDEX> & vlist,
    std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections
  void position_dual_isovertices_centroid
    (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
    const SCALAR_TYPE isovalue,
    const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections
  void position_dual_isovertices_centroid
    (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
    const SCALAR_TYPE isovalue,
    const std::vector<ISO_VERTEX_INDEX> & vlist,
    std::vector<COORD_TYPE> & coord);

  // **************************************************
  // Position using gradients and svd
  // **************************************************

  /// Position dual isosurface vertices using gradients
  void position_dual_isovertices_using_gradients
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   COORD_TYPE * coord,
   SHARPISO_INFO & sharp_info);

  /// Position dual isosurface vertices using gradients
  void position_dual_isovertices_using_gradients
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   std::vector<COORD_TYPE> & coord,
   SHARPISO_INFO & sharp_info);

  /// Position dual isosurface vertices using gradients
  void position_dual_isovertices_using_gradients
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_facet,
   COORD_TYPE * sharp_coord,
   SHARPISO_INFO & sharp_info);

  /// Position dual isosurface vertices using gradients
  void position_dual_isovertices_using_gradients
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_facet,
   std::vector<COORD_TYPE> & sharp_coord,
   SHARPISO_INFO & sharp_info);

  // ********************************************************
  // Position using gradients interpolated on grid edges.
  // ********************************************************

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Approximate gradients using linear interpolation on the edges.
  void position_dual_isovertices_edgeI_interpolate_gradients
  (const ISODUAL_SCALAR_GRID_BASE & grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const ISODUAL_PARAM & isodual_param, COORD_TYPE * coord);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Approximate gradients using linear interpolation on the edges.
  void position_dual_isovertices_edgeI_interpolate_gradients
  (const ISODUAL_SCALAR_GRID_BASE & grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const ISODUAL_PARAM & isodual_param,
   std::vector<COORD_TYPE> & coord);

  // *******************************************************************
  // Position using gradients determining edge-isosurface intersections.
  // *******************************************************************

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Select endpoint gradient which determines edge-isosurface intersection.
  void position_dual_isovertices_edgeI_sharp_gradients
  (const ISODUAL_SCALAR_GRID_BASE & grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const ISODUAL_PARAM & isodual_param, COORD_TYPE * coord);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Select endpoint gradient which determines edge-isosurface intersection.
  void position_dual_isovertices_edgeI_sharp_gradients
  (const ISODUAL_SCALAR_GRID_BASE & grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const ISODUAL_PARAM & isodual_param,
   std::vector<COORD_TYPE> & coord);

  // **************************************************
  // Reposition routine
  // **************************************************

  /// Reposition to separate isosurface vertices
  void reposition_dual_isovertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   COORD_TYPE * isovert_coord,
   SHARPISO_INFO & sharp_info);

  // **************************************************
  // Compute routines
  // **************************************************

  /// Compute centroid of intersections of isosurface and grid edges
  void compute_isosurface_grid_edge_centroid
    (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
    const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
    COORD_TYPE * coord);
};

#endif

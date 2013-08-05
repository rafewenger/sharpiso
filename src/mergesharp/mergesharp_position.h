/// \file mergesharp_position.h
/// Position dual isosurface vertices

/*
Copyright (C) 2011-2013 Rephael Wenger

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


#ifndef _MERGESHARP_POSITION_
#define _MERGESHARP_POSITION_

#include <vector>

#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"


#include "ijkdualtable.h"

namespace MERGESHARP {


  // **************************************************
  // Position at cube centers
  // **************************************************

  /// Position dual isosurface vertices in cube centers
  void position_dual_isovertices_cube_center
    (const SHARPISO_GRID & grid,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in cube centers
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_cube_center
    (const SHARPISO_GRID & grid,
    const std::vector<ISO_VERTEX_INDEX> & vlist,
    std::vector<COORD_TYPE> & coord);

  // **************************************************
  // Position at centroid of isosurface-edge intersections
  // **************************************************

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections
  void position_dual_isovertices_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections.
  /// Version allowing multiple vertices in a grid cube.
  /// @param isodual_table Determines edges between isosurface vertices.
  void position_dual_isovertices_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections.
  /// Version allowing multiple vertices in a grid cube.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
   std::vector<COORD_TYPE> & coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections.
  /// Allow multiple vertices in a grid cube.
  void position_dual_isovertices_centroid_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   COORD_TYPE * coord);

  /// Position dual isosurface vertices in centroid
  ///   of isosurface-edge intersections.
  /// Allowing multiple vertices in a grid cube.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_centroid_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<COORD_TYPE> & coord);

  // **************************************************
  // Position using isovert information
  // **************************************************

  /// Position merged dual isosurface vertices using isovert information.
  /// Allows multiple vertices in a grid cube.
  /// Recompute positions in cubes which are not selected and not smooth.
  void position_merged_dual_isovertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   COORD_TYPE * isov_coord);

  /// Position merged dual isosurface vertices using isovert information.
  /// Allows multiple vertices in a grid cube.
  /// Version using std::vector for array coord[].
  void position_merged_dual_isovertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<COORD_TYPE> & isov_coord);

  /// Position dual isosurface vertices using isovert information.
  /// Allows multiple vertices in a grid cube.
  void position_dual_isovertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   COORD_TYPE * isov_coord);

  /// Position dual isosurface vertices using isovert information.
  /// Allows multiple vertices in a grid cube.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<COORD_TYPE> & isov_coord);

  // ********************************************************
  // Position vertices using SVD on grid edge-isosurface intersections.
  // ********************************************************

  /// Position vertices using SVD on grid edge-isosurface intersections.
  void position_dual_isovertices_edgeI
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const VERTEX_POSITION_METHOD position_method,
   COORD_TYPE * sharp_coord,
   SHARPISO_INFO & sharp_info);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_edgeI
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<ISO_VERTEX_INDEX> & vlist,
   const VERTEX_POSITION_METHOD position_method,
   std::vector<COORD_TYPE> & coord,
   SHARPISO_INFO & sharp_info);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Allows multiple vertices in a grid cube.
  /// Use isodual_table to determine connections.
  void position_dual_isovertices_edgeI_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const VERTEX_POSITION_METHOD position_method,
   COORD_TYPE * coord,
   SHARPISO_INFO & sharp_info);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Allows multiple vertices in a grid cube.
  /// Use isodual_table to determine connections.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_edgeI_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   const VERTEX_POSITION_METHOD position_method,
   std::vector<COORD_TYPE> & coord,
   SHARPISO_INFO & sharp_info);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Version resolving ambiguous cubes.
  /// @param iso_vlist_cube_ambig[] Determines how ambiguities are resolved.
  void position_dual_isovertices_edgeI
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
   const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
   const VERTEX_POSITION_METHOD position_method,
   COORD_TYPE * coord,
   SHARPISO_INFO & sharp_info);

  /// Position vertices using SVD on grid edge-isosurface intersections.
  /// Version allowing multiple vertices in a grid cube.
  /// Version resolving ambiguous cubes.
  /// Version using std::vector for array coord[].
  void position_dual_isovertices_edgeI
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const MERGESHARP_PARAM & mergesharp_param,
   const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
   const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
   const VERTEX_POSITION_METHOD position_method,
   std::vector<COORD_TYPE> & coord,
   SHARPISO_INFO & sharp_info);


  // **************************************************
  // Split dual isosurface vertices.
  // **************************************************

  /// Split dual isosurface vertices.
  /// @param isodual_table Dual isosurface lookup table.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist_cube[] iso_vlist_cube[i] is the cube containing
  ///    isosurface vertex i.
  /// @param iso_vlist_patch[] iso_vlist_patch[j] is the index of the
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  void split_dual_isovert_ambig
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & cube_list,
   const std::vector<AMBIGUITY_TYPE> & cube_ambig,
   const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,
   const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
   std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
   std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   VERTEX_INDEX & num_split);

  /// Split dual isosurface vertices.
  /// @param isodual_table Dual isosurface lookup table.
  /// Version using vector<DUAL_ISOVERT>.
  /// @param[out] iso_vlist[] List of isosurface vertices.
  void split_dual_isovert_ambig
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & cube_list,
   const std::vector<AMBIGUITY_TYPE> & cube_ambig,
   const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
   const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   VERTEX_INDEX & num_split);

  /// Split dual isosurface vertices.
  /// Split all isosurface vertices as determined by the lookup table
  ///   of the non-manifold edges (if flag_split_non_manifold is true.)
  /// Version which returns cube_isovert_data.
  /// @param isodual_table Dual isosurface lookup table.
  void full_split_dual_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
   const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   const MERGESHARP_PARAM & mergesharp_param,
   std::vector<IJKDUALTABLE::TABLE_INDEX> & table_index,
   std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   SHARPISO_INFO & sharp_info);

  /// Split dual isosurface vertices.
  /// Split all isosurface vertices as determined by the lookup table
  ///   of the non-manifold edges (if flag_split_non_manifold is true.)
  /// @param isodual_table Dual isosurface lookup table.
  void full_split_dual_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
   const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   const MERGESHARP_PARAM & mergesharp_param,
   std::vector<DUAL_ISOVERT> & iso_vlist,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   SHARPISO_INFO & sharp_info);

  // **************************************************
  // Compute routines
  // **************************************************

  /// Compute centroid of intersections of isosurface and cube edges.
  void compute_isosurface_grid_edge_centroid
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
    const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
    COORD_TYPE * coord);

  /// Compute centroid of intersections of isosurface and cube edges.
  /// Use only grid cube edges associated with isosurface patch ipatch.
  void compute_isosurface_grid_edge_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX icube,
   const FACET_VERTEX_INDEX ipatch, const IJKDUALTABLE::TABLE_INDEX it,
   const MERGESHARP_CUBE_FACE_INFO & cube,
   COORD_TYPE * coord);

  // **************************************************
  // Copy routine
  // **************************************************

  /// Copy isovert position from data structure isovert
  void copy_isovert_positions
    (const std::vector<GRID_CUBE> & gcube_list,
     COORD_ARRAY & vertex_coord);

};

#endif

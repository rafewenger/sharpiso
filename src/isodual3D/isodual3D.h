/// \file isodual3D.h
/// Generate isosurface in arbitrary dimensions

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2013 Rephael Wenger

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

/*!
  \mainpage ISODUAL3D: 3D DUAL CONTOURING

  ISODUAL3D is a program for generating isosurfaces
  using the dual contouring algorithm.  It computes
  isosurfaces from three dimensional volumetric grid data.
*/

#ifndef _ISODUAL3D_
#define _ISODUAL3D_

#include <string>

#include "ijk.txx"
#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"


/// isodual3D classes and routines.
namespace ISODUAL3D {

  // **************************************************
  // DUAL CONTOURING
  // **************************************************

  /// Dual Contouring Algorithm.
  void dual_contouring
    (const ISODUAL_DATA & isodual_data, const SCALAR_TYPE isovalue,
     DUAL_ISOSURFACE & dual_isosurface, ISODUAL_INFO & isodual_info);

  // **************************************************
  // DUAL CONTOURING USING SCALAR DATA
  // **************************************************

  /// Dual Contouring Algorithm using only scalar data.
  /// Represents each grid edge by a single integer.
  /// @param merge_data = Data structure for merging edges.
  /// Requires memory of size(MERGE_INDEX) for each grid edge.
  void dual_contouring
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data, ISODUAL_INFO & isodual_info);

  /// Dual contouring algorithm.
  /// Position isosurface vertices at cube centers.
  /// @param merge_data = Data structure for merging edges.
  /// Requires memory of size(MERGE_INDEX) for each grid edge.
  void dual_contouring_cube_center
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data,
   ISODUAL_INFO & isodual_info);

  /// Dual contouring algorithm.
  /// Position isosurface vertices at centroid of isosurface-edge intersections.
  void dual_contouring_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data,
   ISODUAL_INFO & isodual_info);

  /// Dual contouring algorithm.
  /// Position isosurface vertices at centroid of isosurface-edge intersections.
  /// Allow multiple isosurface vertices in a grid cube.
  void dual_contouring_centroid_multiv
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const bool flag_separate_neg,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data,
   ISODUAL_INFO & isodual_info);


  // **************************************************
  // DUAL CONTOURING USING SCALAR & GRADIENT DATA
  // **************************************************

  /// Extract dual contouring isosurface.
  /// Returns list of isosurface triangle and quad vertices
  ///   and list of isosurface vertex coordinates.
  /// Use gradients to place isosurface vertices on sharp features. 
  void dual_contouring_sharp_from_grad
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   ISODUAL_INFO & isodual_info);

  /// Extract dual contouring isosurface.
  /// Returns list of isosurface quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @pre isovert contains isovert locations.
  void dual_contouring_extract_isopoly
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   ISODUAL_INFO & isodual_info,
   ISOVERT_INFO & isovert_info);

  /// Extract dual contouring isosurface.
  /// Allow multiple isosurface vertices in a grid cube.
  /// Returns list of isosurface quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @param cube_ambig[]  cube_ambig[i] indicates how ambiguities
  ///        should be handled for i'th cube in isovert.gcube_list[i].
  /// @pre isovert contains isovert locations.
  void dual_contouring_extract_isopoly_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   ISODUAL_INFO & isodual_info,
   ISOVERT_INFO & isovert_info);

  /// Extract dual contouring isosurface.
  /// Allow multiple isosurface vertices in a grid cube.
  /// Resolve ambiguous facets.
  /// Returns list of isosurface quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @param cube_ambig[]  cube_ambig[i] indicates how ambiguities
  ///        should be handled for i'th cube in isovert.gcube_list[i].
  /// @pre isovert contains isovert locations.
  void dual_contouring_extract_isopoly_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   const std::vector<AMBIGUITY_TYPE> & cube_ambig,
   ISODUAL_INFO & isodual_info,
   ISOVERT_INFO & isovert_info);

  // **************************************************
  // MERGE SHARP
  // **************************************************

  /// Extract dual contouring isosurface by merging grid cubes
  ///   around sharp vertices.
  /// Dual Contouring algorithm for sharp isosurface features.
  /// Merges grid cubes around sharp vertices.
  /// Return list of isosurface triangle and quad vertices
  ///   and list of isosurface vertex coordinates.
  /// Use gradients to place isosurface vertices on sharp features. 
  void dual_contouring_merge_sharp_from_grad
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const GRADIENT_GRID_BASE & gradient_grid,
     const SCALAR_TYPE isovalue,
     const ISODUAL_PARAM & isodual_param,
     DUAL_ISOSURFACE & dual_isosurface,
     ISOVERT & isovert,
     ISODUAL_INFO & isodual_info);

  /// Extract dual contouring isosurface by merging grid cubes
  ///   around sharp vertices.
  /// Returns list of isosurface triangle and quad vertices
  ///   and list of isosurface vertex coordinates.
  /// Use input edge-isosurface intersections and normals (hermite data)
  ///   to position isosurface vertices on sharp features.
  void dual_contouring_merge_sharp_from_hermite
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<COORD_TYPE> & edgeI_coord,
   const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   ISODUAL_INFO & isodual_info);

  /// Extract dual contouring isosurface by merging grid cubes
  ///   around sharp vertices.
  /// Returns list of isosurface triangle and quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @pre isovert contains isovert locations.
  void dual_contouring_merge_sharp
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_PARAM & isodual_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   ISODUAL_INFO & isodual_info,
   ISOVERT_INFO & isovert_info);
}

#endif

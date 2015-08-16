/// \file shrec.h
/// generate isosurface using dual contouring algorithm
/// Version 0.1.0

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

/*!
  \mainpage SHREC: SHarp REConstruction

  SHREC is a program for generating isosurfaces with sharp 0-dimensional
  and 1-dimensional features.  It computes isosurfaces from three 
  dimensional volumetric grid data.by merging grid cubes around 
  the features and applying dual contouring to build the isosurface.
*/

#ifndef _SHREC_
#define _SHREC_

#include <string>

#include "ijk.txx"
#include "shrec_types.h"
#include "shrec_datastruct.h"


/// shrec classes and routines.
namespace SHREC {

  // ********************************************************
  // VERSION
  // ********************************************************

  inline const char * ProgramVersion()
    { return("0.1.0"); }

  // **************************************************
  // DUAL CONTOURING
  // **************************************************

  /// Dual Contouring Algorithm.
  void dual_contouring
    (const SHREC_DATA & shrec_data, const SCALAR_TYPE isovalue,
     DUAL_ISOSURFACE & dual_isosurface, SHREC_INFO & shrec_info);

  /// Dual Contouring Algorithm.
  void dual_contouring
    (const SHREC_DATA & shrec_data, const SCALAR_TYPE isovalue,
     DUAL_ISOSURFACE & dual_isosurface, ISOVERT & isovert,
     SHREC_INFO & shrec_info);

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
   const SHREC_PARAM & shrec_param,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data, SHREC_INFO & shrec_info);

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
   SHREC_INFO & shrec_info);

  /// Dual contouring algorithm.
  /// Position isosurface vertices at centroid of isosurface-edge intersections.
  void dual_contouring_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   std::vector<VERTEX_INDEX> & isoquad_vert,
   std::vector<COORD_TYPE> & vertex_coord,
   MERGE_DATA & merge_data,
   SHREC_INFO & shrec_info);

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
   SHREC_INFO & shrec_info);


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
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   SHREC_INFO & shrec_info);

  /// Extract dual contouring isosurface.
  /// Returns list of isosurface quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @pre isovert contains isovert locations.
  void dual_contouring_extract_isopoly
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   SHREC_INFO & shrec_info,
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
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   SHREC_INFO & shrec_info,
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
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   const std::vector<AMBIGUITY_TYPE> & cube_ambig,
   SHREC_INFO & shrec_info,
   ISOVERT_INFO & isovert_info);

  // **************************************************
  // SHREC
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
     const SHREC_PARAM & shrec_param,
     DUAL_ISOSURFACE & dual_isosurface,
     ISOVERT & isovert,
     SHREC_INFO & shrec_info);

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
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   SHREC_INFO & shrec_info);

  /// Extract dual contouring isosurface by merging grid cubes
  ///   around sharp vertices.
  /// Returns list of isosurface triangle and quad vertices
  ///   and list of isosurface vertex coordinates.
  /// @pre isovert contains isovert locations.
  void dual_contouring_merge_sharp
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const SHREC_PARAM & shrec_param,
   DUAL_ISOSURFACE & dual_isosurface,
   ISOVERT & isovert,
   SHREC_INFO & shrec_info,
   ISOVERT_INFO & isovert_info);

  // **************************************************
  // DEFAULT PARAMETERS
  // **************************************************

  class SHREC_DEFAULTS {

  public:

    /// Isosurface vertex position method.
    static const VERTEX_POSITION_METHOD vertex_position_method
      = GRADIENT_POSITIONING; 

    /// Gradient selection.
    static const GRAD_SELECTION_METHOD grad_selection_method = GRAD_5x5x5;

    /// Extend mapping beyond 3x3x3 region.
    static const bool flag_map_extended = true;

    /// Use mod 6 cube selection.
    static const bool flag_select_mod6 = true;

    /// Maximum distance from cube boundary to isosurface vertex outside cube.
    static const COORD_TYPE max_dist = 1.0;

    /// Maximum small eigenvalue.
    /// Eigenvalues below this value are set to zero.
    static const EIGENVALUE_TYPE max_small_eigenvalue = 0.1;

    /// Offset for cube boundary for gradient selection.
    ///   Only gradients whose planes intersect the offset cube
    ///   are selected.
    static const SIGNED_COORD_TYPE grad_selection_cube_offset = 0.5;
  };

}

#endif

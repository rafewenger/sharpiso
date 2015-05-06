/// \file mergesharp_merge.h
/// Merge cubes containing sharp vertices.

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


#ifndef _MERGESHARP_MERGE_
#define _MERGESHARP_MERGE_

#include <unordered_map>
#include <vector>

#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

namespace MERGESHARP {


  // **************************************************
  // Type definitions
  // **************************************************

  class HASH_VERTEX_PAIR;
  class CYCLE_VERTEX;

  typedef std::unordered_map<VERTEX_PAIR, NUM_TYPE, HASH_VERTEX_PAIR> 
  EDGE_HASH_TABLE;

  // **************************************************
  // Merge some isosurface vertices
  // **************************************************

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & isoquad_cube,
   std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info);

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & isoquad_cube,
   SHARPISO_INFO & sharpiso_info);

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  /// Allows multiple isosurface vertices per cube.
  void merge_sharp_iso_vertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   ISOVERT & isovert,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & poly_vert,
   std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info);

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  /// Allows multiple isosurface vertices per cube.
  void merge_sharp_iso_vertices_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<DUAL_ISOVERT> & iso_vlist,
   ISOVERT & isovert,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & poly_vert,
   SHARPISO_INFO & sharpiso_info);


  // **************************************************
  // ROUTINE: is_isopatch_disk3D
  // **************************************************

  /// Return true if isopatch incident on vertex is a disk.
  /// @param tri_vert Triangle vertices.
  /// @param quad_vert Quadrilateral vertices in order around quadrilateral.
  /// @pre Assumes the boundary of the isopatch is the link of some vertex.
  bool is_isopatch_disk3D
  (const std::vector<ISO_VERTEX_INDEX> & tri_vert,
   const std::vector<ISO_VERTEX_INDEX> & quad_vert);

  /// Get list of cubes merged with icube.
  void get_merged_cubes
  (const SHARPISO_GRID & grid,
   const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const AXIS_SIZE_TYPE dist2cube,
   std::vector<VERTEX_INDEX> & merged_cube_list);

  /// Get edges on boundary of merged cubes.
  void get_merged_boundary_edges
  (const SHARPISO_GRID & grid,
   const std::vector<VERTEX_INDEX> & merged_cube_list,
   std::vector<EDGE_INDEX> & boundary_edge_list);

  /// Extract dual isosurface patch with vertex in merged cube.
  /// Returns list of (possibly degenerate) quadrilaterals.
  void extract_dual_quad_isopatch_incident_on
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const AXIS_SIZE_TYPE dist2cube,
   std::vector<ISO_VERTEX_INDEX> & isoquad_gcube,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex);

  /// Extract dual isosurface patch with vertex in merged cube.
  /// Returns non-degenerate triangles and quadrilaterals.
  void extract_dual_isopatch_incident_on
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const AXIS_SIZE_TYPE dist2cube,
   std::vector<ISO_VERTEX_INDEX> & tri_vert,
   std::vector<ISO_VERTEX_INDEX> & quad_vert);

  /// Extract dual isosurface patch with vertex in merged cube.
  /// Allow multiple isosurface vertices in each cube.
  void extract_dual_isopatch_incident_on_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const VERTEX_INDEX cube_index,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const AXIS_SIZE_TYPE dist2cube,
   std::vector<ISO_VERTEX_INDEX> & tri_vert,
   std::vector<ISO_VERTEX_INDEX> & quad_vert);

  /// Extract dual isosurface patch with vertex in merged cube.
  /// Allow multiple isosurface vertices in each cube.
  /// Version returning iso_vlist.
  void extract_dual_isopatch_incident_on_multi
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const ISOVERT & isovert,
   const VERTEX_INDEX cube_index,
   const std::vector<VERTEX_INDEX> & gcube_map,
   const AXIS_SIZE_TYPE dist2cube,
   std::vector<ISO_VERTEX_INDEX> & tri_vert,
   std::vector<ISO_VERTEX_INDEX> & quad_vert,
   std::vector<DUAL_ISOVERT> & iso_vlist);

  /// Insert triangle and quadrilateral edges into edge hash table.
  void insert_tri_quad_edges
  (const std::vector<ISO_VERTEX_INDEX> & tri_vert,
   const std::vector<ISO_VERTEX_INDEX> & quad_vert,
   EDGE_HASH_TABLE & edge_hash);

  /// Construct list of boundary cycle vertices.
  void construct_boundary_cycle
  (const EDGE_HASH_TABLE & edge_hash,
   std::vector<CYCLE_VERTEX> & cycle_vertex);

  // **************************************************
  // EDGE HASH TABLE
  // **************************************************

  struct HASH_VERTEX_PAIR {

    std::hash<VERTEX_INDEX> hash_func;

    HASH_VERTEX_PAIR(){};

    size_t operator() (const VERTEX_PAIR & key) const
    {
      return(hash_func(key.first+key.second));
    };
  };


  class CYCLE_VERTEX {

  public:
    VERTEX_INDEX adjacent[2];
    NUM_TYPE num_adjacent;
    bool is_visited;

    CYCLE_VERTEX() {
      num_adjacent = 0; 
      is_visited = false;
    };
  };

};


#endif

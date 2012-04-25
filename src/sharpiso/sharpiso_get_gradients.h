/// \file sharpiso_get_gradients.h
/// Get gradients in cube and cube neighbors.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2012 Rephael Wenger
 
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

#ifndef _SHARPISO_GET_GRADIENTS_
#define _SHARPISO_GET_GRADIENTS_

#include "sharpiso_cubes.h"
#include "sharpiso_grids.h"
#include "sharpiso_types.h"

#include "ijkgrid.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"


/// Definitions
namespace SHARPISO {
  
  // **************************************************
  // C++ CLASSES
  // **************************************************
  
  class OFFSET_CUBE_111;
  class GET_GRADIENTS_PARAM;
  
  // **************************************************
  // ROUTINES TO GET GRADIENTS
  // **************************************************

  /// Get selected grid vertex gradients.
  /// @param vertex_list[] List of vertices.
  /// @param vertex_flag[] Boolean array. If vertex_flag[i] is true, then
  ///        vertex vertex_list[i] is selected.
  /// @pre Size of vertex_flag[] is at least size of vertex_list[].
  void get_selected_vertex_gradients
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const GRADIENT_GRID_BASE & gradient_grid,
     const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
     const bool vertex_flag[],
     std::vector<COORD_TYPE> & point_coord,
     std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
     std::vector<SCALAR_TYPE> & scalar,
     NUM_TYPE & num_gradients);

  /// Get selected grid vertex gradients.
  /// std::vector variation.
  /// @pre Size of vertex_flag[] is at least size of vertex_list[].
  void get_selected_vertex_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const std::vector<VERTEX_INDEX> & vertex_list,
   const bool vertex_flag[],
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);

  /// Get grid vertex gradients.
  /// @pre point_coord[] is preallocated to size at least num_vertices*DIM3.
  /// @pre gradient_coord[] is preallocated to size at least num_vertices*DIM3.
  /// @pre scalar is preallocated to size at least num_vertices.
  void get_vertex_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
   COORD_TYPE point_coord[], GRADIENT_COORD_TYPE gradient_coord[],
   SCALAR_TYPE scalar[]);

  /// Get all 8 cube gradients
  void get_cube_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   std::vector<COORD_TYPE> & point_coord,
   GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3],
   SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D]);
  
  /// Get gradients in cube and neighboring cubes.
  /// @param sharp_isovert_param Parameters to determine 
  ///   which gradients are selected.
  void get_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GET_GRADIENTS_PARAM & get_gradients_param,
   const OFFSET_CUBE_111 & cube_111,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);
  
  /// Get large gradients at cube vertices.
  void get_large_cube_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_grad,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);
  
  /// Get large gradients at cube and neighboring cube vertices.
  /// @pre cube_index is the index of the lowest/leftmost cube vertex.
  ///      0 <= cube_index < number of grid vertices and
  ///      the vertex with index cube_index is not on the right/top of grid.
  /// @param cube_111 Cube with origin at (1-offset,1-offset,1-offset).
  void get_large_cube_neighbor_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_grad,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);
  
  /// Get selected gradients at cube and neighboring cube vertices.
  /// Selected gradient have magnitudes at least max_small_grad.
  /// Isosurfaces from selected neighboring gradients must intersect cube.
  /// @param cube_index Index of the lowest/leftmost cube vertex.
  ///      0 <= cube_index < number of grid vertices and
  ///      the vertex with index cube_index is not on the right/top of grid.
  /// @param cube_111 Cube with origin at (1-offset,1-offset,1-offset).
  ///      Data structure for processing cube-isosurface intersections.
  void get_selected_cube_neighbor_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients,
   const OFFSET_CUBE_111 & cube_111);

  /// Get gradients at endpoints of grid edges which intersect the isosurface.
  void get_intersected_edge_endpoint_gradients
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const GRADIENT_GRID_BASE & gradient_grid,
     const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
     const SCALAR_TYPE isovalue,
     std::vector<COORD_TYPE> & point_coord,
     std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
     std::vector<SCALAR_TYPE> & scalar,
     NUM_TYPE & num_gradients);

  /// Get gradients of vertices which determine edge isosurface intersections.
  /// @param zero_tolerance No division by numbers less than or equal 
  ///        to zero_tolerance.
  /// @pre zero_tolerance must be non-negative.
  void get_gradients_determining_edge_intersections
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE zero_tolerance,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);

  /// Get gradients at endpoints of neighboring grid edges
  ///   which intersect the isosurface.
  void get_intersected_neighbor_edge_endpoint_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);

  // **************************************************
  // GET VERTICES
  // **************************************************

  /// Get all 8 cube vertices.
  void get_cube_vertices
    (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
     VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D]);

  /// Get vertices at endpoints of cube edges which intersect the isosurface.
  void get_intersected_cube_edge_endpoints
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const VERTEX_INDEX cube_index, const SCALAR_TYPE isovalue,
     VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices);

  /// Get cube vertices determining the intersection of isosurface and edges.
  void get_cube_vertices_determining_edge_intersections
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
     const SCALAR_TYPE isovalue, const GRADIENT_COORD_TYPE zero_tolerance,
     VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices);

  /// Get cube vertices determining the intersection of isosurface and edges.
  /// Allow duplicate gradients.
  void get_cube_vertices_determining_edgeI_allow_duplicates
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
     const SCALAR_TYPE isovalue, const GRADIENT_COORD_TYPE zero_tolerance,
     std::vector<VERTEX_INDEX> & vertex_list);

  /// Get vertices of cube and cube neighbors.
  void get_cube_neighbor_vertices
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
   std::vector<VERTEX_INDEX> & vertex_list);

  /// Get vertices at endpoints of cube and neighbor edges 
  ///   which intersect the isosurface.
  void get_intersected_cube_neighbor_edge_endpoints
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & vertex_list);

  /// Get selected vertices (vertices where vertex_flag[] is true.)
  /// Reorder vertex list so that selected vertices are first.
  /// @pre Size of vertex_flag[] is at least size of vertex_list[].
  void get_selected_vertices
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const bool vertex_flag[],
   const NUM_TYPE num_vertices,
   VERTEX_INDEX vertex_list[], NUM_TYPE & num_selected);


  // **************************************************
  // SELECTION FUNCTIONS
  // **************************************************

  /// Select gradients at cube vertices.
  /// Select large gradients which give a level set intersecting the cube.
  /// @param cube_111 Cube with origin at (1-offset,1-offset,1-offset).
  void select_cube_gradients_based_on_isoplanes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_grad,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients,
   const OFFSET_CUBE_111 & cube_111);

  /// Set vertex_flag[i] to false for any vertex_list[i] 
  ///   with small gradient magnitude.
  /// @pre Array vertex_flag[] is preallocated with size 
  ///      at least vertex_list.size().
  void deselect_vertices_with_small_gradients
  (const GRADIENT_GRID_BASE & gradient_grid, 
   const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
   const GRADIENT_COORD_TYPE max_small_grad,
   bool vertex_flag[]);

  /// Set to false vertex_flag[i] for any vertex_list[i] 
  ///   with small gradient magnitude.
  /// std::vector variation.
  void deselect_vertices_with_small_gradients
  (const GRADIENT_GRID_BASE & gradient_grid, 
   const std::vector<VERTEX_INDEX> & vertex_list,
   const GRADIENT_COORD_TYPE max_small_mag_squared,
   bool vertex_flag[]);

  /// Set vertex_flag[i] to false for any vertex_list[i] 
  ///   determining an isoplane which does not intersect the cube.
  /// @pre Array vertex_flag[] is preallocated with size 
  ///      at least vertex_list.size().
  void deselect_vertices_based_on_isoplanes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid, 
   const GRID_COORD_TYPE * cube_coord, const OFFSET_CUBE_111 & cube_111,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
   bool vertex_flag[]);

  /// Set to false vertex_flag[i] for any vertex_list[i] 
  ///   determining an isoplane which does not intersect the cube.
  /// std::vector variation.
  void deselect_vertices_based_on_isoplanes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid, 
   const GRID_COORD_TYPE * cube_coord, const OFFSET_CUBE_111 & cube_111,
   const SCALAR_TYPE isovalue,
   const std::vector<VERTEX_INDEX> & vertex_list,
   bool vertex_flag[]);


  // **************************************************
  // GET_GRADIENTS_PARAM
  // **************************************************

  /// Parameters for getting gradients in cube or neighboring cubes.
  class GET_GRADIENTS_PARAM {
  protected:
    void Init();
    
  public:
    bool use_only_cube_gradients;
    bool use_selected_gradients;
    bool use_intersected_edge_endpoint_gradients;
    bool use_gradients_determining_edge_intersections;
    bool allow_duplicates;
    SIGNED_COORD_TYPE grad_selection_cube_offset;
    GRADIENT_COORD_TYPE zero_tolerance;

    /// Gradients with magnitude less than max_small_magnitude 
    ///   are treated as zero gradients.
    GRADIENT_COORD_TYPE max_small_magnitude;
    
    /// Constructor
    GET_GRADIENTS_PARAM() { Init(); };
  };

  // **************************************************
  // OFFSET_CUBE_111
  // **************************************************
  
  /// Cube with origin at (1-offset,1-offset,1-offset)
  ///   and edge length (1+2*offset).
  class OFFSET_CUBE_111:public SHARPISO_CUBE {
    
  protected:
    COORD_TYPE offset;
    
  public:
    OFFSET_CUBE_111(const SIGNED_COORD_TYPE offset);
    
    // get functions
    
    COORD_TYPE Offset() const       ///< Cube edge offset.
    { return(offset); }
    
    // set functions
    
    /// Set cube edge offset.
    void SetOffset(const SIGNED_COORD_TYPE offset);
    
    /// Undefine SetVertexCoord
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
  };
  
};

#endif

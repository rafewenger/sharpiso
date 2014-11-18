/// \file sharpiso_get_gradients.h
/// Get gradients in cube and cube neighbors.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2012-2014 Rephael Wenger
 
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

#include <string>

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

  class OFFSET_VOXEL;
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
  /// @pre scalar[] is preallocated to size at least num_vertices.
  void get_vertex_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
   COORD_TYPE point_coord[], GRADIENT_COORD_TYPE gradient_coord[],
   SCALAR_TYPE scalar[]);

  /// Get grid vertex gradients.
  /// point_coord[], gradient_coord[] and scalar[] are class std::vector.
  void get_vertex_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE>  & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar);

  /// Get grid vertex gradients.
  void get_vertex_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const std::vector<VERTEX_INDEX> & vertex_list,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE>  & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar);

  /// Get all 8 cube gradients
  void get_cube_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3],
   SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D]);
  
  /// Get gradients in cube and neighboring cubes.
  /// @param sharpiso_param Determines which gradients are selected.
  /// @param flag_sort If true, sort gradients.  
  ///        Overrides flag_sort in sharpiso_param.
  void get_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GET_GRADIENTS_PARAM & sharpiso_param,
   const OFFSET_VOXEL & voxel,
   const bool flag_sort_gradients,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);

  /// Get gradients from two cubes sharing a facet.
  /// Used in getting gradients around a facet.
  /// @param sharpiso_param Determines which gradients are selected.
  void get_two_cube_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX facet_v0, const NUM_TYPE orth_dir,
   const SCALAR_TYPE isovalue,
   const GET_GRADIENTS_PARAM & sharpiso_param,
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
  /// @param voxel Voxel with dimensions spacing[d]*(1+2*offset_factor).
  ///      Data structure for processing voxel-isosurface intersections.
  void get_selected_cube_neighbor_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients,
   const OFFSET_VOXEL & voxel);

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
  void get_gradients_determining_edge_intersections
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
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

  /// Get gradients from list of edge-isosurface intersections.
  /// @param sharpiso_param Determines which gradients are selected.
  /// @param flag_sort_gradients If true, sort gradients.  
  ///        Overrides flag_sort_gradients in sharpiso_param.
  void get_gradients_from_list
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<COORD_TYPE> & edgeI_coord,
   const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
   const SHARPISO_EDGE_INDEX_GRID & edge_index,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GET_GRADIENTS_PARAM & sharpiso_param,
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

  /// Get vertices at endpoints of cube edges which intersect the isosurface.
  /// Select endpoints whose isoplane intersects the ray 
  ///   pointing to the opposing endpoint.
  void get_intersected_cube_edge_endpoints_select
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index, const SCALAR_TYPE isovalue, 
   VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices);

  /// Get cube vertices determining the intersection of isosurface and edges.
  void get_cube_vertices_determining_edge_intersections
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
     const SCALAR_TYPE isovalue, 
     VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices);

  /// Get cube vertices determining the intersection of isosurface and edges.
  /// Allow duplicate gradients.
  void get_cube_vertices_determining_edgeI_allow_duplicates
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
     const SCALAR_TYPE isovalue, 
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

  /// Get vertex whose gradient determines intersection of isosurface 
  ///    and edge (iv0,iv1)
  void get_vertex_determining_edge_intersection
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
   VERTEX_INDEX & iv2);

  /// Compute point on edge where gradient changes.
  /// @param zero_tolerance No division by numbers less than or equal 
  ///        to zero_tolerance.
  /// @pre zero_tolerance must be non-negative.
  void compute_gradient_change_on_edge
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
   const GRADIENT_COORD_TYPE zero_tolerance,
   bool & flag_no_split,
   SCALAR_TYPE & t_split, 
   SCALAR_TYPE & s_split);

  /// Get vertices at endpoints of cube edges which intersect the isosurface.
  /// Two cubes share a facet.
  /// @param facet_v0 Index of primary (lowest-left) facet vertex.
  /// @param orth_dir Direction orthogonal to facet.
  /// @pre Facet is in interior of grid.
  void get_intersected_two_cube_edge_endpoints
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const VERTEX_INDEX facet_v0, const int orth_dir, 
     const SCALAR_TYPE isovalue, 
     VERTEX_INDEX vertex_list[NUM_TWO_CUBE_VERTICES3D], 
     NUM_TYPE & num_vertices);

  /// Get cube vertices with large gradients.
  void get_cube_vertices_with_large_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_magnitude,
   std::vector<VERTEX_INDEX> & vertex_list);

  /// Get vertices with large gradients magnitudes.
  void get_vertices_with_large_gradient_magnitudes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_magnitude,
   const NUM_TYPE max_dist,
   std::vector<VERTEX_INDEX> & vertex_list);

  /// Get vertices with large gradients magnitudes.
  void get_vertices_with_large_gradient_magnitudes
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const GRADIENT_GRID_BASE & gradient_grid,
     const VERTEX_INDEX cube_index,
     const GET_GRADIENTS_PARAM & gradient_param,
     std::vector<VERTEX_INDEX> & vertex_list);

  /// Get intersected edge endpoints in large neighborhood.
  /// Old, deprecated version.
  void get_ie_endpoints_in_large_neighborhood_old_version
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const GRADIENT_GRID_BASE & gradient_grid,
     const SCALAR_TYPE isovalue,
     const VERTEX_INDEX cube_index,
     const GET_GRADIENTS_PARAM & gradient_param,
     std::vector<VERTEX_INDEX> & vertex_list);

  /// Get intersected edge endpoints in large neighborhood.
  void get_ie_endpoints_in_large_neighborhood
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube_index,
   const GET_GRADIENTS_PARAM & gradient_param,
   std::vector<VERTEX_INDEX> & vertex_list);

  // **************************************************
  // SORT VERTICES
  // **************************************************

  /// Sort vertices based on the distance of the isoplane to point pcoord[].
  void sort_vertices_by_isoplane
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const COORD_TYPE pcoord[DIM3],
   VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices);

  /// Sort vertices based on the distance of the isoplane to point pcoord[].
  void sort_vertices_by_isoplane
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const COORD_TYPE pcoord[DIM3],
   std::vector<VERTEX_INDEX> & vertex_list);

  /// Sort vertices based on the distance of the isoplane distance 
  ///   to cube center.
  void sort_vertices_by_isoplane_dist2cc
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube_index,
   VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices);

  /// Sort vertices based on the distance of the isoplane distance 
  ///   to cube center.
  void sort_vertices_by_isoplane_dist2cc
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX cube_index,
   std::vector<VERTEX_INDEX> & vertex_list);

  // **************************************************
  // SELECTION FUNCTIONS
  // **************************************************

  /// Select gradients at cube vertices.
  /// Select large gradients which give a level set intersecting the cube.
  /// @param voxel Voxel with dimensions spacing[d]*(1+2*offset_factor).
  ///      Data structure for processing voxel-isosurface intersections.
  void select_cube_gradients_based_on_isoplanes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients,
   const OFFSET_VOXEL & voxel);

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
   const GRID_COORD_TYPE * cube_coord, const OFFSET_VOXEL & voxel,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX * vertex_list, const NUM_TYPE num_vertices,
   bool vertex_flag[]);

  /// Set to false vertex_flag[i] for any vertex_list[i] 
  ///   determining an isoplane which does not intersect the cube.
  /// std::vector variation.
  void deselect_vertices_based_on_isoplanes
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid, 
   const GRID_COORD_TYPE * cube_coord, const OFFSET_VOXEL & voxel,
   const SCALAR_TYPE isovalue,
   const std::vector<VERTEX_INDEX> & vertex_list,
   bool vertex_flag[]);

  // **************************************************
  // MAP GRAD_SELECTION_METHOD TO/FROM C++ string
  // **************************************************

  GRAD_SELECTION_METHOD get_grad_selection_method(const std::string & s);
  void get_grad_selection_string
  (const GRAD_SELECTION_METHOD grad_sel, std::string & s);

  // **************************************************
  // GET_GRADIENTS_PARAM
  // **************************************************

  /// Parameters for getting gradients in cube or neighboring cubes.
  class GET_GRADIENTS_PARAM {
  protected:
    void Init();

    GRAD_SELECTION_METHOD grad_selection_method;

  public:
    bool use_only_cube_gradients;
    bool use_large_neighborhood;
    bool use_zero_grad_boundary;
    bool use_diagonal_neighbors;
    bool use_selected_gradients;
    bool select_based_on_grad_dir;
    bool use_intersected_edge_endpoint_gradients;
    bool use_gradients_determining_edge_intersections;
    bool allow_duplicates;
    bool flag_sort_gradients;               ///< Sort gradients.
    bool use_new_version;                  ///< If true, use new version.
    SIGNED_COORD_TYPE grad_selection_cube_offset;
    GRADIENT_COORD_TYPE zero_tolerance;

    /// Gradients with magnitude less than max_small_magnitude 
    ///   are treated as zero gradients.
    GRADIENT_COORD_TYPE max_small_magnitude;

    /// Maximum distance (in edges) of neighborhood vertices from cube.
    NUM_TYPE max_grad_dist;

    /// Constructor
    GET_GRADIENTS_PARAM() { Init(); };

    // Get functions.
    GRAD_SELECTION_METHOD GradSelectionMethod() const
    { return(grad_selection_method); };

    // Set functions.
    void SetGradSelectionMethod
      (const GRAD_SELECTION_METHOD grad_selection_method);
  };

  // **************************************************
  // VOXEL
  // **************************************************
  
  /// 3D grid voxel.
  class VOXEL : public IJK::CUBE_INFO<int,NUM_TYPE> {

  public:
    typedef SHARPISO::COORD_TYPE COORD_TYPE;
    
  protected:
    /// vertex_coord[dimension*k+j] = 
    ///   j'th coordinate of k'th vertex of voxel
    COORD_TYPE * vertex_coord;

    /// max_vertex_index.  Stored for faster processing.
    NUM_TYPE max_vertex_index;

    void InitLocal();
    void FreeLocal();

    void SetMaxVertexIndex();       ///< Set max vertex index.

  public:

    VOXEL();
    ~VOXEL();

    // get functions

    NUM_TYPE MaxVertexIndex() const    ///< Maximum cube vertex index.
    { return(max_vertex_index); }

    NUM_TYPE OppositeVertex            ///< Index of vertex opposite iv
    (const NUM_TYPE iv) const
    { return(MaxVertexIndex()-iv); }

    /// Return pointer to vertex coordinates
    const COORD_TYPE * VertexCoord() const
    { return(vertex_coord); }

    /// Return pointer to coordinates of k'th cube vertex
    const COORD_TYPE * VertexCoord(const NUM_TYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return j'th coordinate of k'th vertex
    const COORD_TYPE VertexCoord         
    (const NUM_TYPE k, const NUM_TYPE j) const
    { return(vertex_coord[this->Dimension()*k+j]); }

    /// Return pointer to coordinates of endpoint 0 of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const COORD_TYPE * DiagonalEnd0Coord(const NUM_TYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return pointer to coordinates of endpoint 1 of diagonal k.
    /// Corner opposite cube corner k is endpoint 1 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const COORD_TYPE * DiagonalEnd1Coord(const NUM_TYPE k) const
    { return(vertex_coord+this->Dimension()*OppositeVertex(k)); }


    // set routines

    /// Redefine SetVertexCoord().
    void SetVertexCoord             ///< Set voxel vertex coordinates.
    (const COORD_TYPE min_coord[DIM3], const COORD_TYPE max_coord[DIM3]);

    /// Undefine SetDimension().
    void SetDimension(const int dimension);
  };

  /// 3D grid voxel and offset.
  class OFFSET_VOXEL : public VOXEL {

  protected:
    SIGNED_COORD_TYPE offset_factor;

  public:
    OFFSET_VOXEL() {};

    /// Return offset factor.
    SIGNED_COORD_TYPE OffsetFactor() const
    { return(offset_factor); }

    /// Set coordinates of voxel plus offset.
    void SetVertexCoord
    (const COORD_TYPE spacing[DIM3], 
     const SIGNED_COORD_TYPE offset_factor);
  };

};

#endif

  /// \file sharpiso_feature.h
  /// Compute sharp isosurface vertices and edges.
  /// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
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

#ifndef _SHARPISO_FEATURES_
#define _SHARPISO_FEATURES_

#include "sharpiso_cubes.h"
#include "sharpiso_eigen.h"
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
  class SVD_INFO {
  public:
	GRADIENT_COORD_TYPE  ray_direction[DIM3];     //ray direction 
	SCALAR_TYPE   ray_initial_point[DIM3]; //ray initial point
	LOC_TYPE location; //defined in sharpiso_types typedef enum
	//{ NONE, CENTROID, CUBE_CENTER } LOC_TYPE ;
	bool ray_intersect_cube;
	};
  
    // **************************************************
    // ROUTINES TO COMPUTE SHARP VERTEX/EDGE
    // **************************************************
  
    /// Compute sharp isosurface vertex using singular valued decomposition.
  void svd_compute_sharp_vertex_in_cube
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_zero_mag,
   const EIGENVALUE_TYPE eigenvalue_tolerance,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_nonzero_eigenvalues,
   SVD_INFO &svd_debug_info);
  
  /// Compute sharp isosurface vertex using simple interpolation on the edges.
   void svd_compute_sharp_vertex_in_cube_edge_based_simple
(const SHARPISO_SCALAR_GRID & scalar_grid, 
 const GRADIENT_GRID & gradient_grid, 
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE max_small_mag,
 const EIGENVALUE_TYPE max_small_eigenvalue,
 COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO &svd_debug_info);
 
 /// Compute sharp isosurface vertex using complex interpolation on the edges
 void svd_compute_sharp_vertex_in_cube_edge_based_cmplx
(const SHARPISO_SCALAR_GRID & scalar_grid, 
 const GRADIENT_GRID & gradient_grid, 
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE max_small_mag,
 const EIGENVALUE_TYPE max_small_eigenvalue,
 COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO &svd_debug_info);
  
    /// Compute sharp isosurface vertex using subgrid sampling of cube.
    /// Use only cube vertex gradients.
  void subgrid_compute_sharp_vertex_in_cube
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
   const VERTEX_INDEX cube_index, const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag, const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3], 
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);
  
    /// Compute sharp isosurface vertex using subgrid sampling of cube.
    /// Use gradients from neighborhood around cube.
  void subgrid_compute_sharp_vertex_neighborhood
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3], 
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error,
   const OFFSET_CUBE_111 & cube_111);
  
    /// Calculate isosurface vertex using regular subgrid of the cube.
  void subgrid_calculate_iso_vertex_in_cube
  (const COORD_TYPE * point_coord, const GRADIENT_COORD_TYPE * gradient_coord,
   const SCALAR_TYPE * scalar, const NUM_TYPE num_points, 
   const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue, 
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3], 
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);
  
    /// Calculate isosurface vertex using regular subgrid of the cube.
  void subgrid_calculate_iso_vertex_in_cube
  (const std::vector<COORD_TYPE> & point_coord, 
   const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   const std::vector<SCALAR_TYPE> & scalar, 
   const NUM_TYPE num_points, 
   const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue, 
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3], 
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);
    //debug/
  
    /// Compute centroid of intersections of isosurface and grid edges
  void compute_isosurface_grid_edge_centroid
  (const SHARPISO_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX iv, 
   COORD_TYPE * coord);
  
    // **************************************************
    // ROUTINES TO GET GRADIENTS
    // **************************************************

  // Get all 8 cube gradients
  void get_cube_gradients
    (const SHARPISO_SCALAR_GRID & scalar_grid, 
     const GRADIENT_GRID & gradient_grid, 
     const VERTEX_INDEX cube_index,
     std::vector<COORD_TYPE> & point_coord,
     GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3],
     SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D]);
  
    /// Get large gradients at cube vertices.
  void get_large_cube_gradients
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
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
  void get_large_cube_neighbor_gradients
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
   const VERTEX_INDEX cube_index,
   const GRADIENT_COORD_TYPE max_small_grad,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients);
  
    /// Get selected gradients at cube and neighboring cube vertices.
    /// Selected gradient have magnitudes at least max_small_grad.
    /// Isosurfaces from selected neighboring gradients must intersect cube.
    /// @param cube_111 Data structure for processing 
    ///       cube-isosurface intersections.
    /// @pre cube_index is the index of the lowest/leftmost cube vertex.
    ///      0 <= cube_index < number of grid vertices and
    ///      the vertex with index cube_index is not on the right/top of grid.
  void get_selected_cube_neighbor_gradients
  (const SHARPISO_SCALAR_GRID & scalar_grid, 
   const GRADIENT_GRID & gradient_grid, 
   const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
   const SCALAR_TYPE isovalue,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients,
   const OFFSET_CUBE_111 & cube_111);
  
    // **************************************************
    // OFFSET_CUBE_111
    // **************************************************
  
    /// Unit cube with origin at (1-offset,1-offset,1-offset)
    ///   and edge length (1+2*offset).
  class OFFSET_CUBE_111:public SHARPISO_CUBE {
    
  protected:
    COORD_TYPE offset;
    
  public:
    OFFSET_CUBE_111(const COORD_TYPE offset);
    
      // get functions
    
    COORD_TYPE Offset() const       ///< Cube edge offset.
    { return(offset); }
    
      // set functions
    
      /// Undefine SetVertexCoord
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
    
    
  };
  
};

#endif

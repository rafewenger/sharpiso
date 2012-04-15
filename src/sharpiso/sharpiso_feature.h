/// \file sharpiso_feature.h
/// Compute sharp isosurface vertices and edges.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2011,2012 Rephael Wenger
 
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

#ifndef _SHARPISO_FEATURE_
#define _SHARPISO_FEATURE_

#include "sharpiso_cubes.h"
#include "sharpiso_eigen.h"
#include "sharpiso_get_gradients.h"
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
  class SHARP_ISOVERT_PARAM;
  class SVD_INFO;
  
  // **************************************************
  // SVD ROUTINES TO COMPUTE SHARP VERTEX/EDGE
  // **************************************************
  
  /// Compute sharp isosurface vertex using singular valued decomposition.
  void svd_compute_sharp_vertex_for_cube
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   const OFFSET_CUBE_111 & cube_111,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);
  
  /// Compute sharp isosurface vertex on the ray
  void compute_vertex_on_ray
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   const COORD_TYPE ray_origin[DIM3],
   const COORD_TYPE ray_direction[DIM3],
   COORD_TYPE sharp_coord[DIM3],
   bool & flag_use_centroid,
   SVD_INFO & svd_info);
  
  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use only cube vertex gradients.
  void svd_compute_sharp_vertex_in_cube
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_zero_mag,
   const EIGENVALUE_TYPE eigenvalue_tolerance,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);
  
  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use selected cube vertex gradients.
  void svd_compute_sharp_vertex_in_cube_S
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_zero_mag,
   const EIGENVALUE_TYPE eigenvalue_tolerance,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info,
   const OFFSET_CUBE_111 & cube_111);
  
  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use gradients from cube and neighboring cubes.
  /// @param cube_111 Cube with origin at (1-offset,1-offset,1-offset).
  void svd_compute_sharp_vertex_neighborhood
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_zero_mag,
   const EIGENVALUE_TYPE eigenvalue_tolerance,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);
  
  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use gradients from cube and neighboring cubes.
  /// @param cube_111 Cube with origin at (1-offset,1-offset,1-offset).
  void svd_compute_sharp_vertex_neighborhood_S
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_zero_mag,
   const EIGENVALUE_TYPE eigenvalue_tolerance,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info,
   const OFFSET_CUBE_111 & cube_111);
  
  /// Compute sharp isosurface vertex using simple interpolation on the edges.
  // with sharp isovert param
  void svd_compute_sharp_vertex_in_cube_edge_based_simple
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const EIGENVALUE_TYPE max_small_eigenvalue,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   SVD_INFO & svd_info);
   // with out the sharp isovert param
  void svd_compute_sharp_vertex_in_cube_edge_based_simple
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const EIGENVALUE_TYPE max_small_eigenvalue,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);
  
  /// Compute sharp isosurface vertex using edge intersection.
  /// Compute edge intersections using endpoint gradients.
  void svd_compute_sharp_vertex_in_cube_edge_based_cmplx
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const EIGENVALUE_TYPE max_small_eigenvalue,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);
  
   void svd_compute_sharp_vertex_in_cube_edge_based_cmplx
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const EIGENVALUE_TYPE max_small_eigenvalue,
   const SCALAR_TYPE cube_offset2,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   SVD_INFO & svd_info);


  // **************************************************
  // SUBGRID ROUTINES TO COMPUTE SHARP VERTEX/EDGE
  // **************************************************
  
  /// Compute sharp isosurface vertex using subgrid sampling of cube.
  /// Use subgrid sampling to locate isosurface vertex on sharp edge/corner.
  void subgrid_compute_sharp_vertex_in_cube
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GET_GRADIENTS_PARAM & get_gradients_param,
   const OFFSET_CUBE_111 & cube_111,
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3],
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);

  /// Calculate isosurface vertex using regular subgrid of the cube.
  void subgrid_calculate_iso_vertex_in_cube
  (const COORD_TYPE * point_coord, const GRADIENT_COORD_TYPE * gradient_coord,
   const SCALAR_TYPE * scalar, const NUM_TYPE num_points,
   const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3],
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);
  
  /// Calculate isosurface vertex using regular subgrid of the cube.
  /// std::vector variation.
  void subgrid_calculate_iso_vertex_in_cube
  (const std::vector<COORD_TYPE> & point_coord,
   const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   const std::vector<SCALAR_TYPE> & scalar,
   const NUM_TYPE num_points,
   const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
   const NUM_TYPE subgrid_axis_size,
   COORD_TYPE sharp_coord[DIM3],
   SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error);
  
  // **************************************************
  // COMPUTE ISO VERTEX AT CENTROID
  // **************************************************
  
  /// Compute centroid of intersections of isosurface and grid edges
  void compute_isosurface_grid_edge_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
   COORD_TYPE * coord);

  // **************************************************
  // CLAMP POINTS TO threshold_cube_offset
  // **************************************************

  // clamp point : shpoint in the global coord 
  void  clamp_point
	  (const float threshold_cube_offset, COORD_TYPE * cube_coord,
      COORD_TYPE * shpoint);
	// clamp point: shpoint in the local coord
	void  clamp_point
	  (const float threshold_cube_offset, COORD_TYPE * shpoint);

  // **************************************************
  // SHARP_ISOVERT_PARAM
  // **************************************************
  
  /// Parameters for computing sharp isosurface vertices using svd.
  class SHARP_ISOVERT_PARAM:public GET_GRADIENTS_PARAM {

  protected:
    void Init();
    
  public:
    bool use_lindstrom;
    SIGNED_COORD_TYPE ray_intersection_cube_offset;

    /// Maximum (Linf) distance from cube to isosurface vertex
    SIGNED_COORD_TYPE max_dist;
    
    /// Normalized eigenvalues with value less than max_small_eigenvalue
    ///   are set to zero.
    EIGENVALUE_TYPE max_small_eigenvalue;
    
    /// Constructor
    SHARP_ISOVERT_PARAM() { Init(); };
  };
  
  // **************************************************
  // SVD_INFO
  // **************************************************
  
  /// Information returned from svd calculations.
  class SVD_INFO {
  public:
    LOC_TYPE location;                       // location type
    GRADIENT_COORD_TYPE ray_direction[DIM3]; // ray direction
    COORD_TYPE ray_initial_point[DIM3];     // point on ray
    COORD_TYPE ray_cube_intersection[DIM3];
    bool ray_intersect_cube;                 // true if ray intersects cube
    bool is_svd_point_in_cube;               // true if svd point is in cube
    
    /// Set ray information.
    void SetRayInfo
    (const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
     const COORD_TYPE intersect[DIM3], const bool flag_intersects_cube);
	};
  
  
};

#endif

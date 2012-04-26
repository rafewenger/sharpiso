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
  
  /// Compute sharp isosurface vertex on the line
  /// @param flag_conflict True if sharp_coord conflicts with other cube.
  void compute_vertex_on_line
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const GRID_COORD_TYPE cube_coord[DIM3],
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   const COORD_TYPE line_origin[DIM3],
   const COORD_TYPE line_direction[DIM3],
   COORD_TYPE sharp_coord[DIM3],
   bool & flag_conflict,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using simple interpolation on the edges.
  void svd_compute_sharp_vertex_in_cube_edge_based_simple
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using edge intersection.
  /// Compute edge intersections using endpoint gradients.
  // *** REDUNDANT PARAMETERS: max_small_mag, max_small_eigenvalue
  // *** MOVE SHARP_ISOVERT_PARAM before output parameters
  void svd_compute_sharp_vertex_in_cube_edge_based_cmplx
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE max_small_mag,
   const EIGENVALUE_TYPE max_small_eigenvalue,
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
  // ROUTINES TO MOVE POINTS
  // **************************************************

  /// Clamp to cube cube_coord[]
  void  clamp_point
  (const float offset,  const GRID_COORD_TYPE cube_coord[DIM3], 
   COORD_TYPE point[DIM3]);
  /// Clamp to unit cube (0,0,0) to (1,1,1).
	void  clamp_point(const float offset, COORD_TYPE point[DIM3]);

  /// Move point which lies in an occupied grid cube.
  void process_conflict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const GRID_COORD_TYPE cube_coord[DIM3],
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE iso_coord[DIM3],
   SVD_INFO & svd_info);

  void process_far_point
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const VERTEX_INDEX cube_index,
     const GRID_COORD_TYPE cube_coord[DIM3],
     const SCALAR_TYPE isovalue,
     const SHARP_ISOVERT_PARAM & sharpiso_param,
     COORD_TYPE iso_coord[DIM3],
     SVD_INFO & svd_info);

  // **************************************************
  // SHARP_ISOVERT_PARAM
  // **************************************************
  
  /// Parameters for computing sharp isosurface vertices using svd.
  class SHARP_ISOVERT_PARAM:public GET_GRADIENTS_PARAM {

  protected:
    void Init();
    
  public:
    bool use_lindstrom;            ///< If true, use Lindstrom formula
    bool flag_allow_conflict;      ///< If true, allow conflicts
    bool flag_clamp_conflict;      ///< If true, clamp conflicts to cube.
    bool flag_clamp_far;           ///< If true, clamp far points.
    bool flag_recompute_eigen2;    ///< If true, recompute with 2 eigenvalues.
    bool flag_round;               ///< Round output coordinates.
    bool use_Linf_dist;            ///< If true, use Linf dist.

    /// Maximum (Linf) distance from cube to isosurface vertex
    SIGNED_COORD_TYPE max_dist;

    /// Snap points within snap distance to cube.
    COORD_TYPE snap_dist;
    
    /// Normalized eigenvalues with value less than max_small_eigenvalue
    ///   are set to zero.
    EIGENVALUE_TYPE max_small_eigenvalue;

    /// Minimum separation distance between isosurface vertices.
    COORD_TYPE separation_distance;

    /// Round to nearest 1/round_denominator
    int round_denominator;

    /// Constructor
    SHARP_ISOVERT_PARAM() { Init(); };

    /// Set
    void Set(const SHARP_ISOVERT_PARAM & param);
  };
  
  // **************************************************
  // SVD_INFO
  // **************************************************
  
  /// Information returned from svd calculations.
  class SVD_INFO {
  public:
    LOC_TYPE location;                       // location type
    GRADIENT_COORD_TYPE ray_direction[DIM3]; // ray direction
    COORD_TYPE ray_initial_point[DIM3];      // point on ray
    COORD_TYPE ray_cube_intersection[DIM3];
    bool flag_conflict;
    
    /// Set ray information.
    void SetRayInfo
    (const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
     const COORD_TYPE intersect[DIM3]);
	};
  
  
};

#endif

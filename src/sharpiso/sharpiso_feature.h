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

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use Lindstrom's formula.
  void svd_compute_sharp_vertex_for_cube_lindstrom
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   const OFFSET_CUBE_111 & cube_111,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use line-cube intersection to compute sharp isosurface vertex
  ///    when number of eigenvalues is 2.
  /// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
  void svd_compute_sharp_vertex_for_cube_lc_intersection
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   const OFFSET_CUBE_111 & cube_111,
   COORD_TYPE sharp_coord[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
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
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using edge-isosurface intersections.
  /// Approximate gradients using linear interpolation on the grid edge.
  void svd_compute_sharp_vertex_edgeI_interpolate_gradients
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using edge-isosurface intersections.
  /// Use sharp formula for computing gradient at intersection.
  void svd_compute_sharp_vertex_edgeI_sharp_gradient
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex near facet 
  ///    using singular valued decomposition.
  /// @param facet_v0 Index of primary (lowest-left) facet vertex.
  /// @param orth_dir Direction orthogonal to facet.
  /// @pre Facet is in interior of grid.
  /// @param[out] sharp_vertex_location  -1, 0, or 1.
  ///    -1: Sharp vertex is below/left of facet.
  ///    0: Vertex not sharp or relative location undetermined.
  void svd_compute_sharp_vertex_near_facet
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX facet_v0,
   const NUM_TYPE facet_orth_dir,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   NUM_TYPE & sharp_vertex_location,
   COORD_TYPE sharp_coord[DIM3],
   GRADIENT_COORD_TYPE line_direction[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues);

  /// Compute sharp isosurface vertex near facet 
  ///    using singular valued decomposition.
  /// Short parameter list.
  void svd_compute_sharp_vertex_near_facet
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX facet_v0,
   const NUM_TYPE facet_orth_dir,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   NUM_TYPE & sharp_vertex_location);


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

  /// Postprocess isosurface vertex coordinates.
  /// Depending on flags in sharpiso_param, move far points,
  ///   move points in occupied grid cubes, and round coordinates.
  void postprocess_isovert_location
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const GRID_COORD_TYPE cube_coord[DIM3],
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE iso_coord[DIM3],
   SVD_INFO & svd_info);

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
    bool flag_allow_conflict;      ///< If true, allow conflicts
    bool flag_clamp_conflict;      ///< If true, clamp conflicts to cube.
    bool flag_clamp_far;           ///< If true, clamp far points.
    bool flag_recompute_eigen2;    ///< If true, recompute with 2 eigenvalues.
    bool flag_round;               ///< Round output coordinates.
    bool flag_remove_gradients;    ///< Remove gradients to resolve conflicts
    bool flag_reselect_gradients;  ///< Select fewer grads to resolve conflicts
    bool flag_centroid_eigen1;     ///< Use centroid for 1 or fewer eigenvalues
    bool use_lindstrom;            ///< If true, use Lindstrom formula
    bool use_lindstrom2;
    bool use_Linf_dist;            ///< If true, use Linf dist.


    /// If true, use distance to centroid of isosurface-edge intersections.
    /// Otherwise, use distance to cube center.
    bool flag_dist2centroid;       

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

  protected:
    void Init();

  public:
    LOC_TYPE location;                       // location type
    GRADIENT_COORD_TYPE ray_direction[DIM3]; // ray direction
    COORD_TYPE ray_initial_point[DIM3];      // point on ray
    COORD_TYPE ray_cube_intersection[DIM3];
    bool flag_conflict;
    bool flag_Linf_iso_vertex_location;
    
    /// Set ray information.
    void SetRayInfo
    (const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
     const COORD_TYPE intersect[DIM3]);

    SVD_INFO() { Init(); };
	};
  
  
};

#endif

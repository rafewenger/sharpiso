/// \file sharpiso_feature.h
/// Compute sharp isosurface vertices and edges.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2012-2014 Arindam Bhattacharya and Rephael Wenger
 
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

  class SHARP_ISOVERT_PARAM;
  class SVD_INFO;
  
  // **************************************************
  // SVD ROUTINES TO COMPUTE SHARP VERTEX/EDGE
  // **************************************************

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use Lindstrom's formula.
  /// @param pointX Compute vertex closest to pointX.
  void svd_compute_sharp_vertex_for_cube_lindstrom
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   const OFFSET_VOXEL & voxel,
   const COORD_TYPE pointX[DIM3],
   COORD_TYPE sharp_coord[DIM3],
   COORD_TYPE edge_direction[DIM3],
   COORD_TYPE orth_direction[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use Lindstrom's formula.
  /// Compute vertex closest to cube center or centroid 
  ///   of interpolated grid edge/isosurface intersections.
  void svd_compute_sharp_vertex_for_cube_lindstrom
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   const OFFSET_VOXEL & voxel,
   COORD_TYPE sharp_coord[DIM3],
   COORD_TYPE edge_direction[DIM3],
   COORD_TYPE orth_direction[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use Lindstrom's formula.
  /// If num_large_eigenvalues == 2, position vertex on plane.
  /// @param pointX Compute vertex closest to pointX.
  /// Also post processes vertices.
  void svd_compute_sharp_vertex_on_plane_lindstrom
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   const OFFSET_VOXEL & voxel,
   const COORD_TYPE pointX[DIM3],
   const COORD_TYPE plane_normal[DIM3],
   COORD_TYPE sharp_coord[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);

  /// Compute sharp isosurface vertex using singular valued decomposition.
  /// Use input edge-isosurface intersections and normals
  ///   to position isosurface vertices on sharp features.
  void svd_compute_sharp_vertex_for_cube_hermite
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<COORD_TYPE> & edgeI_coord,
   const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
   const SHARPISO_EDGE_INDEX_GRID & edge_index,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE sharp_coord[DIM3],
   EIGENVALUE_TYPE eigenvalues[DIM3],
   NUM_TYPE & num_large_eigenvalues,
   SVD_INFO & svd_info);


  // ********************************************************************
  // COMPUTE SHARP VERTEX/EDGE USING SVD & EDGE-ISOSURFACE INTERSECTIONS
  // ********************************************************************

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


  // ********************************************************************
  // COMPUTE SHARP ISOSURFACE VERTEX NEAR FACET
  // ********************************************************************

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
   const OFFSET_VOXEL & offset_voxel,
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

  /// Compute centroid of intersections of isosurface and grid edges.
  /// @param use_sharp_edgeI If true, use sharp formula for
  ///          intersection of isosurface and grid edges.
  void compute_edgeI_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
   const bool use_sharp_edgeI, COORD_TYPE * coord);

  /// Compute centroid of intersections of isosurface and grid edges.
  void compute_edgeI_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
   COORD_TYPE * coord);

  /// Compute centroid of intersections of isosurface and grid edges.
  /// Use sharp formula for computing intersection of isosurface and grid edge.
  void compute_edgeI_sharp_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
   COORD_TYPE * coord);

  /// Compute centroid of intersections of isosurface and grid edges.
  ///   Edge-isosurface intersections are given in an input list.
  void compute_edgeI_centroid
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<COORD_TYPE> & edgeI_coord,
   const SHARPISO_EDGE_INDEX_GRID & edge_index,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX cube_index,
   COORD_TYPE * coord);


  // **************************************************
  // POINT LOCATION/CONFLICT ROUTINES
  // **************************************************

  /// Return true if cube contains point.
  bool cube_contains_point
  (const COORD_TYPE cube_coord[DIM3], const COORD_TYPE point_coord[DIM3],
   const COORD_TYPE spacing[DIM3]);

  /// Return true if cube contains point.
  bool cube_contains_point
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
   const COORD_TYPE point_coord[DIM3]);

  /// Return index of cube containing point.
  /// Set flag_boundary to true if point is on cube boundary.
  /// @pre Point is contained in grid.
  /// @pre grid.AxisSize(d) > 0 for every axis d.
  void get_cube_containing_point
  (const SHARPISO_GRID & grid, const COORD_TYPE * coord,
   VERTEX_INDEX & cube_index, bool & flag_boundary);

  /// Return index of cube containing point.
  /// Set flag_boundary to true if point is on cube boundary.
  /// Set flag_active to true if cube is active.
  /// @pre Point is contained in grid.
  /// @pre grid.AxisSize(d) > 0 for every axis d.
  void get_cube_containing_point
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, const COORD_TYPE * coord,
   VERTEX_INDEX & cube_index, bool & flag_boundary,
   bool & flag_active);

  /// Return true if point lies in an occupied cube
  ///   other than the one given by cube_coord[].
  /// @param[out] conflicting_cube Index of conflicting cube.
  bool check_conflict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const COORD_TYPE cube_coord[DIM3],
   const COORD_TYPE point_coord[DIM3],
   VERTEX_INDEX & conflicting_cube);

  /// Return true if point lies in an occupied cube
  ///   other than the one given by cube_coord[].
  bool check_conflict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const COORD_TYPE cube_coord[DIM3],
   const COORD_TYPE point_coord[DIM3]);



  // **************************************************
  // ROUTINES TO MOVE POINTS
  // **************************************************

  /// Clamp to cube with lower coordinates cube_coord[].
  /// @param cube_coord[] Lower coordinates of cube.
  /// @param spacing[] Spacing along the grid axis.
  /// @param cube_offset Offset scale. 
  ///   Cube is enlarged by cube_offset*spacing[d] in direction d.
  /// @param point[] Clamp point.
  void clamp_point
    (const COORD_TYPE cube_coord[DIM3],
     const COORD_TYPE spacing[DIM3],
     const COORD_TYPE cube_offset,
     COORD_TYPE point[DIM3]);

  /// Postprocess isosurface vertex coordinates.
  /// Depending on flags in sharpiso_param, move far points,
  ///   move points in occupied grid cubes, and round coordinates.
  void postprocess_isovert_location
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const COORD_TYPE cube_coord[DIM3],
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE iso_coord[DIM3],
   SVD_INFO & svd_info);

  /// Move point which lies in an occupied grid cube.
  void process_conflict
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const COORD_TYPE cube_coord[DIM3],
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & sharpiso_param,
   COORD_TYPE iso_coord[DIM3],
   SVD_INFO & svd_info);

  void process_far_point
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const VERTEX_INDEX cube_index,
     const COORD_TYPE cube_coord[DIM3],
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
    bool flag_recompute_isovert;   ///< If true then recompute isovert
    bool flag_round;               ///< Round output coordinates.
    bool use_lindstrom;            ///< If true, use Lindstrom formula
    bool use_lindstrom2;
    bool use_lindstrom_fast;       /// use the garland heckbert approach
    bool use_Linf_dist;            ///< If true, use Linf dist.
	bool flag_map_extended;        ///< If true, then the extended version of mapping used.


    /// Merge sharp vertices which are within this threshold distance.
    COORD_TYPE linf_dist_thresh_merge_sharp;

    /// If true, use sharp formula to calculate intersections
    ///   of isosurface and grid edges.
    bool use_sharp_edgeI;

    /// If true, use distance to centroid of isosurface-edge intersections.
    /// Otherwise, use distance to cube center.
    bool flag_dist2centroid;       

    /// If true, check that merged vertices form a disk.
    bool flag_check_disk;

    /// If true, then check for large triangle angles.
    bool flag_check_triangle_angle;

    /// Maximum (Linf) distance from cube to isosurface vertex
    SIGNED_COORD_TYPE max_dist;

    /// Maximum (Linf) distance from cube center to isosurface vertex
    ///   used in setting other active cube.
    COORD_TYPE max_dist_to_set_other;

    /// Maximum (L2) distance from sharp vertex to sharp edge 
    ///   determined by nearby cube.
    COORD_TYPE max_dist_to_sharp_edge;

    /// Min significant distance in isosurface vertex locations.
    /// Differences in locations under this distance are insignificant.
    /// Should be some positive value.
    /// Suggested value: (0.01)*(minimum grid spacing).
    COORD_TYPE min_significant_distance;

    /// Snap points within snap distance to cube.
    COORD_TYPE snap_dist;
    
    /// Normalized eigenvalues with value less than max_small_eigenvalue
    ///   are set to zero.
    EIGENVALUE_TYPE max_small_eigenvalue;

    /// Maximum value of a small gradient coordinate.
    /// Used in determining whether to use Linf distance.
    COORD_TYPE max_small_grad_coord_Linf;

    /// Width of bin in BIN_GRID.
    int bin_width;

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
    LOC_TYPE location;                       ///< location type
    GRADIENT_COORD_TYPE ray_direction[DIM3]; ///< ray direction
    COORD_TYPE ray_initial_point[DIM3];      ///< point on ray
    COORD_TYPE ray_cube_intersection[DIM3];
    bool flag_conflict;
    VERTEX_INDEX cube_containing_coord;
    bool flag_Linf_iso_vertex_location;

    /// If true, svd point is farther than max_dist.
    bool flag_far;

    /// If true, svd point is on plane.
    bool flag_coord_on_plane;

    /// Compute distance to this point.
    COORD_TYPE central_point[DIM3];
    
    /// Set ray information.
    void SetRayInfo
    (const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
     const COORD_TYPE intersect[DIM3]);

    SVD_INFO() { Init(); };
  };
  
  
};

#endif

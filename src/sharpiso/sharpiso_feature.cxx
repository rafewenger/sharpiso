/// \file sharpiso_feature.cxx
/// Compute sharp isosurface vertices and edges.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
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


#include "sharpiso_feature.h"
#include "sharpiso_svd.h"
#include "sharpiso_findIntersect.h"
#include "sh_point_find.h"

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_scalar.txx"

// Local functions.
bool is_dist_to_cube_le
(const COORD_TYPE coord[DIM3], const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE max_dist);
VERTEX_INDEX get_cube_containing_point
(const SHARPISO_GRID & grid, const COORD_TYPE * coord);
void snap_to_cube
(const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE snap_dist, COORD_TYPE coord[DIM3]);
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const GRID_COORD_TYPE cube_coord[DIM3], 
 const COORD_TYPE point_coord[DIM3]);


// **************************************************
// SVD ROUTINES TO COMPUTE SHARP VERTEX/EDGE
// **************************************************

/// Compute sharp isosurface vertex using singular valued decomposition.
void SHARPISO::svd_compute_sharp_vertex_for_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const OFFSET_CUBE_111 & cube_111,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  if (sharpiso_param.use_lindstrom) {
    svd_compute_sharp_vertex_for_cube_lindstrom
      (scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
       cube_111, sharp_coord, eigenvalues, num_large_eigenvalues,
       svd_info);
  }
  else {

    svd_compute_sharp_vertex_for_cube_lc_intersection
      (scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
       cube_111, sharp_coord, eigenvalues, num_large_eigenvalues,
       svd_info);
  }
}


/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use Lindstrom's formula.
void SHARPISO::svd_compute_sharp_vertex_for_cube_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const OFFSET_CUBE_111 & cube_111,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  const EIGENVALUE_TYPE max_small_eigenvalue =
    sharpiso_param.max_small_eigenvalue;

  NUM_TYPE num_gradients = 0;
  std::vector<COORD_TYPE> point_coord;
  std::vector<GRADIENT_COORD_TYPE> gradient_coord;
  std::vector<SCALAR_TYPE> scalar;

  // Compute coord of the cube.
  GRID_COORD_TYPE cube_coord[DIM3];
  scalar_grid.ComputeCoord(cube_index, cube_coord);

  get_gradients
    (scalar_grid, gradient_grid, cube_index, isovalue,
     sharpiso_param, cube_111, sharpiso_param.flag_sort_gradients,
     point_coord, gradient_coord, scalar, num_gradients);

  svd_info.location = LOC_SVD;
  svd_info.flag_conflict = false;

  // svd_calculate_sharpiso vertex using lindstrom
  COORD_TYPE default_center[DIM3] = {0.5,0.5,0.5};
  COORD_TYPE cube_center[DIM3];
  IJK::add_coord_3D(cube_coord, default_center, cube_center);

  svd_calculate_sharpiso_vertex_using_lindstrom
    (&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
     num_gradients, isovalue, max_small_eigenvalue,
     num_large_eigenvalues, eigenvalues, cube_center, sharp_coord);

  if (!sharpiso_param.flag_allow_conflict && 
      sharpiso_param.flag_remove_gradients ) {
    snap_to_cube(cube_coord, sharpiso_param.snap_dist, sharp_coord);

    if (check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord)) {
      svd_info.flag_conflict = true;

      point_coord.clear();
      gradient_coord.clear();
      scalar.clear();

      // get gradients sorted by distance to cube center.
      get_gradients
        (scalar_grid, gradient_grid, cube_index, isovalue,
         sharpiso_param, cube_111, true,
         point_coord, gradient_coord, scalar, num_gradients);

      // Decrease number of gradients.
      NUM_TYPE numg2 = num_gradients - 1;
      while (numg2 > 0 && 
             check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord)) {

        svd_calculate_sharpiso_vertex_using_lindstrom
          (&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
           numg2, isovalue, max_small_eigenvalue,
           num_large_eigenvalues, eigenvalues, cube_center, sharp_coord);
        numg2--;
      }
    }
  }

  if (sharpiso_param.flag_centroid_eigen1) {
    if (num_large_eigenvalues <= 1) {
      compute_isosurface_grid_edge_centroid
        (scalar_grid, isovalue, cube_index, sharp_coord);
      svd_info.location = CENTROID;
    }
  }

  postprocess_isovert_location
    (scalar_grid, cube_index, cube_coord, isovalue, sharpiso_param, 
     sharp_coord, svd_info);
}

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
///    when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void local_svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_gradients,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info);

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
///    when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void SHARPISO::svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const OFFSET_CUBE_111 & cube_111,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  const EIGENVALUE_TYPE max_small_eigenvalue =
    sharpiso_param.max_small_eigenvalue;
  const COORD_TYPE max_dist = sharpiso_param.max_dist;

  NUM_TYPE num_gradients = 0;
  std::vector<COORD_TYPE> point_coord;
  std::vector<GRADIENT_COORD_TYPE> gradient_coord;
  std::vector<SCALAR_TYPE> scalar;
  GRADIENT_COORD_TYPE line_direction[DIM3];

  // Compute coord of the cube.
  GRID_COORD_TYPE cube_coord[DIM3];
  scalar_grid.ComputeCoord(cube_index, cube_coord);

  get_gradients
    (scalar_grid, gradient_grid, cube_index, isovalue,
     sharpiso_param, cube_111, sharpiso_param.flag_sort_gradients,
     point_coord, gradient_coord, scalar, num_gradients);

  svd_info.location = LOC_SVD;
  svd_info.flag_conflict = false;

  local_svd_compute_sharp_vertex_for_cube_lc_intersection
    (scalar_grid, gradient_grid, cube_index, cube_coord, isovalue, 
     sharpiso_param, point_coord, gradient_coord, scalar, num_gradients,
     sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);

  if (!sharpiso_param.flag_allow_conflict) {
    if (sharpiso_param.flag_remove_gradients ) {
      snap_to_cube(cube_coord, sharpiso_param.snap_dist, sharp_coord);

      if (check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord)) {
        svd_info.flag_conflict = true;

        point_coord.clear();
        gradient_coord.clear();
        scalar.clear();

        // get gradients sorted by distance to cube center.
        get_gradients
          (scalar_grid, gradient_grid, cube_index, isovalue,
           sharpiso_param, cube_111, true,
           point_coord, gradient_coord, scalar, num_gradients);

        // Decrease number of gradients.
        NUM_TYPE numg2 = num_gradients - 1;
        while (numg2 > 0 && 
               check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord) &&
               num_large_eigenvalues > 1) {

          local_svd_compute_sharp_vertex_for_cube_lc_intersection
            (scalar_grid, gradient_grid, cube_index, cube_coord, isovalue, 
             sharpiso_param, point_coord, gradient_coord, scalar, numg2,
             sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);

          numg2--;
        }
      }
    }
    else if (sharpiso_param.flag_reselect_gradients) {
      if (!sharpiso_param.use_only_cube_gradients) {
        SHARP_ISOVERT_PARAM param2 = sharpiso_param;
        param2.use_only_cube_gradients = true;
        param2.use_gradients_determining_edge_intersections = true;

        point_coord.clear();
        gradient_coord.clear();
        scalar.clear();

        // get only cube gradients
        get_gradients
          (scalar_grid, gradient_grid, cube_index, isovalue,
           param2, cube_111, sharpiso_param.flag_sort_gradients,
           point_coord, gradient_coord, scalar, num_gradients);

        local_svd_compute_sharp_vertex_for_cube_lc_intersection
          (scalar_grid, gradient_grid, cube_index, cube_coord, isovalue, 
           sharpiso_param, point_coord, gradient_coord, scalar, num_gradients,
           sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
      }
    }
  }

  postprocess_isovert_location
    (scalar_grid, cube_index, cube_coord, isovalue, sharpiso_param, 
     sharp_coord, svd_info);
}

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
///    when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void local_svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_gradients,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  const EIGENVALUE_TYPE max_eigen = sharpiso_param.max_small_eigenvalue;
  const COORD_TYPE max_dist = sharpiso_param.max_dist;
  GRADIENT_COORD_TYPE line_direction[DIM3];

  svd_calculate_sharpiso_vertex_unit_normals
    (&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]), 
     num_gradients, isovalue, max_eigen, num_large_eigenvalues, 
     eigenvalues, sharp_coord, line_direction);
    
  if (num_large_eigenvalues == 3) {

    // Check for far point or conflict and recompute using 2 eigenvalues.
    if (sharpiso_param.flag_recompute_eigen2) {
      snap_to_cube(cube_coord, sharpiso_param.snap_dist, sharp_coord);
      if (!is_dist_to_cube_le(sharp_coord, cube_coord, max_dist) ||
          check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord)) {

        svd_calculate_sharpiso_vertex_2_svals_unit_normals
          (&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]), 
           num_gradients, isovalue, max_eigen, num_large_eigenvalues,
           eigenvalues, sharp_coord, line_direction);
      }
    }
  }

  if (num_large_eigenvalues == 2) {
    bool flag_conflict;
    COORD_TYPE line_origin[DIM3];
    // if there are 2 sing vals then the coord acts as the line origin
    IJK::copy_coord_3D(sharp_coord, line_origin);

    compute_vertex_on_line
      (scalar_grid, gradient_grid, cube_coord, isovalue, sharpiso_param,
       line_origin, line_direction, sharp_coord, flag_conflict, svd_info);
    svd_info.SetRayInfo(line_origin, line_direction, sharp_coord);
  }
  else if (num_large_eigenvalues  == 1) {
    compute_isosurface_grid_edge_centroid
      (scalar_grid, isovalue, cube_index, sharp_coord);
    svd_info.location = CENTROID;
  }
  else if (num_large_eigenvalues == 0 ) {
    COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};
    IJK::add_coord_3D(cube_coord, cube_center, sharp_coord);
    svd_info.location = CUBE_CENTER;
  }

}

/// Compute sharp isosurface vertex on the line
/// @param flag_conflict True if sharp_coord conflicts with other cube.
void SHARPISO::compute_vertex_on_line
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const COORD_TYPE line_origin[DIM3],
 const COORD_TYPE line_direction[DIM3],
 COORD_TYPE sharp_coord[DIM3],
 bool & flag_conflict,
 SVD_INFO & svd_info)
{
  const SIGNED_COORD_TYPE max_dist = sharpiso_param.max_dist;
  static COORD_TYPE Linf_coord[DIM3];

  // default
  flag_conflict = false;

  // Compute the closest point (L2 distance) on the line to sharp_coord.
  compute_closest_point_to_cube_center
    (cube_coord, line_origin, line_direction, sharp_coord);

  if (is_dist_to_cube_le(sharp_coord, cube_coord, max_dist)) {

    snap_to_cube(cube_coord, sharpiso_param.snap_dist, sharp_coord);
    if (check_conflict(scalar_grid, isovalue, cube_coord, sharp_coord)) {
  
      if (sharpiso_param.use_Linf_dist) {
        // Use Linf distance instead of L1 distance.
        compute_closest_point_to_cube_center_linf
          (cube_coord, line_origin, line_direction, Linf_coord);

        snap_to_cube(cube_coord, sharpiso_param.snap_dist, Linf_coord);
        if (check_conflict(scalar_grid, isovalue, cube_coord, Linf_coord)) 
          { flag_conflict = true; }
        else {
          // No conflict.  Use Linf coord.
          IJK::copy_coord_3D(Linf_coord, sharp_coord);
        }
      }
      else {
        flag_conflict = true;
      }
    }
  }
  else {
    flag_conflict = true;
  }
}


/// Compute sharp isosurface vertex using edge-isosurface intersections.
/// Approximate gradients using linear interpolation on the edges.
void SHARPISO::svd_compute_sharp_vertex_edgeI_interpolate_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  const EIGENVALUE_TYPE max_small_eigenvalue =
    sharpiso_param.max_small_eigenvalue;
  NUM_TYPE num_gradients(0);
  std::vector<COORD_TYPE> point_coord;
  GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3];
  SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D];

  get_cube_gradients
    (scalar_grid, gradient_grid, cube_index,
     point_coord, gradient_coord, scalar);

  // Ray Direction to calculate intersection if there are 2 singular values.
  GRADIENT_COORD_TYPE ray_direction[3]={0.0};

  //tobe added as a parameters
  bool use_cmplx_interp = false;

  bool cube_create = shFindPoint
    (&(gradient_coord[0]), &(scalar[0]), isovalue, use_cmplx_interp,
     max_small_eigenvalue, eigenvalues, num_large_eigenvalues,
     svd_info, coord);

  COORD_TYPE cube_coord[DIM3];
  COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};

  scalar_grid.ComputeCoord(cube_index, cube_coord);
  //check if cube creation failed.
  if(cube_create){
    IJK::add_coord(DIM3, cube_coord, coord, coord);
    svd_info.location = LOC_SVD;
  }
  else{
    IJK::add_coord(DIM3, cube_coord, cube_center, coord);
    svd_info.location = CUBE_CENTER;
  }

  if (sharpiso_param.flag_round)
    { IJK::round_coord
        (sharpiso_param.round_denominator, DIM3, coord, coord);
    }

}

/// Compute sharp isosurface vertex using edge-isosurface intersections.
/// Select endpoint gradient which determines edge-isosurface intersection.
void SHARPISO::svd_compute_sharp_vertex_edgeI_select_gradient
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  const EIGENVALUE_TYPE max_small_eigenvalue =
    sharpiso_param.max_small_eigenvalue;
  NUM_TYPE num_gradients(0);
  std::vector<COORD_TYPE> point_coord;
  GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3];
  SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D];

  get_cube_gradients
    (scalar_grid, gradient_grid, cube_index,
     point_coord, gradient_coord, scalar);

  // Ray Direction to calculate intersection if there are 2 singular values.
  GRADIENT_COORD_TYPE ray_direction[3]={0.0};

  //tobe added as a parameters
  bool use_cmplx_interp = true;
  bool cube_create = shFindPoint
    (&(gradient_coord[0]), &(scalar[0]), isovalue, use_cmplx_interp,
     max_small_eigenvalue, eigenvalues, num_large_eigenvalues,
     svd_info, coord);


  COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};
  GRID_COORD_TYPE cube_coord[DIM3];
  scalar_grid.ComputeCoord(cube_index, cube_coord);

  //check if cube creation failed.
  if(cube_create) {
    IJK::add_coord(DIM3, cube_coord, coord, coord);
    svd_info.location = LOC_SVD;
  }
  else {
    IJK::add_coord(DIM3, cube_coord, cube_center, coord);
    svd_info.location = CUBE_CENTER;
  }

  if (is_dist_to_cube_le(coord, cube_coord, sharpiso_param.max_dist)) {

    if (!sharpiso_param.flag_allow_conflict) {

      snap_to_cube(cube_coord, sharpiso_param.snap_dist, coord);
      bool flag_conflict = 
        check_conflict(scalar_grid, isovalue, cube_coord, coord);

      if (flag_conflict) {
        process_conflict(scalar_grid, cube_index, cube_coord, isovalue,
                         sharpiso_param, coord, svd_info);
      }
    }
  }
  else {
    process_far_point(scalar_grid, cube_index, cube_coord, isovalue,
                      sharpiso_param, coord, svd_info);
  }


  // Clamp point to cube + max_dist.
  clamp_point(sharpiso_param.max_dist, cube_coord, coord);

  if (sharpiso_param.flag_round)
    { IJK::round_coord
        (sharpiso_param.round_denominator, DIM3, coord, coord);
    }
}


// **************************************************
// SUBGRID ROUTINES TO COMPUTE SHARP VERTEX/EDGE
// **************************************************


/// Compute sharp vertex.
/// Use subgrid sampling to locate isosurface vertex on sharp edge/corner.
void SHARPISO::subgrid_compute_sharp_vertex_in_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GET_GRADIENTS_PARAM & get_gradients_param,
 const OFFSET_CUBE_111 & cube_111,
 const NUM_TYPE subgrid_axis_size,
 COORD_TYPE sharp_coord[DIM3],
 SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
  NUM_TYPE num_gradients = 0;
  std::vector<COORD_TYPE> point_coord;
  std::vector<GRADIENT_COORD_TYPE> gradient_coord;
  std::vector<SCALAR_TYPE> scalar;

  get_gradients
    (scalar_grid, gradient_grid, cube_index, isovalue,
     get_gradients_param, cube_111, false,
     point_coord, gradient_coord, scalar, num_gradients);

  IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
  scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

  subgrid_calculate_iso_vertex_in_cube
    (point_coord, gradient_coord, scalar,
     num_gradients, cube_coord.PtrConst(), isovalue, subgrid_axis_size,
     sharp_coord, scalar_stdev, max_abs_scalar_error);
}

/// Calculate isosurface vertex using regular subgrid of the cube.
void SHARPISO::subgrid_calculate_iso_vertex_in_cube
(const COORD_TYPE * point_coord, const GRADIENT_COORD_TYPE * gradient_coord,
 const SCALAR_TYPE * scalar, const NUM_TYPE num_points,
 const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
 const NUM_TYPE subgrid_axis_size,
 COORD_TYPE sharp_coord[DIM3],
 SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
  COORD_TYPE coord[DIM3];
  COORD_TYPE center_coord[DIM3];
  SCALAR_TYPE sharp_stdev_squared(0);
  SCALAR_TYPE sharp_max_abs_error(0);
  COORD_TYPE sharp_dist2center_squared(0);
  IJK::PROCEDURE_ERROR error("subgrid_calculate_iso_vertex_in_cube");

  if (subgrid_axis_size < 1) {
    error.AddMessage
      ("Programming error. Subgrid axis size must be at least 1.");
    error.AddMessage("  Subgrid axis size = ", subgrid_axis_size, ".");
    throw error;
  }

  // Compute center coordinate
  for (NUM_TYPE d = 0; d < DIM3; d++)
    { center_coord[d] = cube_coord[d] + 0.5; }

  const COORD_TYPE h = 1.0/(subgrid_axis_size+1);

  bool flag_set_sharp(false);
  for (NUM_TYPE ix = 0; ix < subgrid_axis_size; ix++) {
    coord[0] = cube_coord[0] + (ix+1)*h;
    for (NUM_TYPE iy = 0; iy < subgrid_axis_size; iy++) {
      coord[1] = cube_coord[1] + (iy+1)*h;
      for (NUM_TYPE iz = 0; iz < subgrid_axis_size; iz++) {

        coord[2] = cube_coord[2] + (iz+1)*h;
        SCALAR_TYPE s, stdev_squared, max_abs_error;

        compute_gradient_based_scalar_diff
          (coord, isovalue, point_coord, gradient_coord, scalar, num_points,
           stdev_squared, max_abs_error);

        if (!flag_set_sharp ||
            stdev_squared < sharp_stdev_squared) {

          IJK::copy_coord(DIM3, coord, sharp_coord);
          sharp_stdev_squared = stdev_squared;
          sharp_max_abs_error = max_abs_error;
          IJK::compute_distance_squared
            (DIM3, coord, center_coord, sharp_dist2center_squared);
          flag_set_sharp = true;
        }
        else if (stdev_squared == sharp_stdev_squared) {
          COORD_TYPE dist2center_squared;
          IJK::compute_distance_squared
            (DIM3, coord, center_coord, dist2center_squared);
          if (dist2center_squared < sharp_dist2center_squared) {
            IJK::copy_coord(DIM3, coord, sharp_coord);
            sharp_stdev_squared = stdev_squared;
            sharp_max_abs_error = max_abs_error;
            sharp_dist2center_squared = dist2center_squared;
            flag_set_sharp = true;
          }
        }

      }
    }
  }

  scalar_stdev = std::sqrt(sharp_stdev_squared);
  max_abs_scalar_error = sharp_max_abs_error;
}


// **************************************************
// COMPUTE ISO VERTEX AT CENTROID
// **************************************************

/// Compute centroid of intersections of isosurface and grid edges
void SHARPISO::compute_isosurface_grid_edge_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  GRID_COORD_TYPE grid_coord[dimension];
  COORD_TYPE vcoord[dimension];
  COORD_TYPE coord0[dimension];
  COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];

  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
    for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
      VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
      VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

      SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
      bool is_end0_positive = true;
      if (s0 < isovalue)
        { is_end0_positive = false; };

      SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);
      bool is_end1_positive = true;
      if (s1 < isovalue)
        { is_end1_positive = false; };

      if (is_end0_positive != is_end1_positive) {

        scalar_grid.ComputeCoord(iend0, coord0);
        scalar_grid.ComputeCoord(iend1, coord1);

        IJK::linear_interpolate_coord
          (dimension, s0, coord0, s1, coord1, isovalue, coord2);

        IJK::add_coord(dimension, vcoord, coord2, vcoord);

        num_intersected_edges++;
      }
    }

  if (num_intersected_edges > 0) {
    IJK::multiply_coord
      (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
  }
  else {
    scalar_grid.ComputeCoord(iv, vcoord);
    for (int d = 0; d < dimension; d++)
      { vcoord[d] += 0.5; };
  }

  IJK::copy_coord(dimension, vcoord, coord);
}



/// Calculate isosurface vertex using regular subgrid of the cube.
void SHARPISO::subgrid_calculate_iso_vertex_in_cube
(const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_points,
 const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
 const NUM_TYPE subgrid_axis_size,
 COORD_TYPE sharp_coord[DIM3],
 SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
  subgrid_calculate_iso_vertex_in_cube
    (&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
     num_points, cube_coord, isovalue, subgrid_axis_size,
     sharp_coord, scalar_stdev, max_abs_scalar_error);
}


// **************************************************
// ROUTINES TO MOVE POINTS
// **************************************************

// Clamp to cube cube_coord[]
void SHARPISO::clamp_point
(const float cube_offset,
 const GRID_COORD_TYPE cube_coord[DIM3],
 COORD_TYPE point[DIM3])
{
  for (int i=0; i<DIM3; i++) {
    float p = point[i] - cube_coord[i];
    if (p < (-cube_offset))
      { point[i] = cube_coord[i] - cube_offset; }
    if (p > 1+cube_offset)
      { point[i]=  cube_coord[i] + 1.0 + cube_offset; }
  }
}

// Clamp to unit cube (0,0,0) to (1,1,1).
void SHARPISO::clamp_point
(const float cube_offset,
 COORD_TYPE point[DIM3])
{
  for (int i=0; i<DIM3; i++) {
    if (point[i]< (-cube_offset))
      { point[i] = -cube_offset; }
    if (point[i] > 1.0 + cube_offset)
      { point[i]=  1.0 + cube_offset; }
  }
}

void SHARPISO::postprocess_isovert_location
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 COORD_TYPE iso_coord[DIM3],
 SVD_INFO & svd_info)
{
  const COORD_TYPE max_dist = sharpiso_param.max_dist;

  if (is_dist_to_cube_le(iso_coord, cube_coord, max_dist)) {

    if (!sharpiso_param.flag_allow_conflict) {
      bool flag_conflict;
      snap_to_cube(cube_coord, sharpiso_param.snap_dist, iso_coord);
      flag_conflict = 
        check_conflict(scalar_grid, isovalue, cube_coord, iso_coord);

      if (flag_conflict) {
        process_conflict(scalar_grid, cube_index, cube_coord, isovalue,
                         sharpiso_param, iso_coord, svd_info);
        svd_info.flag_conflict = true;
      }
    }
  }
  else {
    process_far_point(scalar_grid, cube_index, cube_coord, isovalue,
                      sharpiso_param, iso_coord, svd_info);
  }

  if (sharpiso_param.flag_round)
    { IJK::round_coord
        (sharpiso_param.round_denominator, DIM3, iso_coord, iso_coord); }
}

void SHARPISO::process_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 COORD_TYPE iso_coord[DIM3],
 SVD_INFO & svd_info)
{
  if (!sharpiso_param.flag_allow_conflict) {

    if (sharpiso_param.flag_clamp_conflict) {
      // Clamp point to the cube.
      clamp_point(0, cube_coord, iso_coord);
    }
    else {
      // Use the centroid.
      compute_isosurface_grid_edge_centroid
        (scalar_grid, isovalue, cube_index, iso_coord);
      svd_info.location = CENTROID;
    }
  }

}

void SHARPISO::process_far_point
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 COORD_TYPE iso_coord[DIM3],
 SVD_INFO & svd_info)
{
  if (sharpiso_param.flag_clamp_far) {
    // Clamp point to cube + max_dist.
    clamp_point(sharpiso_param.max_dist, cube_coord, iso_coord);
  }
  else {
    // Use the centroid.
    compute_isosurface_grid_edge_centroid
      (scalar_grid, isovalue, cube_index, iso_coord);
    svd_info.location = CENTROID;
  }
}

// **************************************************
// MISC ROUTINES
// **************************************************

/// Return true if (Linf) distance from coord to cube is at most
bool is_dist_to_cube_le
(const COORD_TYPE coord[DIM3], const GRID_COORD_TYPE cube_coord[DIM3],
 const SCALAR_TYPE max_dist)
{
  for (int d=0; d<DIM3; d++) {
    if (coord[d]+max_dist < cube_coord[d]) { return(false); }
    if (coord[d] > cube_coord[d]+1+max_dist) { return(false); }
  }

  return(true);
}

/// Return index of cube containing point.
/// @pre Point is contained in grid.
/// @pre grid.AxisSize(d) > 0 for every axis d.
VERTEX_INDEX get_cube_containing_point
(const SHARPISO_GRID & grid, const COORD_TYPE * coord)
{
  COORD_TYPE coord2[DIM3];
  VERTEX_INDEX cube_index;

  for (int d = 0; d < DIM3; d++) {
    coord2[d] = floor(coord[d]);
    if (coord2[d] >= grid.AxisSize(d))
      { coord2[d] = grid.AxisSize(d)-1; }
  }

  cube_index = grid.ComputeVertexIndex(coord2);

  return(cube_index);
}

// Snap coordinate to cube.
void snap_to_cube
(const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE snap_dist, COORD_TYPE coord[DIM3])
{
  for (int d = 0; d < DIM3; d++) {
    if (coord[d] < cube_coord[d]) {
      if (coord[d] + snap_dist >= cube_coord[d]) 
        { coord[d] = cube_coord[d]; }
    }
    else {
      GRID_COORD_TYPE right_coord = cube_coord[d]+1;
      if (coord[d] > right_coord) {
        if (coord[d] <= right_coord+snap_dist) 
          { coord[d] = right_coord; }
      }
    }
  }
}


// **************************************************
// Check for conflict
// **************************************************

/// Return true if cube contains point.
/// *** MOVE TO IJKGRID. ***
bool cube_contains_point
(const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE point_coord[DIM3])
{
  for (int d = 0; d < DIM3; d++) {
    if (point_coord[d] < cube_coord[d])
      { return(false); }
    if (point_coord[d] > cube_coord[d]+1)
      { return(false); }
  }

  return(true);
}

/// Return true if point lies in an occupied cube
///   other than the one given by cube_coord[].
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const GRID_COORD_TYPE cube_coord[DIM3], 
 const COORD_TYPE point_coord[DIM3]) 
{
  if (cube_contains_point(cube_coord, point_coord)) {
    // Point is in cube.  No conflict.
    return(false);
  }

  if (scalar_grid.ContainsPoint(point_coord)) {
    // Point is inside grid.
    VERTEX_INDEX icube2 = get_cube_containing_point(scalar_grid, point_coord);

    if (IJK::is_gt_cube_min_le_cube_max(scalar_grid, icube2, isovalue)) {
      // Isosurface intersects icube2. Conflict.
      return(true);
    }
  }

  // No conflict.
  return(false);
}

// **************************************************
// SHARP_ISOVERT_PARAM
// **************************************************

/// Initialize SHARP_ISOVERT_PARAM
void SHARPISO::SHARP_ISOVERT_PARAM::Init()
{
  use_lindstrom = false;
  flag_allow_conflict = false;
  flag_clamp_conflict = true;
  flag_clamp_far = false;
  flag_recompute_eigen2 = true;
  flag_round = false;
  flag_remove_gradients = false;
  flag_reselect_gradients = false;
  flag_centroid_eigen1 = false;
  use_Linf_dist = true;
  max_dist = 1.0;
  snap_dist = 1.0/16.0;
  max_small_eigenvalue = 0.1;
  separation_distance = 0.1;
  round_denominator = 16;
}

// **************************************************
// Class SHARP_ISOVERT_PARAM member functions
// **************************************************

/// Set
void SHARP_ISOVERT_PARAM::Set(const SHARP_ISOVERT_PARAM & param)
{
  *this = param;
}


// **************************************************
// SVD_INFO
// **************************************************

/// Initialize.
void SVD_INFO::Init()
{
  flag_conflict = false;
  location = LOC_SVD;
}

/// Set ray information.
void SVD_INFO::SetRayInfo
(const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
 const COORD_TYPE intersect[DIM3])

{
  IJK::copy_coord_3D(origin, ray_initial_point);
  IJK::copy_coord_3D(direction, ray_direction);
  IJK::copy_coord_3D(intersect, ray_cube_intersection);
}

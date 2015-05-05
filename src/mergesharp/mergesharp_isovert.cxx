// \file mergesharp_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.

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

#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iomanip>  
#include <string>  
#include <stdio.h>

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkgrid_macros.h"
#include "ijkscalar_grid.txx"


#include "sharpiso_array.txx"
#include "mergesharp_types.h"
#include "mergesharp_isovert.h"
#include "mergesharp_position.h"

#include "sharpiso_closest.h"

#include "mergesharp_debug.h"

using namespace std;
using namespace SHARPISO;
using namespace MERGESHARP;
using namespace IJK;


// *** DEBUG ***
#include "ijkprint.txx"
bool flag_debug(false);


namespace {

  void compute_Linf_distance_from_cube_center
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE isovert_coord[DIM3], COORD_TYPE & linf_dist);

  void compute_Linf_distance_from_cube_center
  (const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert);

  void compute_Linf_distance_from_grid_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE coord[DIM3], COORD_TYPE & linf_dist, int & axis);

  void compute_Linf_distance_from_grid_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE coord[DIM3], COORD_TYPE & linf_dist);

  void compute_rescaled_Linf_distance
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE spacing[DIM3], COORD_TYPE & linf_dist, int & axis);

  void compute_rescaled_Linf_distance
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE spacing[DIM3], COORD_TYPE & linf_dist);

  void compute_rescaled_distance_to_line_3D
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE direction[DIM3], const COORD_TYPE spacing[DIM3], 
   COORD_TYPE & distance);

  bool are_connected
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX & cube_index1,
   const VERTEX_INDEX & cube_index2, const SCALAR_TYPE isovalue);

  bool is_angle_large
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const ISOVERT & isovertData, const VERTEX_INDEX iv,
   const SCALAR_TYPE threshold, const VERTEX_INDEX v1,
   const VERTEX_INDEX v2);

  bool is_facet_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, VERTEX_INDEX & cube_index1);

  bool is_facet_or_edge_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   VERTEX_INDEX & cube_index2);

  bool is_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   const VERTEX_INDEX & cube_index2);

  bool is_cube_or_neighbor_covered_by
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1);

  void store_svd_info
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const NUM_TYPE gcube_index,
   const NUM_TYPE num_large_eigenvalues,
   const SVD_INFO svd_info,
   ISOVERT & isovert);

  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, ISOVERT & isovert);

  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, const NUM_TYPE gcube_index,
   ISOVERT & isovert);

  void snap2cube
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index, 
   const COORD_TYPE pointA[DIM3], COORD_TYPE * pointB);

  void snap2cube_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index, 
   const COORD_TYPE pointA[DIM3], COORD_TYPE * pointB);

  // *** MOVED TO mergesharp_select ***
  bool check_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   const ISOVERT & isovert,
   const VERTEX_INDEX & gcube0_index);

  // *** MOVED TO mergesharp_select ***
  void check_and_set_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube_index);

  void check_covered_and_substitute
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SHARPISO_BOOL_GRID & covered_grid,
   const SCALAR_TYPE isovalue,
   const NUM_TYPE gcube_index,
   ISOVERT & isovert);

  void check_not_contained_and_substitute
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const NUM_TYPE gcube_index,
   ISOVERT & isovert);

  void set_covered_grid
  (const ISOVERT & isovert, SHARPISO_BOOL_GRID & covered_grid);
}

void replace_with_substitute_coord
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert);

void set_isovert_position
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, 
 const COORD_TYPE * new_coord, const NUM_TYPE num_eigenvalues,
 ISOVERT & isovert);

void set_isovert_position_from_other_cube
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, 
 const COORD_TYPE * new_coord, const NUM_TYPE num_eigenvalues,
 ISOVERT & isovert);


// **************************************************
// ISOVERT ROUTINES
// **************************************************

/// Compute isosurface vertex positions from gradients 
///   using lindstrom algorithm.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
void compute_isovert_position_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const VERTEX_INDEX & cube_index,
 ISOVERT & isovert)
{
  const INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);
  const GRADIENT_COORD_TYPE max_small_magnitude =
    isovert_param.max_small_magnitude;
	SVD_INFO svd_info;

  if (gcube_index != ISOVERT::NO_INDEX) {

    // this is an active cube
    isovert.gcube_list[gcube_index].flag_centroid_location = false;

    // compute the sharp vertex for this cube
    EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
    COORD_TYPE edge_dir[DIM3], orth_dir[DIM3];
    NUM_TYPE num_large_eigenvalues;

    svd_compute_sharp_vertex_for_cube_lindstrom
      (scalar_grid, gradient_grid, cube_index, isovalue, isovert_param, 
       voxel, isovert.gcube_list[gcube_index].isovert_coord,
       edge_dir, orth_dir, eigenvalues, num_large_eigenvalues, svd_info);

    if (num_large_eigenvalues == 1) {
      IJK::copy_coord_3D(orth_dir, isovert.gcube_list[gcube_index].direction);
    }
    else if (num_large_eigenvalues == 2) {
      IJK::copy_coord_3D(edge_dir, isovert.gcube_list[gcube_index].direction);
    }
    else {
      IJK::set_coord_3D(0, isovert.gcube_list[gcube_index].direction);
    }

    const COORD_TYPE * isovert_coord =
      isovert.gcube_list[gcube_index].isovert_coord;
    const COORD_TYPE * direction =
      isovert.gcube_list[gcube_index].direction;

    if (num_large_eigenvalues > 2 ||
        cube_contains_point(scalar_grid, cube_index, isovert_coord)) {

      IJK::copy_coord_3D
        (isovert_coord, isovert.gcube_list[gcube_index].isovert_coordB);
    }
    else if (num_large_eigenvalues == 2) {
      COORD_TYPE cube_center_coord[DIM3];
      COORD_TYPE closest_point[DIM3];

      scalar_grid.ComputeCubeCenterScaledCoord
        (cube_index, cube_center_coord);

      // Compute point on line closest (L2) to cube center.
      compute_closest_point_on_line
        (cube_center_coord, isovert_coord, direction, max_small_magnitude, 
         closest_point);

      if (cube_contains_point(scalar_grid, cube_index, isovert_coord)) {
        IJK::copy_coord_3D
          (closest_point, isovert.gcube_list[gcube_index].isovert_coordB);
      }
      else {

        // Compute point on line closest (Linf) to cube center.
        compute_closest_point_on_line_unscaled_linf
          (cube_center_coord, isovert_coord,  direction, max_small_magnitude,
           scalar_grid.SpacingPtrConst(),
           isovert.gcube_list[gcube_index].isovert_coordB);
      }
    }
    else if (num_large_eigenvalues == 1) {
      COORD_TYPE cube_center_coord[DIM3];

      scalar_grid.ComputeCubeCenterScaledCoord
        (cube_index, cube_center_coord);

      // Compute point on plane closest (L2) to cube center.
      compute_closest_point_on_plane
        (cube_center_coord, isovert_coord, direction, 
         isovert.gcube_list[gcube_index].isovert_coordB);
    }

    store_svd_info(scalar_grid, cube_index, gcube_index, 
                   num_large_eigenvalues, svd_info, isovert);
  }
}

/// Compute isosurface vertex positions from gradients.
void compute_all_isovert_positions 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;
  SVD_INFO svd_info;

  voxel.SetVertexCoord
    (scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);

  IJK_FOR_EACH_GRID_CUBE(cube_index, scalar_grid, VERTEX_INDEX) {

    compute_isovert_position_lindstrom
      (scalar_grid, gradient_grid, isovalue, isovert_param, voxel,
       cube_index, isovert);
  }

  store_boundary_bits(scalar_grid, isovert.gcube_list);
  set_cube_containing_isovert(scalar_grid, isovalue, isovert);
}

// Compute isosurface vertex positions using isosurface-edge intersections.
void compute_all_isovert_positions_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_POSITION_METHOD vertex_position_method,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  SVD_INFO svd_info;

  if (vertex_position_method == EDGEI_GRADIENT) {

    IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
      NUM_TYPE index = isovert.GCubeIndex(iv);

      if (index != ISOVERT::NO_INDEX) {
        // this is an active cube
        isovert.gcube_list[index].flag_centroid_location = false;

        // compute the sharp vertex for this cube
        EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
        NUM_TYPE num_large_eigenvalues;

        svd_compute_sharp_vertex_edgeI_sharp_gradient
          (scalar_grid, gradient_grid, iv, isovalue, isovert_param,
           isovert.gcube_list[index].isovert_coord,
           eigenvalues, num_large_eigenvalues, svd_info);


        store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
                       svd_info, isovert);
      }
    }
  }
  else {
    // vertex_position_method == EDGEI_INTERPOLATE

    IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
      NUM_TYPE index = isovert.GCubeIndex(iv);

      if (index != ISOVERT::NO_INDEX) {
        // this is an active cube
        isovert.gcube_list[index].flag_centroid_location = false;

        // compute the sharp vertex for this cube
        EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
        NUM_TYPE num_large_eigenvalues;

        svd_compute_sharp_vertex_edgeI_interpolate_gradients
          (scalar_grid, gradient_grid, iv, isovalue, isovert_param,
           isovert.gcube_list[index].isovert_coord,
           eigenvalues, num_large_eigenvalues, svd_info);

        store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
                       svd_info, isovert);

      }
    }
  }

  store_boundary_bits(scalar_grid, isovert.gcube_list);
}


// **************************************************
// RECOMPUTE ROUTINES
// **************************************************

/// Recompute isovert position using lindstrom.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 const NUM_TYPE gcube_index,
 ISOVERT & isovert)
{
  VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  MSDEBUG();
  if (flag_debug) {
    using namespace std;
    cerr << "  Recomputing isovert coord for: " << cube_index << " ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << endl;
    cerr << "    Old coord: ";
    IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coord);
    cerr << " num eigenvalues: " 
         << int(isovert.gcube_list[gcube_index].num_eigenvalues);
    cerr << "  flag: "
         << int(isovert.gcube_list[gcube_index].flag);
    cerr << endl;
    if (isovert.gcube_list[gcube_index].flag_conflict) {
      scalar_grid.PrintIndexAndCoord
        (cerr, "    Isovert coord conflicts with cube ",
         isovert.gcube_list[gcube_index].cube_containing_isovert, "\n");
    }
  }

  compute_isovert_position_lindstrom
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel,
     cube_index, isovert);

  isovert.gcube_list[gcube_index].flag_recomputed_coord = true;
  isovert.gcube_list[gcube_index].flag_recomputed_coord_min_offset = 
    flag_min_offset;

  set_cube_containing_isovert(scalar_grid, isovalue, gcube_index, isovert);

  MSDEBUG();
  if (flag_debug) {
    using namespace std;
    cerr << "    New coord: ";
    IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coord);
    cerr << "  B: ";
    IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coordB);
    cerr << " num eigenvalues: " 
         << int(isovert.gcube_list[gcube_index].num_eigenvalues);
    cerr << "  flag: "
         << int(isovert.gcube_list[gcube_index].flag);
    cerr << endl;
  }
}

/// Recompute isovert positions for cubes containing far points.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_far_points
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

    GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;
    bool flag_recomputed_coord_min_offset =
      isovert.gcube_list[i].flag_recomputed_coord_min_offset;
    if (isovert.gcube_list[i].flag_far) {
      if (!flag_recomputed_coord_min_offset) {

        recompute_isovert_position_lindstrom
          (scalar_grid, gradient_grid, isovalue,
           isovert_param, voxel, flag_min_offset, i, isovert);
      }
    }
  }

}

/// Recompute isovert positions for cubes with far points.
void recompute_far_points
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  if (grad_selection_cube_offset > 0.5) {

    if (isovert_param.min_grad_selection_cube_offset <= 0.5) {

      voxel.SetVertexCoord(scalar_grid.SpacingPtrConst(), 0.5);

      recompute_far_points
        (scalar_grid, gradient_grid, isovalue, isovert_param, 
         voxel, false, isovert);
    }
  }

  if (isovert_param.min_grad_selection_cube_offset <= 0) {

    voxel.SetVertexCoord(scalar_grid.SpacingPtrConst(), 0.0);

    recompute_far_points
      (scalar_grid, gradient_grid, isovalue, isovert_param, 
       voxel, true, isovert);
  }
}


/// Recompute isovert positions for cubes containing covered points.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void MERGESHARP::recompute_covered_point_positions
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

    GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;
    bool flag_recomputed_coord_min_offset =
      isovert.gcube_list[i].flag_recomputed_coord_min_offset;
    if (cube_flag == COVERED_POINT) {
      if (!flag_recomputed_coord_min_offset) {

        recompute_isovert_position_lindstrom
          (scalar_grid, gradient_grid, isovalue,
           isovert_param, voxel, flag_min_offset, i, isovert);

        check_covered_and_substitute
          (scalar_grid, covered_grid, isovalue, i, isovert);
      }
    }
  }

}


/// Recompute isovert positions for cubes containing covered points.
void MERGESHARP::recompute_covered_point_positions
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  if (grad_selection_cube_offset > 0.5) {
    if (isovert_param.min_grad_selection_cube_offset <= 0.5) {

      voxel.SetVertexCoord
        (scalar_grid.SpacingPtrConst(), 0.5);

      recompute_covered_point_positions
        (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param, 
         voxel, false, isovert);
    }
  }

  if (isovert_param.min_grad_selection_cube_offset <= 0) {

    voxel.SetVertexCoord(scalar_grid.SpacingPtrConst(), 0.0);

    recompute_covered_point_positions
      (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param, 
       voxel, true, isovert);
  }
}


/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered using lindstrom
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
void recompute_unselected_uncovered_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

    if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE) ||
        (cube_flag == COVERED_POINT)) 
      {
        const VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
        const COORD_TYPE * isovert_coord = isovert.gcube_list[i].isovert_coord;

        if (!cube_contains_point(scalar_grid, cube_index, isovert_coord)) {

          recompute_isovert_position_lindstrom
            (scalar_grid, gradient_grid, isovalue, isovert_param,
             voxel, flag_min_offset, i, isovert);

          check_not_contained_and_substitute(scalar_grid, i, isovert);

          if (check_covered_point(covered_grid, isovert, i)) 
            { isovert.gcube_list[i].flag = COVERED_POINT; }
        }
      }
  }
}

/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered using lindstrom
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
void recompute_unselected_uncovered_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  if (grad_selection_cube_offset > 0.5) {

    voxel.SetVertexCoord(scalar_grid.SpacingPtrConst(), 0.5);

    recompute_unselected_uncovered_lindstrom
      (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param, 
       voxel, false, isovert);
  }

  voxel.SetVertexCoord
    (scalar_grid.SpacingPtrConst(), 0.0);

  recompute_unselected_uncovered_lindstrom
    (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param, 
     voxel, true, isovert);
}

/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered.
void MERGESHARP::recompute_isovert_positions 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;

  // covered_grid.Scalar(iv) == true, if iv is covered (or selected).
  SHARPISO_BOOL_GRID covered_grid;
  covered_grid.SetSize(scalar_grid);
  set_covered_grid(isovert, covered_grid);

  if (isovert_param.use_lindstrom) {

    recompute_unselected_uncovered_lindstrom
      (scalar_grid, gradient_grid, covered_grid, isovalue, isovert_param, 
       isovert);
  }

  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

    if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE) ||
        (cube_flag == COVERED_POINT)) 
      {
        const COORD_TYPE * isovert_coord = isovert.gcube_list[i].isovert_coord;
        const VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

        if (!cube_contains_point(scalar_grid, cube_index, isovert_coord)) {

          const VERTEX_INDEX cube2_index = 
            isovert.gcube_list[i].cube_containing_isovert;

          if (isovert.GCubeIndex(cube2_index) != ISOVERT::NO_INDEX ||
              isovert.gcube_list[i].linf_dist > linf_dist_threshold) {

            compute_edgeI_centroid
              (scalar_grid, gradient_grid, isovalue, cube_index,
               isovert_param.use_sharp_edgeI, 
               isovert.gcube_list[i].isovert_coord);
            isovert.gcube_list[i].flag_centroid_location = true;

            MSDEBUG();
            if (flag_debug) {
              cerr << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
              scalar_grid.PrintIndexAndCoord
                (cerr, "*** Replace by centroid in cube ", cube_index, "\n");
              scalar_grid.PrintIndexAndCoord
                (cerr, "    Original vertex is in cube ", cube2_index, "\n");
            }
          }
        }
      }
  }

}


/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered.
/// This is not called by the sharp computations
void MERGESHARP::recompute_isovert_positions 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  SHARPISO_EDGE_INDEX_GRID edge_index(DIM3, scalar_grid.AxisSize(), DIM3);

  set_edge_index(edgeI_coord, edge_index);

  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

    if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE)) {

      VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
      compute_edgeI_centroid
        (scalar_grid, edgeI_coord, edge_index, isovalue, cube_index,
         isovert.gcube_list[i].isovert_coord);

      isovert.gcube_list[i].flag_centroid_location = true;
    }
  }

}


inline bool is_uncovered_sharp(const GRID_CUBE_DATA & gcube_data)
{
  GRID_CUBE_FLAG gcube_flag = gcube_data.flag;
  if (gcube_flag == COVERED_POINT || gcube_flag == UNAVAILABLE_GCUBE) {

    if (gcube_data.cube_index == gcube_data.maps_to_cube)
      { return(true); }
  }

  return(false);
}

void MERGESHARP::recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

    if (is_uncovered_sharp(isovert.gcube_list[i])) {
      const VERTEX_INDEX cube_index = isovert.CubeIndex(i);
      bool flag_loc_recomputed;

      recompute_using_adjacent
        (scalar_grid, cube_index, isovert, flag_loc_recomputed);
    }
  }

}

void MERGESHARP::recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const VERTEX_INDEX cube_index, ISOVERT & isovert,
 bool & flag_loc_recomputed)
{
  const NUM_TYPE gcube_index = isovert.GCubeIndex(cube_index);
  const BOUNDARY_BITS_TYPE boundary_bits = 
    isovert.gcube_list[gcube_index].boundary_bits;
  FIXED_ARRAY<NUM_CUBE_NEIGHBORS3D, VERTEX_INDEX, NUM_TYPE>
    neighbor_list;
  COORD_TYPE new_coord[DIM3];


  MSDEBUG();
  if (flag_debug) {
    cerr << "In " << __func__;
    scalar_grid.PrintIndexAndCoord(cerr, "  cube: ", cube_index, "\n");
  }

  flag_loc_recomputed = false;
  if (boundary_bits == 0) {

    for (int d = 0; d < DIM3; d++) {

      for (int iadj = 0; iadj < 2; iadj++) {

        VERTEX_INDEX cube2_index = scalar_grid.AdjacentVertex(cube_index, d, iadj);
        NUM_TYPE gcube2_index = isovert.GCubeIndex(cube2_index);

        if (gcube2_index == ISOVERT::NO_INDEX) { continue; }

        MSDEBUG();
        if (flag_debug) {
          scalar_grid.PrintIndexAndCoord(cerr, "Processing active neighbor ", cube2_index, "\n");
        }

        if (is_uncovered_sharp(isovert.gcube_list[gcube2_index]) == false) {

          // SHOULD CHECK THAT SHARED FACET IS ACTIVE

          VERTEX_INDEX cube3_index = 
            isovert.gcube_list[gcube2_index].maps_to_cube;

          if (neighbor_list.Contains(cube3_index) == false) 
            { neighbor_list.PushBack(cube3_index); }
        }
      }
    }

    MSDEBUG();
    if (flag_debug) {
      cerr << "Neighbor list: ";
      for (int i = 0; i < neighbor_list.NumElements(); i++)
        { cerr << " " << neighbor_list[i]; }
      cerr << endl;
    }

    if (neighbor_list.NumElements() < 3) { return; }

    IJK::set_coord_3D(0, new_coord);
    for (NUM_TYPE i = 0; i < neighbor_list.NumElements(); i++) {
      VERTEX_INDEX cube3_index = neighbor_list[i];
      NUM_TYPE gcube3_index = isovert.GCubeIndex(cube3_index);
      IJK::add_coord_3D(new_coord, isovert.IsoVertCoord(gcube3_index), new_coord);
    }

    IJK::divide_coord
      (DIM3, COORD_TYPE(neighbor_list.NumElements()), new_coord, new_coord);


    MSDEBUG();
    if (flag_debug) {
      scalar_grid.PrintIndexAndCoord(cerr, "*** Recompute: ", cube_index, "\n");
      cerr << "  Old coord: ";
      IJK::print_coord3D(cerr, isovert.IsoVertCoord(gcube_index));
      cerr << "  New: ";
      IJK::print_coord3D(cerr, new_coord);
      cerr << endl;
    }

    IJK::copy_coord_3D(new_coord, isovert.gcube_list[gcube_index].isovert_coord);
    isovert.gcube_list[gcube_index].flag_recomputed_using_adjacent = true;
    isovert.gcube_list[gcube_index].flag_centroid_location = false;
  }
  else {
    // HANDLE BOUNDARY_CUBE
  }

}


/// Set isovert position from grid vertex or grid edge.
void set_isovert_position_from_face
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const COORD_TYPE isovert_coord[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const bool flag_from_vertex,
 ISOVERT & isovert,
 bool & flag_set)
{
  const INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);

  flag_set = false;

  if (gcube_index == ISOVERT::NO_INDEX) { return; }

  MSDEBUG();
  if (flag_debug) {
    cerr << "*** In " << __FUNCTION__ << ".  Cube: " << cube_index << " ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << endl;
    cerr << "    Old coord: ";
    IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coord);
    cerr << "  num eigen: " << isovert.NumEigenvalues(gcube_index) << endl;
    cerr << "    New coord: ";
    IJK::print_coord3D(cerr, isovert_coord);
    cerr << " num eigen: " << num_large_eigenvalues;
    if (flag_from_vertex) 
	  { cerr << "  From vertex."; }
    else
	  { cerr << "  From edge."; }
    cerr << endl;
  }

  set_isovert_position
    (scalar_grid, gcube_index, isovert_coord, num_large_eigenvalues, isovert);

  set_cube_containing_isovert(scalar_grid, isovalue, gcube_index, isovert);

  isovert.gcube_list[gcube_index].flag_using_substitute_coord = false;

  if (flag_from_vertex) 
    { isovert.gcube_list[gcube_index].flag_coord_from_vertex = true; }
  else
    { isovert.gcube_list[gcube_index].flag_coord_from_edge = true; }

  flag_set = true;
}

/// Set isovert position from grid vertex or grid edge.
void set_isovert_position_from_face
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const COORD_TYPE isovert_coord[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const bool flag_from_vertex,
 ISOVERT & isovert)
{
  bool flag_set;
  set_isovert_position_from_face
    (scalar_grid, isovalue, cube_index, isovert_coord, num_large_eigenvalues, 
     flag_from_vertex, isovert, flag_set);
}


/// Return true if sharp edge in cube0 points to cube1.
/// Return false if cube0 has no sharp edge.
bool does_sharp_edge_point_to_cube
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 const VERTEX_INDEX cube0_index, const VERTEX_INDEX cube1_index)
{
  NUM_TYPE gcube0_index = isovert.GCubeIndex(cube0_index);
  COORD_TYPE cube1_center_coord[DIM3];
  COORD_TYPE distance;

  if (isovert.NumEigenvalues(gcube0_index) != 2) 
    { return(false); }

  // Check that sharp edge points to the other cube.
  const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube0_index);
  const COORD_TYPE * edge_dir =
    isovert.gcube_list[gcube0_index].direction;

  grid.ComputeCubeCenterScaledCoord(cube1_index, cube1_center_coord);

  compute_rescaled_distance_to_line_3D
    (cube1_center_coord, isovert_coord, edge_dir, grid.SpacingPtrConst(), 
     distance);

  if (distance > 0.5) { 
    // Sharp edge does not intersect cube 
    //   (or only intersects cube near its boundary)
    return(false);
  }

  return(true);
}

/// Compute isovert coordinate on plane from sharp edges.
/// If no sharp edges or intersections are not near pointX,
///   intersect plane and line containing isovertex coordinates.
/// If isovertices are too close, compute point on plane closest
///   to midpoint of line segment between isovert coordinates.
void compute_isovert_on_plane_from_sharp_edges
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX cube_index[2],
 const COORD_TYPE pointX[DIM3], const COORD_TYPE plane_normal[DIM3],
 COORD_TYPE new_isovert_coord[DIM3])
{
  const ANGLE_TYPE max_angle_LO = 80;
  const COORD_TYPE min_abs_cos_angle_LO = cos(max_angle_LO*M_PI/180.0);
  const COORD_TYPE min_significant_distance =
    isovert_param.min_significant_distance;
  COORD_TYPE intersection_point[DIM3];

  IJK::set_coord_3D(0, new_isovert_coord);
  NUM_TYPE num_added = 0;
  for (int i = 0; i < 2; i++) {
    NUM_TYPE gcube_index = isovert.GCubeIndex(cube_index[i]);

    if (isovert.NumEigenvalues(gcube_index) == 2) {

      const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube_index);
      const COORD_TYPE * edge_dir =
        isovert.gcube_list[gcube_index].direction;

      bool flag_succeeded;
      IJK::intersect_line_plane_3D
        (isovert_coord, edge_dir, pointX, plane_normal, 
         min_abs_cos_angle_LO, intersection_point, flag_succeeded);

      if (flag_succeeded) {
        COORD_TYPE distance;

        compute_rescaled_Linf_distance
          (pointX, intersection_point, grid.SpacingPtrConst(), distance);

        if (distance < 1) {

          IJK::add_coord_3D
            (new_isovert_coord, intersection_point, new_isovert_coord);
          num_added++;
        }
      }
    }
  }

  if (num_added > 0) {
    IJK::divide_coord(DIM3, COORD_TYPE(num_added), new_isovert_coord,
                      new_isovert_coord);
  }
  else {
    // Intersect the plane and the line through the isovert coordinates
    //   in cubes fixed_cube[0] and fixed_cube[1].
    const NUM_TYPE gcube0_index = isovert.GCubeIndex(cube_index[0]);
    const NUM_TYPE gcube1_index = isovert.GCubeIndex(cube_index[1]);
    const COORD_TYPE * isovert_coord0 = isovert.IsoVertCoord(gcube0_index);
    const COORD_TYPE * isovert_coord1 = isovert.IsoVertCoord(gcube1_index);
    COORD_TYPE line_dir[DIM3], magnitude;

    IJK::subtract_coord_3D(isovert_coord1, isovert_coord0, line_dir);
    bool flag_zero;
    IJK::normalize_vector(DIM3, line_dir, min_significant_distance, 
                          line_dir, magnitude, flag_zero);

    if (flag_zero) {
      // Two vertices are close together.
      COORD_TYPE midpoint[DIM3];
      IJK::compute_midpoint(DIM3, isovert_coord0, isovert_coord1, midpoint);

      // Use point on plane closest to midpoint
      compute_closest_point_on_plane
        (midpoint, pointX, plane_normal, new_isovert_coord);
    }
    else {
      bool flag_succeeded;
      IJK::intersect_line_plane_3D
        (isovert_coord0, line_dir, pointX, plane_normal, 
         min_abs_cos_angle_LO, new_isovert_coord, flag_succeeded);

      if (!flag_succeeded) {
        // Just in case...
        COORD_TYPE midpoint[DIM3];
        IJK::compute_midpoint(DIM3, isovert_coord0, isovert_coord1, midpoint);

        // Use point on plane closest to midpoint
        compute_closest_point_on_plane
          (midpoint, pointX, plane_normal, new_isovert_coord);
      }
    }
  }

}

/// Compute isovert coordinate on plane from sharp edges.
/// If no sharp edges or intersections are not near pointX,
///   compute point on plane nearest centroid of isovert coord.
void compute_isovert_on_plane_from_sharp_edges_B
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_INDEX iend0, const int edge_dir,
 const COORD_TYPE pointX[DIM3], const COORD_TYPE plane_normal[DIM3],
 COORD_TYPE new_isovert_coord[DIM3],
 bool & flag_succeeded)
{
  const VERTEX_INDEX * axis_increment = grid.AxisIncrement();
  const ANGLE_TYPE max_angle_LO = 80;
  const COORD_TYPE min_abs_cos_angle_LO = cos(max_angle_LO*M_PI/180.0);
  const COORD_TYPE min_significant_distance =
    isovert_param.min_significant_distance;
  const int d1 = (edge_dir+1)%DIM3;
  const int d2 = (edge_dir+2)%DIM3;
  COORD_TYPE intersection_point[DIM3];
  COORD_TYPE centroid_coord[DIM3];


  flag_succeeded = false;

  IJK::set_coord_3D(0, new_isovert_coord);
  NUM_TYPE num_added = 0;
  for (int j1 = 0; j1 < 2; j1++) {
    for (int j2 = 0; j2 < 2; j2++) {
      VERTEX_INDEX cube_index = 
        iend0 - j1*axis_increment[d1] - j2*axis_increment[d2];
      INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);

      if (gcube_index == ISOVERT::NO_INDEX) { continue; }
      if (isovert.NumEigenvalues(gcube_index) != 2) { continue; }

      const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube_index);
      const COORD_TYPE * edge_dir = isovert.gcube_list[gcube_index].direction;

      bool flag_succeeded_B;
      IJK::intersect_line_plane_3D
        (isovert_coord, edge_dir, pointX, plane_normal, 
         min_abs_cos_angle_LO, intersection_point, flag_succeeded_B);

      MSDEBUG();
      if (flag_debug) {
        cerr << endl;
        grid.PrintIndexAndCoord
          (cerr, "     Intersecting cube: ", cube_index, " with plane.\n");
        cerr << "     Point: ";
        print_coord3D(cerr, isovert_coord);
        cerr << "     Edge dir: ";
        print_coord3D(cerr, edge_dir);
        cerr << endl;
        cerr << "     Point on plane: ";
        print_coord3D(cerr, pointX);
        cerr << "     Plane normal: ";
        print_coord3D(cerr, plane_normal);
        cerr << endl;
      }

      if (flag_succeeded_B) {
        COORD_TYPE distance;

        compute_rescaled_Linf_distance
          (pointX, intersection_point, grid.SpacingPtrConst(), distance);

        MSDEBUG();
        if (flag_debug) {
          cerr << "     Intersection point: ";
          print_coord3D(cerr, intersection_point);
          cerr << " distance: " << distance << endl;
        }

        if (distance < 1) {

          IJK::add_coord_3D
            (new_isovert_coord, intersection_point, new_isovert_coord);
          num_added++;
        }
      }
    }
  }

  if (num_added > 0) {
    IJK::divide_coord(DIM3, COORD_TYPE(num_added), new_isovert_coord,
                      new_isovert_coord);

    flag_succeeded = true;
	MSDEBUG();
    if (flag_debug) {
      cerr << "   *** Computed intersection of plane and edge dir." << endl;
    }
  }
  else {

    if (flag_debug) {
      cerr << "  --- Intersection point NOT found." << endl;
    }

    flag_succeeded = false;
  }

}

/// Recompute isovert position around vertex.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_vertex
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 const VERTEX_INDEX iv,
 ISOVERT & isovert)
{
  const NUM_TYPE index_lastv = NUM_CUBE_VERTICES3D - 1;
  const COORD_TYPE ONE_DIV_SQRT3 = 1.0/std::sqrt(3);
  const VERTEX_INDEX * axis_increment = scalar_grid.AxisIncrement();
  const COORD_TYPE * spacing = scalar_grid.SpacingPtrConst();
  const ANGLE_TYPE max_angle_LO = 80;
  const COORD_TYPE min_abs_cos_angle_LO = cos(max_angle_LO*M_PI/180.0);
  const COORD_TYPE min_significant_distance =
    isovert_param.min_significant_distance;
  BOUNDARY_BITS_TYPE boundary_bits;
  VERTEX_INDEX fixed_cube[2];
  COORD_TYPE pointX[DIM3];
  GRID_COORD_TYPE cube_coord[2][DIM3];
  COORD_TYPE plane_normal[DIM3];
  COORD_TYPE new_isovert_coord[DIM3];
  COORD_TYPE intersection_point[DIM3];
  SVD_INFO svd_info;

  // Note: Compute vertex (not cube) boundary bits.
  scalar_grid.ComputeBoundaryBits(iv, boundary_bits);

  if (boundary_bits != 0) { return; }

  NUM_TYPE num_active_cubes = 0;

  // Number of cubes containing their isovertCoord.
  NUM_TYPE num_fixed_cubes = 0;

  // Number of fixed cubes with exactly two large eigenvalues.
  NUM_TYPE num_fixed_cubes_with_eigenvalue_two = 0;

  VERTEX_INDEX cube0_index = 
    iv - scalar_grid.CubeVertexIncrement(index_lastv);

  for (NUM_TYPE i = 0; i < NUM_CUBE_VERTICES3D; i++) {
    VERTEX_INDEX cube1_index = scalar_grid.CubeVertex(cube0_index, i);
    INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);

    if (gcube1_index == ISOVERT::NO_INDEX) { continue; }

    num_active_cubes++;

    const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube1_index);
    if (cube_contains_point(scalar_grid, cube1_index, isovert_coord)) {

      NUM_TYPE num_eigenvalues = isovert.NumEigenvalues(gcube1_index);

      if (num_eigenvalues < 2) { 
        // This cube will connect sharp cubes, if they exist.
        return; 
      }
      else if (num_eigenvalues == 2) {
        num_fixed_cubes_with_eigenvalue_two++;
      }

      if (num_fixed_cubes >= 2) { return; }

      fixed_cube[num_fixed_cubes] = cube1_index;
      num_fixed_cubes++;
    }
  }

  if (num_fixed_cubes != 2) { return; }
  if (num_active_cubes <= 2) { return; }

  for (int i = 0; i < 2; i++)
    { scalar_grid.ComputeCoord(fixed_cube[i], cube_coord[i]); }

  for (int d = 0; d < DIM3; d++) {

    if (cube_coord[0][d] == cube_coord[1][d]) { return; }
    else if (cube_coord[0][d] < cube_coord[1][d])
      { plane_normal[d] = ONE_DIV_SQRT3; }
    else
      { plane_normal[d] = -ONE_DIV_SQRT3; }
  }

  for (int i = 0; i < 2; i++) {
    NUM_TYPE gcube_index = isovert.GCubeIndex(fixed_cube[i]);

    if (isovert.NumEigenvalues(gcube_index) == 2) {

      if (!does_sharp_edge_point_to_cube
          (scalar_grid, isovert, fixed_cube[i], fixed_cube[1-i])) {
        // Sharp edge in fixed_cube[i] does not point to fixed_cube[1-i].
        return;
      }
    }
  }

  scalar_grid.ComputeScaledCoord(iv, pointX);

  compute_isovert_on_plane_from_sharp_edges
    (scalar_grid, isovert, isovert_param, fixed_cube, pointX, plane_normal,
     new_isovert_coord);

  MSDEBUG();
  if (flag_debug) {
    cerr << "--- Function: " << __func__ << endl;
    scalar_grid.PrintIndexAndCoord(cerr, "  Vertex: ", iv, "\n");
    scalar_grid.PrintIndexAndCoord
      (cerr, "  Fixed cubes: ", fixed_cube[0], " and ");
    scalar_grid.PrintIndexAndCoord(cerr, "", fixed_cube[1], "\n");
    cerr << "  pointX: ";
    IJK::print_coord3D(cerr, pointX);
    cerr << " plane_normal: ";
    IJK::print_coord3D(cerr, plane_normal);
    cerr << endl;
    cerr << "  new_isovert_coord: ";
    IJK::print_coord3D(cerr, new_isovert_coord);
    cerr << endl;
  }

  COORD_TYPE Linf_distance;
  compute_rescaled_Linf_distance
    (pointX, new_isovert_coord, spacing, Linf_distance);
  if (Linf_distance > 1) { return; }

  bool flag_set = false;
  for (NUM_TYPE i = 0; i < NUM_CUBE_VERTICES3D; i++) {
    VERTEX_INDEX cube1_index = scalar_grid.CubeVertex(cube0_index, i);

    if (cube1_index == fixed_cube[0] || cube1_index == fixed_cube[1]) 
      { continue; }

    if (cube_contains_point(scalar_grid, cube1_index, new_isovert_coord)) {

      bool flag_set1;
      set_isovert_position_from_face
        (scalar_grid, isovalue, cube1_index, new_isovert_coord, 2, true, 
         isovert, flag_set1);

      if (flag_set1) { flag_set = true; }
    }
  }

  if (!flag_set) {
    // Cube containing point is not active.
    // Assign point to all active cubes around vertex 
    //   other than fixed_cube[0] and fixed_cube[1].

    for (NUM_TYPE i = 0; i < NUM_CUBE_VERTICES3D; i++) {
      VERTEX_INDEX cube1_index = scalar_grid.CubeVertex(cube0_index, i);

      if (cube1_index == fixed_cube[0] || cube1_index == fixed_cube[1]) 
        { continue; }

      set_isovert_position_from_face
        (scalar_grid, isovalue, cube1_index, new_isovert_coord, 2, 
         true, isovert);
    }
  }

}

/// Recompute isovert position around vertices.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_vertex
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    VERTEX_INDEX cube_index = isovert.CubeIndex(i);

    for (NUM_TYPE j = 0; j < NUM_CUBE_FACET_VERTICES3D; j++) {

      VERTEX_INDEX iv = scalar_grid.FacetVertex(cube_index, 0, j);
      recompute_isovert_position_around_vertex
        (scalar_grid, gradient_grid, isovalue, 
         isovert_param, voxel, flag_min_offset, iv, isovert);
    }
  }
}

/// Recompute isovert position around vertices.
void recompute_isovert_position_around_vertex
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  voxel.SetVertexCoord
    (scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);

  recompute_isovert_position_around_vertex
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel, 
     false, isovert);
}

/// @pre cube_coord[0][d] != cube_coord[1][d]
void set_plane_normal_coord_around_edge
(const int d, const GRID_COORD_TYPE cube_coord[2][DIM3], 
 COORD_TYPE plane_normal[DIM3])
{
  const COORD_TYPE ONE_DIV_SQRT2 = 1.0/std::sqrt(2);

  if (cube_coord[0][d] < cube_coord[1][d])  
    { plane_normal[d] = ONE_DIV_SQRT2; }
  else
    { plane_normal[d] = -ONE_DIV_SQRT2; }
}

/// Recompute isovert position around edge.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_edge
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 const VERTEX_INDEX iend0,
 const int edge_dir,
 ISOVERT & isovert)
{
  const VERTEX_INDEX * axis_increment = scalar_grid.AxisIncrement();
  const COORD_TYPE * spacing = scalar_grid.SpacingPtrConst();
  const ANGLE_TYPE max_angle_LO = 80;
  const COORD_TYPE min_abs_cos_angle_LO = cos(max_angle_LO*M_PI/180.0);
  const COORD_TYPE min_significant_distance =
    isovert_param.min_significant_distance;
  BOUNDARY_BITS_TYPE boundary_bits;
  VERTEX_INDEX fixed_cube[2];
  COORD_TYPE pointX[DIM3];
  GRID_COORD_TYPE cube_coord[2][DIM3];
  COORD_TYPE plane_normal[DIM3];
  COORD_TYPE new_isovert_coord[DIM3];
  COORD_TYPE intersection_point[DIM3];
  SVD_INFO svd_info;

  // Note: Compute vertex (not cube) boundary bits.
  scalar_grid.ComputeBoundaryBits(iend0, boundary_bits);

  // *** SHOULD CHECK THAT EDGE IS INTERNAL ***
  if (boundary_bits != 0) { return; }

  int d1 = (edge_dir+1)%DIM3;
  int d2 = (edge_dir+2)%DIM3;

  NUM_TYPE num_active_cubes = 0;

  // Number of cubes containing their isovertCoord.
  NUM_TYPE num_fixed_cubes = 0;

  // Number of fixed cubes with exactly two large eigenvalues.
  NUM_TYPE num_fixed_cubes_with_eigenvalue_two = 0;

  for (int j1 = 0; j1 < 2; j1++) {
    for (int j2 = 0; j2 < 2; j2++) {
      VERTEX_INDEX cube_index = 
        iend0 - j1*axis_increment[d1] - j2*axis_increment[d2];
      INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);

      if (gcube_index == ISOVERT::NO_INDEX) { continue; }

      num_active_cubes++;

      const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube_index);
      
      if (cube_contains_point(scalar_grid, cube_index, isovert_coord)) {

        NUM_TYPE num_eigenvalues = isovert.NumEigenvalues(gcube_index);

        if (num_eigenvalues < 2) { 
          // This cube will connect sharp cubes, if they exist.
          return; 
        }
        else if (num_eigenvalues == 2) {
          num_fixed_cubes_with_eigenvalue_two++;
        }

        if (num_fixed_cubes >= 2) { return; }

        fixed_cube[num_fixed_cubes] = cube_index;
        num_fixed_cubes++;
      }
    }
  }

  if (num_fixed_cubes != 2) { return; }
  if (num_active_cubes <= 2) { return; }

  for (int i = 0; i < 2; i++)
    { scalar_grid.ComputeCoord(fixed_cube[i], cube_coord[i]); }

  if (cube_coord[0][d1] == cube_coord[1][d1]) { return; }
  if (cube_coord[0][d2] == cube_coord[1][d2]) { return; }

  for (int i = 0; i < 2; i++) {
    NUM_TYPE gcube_index = isovert.GCubeIndex(fixed_cube[i]);

    if (isovert.NumEigenvalues(gcube_index) == 2) {

      if (!does_sharp_edge_point_to_cube
          (scalar_grid, isovert, fixed_cube[i], fixed_cube[1-i])) {
        // Sharp edge in fixed_cube[i] does not point to fixed_cube[1-i].
        return;
      }
    }
  }

  plane_normal[edge_dir] = 0;
  set_plane_normal_coord_around_edge(d1, cube_coord, plane_normal);
  set_plane_normal_coord_around_edge(d2, cube_coord, plane_normal);
  IJK::normalize_vector(DIM3, plane_normal, 0, plane_normal);

  scalar_grid.ComputeScaledCoord(iend0, pointX);
  pointX[edge_dir] += scalar_grid.Spacing(edge_dir)/2.0;

  compute_isovert_on_plane_from_sharp_edges
    (scalar_grid, isovert, isovert_param, fixed_cube, pointX, plane_normal,
     new_isovert_coord);

  MSDEBUG();
  if (flag_debug) {
    cerr << "--- Function: " << __func__ << endl;
    scalar_grid.PrintIndexAndCoord(cerr, "  Edge endpoint: ", iend0, "");
    cerr << "  edge_dir: " << edge_dir << endl;
    scalar_grid.PrintIndexAndCoord
      (cerr, "  Fixed cubes: ", fixed_cube[0], " and ");
    scalar_grid.PrintIndexAndCoord(cerr, "", fixed_cube[1], "\n");
    cerr << "  pointX: ";
    IJK::print_coord3D(cerr, pointX);
    cerr << " plane_normal: ";
    IJK::print_coord3D(cerr, plane_normal);
    cerr << endl;
    cerr << "  new_isovert_coord: ";
    IJK::print_coord3D(cerr, new_isovert_coord);
    cerr << endl;
  }

  COORD_TYPE Linf_distance;
  compute_rescaled_Linf_distance
    (pointX, new_isovert_coord, spacing, Linf_distance);
  if (Linf_distance > 1) { return; }

  bool flag_set = false;
  for (int j1 = 0; j1 < 2; j1++) {
    for (int j2 = 0; j2 < 2; j2++) {
      VERTEX_INDEX cube_index = 
        iend0 - j1*axis_increment[d1] - j2*axis_increment[d2];

      if (cube_index == fixed_cube[0] || cube_index == fixed_cube[1]) 
        { continue; }

      if (cube_contains_point(scalar_grid, cube_index, new_isovert_coord)) {

        bool flag_set1;
        set_isovert_position_from_face
          (scalar_grid, isovalue, cube_index, new_isovert_coord, 2, false, 
           isovert, flag_set1);

        if (flag_set1) { flag_set = true; }
      }
    }
  }

  if (!flag_set) {
    // Cube containing point is not active.
    // Assign point to all active cubes around vertex 
    //   other than fixed_cube[0] and fixed_cube[1].

    for (int j1 = 0; j1 < 2; j1++) {
      for (int j2 = 0; j2 < 2; j2++) {
        VERTEX_INDEX cube_index = 
          iend0 - j1*axis_increment[d1] - j2*axis_increment[d2];

        if (cube_index == fixed_cube[0] || cube_index == fixed_cube[1]) 
          { continue; }

        set_isovert_position_from_face
          (scalar_grid, isovalue, cube_index, new_isovert_coord, 2, 
           false, isovert);
      }
    }

  }

}

/// Recompute isovert position around edges.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_edge
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

    VERTEX_INDEX cube_index = isovert.CubeIndex(i);

    // Process 6 edges per cube.
    for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {

      VERTEX_INDEX iend0 = cube_index;
      recompute_isovert_position_around_edge
        (scalar_grid, gradient_grid, 
         isovalue, isovert_param, voxel, flag_min_offset, 
         iend0, edge_dir, isovert);

      int d1 = (edge_dir+1)%DIM3;
      iend0 = scalar_grid.NextVertex(iend0, d1);
      recompute_isovert_position_around_edge
        (scalar_grid, gradient_grid, 
         isovalue, isovert_param, voxel, flag_min_offset, 
         iend0, edge_dir, isovert);
    }
  }
}

/// Recompute isovert position around edges.
void recompute_isovert_position_around_edge
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  voxel.SetVertexCoord
    (scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);

  recompute_isovert_position_around_edge
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel, 
     false, isovert);
}

/// Recompute isovert position around edge.
/// Recompute for cubes which do not contain their isovert coordinates.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_edge_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 const VERTEX_INDEX iend0,
 const int edge_dir,
 const bool flag_slope_up,
 ISOVERT & isovert)
{
  const VERTEX_INDEX * axis_increment = scalar_grid.AxisIncrement();
  const COORD_TYPE * spacing = scalar_grid.SpacingPtrConst();
  const ANGLE_TYPE max_angle_LO = 80;
  const COORD_TYPE min_abs_cos_angle_LO = cos(max_angle_LO*M_PI/180.0);
  const COORD_TYPE min_significant_distance =
    isovert_param.min_significant_distance;
  const COORD_TYPE ONE_DIV_SQRT2 = 1.0/std::sqrt(2);
  BOUNDARY_BITS_TYPE boundary_bits;
  VERTEX_INDEX orthogonal_cube[2];
  COORD_TYPE pointX[DIM3];
  COORD_TYPE plane_normal[DIM3];
  COORD_TYPE new_isovert_coord[DIM3];
  COORD_TYPE intersection_point[DIM3];
  SVD_INFO svd_info;

  // Note: Compute vertex (not cube) boundary bits.
  scalar_grid.ComputeBoundaryBits(iend0, boundary_bits);

  // *** SHOULD CHECK THAT EDGE IS INTERNAL ***
  if (boundary_bits != 0) { return; }

  int d1 = (edge_dir+1)%DIM3;
  int d2 = (edge_dir+2)%DIM3;

  if (flag_slope_up) {

    // Cubes "orthogonal" to direction (*,1,1) (or (1,*,1) or (1,1,*))
    orthogonal_cube[0] = iend0 - axis_increment[d1];
    orthogonal_cube[1] = iend0 - axis_increment[d2];

    plane_normal[d1] = ONE_DIV_SQRT2;
    plane_normal[d2] = ONE_DIV_SQRT2;
  }
  else {
    // Cubes "orthogonal" to direction (0,1,-1) (or (1,0,-1) or (1,-1,0))
    orthogonal_cube[0] = iend0;
    orthogonal_cube[1] = iend0 - axis_increment[d1] - axis_increment[d2];

    plane_normal[d1] = ONE_DIV_SQRT2;
    plane_normal[d2] = -ONE_DIV_SQRT2;
  }
  IJK::normalize_vector(DIM3, plane_normal, 0, plane_normal);

  bool flag_recompute = false;
  for (int j = 0; j < 2; j++) {
    VERTEX_INDEX cube_index = orthogonal_cube[j];
    INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);

    if (gcube_index == ISOVERT::NO_INDEX) { continue; }

    const COORD_TYPE * isovert_coord = isovert.IsoVertCoord(gcube_index);
    if (cube_contains_point(scalar_grid, cube_index, isovert_coord)) { 
      // Some orthogonal cube contains an isovert coord.
      return;
    }
    else {

      if (isovert.NumEigenvalues(gcube_index) != 2) { continue; }

      // Check isovert edge direction matches slope.
      const COORD_TYPE * edge_dir = isovert.gcube_list[gcube_index].direction;
      if (flag_slope_up) {
        if (edge_dir[d1]*edge_dir[d2] > 0) { flag_recompute = true; }
      }
      else {
        if (edge_dir[d1]*edge_dir[d2] < 0) { flag_recompute = true; }
      }
    }
  }

  if (!flag_recompute) { return; }

  scalar_grid.ComputeScaledCoord(iend0, pointX);
  pointX[edge_dir] += scalar_grid.Spacing(edge_dir)/2.0;

  bool flag_succeeded;
  compute_isovert_on_plane_from_sharp_edges_B
    (scalar_grid, isovert, isovert_param, iend0, edge_dir, 
     pointX, plane_normal, new_isovert_coord, flag_succeeded);

  if (!flag_succeeded) { return; }

  MSDEBUG();
  if (flag_debug) {
    cerr << "--- Function: " << __func__ << endl;
    scalar_grid.PrintIndexAndCoord(cerr, "  Edge endpoint: ", iend0, "");
    cerr << "  edge_dir: " << edge_dir << endl;
    scalar_grid.PrintIndexAndCoord
      (cerr, "  Orthogonal cubes: ", orthogonal_cube[0], " and ");
    scalar_grid.PrintIndexAndCoord(cerr, "", orthogonal_cube[1], "\n");
    cerr << "  pointX: ";
    IJK::print_coord3D(cerr, pointX);
    cerr << " plane_normal: ";
    IJK::print_coord3D(cerr, plane_normal);
    cerr << endl;
    cerr << "  new_isovert_coord: ";
    IJK::print_coord3D(cerr, new_isovert_coord);
    cerr << endl;
  }

  COORD_TYPE Linf_distance;
  GRID_COORD_TYPE end0_coord[DIM3];
  compute_rescaled_Linf_distance
    (pointX, new_isovert_coord, spacing, Linf_distance);
  if (Linf_distance > 1) { return; }
  scalar_grid.ComputeCoord(iend0, end0_coord);
  if ((new_isovert_coord[edge_dir] < end0_coord[edge_dir]) ||
      (new_isovert_coord[edge_dir] > end0_coord[edge_dir]+1))
    { return; }

  bool flag_set = false;
  for (int j = 0; j < 2; j++) {
    VERTEX_INDEX cube_index = orthogonal_cube[j];
    if (cube_contains_point(scalar_grid, cube_index, new_isovert_coord)) {

      bool flag_set1;
      set_isovert_position_from_face
        (scalar_grid, isovalue, cube_index, new_isovert_coord, 2, false, 
         isovert, flag_set1);

      if (flag_set1) { flag_set = true; }
    }
  }

  if (!flag_set) {
    // Cube containing point is not active.
    // Assign point to all active orthogonal_cubes around vertex 

    MSDEBUG();
    if (flag_debug)
      { cerr << "  Cube containing the point is NOT active." << endl; }

    for (int j = 0; j < 2; j++) {
      VERTEX_INDEX cube_index = orthogonal_cube[j];
      set_isovert_position_from_face
        (scalar_grid, isovalue, cube_index, new_isovert_coord, 2, 
         false, isovert);
    }
  }

}

/// Recompute isovert position around edge.
/// Recompute for cubes which do not contain their isovert coordinates.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_edge_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 const VERTEX_INDEX iend0,
 const int edge_dir,
 ISOVERT & isovert)
{
  // Recompute with flag_slope_up = false.
  recompute_isovert_position_around_edge_B
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel,
     flag_min_offset, iend0, edge_dir, false, isovert);

  // Recompute with flag_slope_up = true.
  recompute_isovert_position_around_edge_B
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel,
     flag_min_offset, iend0, edge_dir, true, isovert);
}

/// Recompute isovert position around edges.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_isovert_position_around_edge_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert)
{
  const VERTEX_INDEX * axis_increment = scalar_grid.AxisIncrement();

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    int d1 = (edge_dir+1)%DIM3;
    int d2 = (edge_dir+2)%DIM3;

    bool flag_sharp = false;
    for (int j1 = 0; j1 < 2; j1++) {
      for (int j2 = 0; j2 < 2; j2++) {
        VERTEX_INDEX cube_index = 
          iend0 - j1*axis_increment[d1] - j2*axis_increment[d2];
        INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);

        if (gcube_index == ISOVERT::NO_INDEX) { continue; }
        if (isovert.NumEigenvalues(gcube_index) >= 2) 
          { flag_sharp = true; }
      }
    }

    if (flag_sharp) {
      recompute_isovert_position_around_edge_B
        (scalar_grid, gradient_grid, isovalue, isovert_param,
         voxel, flag_min_offset, iend0, edge_dir, isovert); 
    }
  }

}

/// Recompute isovert position around edges.
void recompute_isovert_position_around_edge_B
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  OFFSET_VOXEL voxel;

  voxel.SetVertexCoord
    (scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);

  recompute_isovert_position_around_edge_B
    (scalar_grid, gradient_grid, isovalue, isovert_param, voxel, 
     false, isovert);
}

void round_down(const COORD_TYPE * coord, GRID_COORD_TYPE min_coord[DIM3])
{
  for (NUM_TYPE d = 0; d < DIM3; d++)
	{ min_coord[d] = int(std::floor(coord[d])); }
}

void MERGESHARP::set_edge_index
(const std::vector<COORD_TYPE> & edgeI_coord,
 SHARPISO_EDGE_INDEX_GRID & edge_index)
{
  GRID_COORD_TYPE min_coord[DIM3];

  edge_index.SetAllCoord(ISOVERT::NO_INDEX);

  for (NUM_TYPE i = 0; i < edgeI_coord.size()/DIM3; i++) {
    round_down(&(edgeI_coord[i*DIM3]), min_coord);

    int edge_dir = 0;
    for (int d = 1; d < DIM3; d++) {
      if (edgeI_coord[i*DIM3+d] != min_coord[d]) { edge_dir = d; }
    }

    VERTEX_INDEX iv0 = edge_index.ComputeVertexIndex(min_coord);

    edge_index.Set(iv0, edge_dir, i);
  }

  // Set edge index from edgeI_coord[] which are on grid vertices.
  for (NUM_TYPE i = 0; i < edgeI_coord.size()/DIM3; i++) {
    round_down(&(edgeI_coord[i*DIM3]), min_coord);

    if (is_coord_equal_3D(&(edgeI_coord[i*DIM3]), min_coord)) {

      VERTEX_INDEX iv0 = edge_index.ComputeVertexIndex(min_coord);
      for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {
        if (edge_index.Vector(iv0, edge_dir) == ISOVERT::NO_INDEX)
          { edge_index.Set(iv0, edge_dir, i); }

        if (min_coord[edge_dir] > 0) {
          VERTEX_INDEX iv1 = edge_index.PrevVertex(iv0, edge_dir);
          if (edge_index.Vector(iv1, edge_dir) == ISOVERT::NO_INDEX)
            { edge_index.Set(iv1, edge_dir, i); }
        }
      }
    }


  }

}

// Compute isosurface vertex positions using hermite data.
void compute_all_isovert_positions 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const std::vector<COORD_TYPE> & edgeI_normal_coord,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isovert_param.grad_selection_cube_offset;
  SHARPISO_EDGE_INDEX_GRID edge_index(DIM3, scalar_grid.AxisSize(), DIM3);
  SVD_INFO svd_info;

  set_edge_index(edgeI_coord, edge_index);

  IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
    NUM_TYPE index = isovert.GCubeIndex(iv);
    if (index!=ISOVERT::NO_INDEX) {
      // this is an active cube
      // compute the sharp vertex for this cube
      EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
      VERTEX_INDEX num_large_eigenvalues;

      svd_compute_sharp_vertex_for_cube_hermite
        (scalar_grid, edgeI_coord, edgeI_normal_coord, edge_index,
         iv, isovalue, isovert_param,
         isovert.gcube_list[index].isovert_coord,
         eigenvalues, num_large_eigenvalues, svd_info);

      //set num eigen
      isovert.gcube_list[index].num_eigenvalues =  
        (unsigned char)num_large_eigenvalues;
      //set the sharp vertex type to be *AVAILABLE*
      if(num_large_eigenvalues > 1)
        { isovert.gcube_list[index].flag=AVAILABLE_GCUBE; }
      else
        { isovert.gcube_list[index].flag=SMOOTH_GCUBE; }

      // store distance
      compute_Linf_distance_from_cube_center
        (scalar_grid, iv, isovert.gcube_list[index].isovert_coord,
         isovert.gcube_list[index].linf_dist);
    }
  }

  store_boundary_bits(scalar_grid, isovert.gcube_list);
}

/// Set isovert position
/// @pre new_coord[] is in cube gcube_index.
void set_isovert_position
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, 
 const COORD_TYPE * new_coord, const NUM_TYPE num_eigenvalues,
 ISOVERT & isovert)
{
  IJK::copy_coord_3D(new_coord,
                     isovert.gcube_list[gcube_index].isovert_coord);
  isovert.gcube_list[gcube_index].num_eigenvalues = num_eigenvalues;
  isovert.gcube_list[gcube_index].cube_containing_isovert = 
    isovert.gcube_list[gcube_index].cube_index;
  isovert.gcube_list[gcube_index].flag_conflict = false;
  isovert.gcube_list[gcube_index].flag_centroid_location = false;

  if (num_eigenvalues > 1) {
    if (isovert.gcube_list[gcube_index].flag == COVERED_POINT ||
        isovert.gcube_list[gcube_index].flag == SMOOTH_GCUBE)
      { isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE; }
  }
  else {
    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE) 
      { isovert.gcube_list[gcube_index].flag = SMOOTH_GCUBE; }
  }

  // store new Linf distance
  compute_Linf_distance_from_cube_center(grid, gcube_index, isovert);
}

/// Set isovert position
/// @pre new_coord[] is in cube gcube_index.
void set_isovert_position_from_other_cube
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, 
 const COORD_TYPE * new_coord, const NUM_TYPE num_eigenvalues,
 ISOVERT & isovert)
{
  set_isovert_position
    (grid, gcube_index, new_coord, num_eigenvalues, isovert);
  isovert.gcube_list[gcube_index].flag_coord_from_other_cube = true;
}

// Swap isovert positions of gcubeA and gcubeB
void swap_isovert_positions
(const SHARPISO_GRID & grid, 
 const NUM_TYPE gcubeA_index, const NUM_TYPE gcubeB_index,
 ISOVERT & isovert)
{
  COORD_TYPE temp_coord[DIM3];
  NUM_TYPE num_eigenvaluesA, num_eigenvaluesB;

  num_eigenvaluesA = isovert.gcube_list[gcubeA_index].num_eigenvalues;
  num_eigenvaluesB = isovert.gcube_list[gcubeB_index].num_eigenvalues;

  const COORD_TYPE * isovert_coordA =
    isovert.gcube_list[gcubeA_index].isovert_coord;
  const COORD_TYPE * isovert_coordB =
    isovert.gcube_list[gcubeB_index].isovert_coord;

  IJK::copy_coord_3D(isovert_coordA, temp_coord);

  set_isovert_position_from_other_cube
    (grid, gcubeA_index, isovert_coordB, num_eigenvaluesB, isovert);
  set_isovert_position_from_other_cube
    (grid, gcubeB_index, temp_coord, num_eigenvaluesA, isovert);
}

// Copy isovert position from from_gcube to to_gcube.
void MERGESHARP::copy_isovert_position
(const SHARPISO_GRID & grid, 
 const NUM_TYPE from_gcube, const NUM_TYPE to_gcube,
 ISOVERT & isovert)
{
  const COORD_TYPE * new_coord =
    isovert.gcube_list[from_gcube].isovert_coord;
  const NUM_TYPE num_eigenvalues =
    isovert.gcube_list[from_gcube].num_eigenvalues;

  set_isovert_position_from_other_cube
    (grid, to_gcube, new_coord, num_eigenvalues, isovert);
}


// Change isovert positions in cubes with conflicts.
void swap_isovert_positions
(const SHARPISO_GRID & grid, 
 const COORD_TYPE max_dist_to_set_other, ISOVERT & isovert)
{
  std::vector<NUM_TYPE> gcube_sharp_list;

  get_corner_or_edge_cubes(isovert.gcube_list, gcube_sharp_list);

  for (NUM_TYPE i = 0; i < gcube_sharp_list.size(); i++) {

    NUM_TYPE gcubeA_index = gcube_sharp_list[i];
    if (isovert.gcube_list[gcubeA_index].flag_conflict) {

      if (isovert.gcube_list[gcubeA_index].linf_dist > 
          max_dist_to_set_other) { continue; }

      VERTEX_INDEX conflicting_cube = 
        isovert.gcube_list[gcubeA_index].cube_containing_isovert;
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(conflicting_cube);

      if (gcubeB_index != ISOVERT::NO_INDEX) {

        if (isovert.gcube_list[gcubeB_index].flag_conflict ||
            isovert.gcube_list[gcubeB_index].flag_centroid_location) {

          VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
          if (isovert.gcube_list[gcubeB_index].cube_containing_isovert == 
              cubeA_index) {

            // swap isovert_coord for gcubeA_index and gcubeB_index
            swap_isovert_positions
              (grid, gcubeA_index, gcubeB_index, isovert);
          }
          else {

            MSDEBUG();
            if (flag_debug) {
              grid.PrintIndexAndCoord
                (cerr, "*** Setting cube: ", conflicting_cube, "");
              grid.PrintIndexAndCoord
                (cerr, " from ", isovert.CubeIndex(gcubeA_index), "\n");
              cerr << "    Cube: " << conflicting_cube;
              cerr << "  isovert_coord: ";
              IJK::print_coord3D(cerr, isovert.gcube_list[gcubeB_index].isovert_coord, "\n");
              cerr << "    flag_far: "
                   << int(isovert.gcube_list[gcubeB_index].flag_far);
              cerr << " flag_conflict: "
                   << int(isovert.gcube_list[gcubeB_index].flag_conflict);
              cerr << " flag_conflict: "
                   << int(isovert.gcube_list[gcubeB_index].flag_centroid_location);
              cerr << endl;
              cerr << "    Cube: " << isovert.CubeIndex(gcubeA_index);
              cerr << "  flag_far: "
                   << int(isovert.gcube_list[gcubeA_index].flag_far);
              cerr << "  flag_centroid: "
                   << int(isovert.gcube_list[gcubeA_index].flag_centroid_location)
                   << endl;
            }

            copy_isovert_position
              (grid, gcubeA_index, gcubeB_index, isovert);
          }
        }
      }
      else {
        // Should never happen, but...
        IJK::PROCEDURE_ERROR error("swap_isovert_positions");
        error.AddMessage("Programming error.  Illegal conflicting cube index: ",
                         conflicting_cube, " for cube ",
                         isovert.CubeIndex(gcubeA_index), ".");
        throw error;

      }
    }
  }

}

// Replace with substitute isovert coord
void replace_with_substitute_coord
(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert)
{
  const COORD_TYPE * isovert_coordB =
    isovert.gcube_list[gcube_index].isovert_coordB;

  MSDEBUG();
  if (flag_debug) {
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    const COORD_TYPE * isovert_coordA =
      isovert.gcube_list[gcube_index].isovert_coord;
    cerr << "    Using substitute isovert_coord for ";
    grid.PrintIndexAndCoord(cerr, "", cube_index, "\n");
    cerr << "      Old coord: ";
    IJK::print_coord3D(cerr, isovert_coordA);
    cerr << ".  New coord: ";
    IJK::print_coord3D(cerr, isovert_coordB);
    cerr << " num eigenvalues: " 
         << int(isovert.gcube_list[gcube_index].num_eigenvalues);
    cerr << endl;
  }

  NUM_TYPE num_eigenvalues =
    isovert.gcube_list[gcube_index].num_eigenvalues;
  set_isovert_position
    (grid, gcube_index, isovert_coordB, num_eigenvalues, isovert);
  isovert.gcube_list[gcube_index].flag_using_substitute_coord = true;
}

// Apply secondary isovert positions
void apply_secondary_isovert_positions
(const SHARPISO_GRID & grid, 
 const COORD_TYPE max_dist_to_set_other, ISOVERT & isovert)
{
  std::vector<NUM_TYPE> gcube_sharp_list;

  get_corner_or_edge_cubes(isovert.gcube_list, gcube_sharp_list);

  for (NUM_TYPE i = 0; i < gcube_sharp_list.size(); i++) {

    NUM_TYPE gcubeA_index = gcube_sharp_list[i];
    if (isovert.gcube_list[gcubeA_index].flag_conflict) {

      const COORD_TYPE * isovert_coordB =
        isovert.gcube_list[gcubeA_index].isovert_coordB;
      VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);

      if (cube_contains_point(grid, cubeA_index, isovert_coordB)) {

        replace_with_substitute_coord(grid, gcubeA_index, isovert);
      }
      else if (grid.ContainsPoint(isovert_coordB)) {

        if (isovert.gcube_list[gcubeA_index].linf_dist > 
            max_dist_to_set_other) { continue; }

        VERTEX_INDEX cubeB_index;
        bool flag_boundary;
        get_cube_containing_point
          (grid, isovert_coordB, cubeB_index, flag_boundary);

        INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

        if (gcubeB_index != ISOVERT::NO_INDEX) {

          if (isovert.gcube_list[gcubeB_index].flag_conflict) {

            MSDEBUG();
            if (flag_debug) {
              const COORD_TYPE * isovert_coordA =
                isovert.gcube_list[gcubeA_index].isovert_coord;
              const COORD_TYPE * directionA =
                isovert.gcube_list[gcubeA_index].direction;
              const COORD_TYPE * directionB =
                isovert.gcube_list[gcubeB_index].direction;
              cerr << "*** Copying substitute isovert_coord from " 
                   << cubeA_index << " ";
              ijkgrid_output_vertex_coord(cerr, grid, cubeA_index);
              cerr << "  to cube " << cubeB_index << " ";
              ijkgrid_output_vertex_coord(cerr, grid, cubeB_index);
              cerr << endl;
              cerr << "  Old coord: ";
              IJK::print_coord3D(cerr, isovert_coordA);
              cerr << "  Old direction: ";
              IJK::print_coord3D(cerr, directionA);
              cerr << endl;
              cerr << "  New coord: ";
              IJK::print_coord3D(cerr, isovert_coordB);
              cerr << endl;
            }

            NUM_TYPE num_eigenvalues =
              isovert.gcube_list[gcubeA_index].num_eigenvalues;

            set_isovert_position_from_other_cube
              (grid, gcubeB_index, isovert_coordB, num_eigenvalues, isovert);
            isovert.gcube_list[gcubeB_index].flag_using_substitute_coord = true;
          }
        }
      }
    }
  }

}

/// Get cubes with more than one eigenvalue.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to sharp cubes
///    sorted by number of large eigenvalues and by distance 
///    of sharp coord from cube center.
void MERGESHARP::get_corner_or_edge_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list)
{
  GCUBE_COMPARE gcube_compare(gcube_list);

  for (int i=0;i<gcube_list.size();i++)
    {
      // *** SHOULD CHECK THAT POINT IS NOT COMPUTED USING CENTROID ***
      if (gcube_list[i].num_eigenvalues > 1)
        { gcube_index_list.push_back(i); }
    }

  sort(gcube_index_list.begin(), gcube_index_list.end(), gcube_compare);
}

/// Get selected cubes.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to selected cubes
///    sorted by number of large eigenvalues and by distance 
///    of sharp coord from cube center.
void MERGESHARP::get_selected_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list)
{
  GCUBE_COMPARE gcube_compare(gcube_list);

  for (int i=0;i<gcube_list.size();i++)
    {
      if (gcube_list[i].flag == SELECTED_GCUBE)
        { gcube_index_list.push_back(i); }
    }

  sort(gcube_index_list.begin(), gcube_index_list.end(), gcube_compare);
}

/// Get selected corner cubes.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to selected corner cubes
///    sorted by distance of sharp coord from cube center.
void MERGESHARP::get_selected_corner_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list)
{
  GCUBE_COMPARE gcube_compare(gcube_list);

  for (int i=0;i<gcube_list.size();i++) {
    if (gcube_list[i].flag == SELECTED_GCUBE) {
      if (gcube_list[i].num_eigenvalues == 3) 
        { gcube_index_list.push_back(i); }
    }
  }

  sort(gcube_index_list.begin(), gcube_index_list.end(), gcube_compare);
}

void MERGESHARP::store_boundary_bits
(const SHARPISO_GRID & grid, GRID_CUBE_DATA_ARRAY & gcube_list)
{
  for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {
    grid.ComputeBoundaryCubeBits
      (gcube_list[i].cube_index, gcube_list[i].boundary_bits);
  }
}

/// Compute the cube index from the gc index
VERTEX_INDEX cube_ind_frm_gc_ind
(const ISOVERT &isoData, const  NUM_TYPE &gc_ind){
  return isoData.gcube_list[gc_ind].cube_index;
}

/// get the *connected_list* of vertices which are selected and
/// connected to vertex *vert_ind*
void get_connected(
                   const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
                   const SCALAR_TYPE isovalue,
                   const VERTEX_INDEX vert_index,
                   const vector<VERTEX_INDEX> &selected_list,
                   vector<VERTEX_INDEX> &connected_list)
{
  for (int iv=0; iv < selected_list.size();iv++){
    if (are_connected(scalar_grid, vert_index, selected_list[iv], isovalue)){
      connected_list.push_back(selected_list[iv]);
    }
  }
}

// Check if selecting this vertex creates a triangle with a large angle.
// @param check_triangle_angle If true, check it triangle has large angles.
// @param bin_grid Contains the already selected vertices.
// @param[out] v1,v2 vertex indices which form a triangle with iv
bool MERGESHARP::creates_triangle 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const bool check_triangle_angle,
 const ISOVERT & isovert,
 const VERTEX_INDEX iv,
 const SCALAR_TYPE isovalue,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 VERTEX_INDEX & v1,
 VERTEX_INDEX & v2)
{
  const SCALAR_TYPE threshold = cos(140*M_PI/180);
  vector <VERTEX_INDEX> selected_list;
  vector <VERTEX_INDEX> connected_list;

  // get the selected vertices around iv
  get_selected(scalar_grid, iv, bin_grid, bin_width, selected_list);

  // get the list of vertices connected to the vertex iv
  get_connected(scalar_grid, isovalue, iv, selected_list, connected_list);

  MSDEBUG();
  bool flag_debug = false;

  int limit = connected_list.size();
  // for each pair jv1 jv2 in the connected list
  for (int i=0; i< limit-1; ++i) {
    for(int j=i+1; j <= (limit-1); ++j) {
      if (are_connected(scalar_grid, connected_list[i],
                        connected_list[j], isovalue)) {
        v1 = connected_list[i];
        v2 = connected_list[j];

        if (check_triangle_angle){
          // checking angle
          bool flag_large_angle = 
            is_angle_large(scalar_grid, isovert, iv, threshold, v1 , v2);

          if (flag_large_angle)
            { 
              // *** DEBUG ***
              /*
                if (flag_debug) {
                using namespace std;
                cerr << "   Cube " << iv << " ";
                ijkgrid_output_vertex_coord(cerr, scalar_grid, iv);
                cerr << " creates large angle triangle with:" << endl;
                cerr << "      " << v1 << " ";
                ijkgrid_output_vertex_coord(cerr, scalar_grid, v1);
                cerr << " and " << v2 << " ";
                ijkgrid_output_vertex_coord(cerr, scalar_grid, v2);
                cerr << endl;
                }
              */

              return true; 
            }
        }
        else 
          {
            // returning true without checking angles angles
            return true;
          }
      }
    }
  }

  return false;
}

// Check if three 3x3x3 regions pairwise overlap.
// Return true if all three pairs overlap.
bool do_3x3x3_regions_overlap
(const SHARPISO_GRID & grid,
 const VERTEX_INDEX v0, const VERTEX_INDEX v1, const VERTEX_INDEX v2)
{
  GRID_COORD_TYPE dist01, dist12, dist02;

  IJK::compute_Linf_distance_between_grid_vertices
    (grid, v0, v1, dist01);
  IJK::compute_Linf_distance_between_grid_vertices
    (grid, v1, v2, dist12);
  IJK::compute_Linf_distance_between_grid_vertices
    (grid, v0, v2, dist02);

  if (dist01 < 3 && dist12 < 3 && dist02 < 3)
    { return(true); }

  return(false);
}


// Check if selecting this vertex creates a triangle with a large angle.
// @param bin_grid Contains the already selected vertices.
// @param[out] iv1, iv2 vertex indices which form a triangle with iv0
bool MERGESHARP::creates_triangle_new
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const ISOVERT & isovert,
 const VERTEX_INDEX iv0,
 const SCALAR_TYPE isovalue,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 VERTEX_INDEX & iv1,
 VERTEX_INDEX & iv2)
{
  const COORD_TYPE threshold = cos(140*M_PI/180);
  vector <VERTEX_INDEX> selected_list;
  vector <VERTEX_INDEX> connected_list;

  // get the selected vertices around iv0
  get_selected(scalar_grid, iv0, bin_grid, bin_width, selected_list);

  // get the list of vertices connected to the vertex iv0
  get_connected(scalar_grid, isovalue, iv0, selected_list, connected_list);

  // for each pair jv1 jv2 in the connected list
  for (int i=0; i+1 < connected_list.size(); ++i)
    for (int j=i+1; j < connected_list.size(); ++j)
      {
        iv1 = connected_list[i];
        iv2 = connected_list[j];
        if (do_3x3x3_regions_overlap(scalar_grid, iv0, iv1, iv2)) {

          if (is_angle_large
              (scalar_grid, isovert, iv0, threshold, iv1 , iv2))
            { return true; }
        }
      }

  return false;
}


/// Set the index_grid, grid and gcube_list.
/// Traverse the scalar grid.
/// If cube icube is not active, set index_grid[icube] to ISOVERT::NO_INDEX.
/// If cube icube is active, then set index_grid[icube] to index.
/// Push the corresponding *grid_cube_data* into the *gcube_list*.
void create_active_cubes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 ISOVERT &isovert)
{
  NUM_TYPE index = 0;

  // Set the size of index_grid and grid.
  isovert.index_grid.SetSize(scalar_grid);
  isovert.grid.SetSize(scalar_grid);

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {
    if (is_gt_cube_min_le_cube_max(scalar_grid, icube, isovalue)) {
      index = isovert.gcube_list.size();
      isovert.index_grid.Set(icube,index);
      GRID_CUBE_DATA gc;
      gc.cube_index = icube;
      scalar_grid.ComputeCoord(icube, gc.cube_coord);
      isovert.gcube_list.push_back(gc);
		}
    else {
      // if it is not intersected then mark as ISOVERT::NO_INDEX
      isovert.index_grid.Set(icube,ISOVERT::NO_INDEX);
		}
	}
}

/// process edge called from are connected
void process_edge
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX v0,
 const VERTEX_INDEX v1,
 const SCALAR_TYPE isovalue,
 bool &is_intersect)
{
  is_intersect = is_gt_min_le_max(scalar_grid, v0, v1, isovalue);
}

/**
 * Compute dual isovert.
 */
void MERGESHARP::compute_dual_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const VERTEX_POSITION_METHOD vertex_position_method,
 ISOVERT & isovert)
{
  const COORD_TYPE max_dist_to_set_other =
    isovert_param.max_dist_to_set_other;
  IJK::PROCEDURE_ERROR error("compute_dual_isovert");

  if (!gradient_grid.Check
      (scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

  create_active_cubes(scalar_grid, isovalue, isovert);

  MSDEBUG();
  flag_debug = false;

  if (vertex_position_method == GRADIENT_POSITIONING) {

    compute_all_isovert_positions 
      (scalar_grid, gradient_grid, isovalue, isovert_param, 
       isovert);

    recompute_far_points
      (scalar_grid, gradient_grid, isovalue, isovert_param, isovert);
    swap_isovert_positions(scalar_grid, max_dist_to_set_other, isovert);
    apply_secondary_isovert_positions
      (scalar_grid, max_dist_to_set_other, isovert);

    recompute_isovert_position_around_vertex
      (scalar_grid, gradient_grid, isovalue, isovert_param, isovert);

    MSDEBUG();
    recompute_isovert_position_around_edge
      (scalar_grid, gradient_grid, isovalue, isovert_param, isovert);
  }
  else {
    compute_all_isovert_positions_edgeI
      (scalar_grid, gradient_grid, isovalue, isovert_param, 
       vertex_position_method, isovert);
  }
}

void MERGESHARP::compute_dual_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT &isovert)
{
  create_active_cubes(scalar_grid, isovalue, isovert);

  compute_all_isovert_positions 
    (scalar_grid, edgeI_coord, edgeI_normal_coord,
     isovalue, isovert_param, isovert);
}


// Select all grid cubes which are not smooth.
void MERGESHARP::select_non_smooth(ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    if (isovert.gcube_list[i].flag != SMOOTH_GCUBE) 
      { isovert.gcube_list[i].flag = SELECTED_GCUBE; };
  }
}

// Get list of grid cubes from isovert.
void MERGESHARP::get_cube_list
(const ISOVERT & isovert, std::vector<VERTEX_INDEX> & cube_list)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();

  cube_list.resize(num_gcube);
  for (NUM_TYPE i = 0; i < num_gcube; i++) {
    cube_list[i] = isovert.gcube_list[i].cube_index;
  }
}

// Set cover type of all covered cubes.
void MERGESHARP::set_cover_type(ISOVERT & isovert)
{

  for (int i = 0; i < isovert.gcube_list.size(); i++)
    { isovert.gcube_list[i].cover_type = NO_ADJACENT_CUBE; }

  // Set VERTEX_ADJACENT.
  for (int i = 0; i < isovert.gcube_list.size(); i++) {

    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

      const BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[i].boundary_bits;

      if (boundary_bits == 0) {
        VERTEX_INDEX cube0_index = isovert.gcube_list[i].cube_index;

        for (int j = 0; j < isovert.grid.NumCubeNeighborsV(); j++) {
          VERTEX_INDEX cube1_index = isovert.grid.CubeNeighborV(cube0_index, j);
          INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
          if (gcube1_index != ISOVERT::NO_INDEX) {
            isovert.gcube_list[gcube1_index].cover_type = VERTEX_ADJACENT;
          }
        }
      }
      else {
        // Handle boundary case.
      }
    }
  }

  // Set EDGE_ADJACENT.
  for (int i = 0; i < isovert.gcube_list.size(); i++) {

    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

      const BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[i].boundary_bits;

      if (boundary_bits == 0) {
        VERTEX_INDEX cube0_index = isovert.gcube_list[i].cube_index;

        for (int j = 0; j < isovert.grid.NumCubeNeighborsE(); j++) {
          VERTEX_INDEX cube1_index = isovert.grid.CubeNeighborE(cube0_index, j);
          INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
          if (gcube1_index != ISOVERT::NO_INDEX) {
            isovert.gcube_list[gcube1_index].cover_type = EDGE_ADJACENT;
          }
        }
      }
      else {
        // Handle boundary case.
      }
    }
  }

  // Set FACET_ADJACENT.
  for (int i = 0; i < isovert.gcube_list.size(); i++) {

    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

      const BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[i].boundary_bits;

      if (boundary_bits == 0) {
        VERTEX_INDEX cube0_index = isovert.gcube_list[i].cube_index;

        for (int j = 0; j < isovert.grid.NumCubeNeighborsF(); j++) {
          VERTEX_INDEX cube1_index = isovert.grid.CubeNeighborF(cube0_index, j);
          INDEX_DIFF_TYPE gcube1_index = isovert.GCubeIndex(cube1_index);
          if (gcube1_index != ISOVERT::NO_INDEX) {
            isovert.gcube_list[gcube1_index].cover_type = FACET_ADJACENT;
          }
        }
      }
      else {
        // Handle boundary case.
      }
    }
  }
  
}

// Store isosurface lookup table index in gcube_list.
void MERGESHARP::store_table_index
(const std::vector<IJKDUALTABLE::TABLE_INDEX> & table_index,
 GRID_CUBE_DATA_ARRAY & gcube_list)
{
  IJK::PROCEDURE_ERROR error("store_table_index");

  if (table_index.size() != gcube_list.size()) {
    error.AddMessage("Programming error.  Numbers of elements in table_index and gcube_list differ.");
    error.AddMessage("  table_index.size() = ", table_index.size(), ".");
    error.AddMessage("  gcube_list.size() = ", gcube_list.size(), ".");
    throw error;
  }

  for (NUM_TYPE i = 0; i < table_index.size(); i++) 
	{ gcube_list[i].table_index = table_index[i]; }
}

/// Compute the overlap region between two 3x3x3 regions around cubes.
/// @param[out] Overlap dimension.  No overlap has overlap dimension -1.
bool MERGESHARP::find_3x3x3_overlap
(const GRID_COORD_TYPE cubeA_coord[DIM3],
 const GRID_COORD_TYPE cubeB_coord[DIM3],
 GRID_COORD_TYPE rmin[DIM3],
 GRID_COORD_TYPE rmax[DIM3],
 int & overlap_dim)
{
  for (int d=0;d<DIM3;d++){
    rmin[d]=max(cubeA_coord[d]-1, cubeB_coord[d]-1);
    rmax[d]=min(cubeA_coord[d]+2, cubeB_coord[d]+2);
  }

  overlap_dim = 0;

  for (int d=0; d<DIM3; d++) {

    if(rmin[d] > rmax[d]) {
      overlap_dim = -1;
      return(false);
    }

    if (rmin[d] < rmax[d]) 
      { overlap_dim++; }
  }

  return(true);
}

/// Compute the overlap region between two 3x3x3 regions around cubes.
/// @param[out] Overlap dimension.  No overlap has overlap dimension -1.
template <typename GRID_TYPE>
bool MERGESHARP::find_3x3x3_overlap
(const GRID_TYPE & grid,
 const VERTEX_INDEX cubeA_index,
 const VERTEX_INDEX cubeB_index,
 GRID_COORD_TYPE rmin[DIM3],
 GRID_COORD_TYPE rmax[DIM3],
 int & overlap_dim)
{
  GRID_COORD_TYPE cubeA_coord[DIM3], cubeB_coord[DIM3];
  bool is_overlap;

  grid.ComputeCoord(cubeA_index, cubeA_coord);
  grid.ComputeCoord(cubeB_index, cubeB_coord);

  is_overlap = 
    find_3x3x3_overlap(cubeA_coord, cubeB_coord, rmin, rmax, overlap_dim);

  return(is_overlap);
}


// **************************************************
// Convert types to strings
// **************************************************

/// Transform GRID_CUBE_FLAG into a string
void MERGESHARP::convert2string
(const GRID_CUBE_FLAG & flag, std::string & s)
{
  switch(flag) {

  case AVAILABLE_GCUBE:
    s = "Available";
    break;

  case SELECTED_GCUBE:
    s = "Selected";
    break;

  case COVERED_A_GCUBE:
    s = "Covered (A)";
    break;

  case COVERED_B_GCUBE:
    s = "Covered (B)";
    break;

  case COVERED_CORNER_GCUBE:
    s = "Covered by corner";
    break;

  case COVERED_POINT:
    s = "Isovert covered";
    break;

  case UNAVAILABLE_GCUBE:
    s = "Unavailable";
    break;

  case NON_DISK_GCUBE:
    s = "Non-disk patch";
    break;

  case SMOOTH_GCUBE:
    s = "Smooth";
    break;

  default:
    s = "Unidentified";
    break;
  }
}

// Transform CUBE_ADJACENCY_TYPE into a string
void MERGESHARP::convert2string
(const CUBE_ADJACENCY_TYPE & adj_type, std::string & s)
{
  switch(adj_type) {

  case VERTEX_ADJACENT:
    s = "Vertex adjacent";
    break;

  case EDGE_ADJACENT:
    s = "Edge adjacent";
    break;

  case FACET_ADJACENT:
    s = "Facet adjacent";
    break;

  case NO_ADJACENT_CUBE:
    s = "No adjacent";
    break;

  default:
    s = "Unidentified";
    break;
  }

}


// **************************************************
// GRID_CUBE_DATA member functions
// **************************************************

void GRID_CUBE_DATA::Init()
{
  num_eigenvalues = 0;
  flag = AVAILABLE_GCUBE;
  cover_type = NO_ADJACENT_CUBE;
  boundary_bits = 0;
  cube_index = 0;
  linf_dist = 0;
  flag_centroid_location = false;
  flag_conflict = false;
  flag_near_corner = false;
  flag_coord_from_other_cube = false;
  flag_coord_from_vertex = false;
  flag_coord_from_edge = false;
  flag_using_substitute_coord = false;
  flag_recomputed_coord = false;
  flag_recomputed_coord_min_offset = false;
  flag_recomputed_using_adjacent = false;
  flag_far = false;
  cube_containing_isovert = 0;
  table_index = 0;
  covered_by = 0;
}

bool GRID_CUBE_DATA::IsCoveredOrSelected() const
{
  if (flag == COVERED_A_GCUBE || flag == COVERED_B_GCUBE ||
      flag == COVERED_CORNER_GCUBE || flag == SELECTED_GCUBE) 
    { return true; }
  else
    { return false; }
}

// **************************************************
// ISOVERT member functions
// **************************************************

bool ISOVERT::isFlag(const int cube_index,  GRID_CUBE_FLAG _flag) const
{
  if (isActive(cube_index)){
    if (gcube_list[index_grid.Scalar(cube_index)].flag == _flag)
      return true;
    else
      return false;
  }
  else
    return false;
}

bool ISOVERT::isActive(const int cube_index) const
{
  if (index_grid.Scalar(cube_index) != NO_INDEX)
    return true;
  else
    return false;
}

INDEX_DIFF_TYPE ISOVERT::GCubeIndex
(const int cube_index, IJK::ERROR & error) const
{
  const INDEX_DIFF_TYPE gcube_index = this->GCubeIndex(cube_index);

  if (gcube_index == NO_INDEX) {
    error.AddMessage("Programming error.  Cube ", cube_index,
                     " is not active.");
  }
  return(gcube_index);
}


// **************************************************
// ISOVERT_INFO member functions
// **************************************************

void ISOVERT_INFO::Clear()
{
  num_sharp_corners = 0;
  num_sharp_edges = 0;
  num_smooth_vertices = 0;
  num_merged_iso_vertices = 0;
  num_conflicts = 0;
  num_Linf_iso_vertex_locations = 0;
}

// **************************************************
// Set ISOVERT_INFO
// **************************************************

/// Count number of vertices on sharp corners or sharp edges.
/// Count number of smooth vertices.
void MERGESHARP::count_vertices
(const ISOVERT & isovert, ISOVERT_INFO & isovert_info)
{
  isovert_info.num_sharp_corners = 0;
  isovert_info.num_sharp_edges = 0;
  isovert_info.num_smooth_vertices = 0;
  for (VERTEX_INDEX i = 0; i < isovert.gcube_list.size(); i++) {
    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
      if (isovert.gcube_list[i].num_eigenvalues == 2)
        { isovert_info.num_sharp_edges++; }
      else if (isovert.gcube_list[i].num_eigenvalues == 3)
        { isovert_info.num_sharp_corners++; }
    }
    else if (isovert.gcube_list[i].flag == SMOOTH_GCUBE) {
      isovert_info.num_smooth_vertices++;
    }
  }
}


// **************************************************
// BIN_GRID ROUTINES
// **************************************************

// local
namespace {

  inline void divide_coord_3D
  (const AXIS_SIZE_TYPE bin_width, GRID_COORD_TYPE coord[DIM3])
  {
    for (int d = 0; d < DIM3; d++)
      coord[d] = IJK::integer_divide(coord[d], bin_width);
  }

}

/// Initialize bin_grid.
/// @param bin_width = number of cubes along each axis.
void MERGESHARP::init_bin_grid
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 BIN_GRID<VERTEX_INDEX> & bin_grid)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);

  for (int d = 0; d < dimension; d++) {
    axis_size[d] = IJK::compute_subsample_size(grid.AxisSize(d), bin_width);
  }

  //bin_grid.SetSize(dimension, axis_size);
  bin_grid.SetSize(dimension, axis_size.PtrConst());

}


void MERGESHARP::bin_grid_insert
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid)
{
  static GRID_COORD_TYPE coord[DIM3];

  grid.ComputeCoord(cube_index, coord);
  divide_coord_3D(bin_width, coord);
  VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
  bin_grid.Insert(ibin, cube_index);
}

void MERGESHARP::bin_grid_remove
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid)
{
  static GRID_COORD_TYPE coord[DIM3];

  grid.ComputeCoord(cube_index, coord);
  divide_coord_3D(bin_width, coord);
  VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
  bin_grid.Remove(ibin, cube_index);
}

/// Get the selected vertices around iv.
void MERGESHARP::get_selected
(const SHARPISO_GRID & grid,
 const VERTEX_INDEX iv,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 std::vector<VERTEX_INDEX> & selected_list)
{
  const int dimension = grid.Dimension();
  static GRID_COORD_TYPE coord[DIM3];
  static GRID_COORD_TYPE min_coord[DIM3];
  static GRID_COORD_TYPE max_coord[DIM3];

  long boundary_bits;

  grid.ComputeCoord(iv, coord);
  divide_coord_3D(bin_width, coord);
  VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
  bin_grid.ComputeBoundaryBits(ibin, boundary_bits);

  selected_list.resize(bin_grid.ListLength(ibin));
  std::copy(bin_grid(ibin).list.begin(), bin_grid(ibin).list.end(),
            selected_list.begin());

  if (boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < bin_grid.NumVertexNeighborsC(); k++) {
      VERTEX_INDEX jbin = bin_grid.VertexNeighborC(ibin, k);
      for (int i = 0; i < bin_grid.ListLength(jbin); i++) {
        selected_list.push_back(bin_grid.List(jbin, i)); 
      }
    }
  }
  else {

    for (int d = 0; d < dimension; d++) {
      if (coord[d] > 0)
        { min_coord[d] = coord[d] - 1; }
      else
        { min_coord[d] = 0; }

      if (coord[d]+1 < bin_grid.AxisSize(d))
        { max_coord[d] = coord[d] + 1; }
      else
        { max_coord[d] = coord[d]; }
    }

    for (coord[0] = min_coord[0]; coord[0] <= max_coord[0];
         coord[0]++)
      for (coord[1] = min_coord[1]; coord[1] <= max_coord[1];
           coord[1]++)
        for (coord[2] = min_coord[2]; coord[2] <= max_coord[2];
             coord[2]++) {

          VERTEX_INDEX jbin = bin_grid.ComputeVertexIndex(coord);
          if (ibin != jbin) {
            for (int i = 0; i < bin_grid.ListLength(jbin); i++)
              { selected_list.push_back(bin_grid.List(jbin, i)); }
          }
        }
  }
}


// **************************************************
// Local routines
// **************************************************

namespace {

  // Compute rescaled L-infinity distance.
  // Compute L-infinity distance where 
  //   coord[d] is mapped to coord[d]/spacing[d],
  //   i.e., where grid cubes are mapped to unit cubes.
  void compute_rescaled_Linf_distance
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE spacing[DIM3], COORD_TYPE & linf_dist, int & axis)
  {
    linf_dist = 0;

    for (int d=0; d<DIM3; d++) {
      COORD_TYPE x = abs(coord0[d]-coord1[d])/spacing[d];
      if (x > linf_dist) {
        linf_dist = x;
        axis = d;
      }
    }
  }

  // Compute rescaled L-infinity distance.
  // Version which does not return axis.
  //   coord[d] is mapped to coord[d]/spacing[d],
  //   i.e., where grid cubes are mapped to unit cubes.
  void compute_rescaled_Linf_distance
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE spacing[DIM3], COORD_TYPE & linf_dist)
  {
    int axis;
    compute_rescaled_Linf_distance(coord0, coord1, spacing, linf_dist, axis);
  }

  // Compute the scaled L-infinity distance from the cube center.
  void compute_Linf_distance_from_cube_center
  (const SHARPISO_GRID & grid, const VERTEX_INDEX icube,
   const COORD_TYPE coord[DIM3], COORD_TYPE & linf_dist)
  {
    COORD_TYPE cc_coord[DIM3];

    grid.ComputeCubeCenterScaledCoord(icube, cc_coord);

    linf_dist = 0;
    for (int d=0; d<DIM3; d++) {
      COORD_TYPE x = abs(coord[d]-cc_coord[d])/grid.Spacing(d);
      if (x > linf_dist)
        { linf_dist = x; }
    }
  }

  // Compute the scaled L-infinity distance from the cube center.
  // Store in isovert.gcube_list[gcube_index].
  void compute_Linf_distance_from_cube_center
  (const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert)
  {
    const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    compute_Linf_distance_from_cube_center
      (grid, cube_index, isovert.gcube_list[gcube_index].isovert_coord,
       isovert.gcube_list[gcube_index].linf_dist);
  }

  // Compute the scaled L-infinity distance from vertex iv.
  void compute_Linf_distance_from_grid_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE coord[DIM3], COORD_TYPE & linf_dist, int & axis)
  {
    COORD_TYPE vcoord[DIM3];

    grid.ComputeCubeCenterScaledCoord(iv,vcoord);
    linf_dist = 0;

    for (int d=0; d<DIM3; d++) {
      COORD_TYPE x = abs(vcoord[d]-coord[d])/grid.Spacing(d);
      if (x > linf_dist) {
        linf_dist = x;
        axis = d;
      }
    }
  }

  // Compute the scaled L-infinity distance from vertex iv.
  void compute_Linf_distance_from_grid_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE coord[DIM3], COORD_TYPE & linf_dist)
  {
    int axis;
    compute_Linf_distance_from_grid_vertex(grid, iv, coord, linf_dist, axis);
  }

  // Compute rescaled distance to line.
  void compute_rescaled_distance_to_line_3D
  (const COORD_TYPE coord0[DIM3], const COORD_TYPE coord1[DIM3],
   const COORD_TYPE direction[DIM3], const COORD_TYPE spacing[DIM3], 
   COORD_TYPE & distance)
  {
    COORD_TYPE rescaled_coord0[DIM3];
    COORD_TYPE rescaled_coord1[DIM3];
    COORD_TYPE rescaled_direction[DIM3];

    for (int d = 0; d < DIM3; d++) {
      rescaled_coord0[d] = coord0[d]/spacing[d];
      rescaled_coord1[d] = coord1[d]/spacing[d];
      rescaled_direction[d] = direction[d]/spacing[d];
    }

    IJK::normalize_vector(DIM3, rescaled_direction, 0, rescaled_direction);
    IJK::compute_distance_to_line_3D
      (rescaled_coord0, rescaled_coord1, rescaled_direction, distance);
  }

  // Return true if two  cube-indices are connected
  bool are_connected 
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX & cube_index1,
   const VERTEX_INDEX & cube_index2,
   const SCALAR_TYPE isovalue )
  {
    // find the overlap region
    GRID_COORD_TYPE rmin[DIM3], rmax[DIM3];

    int overlap_dim;
    find_3x3x3_overlap
      (scalar_grid, cube_index1, cube_index2, rmin, rmax, overlap_dim);

    if (overlap_dim < 2) { return false; }

    COORD_TYPE vbase = scalar_grid.ComputeVertexIndex(rmin);
    COORD_TYPE base_coord[DIM3];
    VERTEX_INDEX base_index;
    VERTEX_INDEX vnext;
    COORD_TYPE diff[DIM3]={0.0,0.0,0.0};
    //visit each edge along the way
    VERTEX_INDEX v0,v1,v2;

    int d0=0,d1=0,d2=0;
    bool is_intersect=false;

    int num=0;
    // if there is an overlap
    for(int d0=0;d0<DIM3;d0++){
      d1 = (d0+1)%3;
      d2 = (d0+2)%3;
      v1 = vbase;
      for (int i1=rmin[d1]; i1 <=rmax[d1];i1++){
        v2=v1;
        for (int i2=rmin[d2];i2<=rmax[d2];i2++){
          v0=v2;
          for(int i0=rmin[d0];i0<rmax[d0];i0++){
            vnext = scalar_grid.NextVertex(v0,d0);
            process_edge(scalar_grid, v0,vnext,isovalue, is_intersect);
            num++;
            if (is_intersect)
              return true;
            v0=vnext;
          }
          v2=scalar_grid.NextVertex(v2,d2);
        }
        v1=scalar_grid.NextVertex(v1,d1);
      }
    }

    return false;
  }

  // Compute the cosine of the angle between (v2,v1) and (v2,v3)
  void compute_cos_angle
  (const ISOVERT & isovert,
   const VERTEX_INDEX gcube_list_index_v1,
   const VERTEX_INDEX gcube_list_index_v2,
   const VERTEX_INDEX gcube_list_index_v3,
   SCALAR_TYPE & cos_angle)
  {
    COORD_TYPE coord0[DIM3], coord1[DIM3],
      coord2[DIM3],vec12[DIM3], vec32[DIM3];

    subtract_coord(DIM3, isovert.gcube_list[gcube_list_index_v2].isovert_coord,
                   isovert.gcube_list[gcube_list_index_v1].isovert_coord, vec12);
    normalize_vector(DIM3, vec12, 0.001, vec12);

    subtract_coord(DIM3, isovert.gcube_list[gcube_list_index_v2].isovert_coord,
                   isovert.gcube_list[gcube_list_index_v3].isovert_coord, vec32);
    normalize_vector(DIM3, vec32, 0.001, vec32);

    compute_inner_product(DIM3, vec12, vec32, cos_angle);
  }


  // Return true if triangle(iv,v1,v2) has a large angle
  // iv, v1 and v2 are cube_indices.
  bool is_angle_large
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const ISOVERT &isovert,
   const VERTEX_INDEX iv,
   const SCALAR_TYPE threshold,
   const VERTEX_INDEX v1,
   const VERTEX_INDEX v2)
  {

    SCALAR_TYPE cos_angle_iv, cos_angle_v1, cos_angle_v2;
    VERTEX_INDEX gcube_list_index_v1;
    VERTEX_INDEX gcube_list_index_v2;
    VERTEX_INDEX gcube_list_index_v3;

    gcube_list_index_v1=isovert.index_grid.Scalar(v1);
    gcube_list_index_v2=isovert.index_grid.Scalar(v2);
    gcube_list_index_v3=isovert.index_grid.Scalar(iv);


    compute_cos_angle(isovert, gcube_list_index_v1,
                      gcube_list_index_v3,  gcube_list_index_v2, cos_angle_iv);

    compute_cos_angle(isovert, gcube_list_index_v3,
                      gcube_list_index_v1,  gcube_list_index_v2, cos_angle_v1);

    compute_cos_angle(isovert, gcube_list_index_v3,
                      gcube_list_index_v2,  gcube_list_index_v1, cos_angle_v2);

    if ((cos_angle_iv < threshold) || (cos_angle_v1 < threshold)
        || (cos_angle_v2 < threshold))
      {
        return true;
      }
    else
      return false;
  }

  // Return true if some facet neighbor is covered by a sharp cube.
  bool is_facet_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, VERTEX_INDEX & cube_index1)
  {
    cube_index1 = 0;
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborF(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.index_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        if (isovert.gcube_list[index_gcube].covered_by != icube) {
          cube_index1 = isovert.gcube_list[index_gcube].covered_by;
          return(true);
        }
      }
    }
    return(false);
  }

  // Return true if some facet or edge neighbor is covered by a sharp cube
  //   other than cube_index1.
  // @pre cube_index0 is not on grid boundary.
  bool is_facet_or_edge_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   VERTEX_INDEX & cube_index2)
  {
    cube_index2 = 0;

    // Check facet neighbors
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborF(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.index_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1) {
          cube_index2 = covered_by;
          return(true);
        }
      }
    }

    // Check edge neighbors
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborE(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.index_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1) {
          cube_index2 = covered_by;
          return(true);
        }
      }
    }

    return(false);
  }

  // Return true if some facet, edge or vertex neighbor is covered 
  //   by a sharp cube other than cube_index1 or cube_index2
  // @pre cube_index0 is not on grid boundary.
  bool is_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   const VERTEX_INDEX & cube_index2)
  {
    for (NUM_TYPE j = 0; j < grid.NumVertexNeighborsC(); j++) {
      VERTEX_INDEX icube = grid.VertexNeighborC(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.index_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1 && 
            covered_by != cube_index2) {
          return(true);
        }
      }
    }
    return(false);
  }

  // Return true if cube or some cube neighbor of cube_index0 
  //   is covered by cube_index1
  bool is_cube_or_neighbor_covered_by
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1)
  {
    VERTEX_INDEX gcube_index;

    gcube_index = isovert.index_grid.Scalar(cube_index0);
    if (gcube_index != ISOVERT::NO_INDEX) {
      VERTEX_INDEX covered_by = isovert.gcube_list[gcube_index].covered_by;
      if (covered_by == cube_index1) { return(true); }
    }

    if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

      for (NUM_TYPE j = 0; j < grid.NumVertexNeighborsC(); j++) {
        VERTEX_INDEX icube = grid.VertexNeighborC(cube_index0, j);

        gcube_index = isovert.index_grid.Scalar(icube);

        if (gcube_index != ISOVERT::NO_INDEX) {
          VERTEX_INDEX covered_by = isovert.gcube_list[gcube_index].covered_by;
          if (covered_by == cube_index1) { return(true); }
        }
      }
    }
    else {
      // Handle boundary case
    }

    return(false);
  }

  // Store information from singular valued decomposition.
  void store_svd_info
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const NUM_TYPE gcube_index,
   const NUM_TYPE num_large_eigenvalues,
   const SVD_INFO svd_info,
   ISOVERT & isovert)
  {
    // set number of eigenvalues.
    isovert.gcube_list[gcube_index].num_eigenvalues =
      (unsigned char) num_large_eigenvalues;

    // set the sharp vertex type to be *AVAILABLE*.
    if (num_large_eigenvalues > 1 && svd_info.location == LOC_SVD)
      { isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE; }
    else 
      { isovert.gcube_list[gcube_index].flag = SMOOTH_GCUBE; }

    if (svd_info.location == CENTROID) 
      { isovert.gcube_list[gcube_index].flag_centroid_location = true; }
    else
      { isovert.gcube_list[gcube_index].flag_centroid_location = false; }

    isovert.gcube_list[gcube_index].flag_conflict = svd_info.flag_conflict;
    isovert.gcube_list[gcube_index].flag_far = svd_info.flag_far;

    // Initialize
    isovert.gcube_list[gcube_index].flag_using_substitute_coord = false;
    isovert.gcube_list[gcube_index].flag_coord_from_other_cube = false;

    // store distance.
    compute_Linf_distance_from_cube_center
      (scalar_grid, cube_index, isovert.gcube_list[gcube_index].isovert_coord,
       isovert.gcube_list[gcube_index].linf_dist);

    MSDEBUG();
    if (flag_debug) {
      scalar_grid.PrintIndexAndCoord
        (cerr, "*** Storing svd for cube: ", cube_index, "\n");
      cerr << "  isovert_coord: ";
      IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coord);
      cerr << "  isovert_coordB: ";
      IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coordB);
      cerr << "  Linf dist: " << isovert.gcube_list[gcube_index].linf_dist << endl;
    }
  }

  // Set cube_containing_isovert
  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, const NUM_TYPE gcube_index,
   ISOVERT & isovert)
  {
    COORD_TYPE cube_coord[DIM3];
    VERTEX_INDEX cube_containing_point;

    const COORD_TYPE * point_coord = 
      isovert.gcube_list[gcube_index].isovert_coord;
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    scalar_grid.ComputeScaledCoord(cube_index, cube_coord);

    if (cube_contains_point(cube_coord, point_coord, 
                            scalar_grid.SpacingPtrConst()) ||
        // *** SCALING DEPENDENT ***
        !scalar_grid.ContainsPoint(point_coord)) {

      isovert.gcube_list[gcube_index].cube_containing_isovert = cube_index;
      isovert.gcube_list[gcube_index].flag_conflict = false;
    }
    else {

      bool flag_boundary;
      bool flag_active;

      get_cube_containing_point
        (scalar_grid, isovalue, point_coord, cube_containing_point,
         flag_boundary, flag_active);
      
      isovert.gcube_list[gcube_index].cube_containing_isovert = 
        cube_containing_point;

      if (flag_active)
        { isovert.gcube_list[gcube_index].flag_conflict = true; }
      else if (!flag_boundary) 
        { isovert.gcube_list[gcube_index].flag_conflict = false; }
      else {

        // Cube cube_containing_point is not active but
        //   some other cube may also contain the point.
        // Check whether any cube containing the point is active.
        VERTEX_INDEX conflicting_cube;
        if (check_conflict(scalar_grid, isovalue, cube_coord, 
                           isovert.gcube_list[gcube_index].isovert_coord, 
                           conflicting_cube)) {
            
          isovert.gcube_list[gcube_index].cube_containing_isovert = 
            conflicting_cube;
          isovert.gcube_list[gcube_index].flag_conflict = true;
        }
        else {
          isovert.gcube_list[gcube_index].flag_conflict = false;
        }
      }
    }
  }

  // Set cube_containing_isovert
  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, ISOVERT & isovert)
  {
    COORD_TYPE cube_coord[DIM3];
    VERTEX_INDEX cube_containing_point;

    for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
      set_cube_containing_isovert(scalar_grid, isovalue, i, isovert);
    }
  }

  /// Snap pointA to point in cube cube_index.
  void snap2cube
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index, 
   const COORD_TYPE pointA[DIM3], COORD_TYPE * pointB)
  {
    COORD_TYPE cube_coord[DIM3];

    grid.ComputeScaledCoord(cube_index, cube_coord);

    for (int d = 0; d < DIM3; d++) {
      pointB[d] = pointA[d];
      if (pointB[d] < cube_coord[d]) { pointB[d] = cube_coord[d]; }
      else {
        COORD_TYPE x = cube_coord[d] + grid.Spacing(d);
        if (pointB[d] > x) { pointB[d] = x; }
      }
    }
  }

  /// Snap pointA to vertex of cube cube_index.
  void snap2cube_vertex
  (const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index, 
   const COORD_TYPE pointA[DIM3], COORD_TYPE * pointB)
  {
    COORD_TYPE cube_coord[DIM3];

    grid.ComputeScaledCoord(cube_index, cube_coord);

    for (int d = 0; d < DIM3; d++) {

      COORD_TYPE x0 = cube_coord[d];
      COORD_TYPE x1 = cube_coord[d]+grid.Spacing(d);
      COORD_TYPE midx = (x0+x1)/2.0;

      if (pointA[d] <= midx) 
        { pointB[d] = x0; }
      else 
        { pointB[d] = x1; }
    }
  }

  /// If sharp vertex is not contained in the cube,
  ///   replace with substitute coord isovert_coordB.
  void check_not_contained_and_substitute
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const NUM_TYPE gcube_index,
   ISOVERT & isovert)
  {
    const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    const COORD_TYPE * isovert_coord =
      isovert.gcube_list[gcube_index].isovert_coord;
    const COORD_TYPE * isovert_coordB =
      isovert.gcube_list[gcube_index].isovert_coordB;

    if (!cube_contains_point(scalar_grid, cube_index, isovert_coord)) {
      if (cube_contains_point(scalar_grid, cube_index, isovert_coordB)) {
        replace_with_substitute_coord(scalar_grid, gcube_index, isovert);
      }
    }
  }

  /// If sharp vertex is in a cube which is already covered,
  ///   replace with substitute coord isovert_coordB.
  void check_covered_and_substitute
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SHARPISO_BOOL_GRID & covered_grid,
   const SCALAR_TYPE isovalue,
   const NUM_TYPE gcube_index,
   ISOVERT & isovert)
  {
    if (check_covered_point(covered_grid, isovert, gcube_index)) {
      // *** CHECK THAT SUBSTITUTE IS AN IMPROVEMENT ***
      replace_with_substitute_coord(scalar_grid, gcube_index, isovert);
      set_cube_containing_isovert(scalar_grid, isovalue, gcube_index, isovert);
      check_and_set_covered_point(covered_grid, isovert, gcube_index);
    }
  }


  // *** MOVED TO mergesharp_select
  /// Return true if sharp vertex is in a cube which is covered.
  bool check_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   const ISOVERT & isovert,
   const VERTEX_INDEX & gcube0_index)
  {
    const VERTEX_INDEX cube1_index =
      isovert.gcube_list[gcube0_index].cube_containing_isovert;
  
    if (covered_grid.Scalar(cube1_index)) { return true; }

    // *** IS THIS NECESSARY? ***
    // *** SHOULD USE IsCoveredOrSelected() ***
    // if the sharp vertex is not in the same cube. 
    if (isovert.GCubeIndex(cube1_index) != ISOVERT::NO_INDEX) {
      if (isovert.gcube_list[gcube0_index].cube_index != cube1_index){
        if (isovert.isFlag(cube1_index, COVERED_A_GCUBE))
          { return true; }
      }
    }

    return false;
  }


  /// If sharp vertex is in covered cube, set cube to COVERED_POINT.
  void check_and_set_covered_point
  (const SHARPISO_BOOL_GRID & covered_grid,
   ISOVERT & isovert,
   const VERTEX_INDEX & gcube_index)
  {
    if (isovert.gcube_list[gcube_index].flag == AVAILABLE_GCUBE) {

      if (isovert.gcube_list[gcube_index].cube_index !=
          isovert.gcube_list[gcube_index].cube_containing_isovert) {

        if (check_covered_point
            (covered_grid, isovert, gcube_index)) {

          isovert.gcube_list[gcube_index].flag = COVERED_POINT;
        }
      }
    }
  }

  void set_covered_grid
  (const ISOVERT & isovert, SHARPISO_BOOL_GRID & covered_grid)
  {
    SHARPISO_GRID_NEIGHBORS grid;
    BOUNDARY_BITS_TYPE mask[2*DIM3];

    grid.SetSize(covered_grid);

    for (int i = 0; i < 2*DIM3; i++) 
      { mask[i] = (BOUNDARY_BITS_TYPE(1) << i); }

    covered_grid.SetAll(false);
    for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

      const VERTEX_INDEX cube_index = isovert.CubeIndex(i);

      covered_grid.Set(cube_index, true);

      const BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[i].boundary_bits;

      if (boundary_bits == 0) {

        // mark all the neighbors as covered
        for (int i=0; i < grid.NumVertexNeighborsC(); i++) {
          VERTEX_INDEX cube2_index = grid.VertexNeighborC(cube_index, i);

          covered_grid.Set(cube2_index, true);
        }

      }
      else {
        // Handle boundary case.
        for (int j0 = -1; j0 < 2; j0++) {

          if ((boundary_bits & mask[0]) && j0 < 0) { continue; }
          if ((boundary_bits & mask[1]) && j0 > 0) { continue; }

          for (int j1 = -1; j1 < 2; j1++) {

            if ((boundary_bits & mask[2]) && j1 < 0) { continue; }
            if ((boundary_bits & mask[3]) && j1 > 0) { continue; }

            for (int j2 = -1; j2 < 2; j2++) {

              if ((boundary_bits & mask[4]) && j2 < 0) { continue; }
              if ((boundary_bits & mask[5]) && j2 > 0) { continue; }

              VERTEX_INDEX cube2_index =
                cube_index + j0*grid.AxisIncrement(0) +
                j1*grid.AxisIncrement(1) + j2*grid.AxisIncrement(2);

              covered_grid.Set(cube2_index, true);
            }
          }
        }
      }
    }
  }

}

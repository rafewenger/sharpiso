/// \file sharpiso_get_gradients.cxx
/// Get gradients in cube and cube neighbors.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
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


#include "sharpiso_get_gradients.h"

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_scalar.txx"


// **************************************************
// GET GRADIENTS
// **************************************************

// local namespace
namespace {

  using namespace SHARPISO;

  inline void add_gradient
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX iv,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients)
  {
    NUM_TYPE ic = point_coord.size();
    point_coord.resize(ic+DIM3);
    gradient_grid.ComputeCoord(iv, &(point_coord[ic]));

    gradient_coord.resize(ic+DIM3);
    std::copy(gradient_grid.VectorPtrConst(iv),
              gradient_grid.VectorPtrConst(iv)+DIM3,
              &(gradient_coord[ic]));

    scalar.push_back(scalar_grid.Scalar(iv));

    num_gradients++;
  }

  inline void add_large_gradient
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX iv,
   const GRADIENT_COORD_TYPE max_small_mag_squared,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients)
  {
    GRADIENT_COORD_TYPE magnitude_squared =
      gradient_grid.ComputeMagnitudeSquared(iv);

    if (magnitude_squared > max_small_mag_squared) {
      add_gradient(scalar_grid, gradient_grid, iv,
                   point_coord, gradient_coord, scalar, num_gradients);
    }
  }

  inline void add_selected_gradient
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX iv,
   const GRID_COORD_TYPE * cube_coord,
   const GRADIENT_COORD_TYPE max_small_mag_squared,
   const SCALAR_TYPE isovalue,
   const OFFSET_CUBE_111 & cube_111,
   std::vector<COORD_TYPE> & point_coord,
   std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
   std::vector<SCALAR_TYPE> & scalar,
   NUM_TYPE & num_gradients)
  {
    typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

    // static so not reallocated at each call
    static GRID_COORD_TYPE vertex_coord[DIM3];
    static GRID_COORD_TYPE coord[DIM3];

    GRADIENT_COORD_TYPE magnitude_squared =
      gradient_grid.ComputeMagnitudeSquared(iv);

    if (magnitude_squared > max_small_mag_squared) {

      gradient_grid.ComputeCoord(iv, vertex_coord);
      // Add (1,1,1) since cube_111 has origin near (1,1,1)
      // Ensures that coord[] is not negative.
      for (DTYPE d = 0; d < DIM3; d++)
        { coord[d] = (vertex_coord[d]+1) - cube_coord[d]; }
      const GRADIENT_COORD_TYPE * vertex_gradient_coord =
        gradient_grid.VectorPtrConst(iv);
      SCALAR_TYPE s = scalar_grid.Scalar(iv);

      if (iso_intersects_cube
          (cube_111, coord, vertex_gradient_coord, s, isovalue)) {
        add_gradient(scalar_grid, gradient_grid, iv,
                     point_coord, gradient_coord, scalar, num_gradients);
      }
    }
  }

}

// Get selected grid vertex gradients.
/// @pre Size of vertex_flag[] is at least size of vertex_list[].
void SHARPISO::get_selected_vertex_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
 const bool vertex_flag[],
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  num_gradients = 0;

  for (NUM_TYPE i = 0; i < num_vertices; i++) {

    if (vertex_flag[i]) {

      VERTEX_INDEX iv = vertex_list[i];
      NUM_TYPE ic = point_coord.size();
      point_coord.resize(ic+DIM3);
      gradient_grid.ComputeCoord(iv, &(point_coord[ic]));

      gradient_coord.resize(ic+DIM3);
      std::copy(gradient_grid.VectorPtrConst(iv),
                gradient_grid.VectorPtrConst(iv)+DIM3,
                &(gradient_coord[ic]));

      scalar.push_back(scalar_grid.Scalar(iv));

      num_gradients++;
    }
  }
}

/// Get selected grid vertex gradients.
/// std::vector variation.
void SHARPISO::get_selected_vertex_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const std::vector<VERTEX_INDEX> & vertex_list,
 const bool vertex_flag[],
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  const NUM_TYPE num_vertices = vertex_list.size();

  if (num_vertices > 0) {

    get_selected_vertex_gradients
      (scalar_grid, gradient_grid, &(vertex_list[0]), num_vertices,
       vertex_flag, point_coord, gradient_coord, scalar, 
       num_gradients);
  }
  else {
    num_gradients = 0;
  }
}

// Get all 8 cube gradients
void SHARPISO::get_cube_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 std::vector<COORD_TYPE> & point_coord,
 GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3],
 SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D])
{

  for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) {
    VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
    scalar[k] = scalar_grid.Scalar(iv);
    IJK::copy_coord(DIM3, gradient_grid.VectorPtrConst(iv),
                    gradient_coord+k*DIM3);
  }
}

/// Get gradients.
/// @param sharp_isovert_param Parameters to determine
///   which gradients are selected.
void SHARPISO::get_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GET_GRADIENTS_PARAM & sharpiso_param,
 const OFFSET_CUBE_111 & cube_111,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  const GRADIENT_COORD_TYPE max_small_mag = sharpiso_param.max_small_magnitude;
  const GRADIENT_COORD_TYPE max_small_mag_squared = 
    max_small_mag * max_small_mag;
  const GRADIENT_COORD_TYPE zero_tolerance = sharpiso_param.zero_tolerance;

  // static so not reallocated at each call
  static GRID_COORD_TYPE cube_coord[DIM3];

  if (sharpiso_param.use_only_cube_gradients) {

    // NOTE: This vertex_list is an array.
    VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D];
    bool vertex_flag[NUM_CUBE_VERTICES3D];
    NUM_TYPE num_vertices(0);

    // initialize vertex_flag
    for (NUM_TYPE i = 0; i < NUM_CUBE_VERTICES3D; i++)
      { vertex_flag[i] = true; }

    if (sharpiso_param.use_intersected_edge_endpoint_gradients) {
      get_intersected_cube_edge_endpoints
        (scalar_grid, cube_index, isovalue, vertex_list, num_vertices);
    }
    else if (sharpiso_param.use_gradients_determining_edge_intersections) {
      get_cube_vertices_determining_edge_intersections
        (scalar_grid, gradient_grid, cube_index, isovalue, zero_tolerance,
         vertex_list, num_vertices);
    }
    else {
      get_cube_vertices(scalar_grid, cube_index, vertex_list);
      num_vertices = NUM_CUBE_VERTICES3D;
    }

    deselect_vertices_with_small_gradients
      (gradient_grid, vertex_list, num_vertices, max_small_mag_squared,
       vertex_flag);

    if (sharpiso_param.use_selected_gradients) {

      scalar_grid.ComputeCoord(cube_index, cube_coord);

      deselect_vertices_based_on_isoplanes
        (scalar_grid, gradient_grid, cube_coord, cube_111,
         isovalue, vertex_list, num_vertices, vertex_flag);
    }

    get_selected_vertex_gradients
      (scalar_grid, gradient_grid, vertex_list, num_vertices, vertex_flag,
       point_coord, gradient_coord, scalar, num_gradients);
  }
  else {

    // NOTE: This vertex_list is a C++ vector.
    std::vector<VERTEX_INDEX> vertex_list;

    if (sharpiso_param.use_intersected_edge_endpoint_gradients) {
      get_intersected_cube_neighbor_edge_endpoints
        (scalar_grid, cube_index, isovalue, vertex_list);
    }
    else {
      get_cube_neighbor_vertices(scalar_grid, cube_index, vertex_list);
    }

    IJK::ARRAY<bool> vertex_flag(vertex_list.size());

    // initialize vertex_flag
    for (NUM_TYPE i = 0; i < vertex_list.size(); i++)
      { vertex_flag[i] = true; }

    deselect_vertices_with_small_gradients
      (gradient_grid, vertex_list, max_small_mag_squared, 
       vertex_flag.Ptr());

    if (sharpiso_param.use_selected_gradients) {

      scalar_grid.ComputeCoord(cube_index, cube_coord);

      deselect_vertices_based_on_isoplanes
        (scalar_grid, gradient_grid, cube_coord, cube_111,
         isovalue, vertex_list, vertex_flag.Ptr());
    }

    get_selected_vertex_gradients
      (scalar_grid, gradient_grid, vertex_list, vertex_flag.Ptr(),
       point_coord, gradient_coord, scalar, num_gradients);
  }
}

void SHARPISO::get_large_cube_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const GRADIENT_COORD_TYPE max_small_mag,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  IJK::PROCEDURE_ERROR error("get_large_cube_gradients");

  // Initialize num_gradients
  num_gradients = 0;

  for (NUM_TYPE k = 0; k < scalar_grid.NumCubeVertices(); k++) {
    VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
    add_large_gradient
      (scalar_grid, gradient_grid, iv, max_small_mag_squared,
       point_coord, gradient_coord, scalar, num_gradients);
  }

}

void SHARPISO::get_large_cube_neighbor_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const GRADIENT_COORD_TYPE max_small_mag,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
  IJK::PROCEDURE_ERROR error("get_large_cube_neighbor_gradients");

  // Initialize num_gradients
  num_gradients = 0;

  scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

  get_large_cube_gradients
    (scalar_grid, gradient_grid, cube_index, max_small_mag,
     point_coord, gradient_coord, scalar, num_gradients);

  for (DTYPE d = 0; d < DIM3; d++) {

    if (cube_coord[d] > 0) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
        add_large_gradient
          (scalar_grid, gradient_grid, iv0, max_small_mag_squared,
           point_coord, gradient_coord, scalar, num_gradients);
      }

    }

    if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv2 = iv1 + 2*scalar_grid.AxisIncrement(d);
        add_large_gradient
          (scalar_grid, gradient_grid, iv2, max_small_mag_squared,
           point_coord, gradient_coord, scalar, num_gradients);
      }

    }

  }

}

void SHARPISO::get_selected_cube_neighbor_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients,
 const OFFSET_CUBE_111 & cube_111)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  IJK::ARRAY<GRID_COORD_TYPE> offset_111(DIM3, 1);
  IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
  IJK::ARRAY<GRID_COORD_TYPE> vertex_coord(DIM3);
  IJK::ARRAY<COORD_TYPE> coord(DIM3);
  IJK::ARRAY<COORD_TYPE> cube_diagonal_coord(DIM3*NUM_CUBE_VERTICES3D);
  IJK::PROCEDURE_ERROR error("get_large_cube_neighbor_gradients");

  // Initialize num_gradients
  num_gradients = 0;

  scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

  select_cube_gradients_based_on_isoplanes
    (scalar_grid, gradient_grid, cube_index, max_small_mag, isovalue,
     point_coord, gradient_coord, scalar, num_gradients, cube_111);

  for (DTYPE d = 0; d < DIM3; d++) {

    if (cube_coord[d] > 0) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);

        add_selected_gradient
          (scalar_grid, gradient_grid, iv0, cube_coord.PtrConst(),
           max_small_mag_squared, isovalue, cube_111,
           point_coord, gradient_coord, scalar, num_gradients);
      }

    }

    if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv2 = iv1 + 2*scalar_grid.AxisIncrement(d);

        add_selected_gradient
          (scalar_grid, gradient_grid, iv2, cube_coord.PtrConst(),
           max_small_mag_squared, isovalue, cube_111,
           point_coord, gradient_coord, scalar, num_gradients);
      }

    }
  }

}

namespace {
  int axis_size_222[DIM3] = { 2, 2, 2 };
  SHARPISO_GRID grid_222(DIM3, axis_size_222);

  void flag_intersected_cube_edge_endpoints
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   bool corner_flag[NUM_CUBE_VERTICES3D])
  {
    typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    for (VERTEX_INDEX i = 0; i < NUM_CUBE_VERTICES3D; i++)
      { corner_flag[i] = false; }

    for (DTYPE d = 0; d < dimension; d++) {
      for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
        VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {
    
          VERTEX_INDEX icorner0 = grid_222.FacetVertex(0, d, k);
          VERTEX_INDEX icorner1 = grid_222.NextVertex(icorner0, d);

          corner_flag[icorner0] = true;
          corner_flag[icorner1] = true;
        }
      }
    }

  }

}

void SHARPISO::get_intersected_edge_endpoint_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const DTYPE dimension = scalar_grid.Dimension();
  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  bool corner_flag[NUM_CUBE_VERTICES3D];

  flag_intersected_cube_edge_endpoints
    (scalar_grid, cube_index, isovalue, corner_flag);

  for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
    if (corner_flag[j]) {
      // grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
      VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
      VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

      add_large_gradient
        (scalar_grid, gradient_grid, iv, max_small_mag_squared, 
         point_coord, gradient_coord, scalar, num_gradients);
    }
  }

}


namespace {

  /// Get vertex whose gradient determines intersection of isosurface 
///    and edge (iv0,iv1)
/// @param zero_tolerance No division by numbers less than or equal 
///        to zero_tolerance.
/// @pre zero_tolerance must be non-negative.
void get_gradient_determining_edge_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
 const GRADIENT_COORD_TYPE zero_tolerance,
 VERTEX_INDEX & iv2)
{
  const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
  const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

  const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
  const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
  const GRADIENT_COORD_TYPE gdiff = g0 - g1;
    
  iv2 = iv0;
  if (abs(gdiff) <= zero_tolerance) { return; }

  const SCALAR_TYPE s2 = g0*(g1+s1-s0)/gdiff;

  if (s0 <= s1) {
    if (isovalue < s2) { iv2 = iv0; }
    else { iv2 = iv1; }
  }
  else {
    if (isovalue < s2) { iv2 = iv1; }
    else { iv2 = iv0; }
  }
}

  void flag_cube_gradients_determining_edge_intersections
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   const GRADIENT_COORD_TYPE zero_tolerance,
   bool corner_flag[NUM_CUBE_VERTICES3D])
  {
    typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    for (VERTEX_INDEX i = 0; i < NUM_CUBE_VERTICES3D; i++)
      { corner_flag[i] = false; }

    for (DTYPE d = 0; d < dimension; d++) {
      for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
        VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {
    
          VERTEX_INDEX icorner0 = grid_222.FacetVertex(0, d, k);
          VERTEX_INDEX icorner1 = grid_222.NextVertex(icorner0, d);

          VERTEX_INDEX iv2;
          get_gradient_determining_edge_intersection
            (scalar_grid, gradient_grid, isovalue, iv0, iv1, d,
             zero_tolerance, iv2);

          if (iv2 == iv0) { corner_flag[icorner0] = true; }
          else { corner_flag[icorner1] = true; }
        }
      }
    }

  }

}

/// Get gradients of vertices which determine edge isosurface intersections.
/// @param zero_tolerance No division by numbers less than or equal 
///        to zero_tolerance.
/// @pre zero_tolerance must be non-negative.
void SHARPISO::get_gradients_determining_edge_intersections
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE zero_tolerance,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const DTYPE dimension = scalar_grid.Dimension();
  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  bool corner_flag[NUM_CUBE_VERTICES3D];

  flag_cube_gradients_determining_edge_intersections
    (scalar_grid, gradient_grid, cube_index, isovalue, zero_tolerance,
     corner_flag);

  for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
    if (corner_flag[j]) {
      // grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
      VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
      VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

      add_large_gradient
        (scalar_grid, gradient_grid, iv, max_small_mag_squared, 
         point_coord, gradient_coord, scalar, num_gradients);
    }
  }

}

namespace {

  int axis_size_444[DIM3] = { 4, 4, 4 };
  SHARPISO_GRID grid_444(DIM3, axis_size_444);

  void add_cube_flags_to_grid_444_flags
  (const bool cube_corner_flag[NUM_CUBE_VERTICES3D],
   const VERTEX_INDEX iw0,
   bool * vertex_flag)
  {
    for (VERTEX_INDEX k = 0; k < NUM_CUBE_VERTICES3D; k++)
      if (cube_corner_flag[k]) {
        VERTEX_INDEX iw = grid_444.CubeVertex(iw0, k);
        vertex_flag[iw] = true;
      }
  }

  void flag_intersected_neighbor_edge_endpoints
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const SCALAR_TYPE isovalue,
   bool * vertex_flag)
  {
    typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
    bool corner_flag[NUM_CUBE_VERTICES3D];

    // Index of vertex in grid_444 with coordinates (1,1,1).
    const VERTEX_INDEX iw_111 = 21;

    for (VERTEX_INDEX i = 0; i < grid_444.NumVertices(); i++)
      { vertex_flag[i] = false; }

    flag_intersected_cube_edge_endpoints
      (scalar_grid, cube_index, isovalue, corner_flag);
    add_cube_flags_to_grid_444_flags(corner_flag, iw_111, vertex_flag);

    scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

    for (DTYPE d = 0; d < DIM3; d++) {

      if (cube_coord[d] > 0) {

        VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, d);
        VERTEX_INDEX iw0 = grid_444.PrevVertex(iw_111, d);

        // *** Unnecessarily recalculates some edge intersections.
        flag_intersected_cube_edge_endpoints
          (scalar_grid, iv0, isovalue, corner_flag);
        add_cube_flags_to_grid_444_flags(corner_flag, iw0, vertex_flag);
      }

      if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {
        VERTEX_INDEX iv0 = scalar_grid.NextVertex(cube_index, d);
        VERTEX_INDEX iw0 = grid_444.NextVertex(iw_111, d);

        // *** Unnecessarily recalculates some edge intersections.
        flag_intersected_cube_edge_endpoints
          (scalar_grid, iv0, isovalue, corner_flag);
        add_cube_flags_to_grid_444_flags(corner_flag, iw0, vertex_flag);

      }
    }


  }

}

void SHARPISO::get_intersected_neighbor_edge_endpoint_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const DTYPE dimension = scalar_grid.Dimension();
  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  bool vertex_flag[grid_444.NumVertices()];
  GRID_COORD_TYPE cube_coord[DIM3];
  GRID_COORD_TYPE coord_inc[DIM3];
  GRID_COORD_TYPE vcoord[DIM3];
  GRID_COORD_TYPE coord_111[DIM3] = { 1, 1, 1 };

  flag_intersected_neighbor_edge_endpoints
    (scalar_grid, cube_index, isovalue, vertex_flag);

  scalar_grid.ComputeCoord(cube_index, cube_coord);

  for (VERTEX_INDEX j = 0; j < grid_444.NumVertices(); j++) {
    if (vertex_flag[j]) {

      // *** SLOW COMPUTATION ***
      grid_444.ComputeCoord(j, coord_inc);
      IJK::add_coord_3D(cube_coord, coord_inc, vcoord);
      IJK::subtract_coord_3D(vcoord, coord_111, vcoord);
      VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(vcoord);

      add_large_gradient
        (scalar_grid, gradient_grid, iv, max_small_mag_squared, 
         point_coord, gradient_coord, scalar, num_gradients);
    }
  }

}


// **************************************************
// GET VERTICES
// **************************************************

// Get all 8 cube vertices.
void SHARPISO::get_cube_vertices
(const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
 VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D])
{

  for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) 
    { vertex_list[k] = grid.CubeVertex(cube_index, k); }
}


// Get vertices at endpoints of cube edges which intersect the isosurface.
void SHARPISO::get_intersected_cube_edge_endpoints
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const DTYPE dimension = scalar_grid.Dimension();
  bool corner_flag[NUM_CUBE_VERTICES3D];

  // Initialize.
  num_vertices = 0;

  flag_intersected_cube_edge_endpoints
    (scalar_grid, cube_index, isovalue, corner_flag);

  for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
    if (corner_flag[j]) {
      // grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
      VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
      VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

      vertex_list[num_vertices] = iv;
      num_vertices++;
    }
  }

}

// Get cube vertices which determining intersection of isosurface and edges.
void SHARPISO::get_cube_vertices_determining_edge_intersections
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue, const GRADIENT_COORD_TYPE zero_tolerance,
 VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  const DTYPE dimension = scalar_grid.Dimension();
  bool corner_flag[NUM_CUBE_VERTICES3D];

  // Initialize.
  num_vertices = 0;

  flag_cube_gradients_determining_edge_intersections
    (scalar_grid, gradient_grid, cube_index, isovalue, zero_tolerance,
     corner_flag);

  for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
    if (corner_flag[j]) {
      // grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
      VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
      VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

      vertex_list[num_vertices] = iv;
      num_vertices++;
    }
  }

}

/// Get vertices of cube and cube neighbors.
void SHARPISO::get_cube_neighbor_vertices
(const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & vertex_list)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  // static so not reallocated at each call
  static GRID_COORD_TYPE cube_coord[DIM3];

  vertex_list.resize(NUM_CUBE_VERTICES3D);
  get_cube_vertices(grid, cube_index, &vertex_list[0]);

  grid.ComputeCoord(cube_index, cube_coord);

  for (DTYPE d = 0; d < DIM3; d++) {

    if (cube_coord[d] > 0) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv0 = grid.PrevVertex(iv1, d);
        vertex_list.push_back(iv0);
      }

    }

    if (cube_coord[d]+2 < grid.AxisSize(d)) {

      for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
        VERTEX_INDEX iv1 = grid.FacetVertex(cube_index, d, k);
        VERTEX_INDEX iv2 = iv1 + 2*grid.AxisIncrement(d);
        vertex_list.push_back(iv2);
      }

    }

  }

}

// Get vertices at endpoints of cube and neighbor edges 
//   which intersect the isosurface.
void SHARPISO::get_intersected_cube_neighbor_edge_endpoints
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & vertex_list)
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  bool vertex_flag[grid_444.NumVertices()];
  GRID_COORD_TYPE cube_coord[DIM3];
  GRID_COORD_TYPE coord_inc[DIM3];
  GRID_COORD_TYPE vcoord[DIM3];
  GRID_COORD_TYPE coord_111[DIM3] = { 1, 1, 1 };

  flag_intersected_neighbor_edge_endpoints
    (scalar_grid, cube_index, isovalue, vertex_flag);

  scalar_grid.ComputeCoord(cube_index, cube_coord);

  for (VERTEX_INDEX j = 0; j < grid_444.NumVertices(); j++) {
    if (vertex_flag[j]) {

      // *** SLOW COMPUTATION ***
      grid_444.ComputeCoord(j, coord_inc);
      IJK::add_coord_3D(cube_coord, coord_inc, vcoord);
      IJK::subtract_coord_3D(vcoord, coord_111, vcoord);
      VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(vcoord);
      vertex_list.push_back(iv);
    }
  }

}

// **************************************************
// SELECTION FUNCTIONS
// **************************************************

// Select gradients at cube vertices.
// Select large gradients which give a level set intersecting the cube.
void SHARPISO::select_cube_gradients_based_on_isoplanes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const GRADIENT_COORD_TYPE max_small_mag,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients,
 const OFFSET_CUBE_111 & cube_111)
{
  const GRADIENT_COORD_TYPE max_small_mag_squared =
    max_small_mag * max_small_mag;
  IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
  IJK::PROCEDURE_ERROR error("get_large_cube_gradients");

  // Initialize num_gradients
  num_gradients = 0;

  scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

  for (NUM_TYPE k = 0; k < scalar_grid.NumCubeVertices(); k++) {
    VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
    add_selected_gradient
      (scalar_grid, gradient_grid, iv, cube_coord.PtrConst(),
       max_small_mag_squared, isovalue, cube_111,
       point_coord, gradient_coord, scalar, num_gradients);
  }

}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   with small gradient magnitude.
/// @pre Array vertex_flag[] is preallocated with size 
///      at least vertex_list.size().
void SHARPISO::deselect_vertices_with_small_gradients
(const GRADIENT_GRID_BASE & gradient_grid, 
 const VERTEX_INDEX * vertex_list, const NUM_TYPE num_vertices,
 const GRADIENT_COORD_TYPE max_small_mag_squared,
 bool vertex_flag[])
{
  for (NUM_TYPE i = 0; i < num_vertices; i++)
    if (vertex_flag[i]) {

      VERTEX_INDEX iv = vertex_list[i];
      GRADIENT_COORD_TYPE mag_squared = 
        gradient_grid.ComputeMagnitudeSquared(iv);

      if (mag_squared <= max_small_mag_squared) 
        { vertex_flag[i] = false; }
    }
}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   with small gradient magnitude.
/// std::vector variation.
void SHARPISO::deselect_vertices_with_small_gradients
(const GRADIENT_GRID_BASE & gradient_grid, 
 const std::vector<VERTEX_INDEX> & vertex_list,
 const GRADIENT_COORD_TYPE max_small_mag_squared,
 bool vertex_flag[])
{
  const NUM_TYPE num_vertices = vertex_list.size();

  if (num_vertices > 0) {

    deselect_vertices_with_small_gradients
      (gradient_grid,  &(vertex_list.front()), num_vertices, 
       max_small_mag_squared, vertex_flag);
  }
}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   determining an isoplane which does not intersect the cube.
/// @pre Array vertex_flag[] is preallocated with size 
///      at least vertex_list.size().
void SHARPISO::deselect_vertices_based_on_isoplanes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid, 
 const GRID_COORD_TYPE * cube_coord, const OFFSET_CUBE_111 & cube_111,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX * vertex_list, const NUM_TYPE num_vertices,
 bool vertex_flag[])
{
  typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

  // static so not reallocated at each call
  static GRID_COORD_TYPE vertex_coord[DIM3];
  static GRID_COORD_TYPE coord[DIM3];

  for (NUM_TYPE i = 0; i < num_vertices; i++) {

    if (vertex_flag[i]) {

      VERTEX_INDEX iv = vertex_list[i];

      gradient_grid.ComputeCoord(iv, vertex_coord);
      // Add (1,1,1) since cube_111 has origin near (1,1,1)
      // Ensures that coord[] is not negative.
      for (DTYPE d = 0; d < DIM3; d++)
        { coord[d] = (vertex_coord[d]+1) - cube_coord[d]; }
      const GRADIENT_COORD_TYPE * vertex_gradient_coord =
        gradient_grid.VectorPtrConst(iv);
      SCALAR_TYPE s = scalar_grid.Scalar(iv);

      if (!iso_intersects_cube
          (cube_111, coord, vertex_gradient_coord, s, isovalue)) 
        { vertex_flag[i] = false; }
    }
  }

}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   determining an isoplane which does not intersect the cube.
/// std::vector variation.
void SHARPISO::deselect_vertices_based_on_isoplanes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid, 
 const GRID_COORD_TYPE * cube_coord, const OFFSET_CUBE_111 & cube_111,
 const SCALAR_TYPE isovalue,
 const std::vector<VERTEX_INDEX> & vertex_list,
 bool vertex_flag[])
{
  const NUM_TYPE num_vertices = vertex_list.size();

  if (num_vertices > 0) {

    deselect_vertices_based_on_isoplanes
      (scalar_grid, gradient_grid, cube_coord, cube_111, isovalue,
       &(vertex_list.front()), num_vertices, vertex_flag);
  }
}

// **************************************************
// OFFSET_CUBE_111
// **************************************************

SHARPISO::OFFSET_CUBE_111::OFFSET_CUBE_111
(const SIGNED_COORD_TYPE offset)
{
  SetOffset(offset);
}

void SHARPISO::OFFSET_CUBE_111::SetOffset(const SIGNED_COORD_TYPE offset)
{
  IJK::PROCEDURE_ERROR error("OFFSET_CUBE_111::SetOffset");

  this->offset = 0;

  if (offset > 1) {
    error.AddMessage
      ("Programming error.  Offset must be less than or equal to 1.");
    error.AddMessage("  offset = ", offset, ".");
    throw error;
  }

  if (offset <= -1) {
    error.AddMessage
      ("Programming error.  Offset must be greater than -1.");
    error.AddMessage("  offset = ", offset, ".");
    throw error;
  }

  IJK::ARRAY<COORD_TYPE> v0_coord(this->Dimension(), 1-offset);
  SHARPISO_CUBE::SetVertexCoord(v0_coord.Ptr(), 1+2*offset);

  this->offset = offset;
}

// **************************************************
// GET_GRADIENTS_PARAM
// **************************************************

/// Initialize GET_GRADIENTS_PARAM
void SHARPISO::GET_GRADIENTS_PARAM::Init()
{
  use_only_cube_gradients = false;
  use_selected_gradients = true;
  use_intersected_edge_endpoint_gradients = false;
  use_gradients_determining_edge_intersections = false;
  grad_selection_cube_offset = 0;
  max_small_magnitude = 0.0;
  zero_tolerance = 0.0000001;
}

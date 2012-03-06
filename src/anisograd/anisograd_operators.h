/// \file anisograd_operators.h
/// compute gradients from scalar data
/// Version 0.0.1

/*
 IJK: Isosurface Jeneration Kode
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


#ifndef _ANISOGRAD_OPERATORS_
#define _ANISOGRAD_OPERATORS_

#include "sharpiso_grids.h"

// **************************************************
// TYPE DECLARATIONS
// **************************************************
using SHARPISO::SHARPISO_SCALAR_GRID_BASE;
using SHARPISO::GRADIENT_GRID_BASE;
using SHARPISO::GRADIENT_GRID;
using SHARPISO::VERTEX_INDEX;
using SHARPISO::GRADIENT_COORD_TYPE;
using SHARPISO::DIM3;

// **************************************************
// COMPUTE FORWARD DIFFERENCE
// **************************************************

/// Calculate the forward difference in the 'd' direction
/// Calculates for the scalar_grid
void compute_forward_difference_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE &fwd_diff_d);

/// Calculate the forward difference of the gradients in the 'd' direction
void compute_forward_difference_d_normals
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE fwd_diff_d_normals[DIM3]);

/// Calculate the forward difference of the gradients in the 'd' direction
/// for the Normals[index]
void compute_forward_difference_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int direction,
 const int index,
 GRADIENT_COORD_TYPE & fwd_diff_d_normals_index);

// **************************************************
// COMPUTE CENTRAL DIFFERENCE
// **************************************************

// Central difference operators
// compute central difference in the 'd' direction
void compute_central_difference_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, 
 const int d,
 GRADIENT_COORD_TYPE &cntrl_diff_d);

// Compute gradient of normals using central difference
// Used for anisogradinfo.
void compute_gradient_normals
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 GRADIENT_COORD_TYPE  gradient_normals[DIM3*DIM3]);

// Calculates the central  difference in th 'd' direction 
// for the Normals[index]
void compute_central_difference_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int direction,
 const int index,
 GRADIENT_COORD_TYPE & cntrl_diff_d_normals_index);

// **************************************************
// COMPUTE GRADIENTS
// **************************************************

void compute_gradH_d_scalar_grid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE gradientH_d_scalar_grid[DIM3]);

// Compute operator gradH for the direction 'd'
// for the Normal field
void compute_gradH_d_normals
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE  gradientH_d_Normals[DIM3*DIM3]);

// Compute operator gradH for the direction 'd'
// for Norm[index]
void compute_gradH_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
 const int direction,
 const int index,
 const VERTEX_INDEX iv1,
 GRADIENT_COORD_TYPE  gradientH_d_Normals_per_index[DIM3]);

// **************************************************
// COMPUTE VECTORS m, c, w
// **************************************************

// Compute M d for direction 'd' for  vertex iv1
void compute_m_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE m[DIM3]);

// Compute vector c_d
void compute_c_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE c[DIM3]);

// Compute vector w
void compute_w
(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const VERTEX_INDEX iv1, const int flag_aniso,
 const GRADIENT_GRID_BASE & gradient_grid,
 GRADIENT_COORD_TYPE w[DIM3]);

// **************************************************
// COMPUTE CURVATURE AND EXP FUNCTION
// **************************************************

// compute the curvature k for a vertex iv1
void compute_curvature_iv 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 GRADIENT_COORD_TYPE K[DIM3]);

// Compute gx as e^(-x/2*mu^2)
void compute_g_x(const float mu, const float param, 
                 const int flag_aniso, float & result);

// **************************************************
// VECTOR OPERATORS
// **************************************************

// Calculate the sum of squares of all elements in a vector 'vec'
// of size 'num_elements' and return the 'sum'
void vector_sum_of_squares 
(const float *vec, const int num_elements, float &sum);


// Vector dot product
void vector_dot_pdt 
(const float * A, const float *B, const int num_elements, float &res);

// Calculate vector magnitude.
void vector_magnitude (const float * vec, const int num_elements, float & mag);

/// Normalize the vectors.
/// Set small magnitude vectors to zero.
void normalize
(float *vec, float & magnitude, const int num_elements, 
 const float max_small_mag);

/// Normalize the vectors.
/// Set small magnitude vectors to zero.
void normalize
(float *vec, const int num_elements, const float max_small_mag);

#endif

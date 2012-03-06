/// \file anisograd.h
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

#ifndef _ANISOGRAD_
#define _ANISOGRAD_


#include <iostream>
#include "sharpiso_grids.h"
#include "anisograd_operators.h"

// global constants.
const int DIM9 = DIM3*DIM3;
const float EPSILON = 0.001;

using namespace SHARPISO;
using namespace std;
void compute_gradient_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const int icube, GRADIENT_GRID & gradient_grid);

// Calculate the anisotropic diff of the gradients.
void anisotropic_diff
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int num_iter,
 const int flag_aniso,
 const int icube,
 GRADIENT_GRID & gradient_grid);


// Calculate the anisotropic diff per iteration
void anisotropic_diff_iter_k
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int iter_k,
 const int flag_aniso,
 const int icube,
 const int dimension,
 GRADIENT_GRID & gradient_grid);


// Compute C d for direction 'd' for  vertex iv1
// called from the compute_m_d function 

void compute_c_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE c[DIM3]);

// Compute the gradient of the normal vector[index] in the 
// 'd' direction

void compute_gradH_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
 const int direction,
 const int index,
 const VERTEX_INDEX iv1,
 GRADIENT_COORD_TYPE  gradientH_d_Normals_per_index[DIM3]);

// Compute operator gradH for the direction 'd'
// for the Normal field
void compute_gradH_d_normals
(const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE  gradientH_d_Normals[DIM9]);


void compute_grad_H_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE * gradient);

void compute_m_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const int icube,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_COORD_TYPE m[DIM3]);

// Normalize the gradients, and store the magnitudes
// before calling anisotropic diff
// mag_list has the magnitudes of the original gradients.
void normalize_and_store_gradient_magnitudes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float max_small_mag,
 GRADIENT_GRID & gradient_grid,
 vector<GRADIENT_COORD_TYPE> &mag_list);

// Re Normalize the gradients and reset the magnitudes.
// after computing the anisotropic diffusion
void reset_gradient_magnitudes
(const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid,
 const float max_small_mag,
 GRADIENT_GRID & gradient_grid,
 const vector<GRADIENT_COORD_TYPE> mag_list);

#endif

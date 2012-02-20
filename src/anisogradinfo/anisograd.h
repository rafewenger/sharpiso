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

#include "isodual3D_datastruct.h"
#include "ijkscalar_grid.txx"


const float EPSILON = 0.001;
const int DIM3 = 3;
const int DIM9 = 9;


void compute_gradient_central_difference
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const int icube,
 ISODUAL3D::GRADIENT_GRID & gradient_grid);

// Calculate the anisotropic diff of the gradients.
void anisotropic_diff
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int num_iter,
 const int flag_aniso,
 const int icube,
  ISODUAL3D::GRADIENT_GRID & gradient_grid);


void compute_anisotropic_gradient_filtering
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid);

// Normalize the vectors.
void normalize (float *vec, const int num_elements);
// Calculate vector magnitude.
void vector_magnitude (const float * vec, const int num_elements, float & mag);
#endif

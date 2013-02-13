/// \file cgradient.h
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

#include "isodual3D_datastruct.h"
#include "ijkscalar_grid.txx"
#include "cgradient_inputinfo.h"

void compute_gradient_central_difference
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid, const INPUT_INFO & io_info);

void compute_reliable_gradients
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid, INPUT_INFO & io_info);

/*
 * Compute reliable gradients far, computes the central difference on all grid vertex.
 * if the gradient matches its grid neighbors then it is reliable
 * for far we look at 2 vertices one vertex over A-*-B
 */
void compute_reliable_gradients_far
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid, INPUT_INFO & io_info);

/*
 * Compute reliable gradients based on how good the gradients predict the scalars around it.
 * SBP stands for scalar based prediction.
 */
void compute_reliable_gradients_SBP
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid, INPUT_INFO & io_info);



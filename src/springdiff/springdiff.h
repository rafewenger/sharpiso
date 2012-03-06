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


#include "ijkscalar_grid.txx"
#include "sharpiso_grids.h"

void compute_gradient_central_difference
(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 SHARPISO::GRADIENT_GRID & gradient_grid);


void compute_spring_diffusion
(const SHARPISO::SHARPISO_SCALAR_GRID_BASE &scalar_grid,
 const int num_iter,
 const float lambda,
 const float mu,
 SHARPISO::GRADIENT_GRID & gradient_grid);



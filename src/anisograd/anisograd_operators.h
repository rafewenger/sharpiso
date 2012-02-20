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

#include "anisograd.h"
#include "ijkscalar_grid.txx"
#include "isodual3D_datastruct.h"

using namespace ISODUAL3D;

///
/// FORWARD DIFFERENCE CALCUALTORS 
///

// Calculate the FORWARD difference in the 'd' direction
// Calculates for the scalar_grid
void compute_forward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE &fwd_diff_d);

// Calculates the forward difference in the 'd' direction 
// Calculates for the Normals[index_coord]
void compute_forward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE &fwd_diff_d);



// Calculates the forward difference in the 'd' direction
// Calculates for the Normals
// [1X3]
void compute_forward_difference_d_normals
( const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE &fwd_diff_d_normals); 


///
/// BACKWARD DIFFERNCE CALCULATORS 
///

// Calculate the BACKWARD difference in the 'd' direction
// Calculates for the scalar_grid
void compute_backward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d, 
 GRADIENT_TYPE &bkwd_diff_d);

// Calculates the BACKWARD difference in the 'd' direction 
// Calculates for the Normals[index_coord]
void compute_backward_difference_d
( const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE &bkwd_diff_d);


// Compute operator gradH for the direction 'd'
// for the Normal field
void compute_gradH_d_normals
( const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE  gradientH_d_Normals[DIM9]
 );


// Central difference operators
// compute central difference in the 'd' direction
void compute_central_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, 
 const DIRECTION d,
 GRADIENT_TYPE &cntrl_diff_d);


// Compute M d for dierection 'd' for  vertex iv1
// 
void compute_m_d
(
 const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const int d,
 GRADIENT_TYPE m[DIM3]
 );

void compute_c_d
(
 const ISODUAL_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE c[DIM3]
 );

///// other functions


// Calculate the sum of squares of all elements in a vector 'vec'
// of size 'num_elements' and return the 'sum'
void vector_sum_of_squares 
(const float *vec, const int num_elements, float &sum);


// vector dot pdt
void vector_dot_pdt 
(const float * A, const float *B, const int num_elements, float &res);

#endif

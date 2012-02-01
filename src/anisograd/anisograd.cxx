/// \file anisograd.cxx
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

#include "anisograd.h"
#include "ijkscalar_grid.txx"
#include "isodual3D_datastruct.h"

using namespace ISODUAL3D;

// local type definition
namespace {
    
    typedef IJK::BOOL_GRID_BASE<ISODUAL_GRID> BOOL_GRID_BASE;
    typedef IJK::BOOL_GRID<ISODUAL_GRID> BOOL_GRID;
    
};

// local routines
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient);

void compute_central_differnce_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,  const DIRECTION &d, GRADIENT_TYPE cntrl_diff_d);

void compute_forward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, const DIRECTION &d,
 GRADIENT_TYPE fwd_diff_d);

void compute_backward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, const DIRECTION & d, 
 GRADIENT_TYPE bkwd_diff_d);


void compute_boundary_gradient
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient);


// Compute central difference main function 
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 GRADIENT_GRID & gradient_grid)
{
    const int dimension = scalar_grid.Dimension();
    
    gradient_grid.SetSize(scalar_grid, dimension);
    
    BOOL_GRID boundary_grid;
    boundary_grid.SetSize(scalar_grid);
    compute_boundary_grid(boundary_grid);
    
    for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
        if (boundary_grid.Scalar(iv)) {
            compute_boundary_gradient(scalar_grid, iv, gradient_grid.VectorPtr(iv));
        }
        else {
            compute_gradient_central_difference
            (scalar_grid, iv, gradient_grid.VectorPtr(iv));
        }
    }
}


// Compute gradH 
void compute_grad_H_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, const DIRECTION &d, GRADIENT_TYPE * gradient)
{
    GRADIENT_TYPE   temp0, temp1;
    for (int i=0; i<dimension; i++) {
        if (i==d) {
            compute_forward_difference_d(scalar_grid, iv1, i, gradient[i]);
        }
        else{
            compute_central_differnce_d(scalar_grid, iv1, i, temp0);
            VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
            compute_central_differnce_d(scalar_grid, iv2, i, temp1);
            gradient[i] = temp0 + temp1;
        }
    }
}

void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient)
{
    const int dimension = scalar_grid.Dimension();
    for (int d = 0; d < dimension; d++) {
        VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
        VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
        gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2;
    }
}


void compute_central_differnce_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,  const DIRECTION &d, GRADIENT_TYPE cntrl_diff_d)
{
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
    VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
    cntrl_diff_d = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv1))/2.0;
}


// Computes the forward diff in the d th dim
void compute_forward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, const DIRECTION &d, GRADIENT_TYPE fwd_diff_d)
{
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv0 = scalar_grid.NextVertex(iv1, d);
    fwd_diff_d = scalar_grid.Scalar(iv0) - scalar_grid.Scalar(iv1);
};

// Computes the backward diff in the d th dim
void compute_backward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, const DIRECTION & d, 
 GRADIENT_TYPE bkwd_diff_d)
{
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
    bkwd_diff_d = scalar_grid.Scalar(iv1) - scalar_grid.Scalar(iv0);
};



void compute_boundary_gradient
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient)
{
    const int dimension = scalar_grid.Dimension();
    GRID_COORD_TYPE coord[dimension];
    
    scalar_grid.ComputeCoord(iv1, coord);
    
    for (int d = 0; d < dimension; d++) {
        if (coord[d] > 0) {
            VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
            if (coord[d]+1 < scalar_grid.AxisSize(d)) {
                VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
                // use central difference
                gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2;
            }
            else {
                gradient[d] = scalar_grid.Scalar(iv1) - scalar_grid.Scalar(iv0);
            }
        }
        else if (coord[d]+1 < scalar_grid.AxisSize(d)) {
            VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
            gradient[d] = scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv1);
        }
        else {
            gradient[d] = 0;
        }
    }
}


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


using namespace ISODUAL3D;


const float EPSILON = 0.001;
const int DIM3 = 3;
const int DIM9 = 9;

class ANISO_INFO_TYPE{
    public :
    VERTEX_INDEX iv;
    int num_iter;
    GRADIENT_TYPE Normals[DIM3];
    GRADIENT_TYPE Normals2[DIM3];
    GRADIENT_TYPE mX[DIM3];
    GRADIENT_TYPE mX_prev_vert_X[DIM3];
    
    GRADIENT_TYPE mY[DIM3];
    GRADIENT_TYPE mY_prev_vert_Y[DIM3];
    
    GRADIENT_TYPE mZ[DIM3];
    GRADIENT_TYPE mZ_prev_vert_Z[DIM3];
    VERTEX_INDEX prev_vert[DIM3];
    
    SCALAR_TYPE K[DIM3];
    SCALAR_TYPE gK[DIM3];
    
    SCALAR_TYPE Kprev[DIM3];
    SCALAR_TYPE gKprev[DIM3];
    
    GRADIENT_TYPE   gradHNd[DIM9*DIM3];
    GRADIENT_TYPE   gradHNd_prev[DIM9*DIM3];
    GRADIENT_TYPE   c[DIM3*DIM3];
    GRADIENT_TYPE   c_prev[DIM3*DIM3];
    GRADIENT_TYPE   w[DIM3];
    GRADIENT_TYPE   wN ;
    GRADIENT_TYPE   w_dash[DIM3];
    void mycopy(const GRADIENT_TYPE temp[], GRADIENT_TYPE temp1[], const int  &n)
    {
        for(int i=0;i<n;i++)
        {
            temp1[i] = temp[i];
        }
    };

};

//print info
void debug_print(ANISO_INFO_TYPE & aniso_info);

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

// for debug purposes 
void anisotropic_diff_debug
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int num_iter,
 const int flag_aniso,
 const int icube,
 ISODUAL3D::GRADIENT_GRID & gradient_grid,
 ANISO_INFO_TYPE & aniso_infotype);


void compute_anisotropic_gradient_filtering
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 ISODUAL3D::GRADIENT_GRID & gradient_grid);

// Normalize the vectors.
void normalize (float *vec, const int num_elements);
// Calculate vector magnitude.
void vector_magnitude (const float * vec, const int num_elements, float & mag);
#endif

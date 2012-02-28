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
#include "anisograd_operators.h"

using namespace ISODUAL3D;

// local type definition
namespace {
    
    typedef IJK::BOOL_GRID_BASE<ISODUAL_GRID> BOOL_GRID_BASE;
    typedef IJK::BOOL_GRID<ISODUAL_GRID> BOOL_GRID;
    
};

// help routines 


// Normalize the vectors.
void normalize (float *vec, const int num_elements);



// local routines
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 GRADIENT_TYPE  gradient[DIM3]);

// Compute central difference in one direction and return a single value
void compute_central_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE &cntrl_diff_d);





// Computes the central difference in the d th dim 
// this is for the normals
void compute_central_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE &cntrl_diff_d);





void compute_boundary_gradient
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient);

void compute_grad_H_d
( const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1, const DIRECTION d, const int index_coord, GRADIENT_TYPE *gradient);



// Compute central difference main function 
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX icube,
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



// Compute anisotropic filtering of gradients main function 
using namespace std;

// Compute Gd
void compute_grad_H_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE * gradient)
{
    const int dimension = scalar_grid.Dimension();
    GRADIENT_TYPE   temp0, temp1;
    for (int i=0; i<dimension; i++) {
        if (i==d) {
            compute_forward_difference_d(scalar_grid, iv1, i, gradient[i]);
        }
        else{
            compute_central_difference_d(scalar_grid, iv1, i, temp0);
            VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
            compute_central_difference_d(scalar_grid, iv2, i, temp1);
            gradient[i] = (temp0 + temp1)/2.0;
        }
    }
}


// this is for the normals [help function for Gnd]
void compute_grad_H_d
( const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE *gradient)

{
    const int dimension  = gradient_grid.Dimension();
    GRADIENT_TYPE   temp0, temp1;
    for (int i=0; i<dimension; i++) {
        if (i==d) {
            compute_forward_difference_d(scalar_grid, gradient_grid, iv1, i, index_coord, gradient[i]);
        }
        else{
            compute_central_difference_d(scalar_grid, gradient_grid, iv1, i, index_coord, temp0);
            
            VERTEX_INDEX iv2 = gradient_grid.NextVertex(iv1, d);
            compute_central_difference_d(scalar_grid, gradient_grid, iv2, i, index_coord, temp1);
            
            gradient[i] = (temp0 + temp1)/2.0;
        }
        
    }
}


// compute GNd
void compute_grad_H_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 GRADIENT_TYPE * gradient)
{
    // This generates a  3X3 matrix
    const int dimension = scalar_grid.Dimension();
    for (int index_coord=0; index_coord<dimension; index_coord++) {
        
        compute_grad_H_d (scalar_grid, gradient_grid, iv1, d, index_coord, gradient + 3*index_coord);
    }
}

void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 GRADIENT_TYPE gradient[DIM3])
{
    const int dimension = scalar_grid.Dimension();
    for (int d = 0; d < dimension; d++) {
        VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
        VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
        gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2.0;
    }
}







// Computes the backward diff in the d th dim
void compute_backward_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d, 
 GRADIENT_TYPE &bkwd_diff_d)
{
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
    bkwd_diff_d = scalar_grid.Scalar(iv1) - scalar_grid.Scalar(iv0);
};


void compute_central_difference_d
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1, 
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE &cntrl_diff_d)
{
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
    VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
    
    const GRADIENT_TYPE * vertex_gradient_coord0 = gradient_grid.VectorPtrConst(iv0);
    
    const GRADIENT_TYPE * vertex_gradient_coord2 = gradient_grid.VectorPtrConst(iv2);
    
    cntrl_diff_d = (vertex_gradient_coord2[index_coord] - vertex_gradient_coord0[index_coord])/2.0;
}




// Computes the backward diff in the d th dim
void compute_backward_difference_d
( const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1,
 const DIRECTION d,
 const int index_coord,
 GRADIENT_TYPE &bkwd_diff_d)
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


// Compute C for eacn dimension

void compute_C_d 
(GRADIENT_TYPE * GNd,
 GRADIENT_TYPE * Gd,
 const int dimension,
 GRADIENT_TYPE * Cd)
{
    //normalize the Gd
    normalize(Gd, dimension);
    // dot pdt
    for (int i=0; i<dimension; i++)
    {
        GRADIENT_TYPE  temp[3]={0.0};
        
        for (int j=0; j<dimension; j++) 
        {
            temp[i] = GNd[3*j+i];
        }
        vector_dot_pdt (temp, Gd, dimension, Cd[i]);
    }
};

// Compute gx as e^(-x/2*mu^2)
void compute_g_x(const float mu, const float param, const int flag_aniso, float & result )
{
    // the flag_aniso when set to zero calculates isotropic diffusion 
    // when set to 1 calculates the anisotropic diffusion 
    result = exp ((param*param*float(flag_aniso))/(-2.0*mu*mu));
    
};



// new version of anisotropic gradient filtering per vertex
void anisotropic_diff_per_vert
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const VERTEX_INDEX iv1,
 const int flag_aniso,
 const int icube,
 const GRADIENT_GRID & gradient_grid,
 GRADIENT_GRID & temp_gradient_grid)
{
    GRADIENT_TYPE mX[DIM3]={0.0};
    GRADIENT_TYPE mX_prev_vert_X[DIM3]={0.0};
    
    GRADIENT_TYPE mY[DIM3]={0.0};
    GRADIENT_TYPE mY_prev_vert_Y[DIM3]={0.0};
    
    GRADIENT_TYPE mZ[DIM3]={0.0};
    GRADIENT_TYPE mZ_prev_vert_Z[DIM3]={0.0};
    
    // Compute M d for direction 'd' for  vertex iv1
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 0, mX);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 1, mY);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 2, mZ);
    
    // compute prev vertex in 0,1,2 direction 
    VERTEX_INDEX prev_vert[DIM3];
    prev_vert[0] = scalar_grid.PrevVertex(iv1, 0);
    prev_vert[1] = scalar_grid.PrevVertex(iv1, 1);
    prev_vert[2] = scalar_grid.PrevVertex(iv1, 2);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, prev_vert[0], 0, mX_prev_vert_X);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, prev_vert[1], 1, mY_prev_vert_Y);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, prev_vert[2], 2, mZ_prev_vert_Z);
    
    // Compute 'k' for each dimensions , 
    // used for anisotropic gradients
    
     SCALAR_TYPE K[DIM3]={0.0};
     SCALAR_TYPE gK[DIM3]={0.0};
     
     SCALAR_TYPE Kprev[DIM3]={0.0};
     SCALAR_TYPE gKprev[DIM3]={0.0};
     
     for (int d=0; d<DIM3; d++) {
     
     // compute gradHN_d
     GRADIENT_TYPE   gradHN_d[DIM9]={0.0};
     compute_gradH_d_normals (gradient_grid, iv1, d,gradHN_d);
     // compute C_d
     GRADIENT_TYPE c[DIM3]={0.0};
     compute_c_d(scalar_grid, gradient_grid, iv1, d, c);
     
     SCALAR_TYPE sum_gradHNd = 0.0, c_square = 0.0;
     
     vector_sum_of_squares(gradHN_d, DIM9, sum_gradHNd);
     
     
     vector_dot_pdt(c, c, DIM3, c_square);
     
     K[d] = sum_gradHNd  - c_square;

     // *** DEBUG ***
     /*
     if (iv1 == icube){
       cout << "c: ";
       IJK::ijkgrid_output_coord(cout, DIM3, c);
       cout <<" gradHNd ";
       IJK::ijkgrid_output_coord(cout, 9, gradHN_d);
       cout << endl;
       cout <<"\nd "<<d<<" sum_gradHNd " << sum_gradHNd 
            <<" csq "<<c_square<<endl;	
     }
     */

     compute_g_x(mu, K[d], flag_aniso, gK[d]);
     }

    // *** DEBUG ***
    /*
    if (iv1 == icube) {
        cout <<"gk ";
        for (int q=0; q<3; q++) {
            cout <<" "<<gK[q];
        }
        cout <<endl;
        cout <<"K ";
        for (int q=0; q<3; q++) {
            cout <<" "<<K[q];
        }
        cout <<endl;
    }
    */

    for (int d=0; d<DIM3; d++) 
    {
        
        // compute gradHN_d
        GRADIENT_TYPE   gradHN_d[DIM9]={0.0};
        compute_gradH_d_normals (gradient_grid, prev_vert[d], d,gradHN_d);
        // compute C_d
        GRADIENT_TYPE c[DIM3]={0.0};
        compute_c_d(scalar_grid, gradient_grid, prev_vert[d], d, c);
        
        SCALAR_TYPE sum_gradHNd = 0.0, c_square = 0.0;
        
        vector_sum_of_squares(gradHN_d, DIM9, sum_gradHNd);
        
        vector_dot_pdt(c, c, DIM3, c_square); 
        
        K[d] = sum_gradHNd  - c_square;
        
        compute_g_x(mu, K[d],flag_aniso, gKprev[d]);
    }
    //debug
    //
    // *** DEBUG ***
    if (iv1 == icube) {
        cout <<"gkprev ";
        for (int q=0; q<3; q++) {
            cout <<" "<<gKprev[q];
        }
        cout <<endl;
        cout <<"K ";
        for (int q=0; q<3; q++) {
            cout <<" "<<K[q];
        }
        cout <<endl;
    }

    GRADIENT_TYPE    w[DIM3]={0.0};
    /*
    for (int i =0; i<DIM3; i++) {
        w[i] =  mX[i] - mX_prev_vert_X[i] +
        mY[i] - mY_prev_vert_Y[i] +
        mZ[i] - mZ_prev_vert_Z[i] ;
    }
     */
    
     for (int i =0; i<DIM3; i++) {
     w[i] =  gK[i]*mX[i] - gKprev[i]*mX_prev_vert_X[i] +
     gK[i]*mY[i] - gKprev[i]*mY_prev_vert_Y[i] +
     gK[i]*mZ[i] - gKprev[i]*mZ_prev_vert_Z[i] ;
     }
    
    GRADIENT_TYPE  wN = 0.0;
    
    const GRADIENT_TYPE * Normals = gradient_grid.VectorPtrConst(iv1);
    for (int i=0; i<DIM3; i++) {
        wN = wN + Normals[i]*w[i];
    }
    
    for (int i=0; i<DIM3; i++) {
        w[i] = w[i] - wN*Normals[i];
    }
    
    // *** DEBUG ***
    /*
    if (iv1 == icube) {
        cout <<" wN "<< wN<<" No "<<Normals[0]<<" "<<Normals[1]<<" "<<Normals[2]<<endl;;
        cout <<" w: "<< w[0]<<" "<< w[1]<<" "<< w[2] <<endl;
        for(int h=0;h<3;h++){
        cout <<"gK[i]*mX[i]                 "<<h<<" "<< gK[h]*mX[h]<<endl;
        cout <<"gKprev[h]*mX_prev_vert_X[h] "<<h<<" "<<gKprev[h]*mX_prev_vert_X[h]<<endl;
        cout <<"gK[i]*mY[i]                 "<<h<<" "<< gK[h]*mY[h]<<endl ;  
        cout <<"gKprev[i]*mY_prev_vert_Y[i] "<<h<<" "<< gKprev[h]*mY_prev_vert_Y[h]<<endl ; 
        cout <<"gK[i]*mZ[i]                 "<<h<<" "<< gK[h]*mZ[h] <<endl;
        cout <<"gKprev[i]*mZ_prev_vert_Z[i] "<<h<<" "<< gKprev[h]*mZ_prev_vert_Z[h]<<endl ;  
           double  w_test =  gK[h]*mX[h] - gKprev[h]*mX_prev_vert_X[h] +
            gK[h]*mY[h] - gKprev[h]*mY_prev_vert_Y[h] +
            gK[h]*mZ[h] - gKprev[h]*mZ_prev_vert_Z[h] ;  
            cout<<"wh                           "<<h<<" "<<w_test<<endl;

        }
        cout <<"gk ";
        for (int q=0; q<3; q++) {
            cout <<" "<<gK[q];
        }
        cout <<endl;
        cout <<"K ";
        for (int q=0; q<3; q++) {
            cout <<" "<<K[q];
        }
        cout <<endl;
    }
    */

    // *** DEBUG ***
    /* 
    if (iv1 == icube)
    cout <<" w: "<< w[0]<<" "<< w[1]<<" "<< w[2] <<endl;
    */

    GRADIENT_TYPE Normals2[DIM3];
    for (int i=0; i<DIM3; i++) {
        Normals2[i] = Normals[i]  + lambda *w[i];
    }

    // update the temporary gradient grid
    temp_gradient_grid.Set(iv1, Normals2);
}


// Calculate the anisotropic diff of the gradients.
void anisotropic_diff
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int num_iter,
 const int flag_aniso,
 const int icube,
 GRADIENT_GRID & gradient_grid)
{
    
    const int dimension = scalar_grid.Dimension();
    
    gradient_grid.SetSize(scalar_grid, dimension);
    
    BOOL_GRID boundary_grid;
    boundary_grid.SetSize(scalar_grid);
    compute_boundary_grid(boundary_grid);
    
    for (int k=0; k<num_iter; k++) 
    {
        // create a temporary grid for each iteration 
        GRADIENT_GRID temp_gradient_grid;
        temp_gradient_grid.SetSize(scalar_grid, dimension);

        for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++)
        {

            if(!boundary_grid.Scalar(iv))
            {
               // cout <<"v "<<iv<<" k "<<k<<endl;
                // send in 2 seperate grids
                anisotropic_diff_per_vert
                (scalar_grid, mu, lambda, iv, flag_aniso, icube, gradient_grid, temp_gradient_grid );
            }
            else
            {
                 compute_boundary_gradient(scalar_grid, iv, temp_gradient_grid.VectorPtr(iv));
            }
            using namespace std;
            
            //debug
            COORD_TYPE coord[DIM3];
            scalar_grid.ComputeCoord(icube, coord);
            if(iv == icube && flag_aniso == 0 )
            {
                GRADIENT_TYPE * grad = gradient_grid.VectorPtr(icube);
            }
            
            if(iv == icube && flag_aniso != 0 )
            {
                GRADIENT_TYPE * grad = gradient_grid.VectorPtr(icube);
            }
        }

        // Copy temp_gradient_grid to gradient_grid.
        for (int l=0; l<scalar_grid.NumVertices(); l++) {
            GRADIENT_TYPE * grad = temp_gradient_grid.VectorPtr(l);
            gradient_grid.Set(l, grad);
        }

        // Normalize gradient grid vectors.
        for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++)
        {
            GRADIENT_TYPE  * N = gradient_grid.VectorPtr(iv);
            GRADIENT_TYPE   mag = 0.0;
            vector_magnitude (N, DIM3, mag);
            
            if (mag > 0.0001)
            {
               
                normalize (N, DIM3);
                gradient_grid.Set(iv, N);
            }
        }

    }
};


// helper routines 

// Calculate vector magnitude.
void vector_magnitude (const float * vec, const int num_elements, float & mag)
{
    float sum = 0.0;
    mag = 0.0;
    for (int i=0; i<num_elements; i++) {
        sum = sum + vec[i]*vec[i];
    }
    mag = sqrt(sum);
}





// Normalize the vectors.
void normalize (float *vec, const int num_elements)
{
    float mag = 0.0;
    vector_magnitude (vec, num_elements, mag);
    
    for (int i=0; i<num_elements; i++) {
        if (abs(mag-0.0) < EPSILON) {
            vec[i] = 0.0;
        }
        else
            vec[i] = vec[i] / mag;
    }
}


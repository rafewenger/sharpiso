/// \file springdiff.cxx
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
 Foundation, Inc., 89 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <iomanip>
#include "springdiff.h"
#include "ijkscalar_grid.txx"
#include "sharpiso_scalar.txx"

#include "sharpiso_grids.h"

using namespace SHARPISO;
using namespace std;

const float EPSILON = 0.0001;

// local type definition
namespace {
  
  typedef IJK::BOOL_GRID_BASE<SHARPISO_GRID> BOOL_GRID_BASE;
  typedef IJK::BOOL_GRID<SHARPISO_GRID> BOOL_GRID;
  
};

// local routines
void compute_gradient_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient);
void compute_boundary_gradient
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient);
// Calculate vector magnitude.
void vector_magnitude 
(const float * vec, const int num_elements, float & mag);


void compute_gradient_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
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

void compute_gradient_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient)
{
  const int dimension = scalar_grid.Dimension();
  for (int d = 0; d < dimension; d++) {
    VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
    VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
    gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2;
  }
}

void compute_boundary_gradient
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient)
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

//local function 
// computes the distance between the the vertex iv and v;
void compute_dist_for_v
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const VERTEX_INDEX iv, 
 const VERTEX_INDEX u,
 float & distance
 ){
  const GRADIENT_COORD_TYPE * grad_u = gradient_grid.VectorPtrConst(u);
  const GRADIENT_COORD_TYPE * grad_iv = gradient_grid.VectorPtrConst(iv);
  GRID_COORD_TYPE coord_u[DIM3], coord_iv[DIM3];
  scalar_grid.ComputeCoord(u, coord_u);
  scalar_grid.ComputeCoord(iv, coord_iv);
  
  float mag_u={0.0};
  float mag_v={0.0};
  vector_magnitude (grad_u, DIM3,  mag_u);
  vector_magnitude (grad_u, DIM3,  mag_v);
  compute_distance_to_gfield_plane
  (grad_u, coord_u, scalar_grid.Scalar(u), coord_iv, scalar_grid.Scalar(iv), distance);    
  /*
   //debug 
   cout <<" dist     "<<distance<<endl;
   cout <<" grad_u  ("<<grad_u[0]<<","<<grad_u[1]<<","<<grad_u[2]<<")"<<endl;
   cout <<" grad_iv  ("<<grad_iv[0]<<","<<grad_iv[1]<<","<<grad_iv[2]<<")"<<endl;
   cout <<" coord_u ("<<coord_u[0]<<","<<coord_u[1]<<","<<coord_u[2]<<")"<<endl;
   cout <<" coord_iv("<<coord_iv[0]<<","<<coord_iv[1]<<","<<coord_iv[2]<<")"<<endl;
   cout <<" scalar_u("<<scalar_grid.Scalar(u)<<") scalar v("<<scalar_grid.Scalar(iv)<<")"<<endl;
   */
   };


//local function to calculate the weighted distance
void f
(const float lambda,
 const float mu, 
 const float distance,
 float &weighted_dist)
{
  float dist_sq = distance*distance;
  weighted_dist = lambda*exp(-dist_sq/(mu*mu));
  
}
// Local function 
// Upgrade the gdiff based on the calculated distance
void update_gdiff 
(
 const GRADIENT_GRID & gradient_grid,
 const float dist,
 const float lambda,
 const float mu,
 const VERTEX_INDEX u,
 const VERTEX_INDEX iv,
 GRADIENT_COORD_TYPE gdiff[DIM3])
{
  //calculate the weight based on distance 
  float weighted_dist=0.0;
  f(lambda, mu, dist, weighted_dist);
  const GRADIENT_COORD_TYPE * grad_u = gradient_grid.VectorPtrConst(u);
  const GRADIENT_COORD_TYPE * grad_iv = gradient_grid.VectorPtrConst(iv);
  for (int d=0; d<DIM3; d++) {
    gdiff[d] = gdiff[d] + weighted_dist*(grad_u[d] - grad_iv[d]);
  }
}
// Function to update all the gradients from the temporary gradient_grid
// to the original gradient grid.
void update_all_gradients
( GRADIENT_GRID & temp_gradient_grid,
 GRADIENT_GRID & gradient_grid)
{
  // Copy temp_gradient_grid to gradient_grid.
  for (int l=0; l<gradient_grid.NumVertices(); l++) {
    GRADIENT_COORD_TYPE * grad = temp_gradient_grid.VectorPtr(l);
    gradient_grid.Set(l, grad);
  }
}
// main function for computing spring based diffusion
void compute_spring_diffusion
(const SHARPISO_SCALAR_GRID_BASE &scalar_grid,
 const int num_iter,
 const float lambda,
 const float mu,
 GRADIENT_GRID & gradient_grid){
  
  GRADIENT_GRID temp_gradient_grid;
  const int dimension = scalar_grid.Dimension();
  temp_gradient_grid.SetSize(scalar_grid, dimension);
  
  BOOL_GRID boundary_grid;
  boundary_grid.SetSize(scalar_grid);
  compute_boundary_grid(boundary_grid);
  
  float distance=0.0;
  for (int i=0; i<num_iter; i++) {
    for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++){
      
      if (boundary_grid.Scalar(iv)) {
        compute_boundary_gradient(scalar_grid, iv, temp_gradient_grid.VectorPtr(iv));
      }
      else {
        
        COORD_TYPE coord[DIM3]={0.0};
        scalar_grid.ComputeCoord(iv, coord);
        GRADIENT_COORD_TYPE gdiff[DIM3]={0.0};
        GRADIENT_COORD_TYPE * grad_iv = gradient_grid.VectorPtr(iv);
        
          for (int d=0; d<DIM3; d++) {
            // prev vertices to v
            VERTEX_INDEX prev = scalar_grid.PrevVertex(iv, d);
            if(gradient_grid.ComputeMagnitudeSquared(prev) > EPSILON){
            compute_dist_for_v(scalar_grid, gradient_grid, iv, prev, distance);
            update_gdiff(gradient_grid, distance, lambda, mu, prev, iv, gdiff);
            }
            // next vertices to iv

            VERTEX_INDEX next = scalar_grid.NextVertex(iv, d);
            if(gradient_grid.ComputeMagnitudeSquared(next) > EPSILON){
            compute_dist_for_v(scalar_grid, gradient_grid, iv, next, distance);
            update_gdiff(gradient_grid, distance, lambda, mu, next, iv, gdiff);
            }
          }
        
        GRADIENT_COORD_TYPE temp_grad[DIM3]={0.0};
        // update the gradients of iv based on the weights
        for (int d=0; d<DIM3; d++) {
          temp_grad[d] = grad_iv[d] + gdiff[d];
        }
        temp_gradient_grid.Set(iv, temp_grad);
      }
    }
    //update all gradients
    update_all_gradients(temp_gradient_grid, gradient_grid);
  }
};


//local routines 
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



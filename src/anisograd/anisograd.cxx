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

#include "ijkcoord.txx"

// local type definition
namespace {
  typedef IJK::BOOL_GRID_BASE<SHARPISO_GRID> BOOL_GRID_BASE;
  typedef IJK::BOOL_GRID<SHARPISO_GRID> BOOL_GRID;

};

// **************************************************
// COMPUTE GRADIENT USING CENTRAL DIFFERENCE
// **************************************************

// Compute central difference main function
void compute_gradient_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const int icube,
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
      compute_central_difference
        (scalar_grid, iv, gradient_grid.VectorPtr(iv));
    }
  }
}

// **************************************************
// ANISOTROPIC GRADIENT FILTERING
// **************************************************

//////
// anisotropic gradient filtering per vertex
void anisotropic_diff_per_vert
(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const VERTEX_INDEX iv1,
 const int flag_aniso,
 const int icube,
 const GRADIENT_GRID & gradient_grid,
 GRADIENT_GRID & temp_gradient_grid)
{
  const GRADIENT_COORD_TYPE * N = gradient_grid.VectorPtrConst(iv1);
  GRADIENT_COORD_TYPE w[DIM3], w_orth[DIM3], newN[DIM3];

  compute_w(scalar_grid, mu, iv1, flag_aniso, gradient_grid, w);

  IJK::compute_orthogonal_vector(DIM3, w, N, w_orth);

  for (int d=0; d<DIM3; d++) 
    { newN[d] = N[d]  + lambda * w_orth[d]; }

  // update the temporary gradient grid
  temp_gradient_grid.Set(iv1, newN);
}


/// calculate the anisotropic diff of the gradients per iteration
void anisotropic_diff_iter_k
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int iter_k,
 const int flag_aniso,
 const int icube,
 const int dimension,
 GRADIENT_GRID & gradient_grid)
{
  // create a temporary grid for each iteration
  GRADIENT_GRID temp_gradient_grid;
  temp_gradient_grid.SetSize(scalar_grid, dimension);

  for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++)
    {
      GRID_COORD_TYPE coord[DIM3];
      scalar_grid.ComputeCoord(iv, coord);

      if (2 <= coord[0] && coord[0]+2 < scalar_grid.AxisSize(0) &&
          2 <= coord[1] && coord[1]+2 < scalar_grid.AxisSize(1) &&
          2 <= coord[2] && coord[2]+2 < scalar_grid.AxisSize(2)) {

        // Store new gradient values in temp_gradient_grid
        anisotropic_diff_per_vert
          (scalar_grid, mu, lambda, iv, flag_aniso, icube, 
           gradient_grid, temp_gradient_grid );
      }
      else {
        compute_boundary_gradient
          (scalar_grid, iv, temp_gradient_grid.VectorPtr(iv));
      }
    }

  // Copy temp_gradient_grid to gradient_grid.
  for (int l=0; l<scalar_grid.NumVertices(); l++) {
    GRADIENT_COORD_TYPE * grad = temp_gradient_grid.VectorPtr(l);
    // debug normalize
    normalize(grad, DIM3, EPSILON);
    gradient_grid.Set(l, grad);
  }
}

// Calculate the anisotropic diff of the gradients.
void anisotropic_diff
//(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float mu,
 const float lambda,
 const int num_iter,
 const int flag_aniso,
 const int icube,
 GRADIENT_GRID & gradient_grid)
{

  const int dimension = scalar_grid.Dimension();
  gradient_grid.SetSize(scalar_grid, dimension);
  for (int k=0; k<num_iter; k++) {
    anisotropic_diff_iter_k
      (scalar_grid, mu, lambda, k, flag_aniso, icube, dimension, gradient_grid);
  }
};


// Normalize the gradients, and store the magnitudes
// before calling anisotropic diff
void normalize_and_store_gradient_magnitudes
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const float max_small_mag,
 GRADIENT_GRID & gradient_grid,
 vector<GRADIENT_COORD_TYPE> &mag_list)
// mag_list has the magnitudes of the original gradients.
{
  // Normalize the gradients and
  // store the magnitudes so that they can be later added back
  for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++)
    {
      GRADIENT_COORD_TYPE  * N = gradient_grid.VectorPtr(iv);
      GRADIENT_COORD_TYPE   magnitude = 0.0;

      normalize(N, magnitude, DIM3, max_small_mag);
      mag_list.push_back(magnitude);
    }
};

// Re Normalize the gradients and reset the magnitudes.
// after computing the anisotropic diffusion
void reset_gradient_magnitudes
(const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid,
 const float max_small_mag,
 GRADIENT_GRID & gradient_grid,
 const vector<GRADIENT_COORD_TYPE> mag_list)
{
  // reset the gradients to be normalized
  for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
    {
      GRID_COORD_TYPE coord[DIM3];

      GRADIENT_COORD_TYPE  * N = gradient_grid.VectorPtr(iv);
      GRADIENT_COORD_TYPE   mag;

      normalize(N, mag, DIM3, max_small_mag);

      for (int i=0; i<DIM3; i++) 
        { N[i] = N[i]*mag_list[iv]; }

      gradient_grid.Set(iv, N);
    }
}


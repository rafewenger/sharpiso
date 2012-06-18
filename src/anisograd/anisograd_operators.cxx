/// \file anisograd_operators.cxx
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

#include "anisograd_operators.h"

using namespace std;

using namespace SHARPISO;

// **************************************************
// COMPUTE FORWARD DIFFERENCE
// **************************************************


// Calculate the FORWARD difference in the 'd' direction
// Calculates for the scalar_grid
void compute_forward_difference_d(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE &fwd_diff_d
)
{
	VERTEX_INDEX iv0 = scalar_grid.NextVertex(iv1, d);
	fwd_diff_d = scalar_grid.Scalar(iv0) - scalar_grid.Scalar(iv1);
}


///////
// Calculates the forward difference in the 'd' direction
// Calculates for the Normals
void compute_forward_difference_d_normals
(const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE fwd_diff_d_normals[DIM3]) //[3 terms in it]
{
	const int dimension = gradient_grid.Dimension();

	for (int index=0; index<dimension; index++)
		compute_forward_difference_d_normals_per_index
		(gradient_grid, iv1, d, index, fwd_diff_d_normals[index]);
}

// Calculate the forward difference of the gradients in the 'd' direction
//   for i'th component of the gradient.
void compute_forward_difference_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int direction,
		const int index,
		GRADIENT_COORD_TYPE & fwd_diff_d_normals_index)
{
	VERTEX_INDEX next_vert = gradient_grid.NextVertex(iv1, direction);

	const GRADIENT_COORD_TYPE * vertex_grad_iv1 =
			gradient_grid.VectorPtrConst (iv1);
	const GRADIENT_COORD_TYPE * vertex_grad_next_vert =
			gradient_grid.VectorPtrConst (next_vert);

	fwd_diff_d_normals_index =
			vertex_grad_next_vert[index] - vertex_grad_iv1[index];
}


// **************************************************
// COMPUTE CENTRAL DIFFERENCE
// **************************************************

/// Compute central difference on the scalar grid in the 'd' direction
void compute_central_difference_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE &cntrl_diff_d)
{
	VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
	VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
	cntrl_diff_d = 0.5*(scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0));
}

/// Compute central difference on a vertex in the scalar grid
void compute_central_difference
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1,
		GRADIENT_COORD_TYPE gradient[DIM3])
{
	const int dimension = scalar_grid.Dimension();
	for (int d = 0; d < dimension; d++) {
		VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
		VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
		gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2.0;
	}
}

// Calculates the central  difference in th 'd' direction 
// for the Normals[index]
void compute_central_difference_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int direction,
		const int index,
		GRADIENT_COORD_TYPE & cntrl_diff_d_normals_index)
{
	GRID_COORD_TYPE coord[DIM3];
	gradient_grid.ComputeCoord(iv1, coord);
	if (0 < coord[direction] &&
			coord[direction] + 1 < gradient_grid.AxisSize(direction)) {
		VERTEX_INDEX next_vert = gradient_grid.NextVertex(iv1, direction);
		VERTEX_INDEX prev_vert = gradient_grid.PrevVertex(iv1, direction);

		cntrl_diff_d_normals_index =
				(gradient_grid.Vector(next_vert, index) -
						gradient_grid.Vector(prev_vert, index))*0.5;

	}
	else {
		cout <<"error"<<endl;
		cntrl_diff_d_normals_index=0.0;
	}

}

/////////////
// Compute gradient of normals using central difference
// Used for anisogradinfo.
void compute_gradient_normals
(const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		GRADIENT_COORD_TYPE  gradient_normals[DIM3*DIM3])
{
	GRADIENT_COORD_TYPE cdiff;

	// for each index to the normal vector N[0]. N[1]. N[2].
	for (int i=0; i<DIM3; i++) {

		for (int dir=0; dir<DIM3; dir++) {
			compute_central_difference_d_normals_per_index
			(gradient_grid, iv1, dir, i, cdiff);
			gradient_normals[i*DIM3+dir] = cdiff;
		}
	}
}


// **************************************************
// COMPUTE BOUNDARY GRADIENT
// **************************************************

void compute_boundary_gradient
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE gradient[DIM3])
{
	GRID_COORD_TYPE coord[DIM3];

	scalar_grid.ComputeCoord(iv1, coord);

	for (int d = 0; d < DIM3; d++) {
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


// **************************************************
// COMPUTE GRADIENTS
// **************************************************

// Compute Gd
void compute_grad_H_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE * gradient)
{
	const int dimension = scalar_grid.Dimension();
	GRADIENT_COORD_TYPE   temp0, temp1;
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

/////////////
// Compute operator gradH for the direction 'd'
// for the Normal field
void compute_gradH_d_normals
(const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE  gradientH_d_Normals[DIM3*DIM3])
{
	const int dimension = gradient_grid.Dimension();
	for (int index=0; index<dimension; index++) {
		// for each index to the normal vector N[0]. N[1]. N[2].
		compute_gradH_d_normals_per_index
		(gradient_grid, d, index, iv1, &gradientH_d_Normals[DIM3*index]);
	}
}

////////////
// Compute the gradient of the normal vector[index] in the 
// 'd' direction
void compute_gradH_d_normals_per_index
(const GRADIENT_GRID_BASE & gradient_grid,
		const int direction,
		const int index,
		const VERTEX_INDEX iv1,
		GRADIENT_COORD_TYPE  gradientH_d_Normals_per_index[DIM3])
{
	for (int i=0; i<DIM3; i++) {
		if(i==direction) {
			// compute_forward_difference_d_normals_per_index
			compute_forward_difference_d_normals_per_index
			(gradient_grid, iv1, direction, index,
					gradientH_d_Normals_per_index[i]);
		}
		else {
			GRADIENT_COORD_TYPE    temp1=0.0, temp2=0.0;
			compute_central_difference_d_normals_per_index
			(gradient_grid, iv1, i, index, temp1);

			VERTEX_INDEX nextvert = gradient_grid.NextVertex(iv1, direction);

			compute_central_difference_d_normals_per_index
			(gradient_grid, nextvert, i, index, temp2);

			gradientH_d_Normals_per_index[i] = 0.5*(temp1 + temp2);
		}
	}

}


// Compute operator gradH for the direction 'd' for the scalar grid
void compute_gradH_d_scalar_grid(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE gradientH_d_scalar_grid[DIM3]
)
{
	for (int i=0; i<DIM3; i++) {
		if (i==d) {
			compute_forward_difference_d (scalar_grid, iv1, d, gradientH_d_scalar_grid[i]);
		}
		else {
			GRADIENT_COORD_TYPE temp1, temp2;
			compute_central_difference_d ( scalar_grid, iv1, i, temp1);

			VERTEX_INDEX nextvert  = scalar_grid.NextVertex(iv1, d);
			compute_central_difference_d (scalar_grid, nextvert, i, temp2);

			gradientH_d_scalar_grid[i] = 0.5 * (temp1 + temp2);
		}
	}
}

// **************************************************
// COMPUTE CURVATURE AND EXP FUNCTION
// **************************************************

// compute the curvature k for a vertex iv1, direction d
void compute_curvature_iv_d
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE & K)
{
	// compute gradHN_d
	GRADIENT_COORD_TYPE   gradHN_d[DIM3*DIM3]={0.0};
	compute_gradH_d_normals(gradient_grid, iv1, d,gradHN_d);

	// compute C_d
	GRADIENT_COORD_TYPE c[DIM3]={0.0};
	compute_c_d(scalar_grid, gradient_grid, iv1, d, c);

	SCALAR_TYPE sum_gradHNd = 0.0, c_square = 0.0;

	vector_sum_of_squares(gradHN_d, DIM3*DIM3, sum_gradHNd);

	vector_dot_pdt(c, c, DIM3, c_square);

	GRADIENT_COORD_TYPE gr[DIM3];
	compute_grad_H_d(scalar_grid, iv1, d, gr);

	GRADIENT_COORD_TYPE mag=0.0;
	vector_magnitude (gr, DIM3, mag);

	K = sum_gradHNd  - c_square*mag*mag;
}

// compute the curvature k for a vertex iv1
void compute_curvature_iv
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		GRADIENT_COORD_TYPE K[DIM3])
{

	for (int d = 0; d<DIM3; d++)
	{ compute_curvature_iv_d(scalar_grid, gradient_grid, iv1, d, K[d]); }
}

// Compute gx as e^(-x^2/2*mu^2)
void compute_g_x
(const float mu, const float param, const int flag_aniso, float & result )
{
	// the flag_aniso when set to zero calculates isotropic diffusion
	// when set to 1 calculates the anisotropic diffusion
	result = exp ((param*param*float(flag_aniso))/(-2.0*mu*mu));

}

// **************************************************
// COMPUTE VECTORS m, c, w
// **************************************************

// Compute C d for direction 'd' for  vertex iv1
// called from the compute_m_d function 
void compute_c_d
(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE c[DIM3]
)
{
	// compute the gradient of the Normal vector
	GRADIENT_COORD_TYPE gradientN_d[DIM3*DIM3]={0.0};
	compute_gradH_d_normals (gradient_grid, iv1, d, gradientN_d);

	// compute the gradient of the scalar grid
	GRADIENT_COORD_TYPE gradientS_d[DIM3]={0.0};
	compute_gradH_d_scalar_grid (scalar_grid, iv1, d, gradientS_d);
	float gradientS_d_magnitude=0.0;
	vector_magnitude(gradientS_d, DIM3, gradientS_d_magnitude);
	float gradientS_d_magnitude_sq=abs(gradientS_d_magnitude)*abs(gradientS_d_magnitude);

	gradientS_d[0]=gradientS_d[0]/gradientS_d_magnitude_sq;
	gradientS_d[1]=gradientS_d[1]/gradientS_d_magnitude_sq;
	gradientS_d[2]=gradientS_d[2]/gradientS_d_magnitude_sq;

	for (int i=0; i<DIM3; i++)
		vector_dot_pdt(&(gradientN_d[3*i]), gradientS_d, DIM3, c[i]);

	//DEBUG

	/*
	for (int i=0; i<DIM3; i++)
	{ vector_dot_pdt(&(gradientN_d[3*i]), gradientS_d, DIM3, c[i]); }
	GRADIENT_COORD_TYPE mag_gradS_d = 0.0;

	vector_dot_pdt (gradientS_d, gradientS_d, DIM3, mag_gradS_d);
	for (int k=0; k<DIM3; k++)
	{ c[k]=c[k]/mag_gradS_d; }
	 */


}


// Compute M d for direction 'd' for  vertex iv1

void compute_m_d
(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv1,
		const int d,
		GRADIENT_COORD_TYPE m[DIM3]
)
{
	GRADIENT_COORD_TYPE   c[DIM3] = {0.0};
	GRADIENT_COORD_TYPE   fwd_diff_d_normals[DIM3] = {0.0};
	GRADIENT_COORD_TYPE   fwd_diff_d = 0.0;

	// calculate C for direction d
	compute_c_d(scalar_grid, gradient_grid, iv1, d, c);

	//compute forward difference of normals
	compute_forward_difference_d_normals
	( gradient_grid, iv1, d, fwd_diff_d_normals);

	compute_forward_difference_d
	(scalar_grid, iv1, d, fwd_diff_d);

	for (int i=0; i<DIM3; i++)
	{ m[i] = fwd_diff_d_normals[i] - (fwd_diff_d * c[i]); }
}

// Compute w
// 
void compute_w
(
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const float mu,
		const VERTEX_INDEX iv1, const int flag_aniso,
		const GRADIENT_GRID_BASE & gradient_grid,
		GRADIENT_COORD_TYPE w[DIM3]
)
{
	GRADIENT_COORD_TYPE mX[DIM3]={0.0}, mY[DIM3]={0.0}, mZ[DIM3]={0.0};
	GRADIENT_COORD_TYPE mXprevX[DIM3]={0.0};
	GRADIENT_COORD_TYPE mYprevY[DIM3]={0.0};
	GRADIENT_COORD_TYPE mZprevZ[DIM3]={0.0};

	// Compute M d for direction 'd' for  vertex iv1
	compute_m_d(scalar_grid, gradient_grid, iv1, 0, mX);
	compute_m_d(scalar_grid, gradient_grid, iv1, 1, mY);
	compute_m_d(scalar_grid, gradient_grid, iv1, 2, mZ);

	// compute prev vertex in 0,1,2 direction
	VERTEX_INDEX prev_vert[DIM3];
	prev_vert[0] = scalar_grid.PrevVertex(iv1, 0);
	prev_vert[1] = scalar_grid.PrevVertex(iv1, 1);
	prev_vert[2] = scalar_grid.PrevVertex(iv1, 2);

	// Compute m_d for the previous vertices.
	compute_m_d
	(scalar_grid, gradient_grid, prev_vert[0], 0, mXprevX);
	compute_m_d
	(scalar_grid, gradient_grid, prev_vert[1], 1, mYprevY);
	compute_m_d
	(scalar_grid, gradient_grid, prev_vert[2], 2, mZprevZ);

	// Compute 'k' for each dimensions ,
	// used for anisotropic gradients

	SCALAR_TYPE K[DIM3]={0.0};
	SCALAR_TYPE gK[DIM3]={0.0};

	SCALAR_TYPE gKprev[DIM3]={0.0};

	GRADIENT_COORD_TYPE   mag=0.0;

	// compute_curvature for the present vertex
	compute_curvature_iv(scalar_grid, gradient_grid, iv1, K);

	// compute the gx
	for (int d=0; d<DIM3; d++)
		compute_g_x(mu, K[d], flag_aniso, gK[d]);

	// compute k _d and gkd for the previous vertices
	for (int d=0; d<DIM3; d++) {
		compute_curvature_iv_d
		(scalar_grid, gradient_grid, prev_vert[d], d, K[d]);

		compute_g_x(mu, K[d], flag_aniso, gKprev[d]);
	}
   // DEBUG
	/*
	for (int i=0; i<DIM3; i++) {
		w[i] =  gK[i]*mX[i] - gKprev[i]*mXprevX[i] +
				gK[i]*mY[i] - gKprev[i]*mYprevY[i] +
				gK[i]*mZ[i] - gKprev[i]*mZprevZ[i];
	}
	*/
	for (int i=0; i<DIM3; i++) {
			w[i] =  mX[i] - mXprevX[i] +
					mY[i] - mYprevY[i] +
					mZ[i] - mZprevZ[i];
		}
}

// **************************************************
// VECTOR OPERATORS
// **************************************************

// Calculate the sum of squares of all elements in a vector 'vec'
// of size 'num_elements' and return the 'sum'
void vector_sum_of_squares 
(const float *vec, const int num_elements, float &sum)
{
	for (int i=0; i<num_elements; i++)
	{ sum += vec[i]*vec[i]; }
}

// vector dot pdt
void vector_dot_pdt 
(const float * A, const float *B,
		const int num_elements, float &res)
{
	res = 0.0;
	for (int i=0; i<num_elements; i++)
	{ res = res + A[i]*B[i]; }
}

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

/// Normalize the vectors.
/// Set small magnitude vectors to zero.
void normalize
(float *vec, float & magnitude, const int num_elements, 
		const float max_small_mag)
{
	vector_magnitude(vec, num_elements, magnitude);

	if (magnitude <= max_small_mag) {
		magnitude = 0.0;
		for (int i=0; i<num_elements; i++)
		{ vec[i] = 0.0; }
	}
	else {
		for (int i=0; i<num_elements; i++)
		{ vec[i] = vec[i] / magnitude; }
	}

}

/// Normalize the vectors.
/// Set small magnitude vectors to zero.
void normalize
(float *vec, const int num_elements, const float max_small_mag)
{
	float magnitude;
	normalize(vec, magnitude, num_elements, max_small_mag);
}

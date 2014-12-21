/// \file sharpiso_svd.h
/// Compute sharp isosurface vertices and edges 
///   using singular valued decomposition.

/*
 Copyright (C) 2011-2014 Arindam Bhattacharya
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public License
 (LGPL) as published by the Free Software Foundation; either
 version 2.1 of the License, or any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include<iostream>
#include <Eigen/Dense>
#include "sharpiso_types.h"
#include "sharpiso_eigen.h"

using namespace std;
using namespace Eigen;
using namespace SHARPISO;


const double TOLERANCE = 0.0001;

//FUNCTION DEFINITION
/*
 * inputs:
 * grid_vertex_coords
 * gris_vertex_scalars
 * grid_vertex_gradients
 * Number of grid vertex
 * Isovalue
 * EigenValue Tolerance
 *
 */


void svd_calculate_sharpiso_vertex
(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords,
		GRADIENT_COORD_TYPE *ray_direction);



// Calculate the svd based sharp isovertex but force it to have 2 singular values.
void svd_calculate_sharpiso_vertex_2_svals
(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords,
		GRADIENT_COORD_TYPE *ray_direction);

void svd_calculate_sharpiso_vertex_unit_normals
(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isoVertcoords[DIM3],
		GRADIENT_COORD_TYPE ray_direction[DIM3]);

// Calculate the sharp iso vertex using SVD,
// and the lindstrom approach
// this is called from svd_compute_sharp_vertex_for_cube in sharpiso_feature.cxx
void svd_calculate_sharpiso_vertex_using_lindstrom
(
		const bool useLindstrom2,
		const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords);

// Calculate the sharp iso vertex using SVD and lindstrom approach.
// Input is isosurface-edge intersections.
void svd_calculate_sharpiso_vertex_using_lindstrom(
		const bool useLindstrom2,
		const COORD_TYPE * edgeI_coords, 
    const GRADIENT_COORD_TYPE * edgeI_normal_coord,
		const NUM_TYPE num_intersections,
		const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3], 
		NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
    COORD_TYPE isoVertcoords[DIM3]);

void svd_calculate_sharpiso_vertex_using_lindstrom_fast(
		const NUM_TYPE num_vert,
		const EIGENVALUE_TYPE err_tolerance,
		const SCALAR_TYPE isovalue,
		const SCALAR_TYPE * vert_scalars,
		const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isoVertcoords[DIM3]);

// Calculate the svd based sharp isovertex but force it to have 2 singular values.
void svd_calculate_sharpiso_vertex_2_svals_unit_normals
(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords,
		GRADIENT_COORD_TYPE *ray_direction);


void compute_cube_vertex
(const MatrixXf &A, const RowVectorXf &b, MatrixXf &singular_values,
		const float  err_tolerance, int &num_singular_vals,
		const RowVectorXf &centroid, float * sharp_point);

// Compute Cube Vertex using the modified Lindstrom formula
void compute_cube_vertex_lind2(const MatrixXf &A, const RowVectorXf &b,
		MatrixXf &singular_values, const float err_tolerance,
		int & num_large_sval, const RowVectorXf &centroid, float * sharp_point);

// Compute A inverse using svd
/*
MatrixXf compute_A_inverse
(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance,
 MatrixXf &singularValues, NUM_TYPE & num_singular_vals );
 */
void compute_A_inverse
(const MatrixXf &A, const EIGENVALUE_TYPE  err_tolerance,
		MatrixXf &singularValues, NUM_TYPE & num_singular_vals, MatrixXf & );

// Compute X as Ainverse times B
void compute_X(const MatrixXf &Inv_A, RowVectorXf &B,RowVectorXf &);

// FUNCTION compute w
void calculate_w
(const MatrixXf & inA, const MatrixXf & A, const MatrixXf &I, RowVectorXf &w);

// function to normalize an array
void normalize(const GRADIENT_COORD_TYPE *intial, GRADIENT_COORD_TYPE  *normalized);

void svd_calculate_sharpiso_vertex_edge_based
(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars,
		const NUM_TYPE  num_vert,
		const SCALAR_TYPE isovalue,
		const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE * singular_vals,
		COORD_TYPE * isoVertcoords,
		GRADIENT_COORD_TYPE *ray_direction);

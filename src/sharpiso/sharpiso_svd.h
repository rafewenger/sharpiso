/*
 *  SharpIso_findVert.h
 *  SHARPISO
 *
 *  Created by arindam bhattacharya on 11/28/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
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
		COORD_TYPE * isoVertcoords,
		GRADIENT_COORD_TYPE *ray_direction);

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
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * cubecenter,
		COORD_TYPE * isoVertcoords);



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

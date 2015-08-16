/// \file sharpiso_svd.cxx
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
#include<vector>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/SVD>
using namespace Eigen;


#include "sharpiso_svd.h"
#include "ijkcoord.txx"
#include "ijk.txx"


using namespace std;
using namespace SHARPISO;
using namespace IJK;

// FUNCTION PROTOTYPES

// Compute A where A is g_i's
void compute_A(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert,
		MatrixXf &);

// Compute A where A is g_i's
void compute_A_normalize(const GRADIENT_COORD_TYPE * vert_grads,
		const int num_vert, MatrixXf &);

// Compute B where B is isovalue - s_i;
void compute_B(const COORD_TYPE *vert_cooords,
		const GRADIENT_COORD_TYPE *vert_grads, const SCALAR_TYPE *vert_scalars,
		const int num_vert, const SCALAR_TYPE isovalue, RowVectorXf&);

// Compute B where B is isovalue - s_i;
void compute_B_normalize(const COORD_TYPE *vert_cooords,
		const GRADIENT_COORD_TYPE *vert_grads, const SCALAR_TYPE *vert_scalars,
		const int num_vert, const SCALAR_TYPE isovalue, RowVectorXf&);

// Compute A inverse using svd
void compute_A_inverse_top_2(const MatrixXf &A,
		const EIGENVALUE_TYPE err_tolerance, MatrixXf &singularValues,
		NUM_TYPE & num_singular_vals, MatrixXf&);

// Normalize a vector
void normalize(
	const GRADIENT_COORD_TYPE *initial,
		GRADIENT_COORD_TYPE *normalized) {
	double sum(0.0), mag(0.0);
	for (int i = 0; i < DIM3; i++) {
		sum = sum + initial[i] * initial[i];
	}
	if (sum > 0) {
		mag = sqrt(sum);
		for (int j = 0; j < 3; j++) {
			normalized[j] = initial[j] / mag;
		}
	}

}

// FUNCTION compute x which computes the position x
//   based on x - Ainverse times B
inline void compute_X(const MatrixXf &Inv_A, RowVectorXf &B, RowVectorXf &X) 
{	X = Inv_A * B.transpose(); }


/// Compute pseudo inverse sigma matrix.
/// Set small singular values to zero.
/// @pre max_small_singular_value >= 0.
void compute_pseudo_inverse_sigma
(const Vector3f & singVals, const int num_non_zero_singVals,
 const COORD_TYPE max_small_singular_value_ratio,
 COORD_TYPE singular_vals[DIM3],
 Matrix3f & pseudo_inv_sigma,
 int & num_non_zero_in_pseudo_inv_sigma)
{
  const GRADIENT_COORD_TYPE max_small_singular_value =
    max_small_singular_value_ratio * singVals[0];

  num_non_zero_in_pseudo_inv_sigma = 0;

  pseudo_inv_sigma = Matrix3f::Zero();
  IJK::set_coord_3D(0, singular_vals);

	for (int i=0; i < num_non_zero_singVals; i++)
	{
		if (singVals[i] > max_small_singular_value) {
      singular_vals[i] = singVals[i];
			pseudo_inv_sigma(i,i) =  1.0/singVals[i];
			num_non_zero_in_pseudo_inv_sigma++;
		}
	}
}

/// Compute edge direction.
void compute_edge_direction
(const MatrixXf & A_pseudo_inverse, const MatrixXf & A, 
 RowVector3f & edge_dir) 
{
  Matrix3f I = Matrix3f::Identity();
  Matrix3f M;
  Vector3f magnitude_squared;

  M = (I - A_pseudo_inverse*A);
  for (int i = 0; i < DIM3; i++) {
    magnitude_squared[i] = (M.row(i)).squaredNorm();
  }

  int index = 0;
  for (int i = 1; i < DIM3; i++) 
    if (magnitude_squared[i] > magnitude_squared[index])
      { index = i; }

  edge_dir = M.row(index);
}

/// @param pointX Compute sharp point closest to pointX.
void compute_sharp_point_lindstrom_3x3
(		SCALAR_TYPE * A,
		const SCALAR_TYPE _b[3],
		const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isoVertCoords[DIM3])
{
	Map<const Matrix3f> eigenA(A);
	Map<const RowVector3f> eigen_pX(pointX);
	Map<const RowVector3f> b(_b);
	Matrix3f pseudo_inv_sigma;
	JacobiSVD<Matrix3f> svd(eigenA, ComputeFullU | ComputeFullV);

	// compute the singular values
	Vector3f singVals = svd.singularValues();

  compute_pseudo_inverse_sigma
    (singVals, svd.nonzeroSingularValues(), err_tolerance,
     singular_vals, pseudo_inv_sigma, num_singular_vals);

  // FORMULA CORRECTED: 12-11-2014 by R.W.
	Vector3f sharpCoord = eigen_pX.transpose() +
				svd.matrixV()*pseudo_inv_sigma*(svd.matrixU().transpose())*
				(b.transpose() - eigenA*eigen_pX.transpose());

	// set isoVertCoords
	for (int d=0; d<DIM3; d++) 
    { isoVertCoords[d] = sharpCoord[d]; }
}

/// @param pointX Compute sharp point closest to pointX.
/// @param edge_direction[] If num_singular_vals = 2, return sharp edge direction.
/// @param orth_direction[] If num_singular_vals = 1, return direction orthogonal to surface.
void compute_sharp_point_lindstrom_3x3
(		SCALAR_TYPE * A,
		const SCALAR_TYPE _b[3],
		const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isoVertCoords[DIM3],
    COORD_TYPE edge_direction[DIM3],
    COORD_TYPE orth_direction[DIM3])
{
	Map<const Matrix3f> eigenA(A);
	Map<const RowVector3f> eigen_pX(pointX);
	Map<const RowVector3f> b(_b);
	Matrix3f pseudo_inv_sigma;
	Matrix3f A_pseudo_inv;
	JacobiSVD<Matrix3f> svd(eigenA, ComputeFullU | ComputeFullV);

	// compute the singular values
	Vector3f singVals = svd.singularValues();

  compute_pseudo_inverse_sigma
    (singVals, svd.nonzeroSingularValues(), err_tolerance,
     singular_vals, pseudo_inv_sigma, num_singular_vals);

  A_pseudo_inv = svd.matrixV()*pseudo_inv_sigma*(svd.matrixU().transpose());

  // FORMULA CORRECTED: 12-11-2014 by R.W.
	Vector3f sharpCoord = 
    eigen_pX.transpose() + 
    A_pseudo_inv*(b.transpose() - eigenA*eigen_pX.transpose());

	// set isoVertCoords
	for (int d=0; d<DIM3; d++) 
    { isoVertCoords[d] = sharpCoord[d]; }

  if (num_singular_vals == 1) {
    
    for (int d = 0; d < DIM3; d++) 
      { orth_direction[d] = svd.matrixU()(d,0); }
		normalize(orth_direction, orth_direction);

    IJK::set_coord_3D(0, edge_direction);
  }
  else if (num_singular_vals == 2) {
    RowVector3f edge_dir;
    compute_edge_direction(A_pseudo_inv, eigenA, edge_dir);

    for (int d = 0; d < DIM3; d++)
      { edge_direction[d] = edge_dir[d]; }

		normalize(edge_direction, edge_direction);
    IJK::set_coord_3D(0, orth_direction);
  }
  else {
    IJK::set_coord_3D(0, edge_direction);
    IJK::set_coord_3D(0, orth_direction);
  }
}


/// Compute sharp point closest to pointX.
/// If num_singular_vals == 2, compute sharp point on plane
///   through pointX with normal plan_normal
/// @param edge_direction If num_singular_vals = 2, return sharp edge direction.
void compute_sharp_point_on_plane_lindstrom_3x3
(		SCALAR_TYPE * A,
		const SCALAR_TYPE _b[3],
		const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3],
		const COORD_TYPE plane_normal[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isovert_coord[DIM3],
    bool & flag_coord_on_plane)
{
  const ANGLE_TYPE max_angle_normal_and_edge_dir = 80;
  const COORD_TYPE min_cos_normal_and_edge_dir
    = cos(max_angle_normal_and_edge_dir*M_PI/180.0);
  COORD_TYPE edge_direction[DIM3];
  COORD_TYPE orth_direction[DIM3];
  COORD_TYPE diff[DIM3];

  compute_sharp_point_lindstrom_3x3
    (A, _b, err_tolerance, pointX, num_singular_vals, singular_vals,
     isovert_coord, edge_direction, orth_direction);

  flag_coord_on_plane = false;
  if (num_singular_vals == 2) {

    COORD_TYPE a1, a2;
    IJK::compute_inner_product_3D(plane_normal, edge_direction, a1);

    // Ignore plane if angle between line through plane_normal and line through edge_direction 
    //   is near 90.
    if (abs(a1) >= min_cos_normal_and_edge_dir) {
      IJK::subtract_coord_3D(pointX, isovert_coord, diff);
      IJK::compute_inner_product_3D(plane_normal, diff, a2);
      COORD_TYPE t = a2/a1;
      IJK::add_scaled_coord_3D(t, edge_direction, isovert_coord, isovert_coord);
      flag_coord_on_plane = true;
    }
  }
}

template <typename T1, typename T2>
inline T1 M3x3_index(const T1 i, const T2 j)
{
  return(i*DIM3+j);
}

/// Compute 3x3 matrix A and vector B.
void compute_A_B
(const NUM_TYPE num_vert,
 const EIGENVALUE_TYPE err_tolerance,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE * vert_scalars,
 const COORD_TYPE * vert_coords,
 const GRADIENT_COORD_TYPE * vert_grads,
 COORD_TYPE A[DIM3*DIM3],
 COORD_TYPE B[DIM3])
{
	IJK::PROCEDURE_ERROR error("compute_A_B");

	if (num_vert < 1) {
		error.AddMessage("Programming error. Number of gradients is 0.");
		throw error;
	}

  IJK::set_coord_3D(0, B);
  IJK::set_coord(DIM3*DIM3, 0, A);

	for (int i=0; i<num_vert; i++){

		GRADIENT_COORD_TYPE tempN[DIM3];
		normalize(vert_grads+i*DIM3, tempN);

		A[M3x3_index(0,0)] += tempN[0]*tempN[0];
		A[M3x3_index(1,1)] += tempN[1]*tempN[1];
		A[M3x3_index(2,2)] += tempN[2]*tempN[2];

		A[M3x3_index(0,1)] += (tempN[0]*tempN[1]);
		A[M3x3_index(0,2)] += (tempN[0]*tempN[2]);
		A[M3x3_index(1,2)] += (tempN[1]*tempN[2]);

		// compute isovalue - s_i + g_i*p_i;
		SCALAR_TYPE iprod, x;
		compute_inner_product(DIM3, tempN, vert_coords+i*DIM3, iprod);
		SCALAR_TYPE gradient_magnitude;

		IJK::compute_magnitude_3D(vert_grads+i*DIM3, gradient_magnitude);
	
		x = -1.0*iprod + (vert_scalars[i] - isovalue)/gradient_magnitude;

    IJK::add_scaled_coord_3D(-x, tempN, B, B);
	}

	A[M3x3_index(1,0)] = A[M3x3_index(0,1)];
	A[M3x3_index(2,0)] = A[M3x3_index(0,2)];
	A[M3x3_index(2,1)] = A[M3x3_index(1,2)];
}



/// Calculate the sharp vertex using svd and the faster garland heckbert way
/// of storing normals
/// @param pointX Compute vertex closest to pointX.
void svd_calculate_sharpiso_vertex_using_lindstrom_fast
(		const NUM_TYPE num_vert,
		const EIGENVALUE_TYPE err_tolerance,
		const SCALAR_TYPE isovalue,
		const SCALAR_TYPE * vert_scalars,
		const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isovert_coords[DIM3])
{
  COORD_TYPE A[DIM3*DIM3];
  COORD_TYPE B[DIM3];
  compute_A_B(num_vert, err_tolerance, isovalue, vert_scalars, vert_coords,
              vert_grads, A, B);

	compute_sharp_point_lindstrom_3x3
    (A, B, err_tolerance, pointX, num_singular_vals, 
     singular_vals, isovert_coords);
}

/// Calculate the sharp vertex using svd and the faster garland heckbert way
/// of storing normals
/// @param pointX Compute vertex closest to pointX.
void svd_calculate_sharpiso_vertex_using_lindstrom_fast
(		const NUM_TYPE num_vert,
		const EIGENVALUE_TYPE err_tolerance,
		const SCALAR_TYPE isovalue,
		const SCALAR_TYPE * vert_scalars,
		const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const COORD_TYPE pointX[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isovert_coords[DIM3],
    COORD_TYPE edge_direction[DIM3],
    COORD_TYPE orth_direction[DIM3])
{
  COORD_TYPE A[DIM3*DIM3];
  COORD_TYPE B[DIM3];
  compute_A_B(num_vert, err_tolerance, isovalue, vert_scalars, vert_coords,
              vert_grads, A, B);
	compute_sharp_point_lindstrom_3x3
    (A, B, err_tolerance, pointX, num_singular_vals, 
     singular_vals, isovert_coords, edge_direction, orth_direction);
}

/// Calculate the sharp vertex using svd and the faster garland heckbert way
/// of storing normals.
/// If num_singular_vals is 2, position isovert_coords on plane.
/// @param pointX Compute vertex closest to pointX.
void svd_calculate_sharpiso_vertex_on_plane_using_lindstrom_fast
(		const NUM_TYPE num_vert,
		const EIGENVALUE_TYPE err_tolerance,
		const SCALAR_TYPE isovalue,
		const SCALAR_TYPE * vert_scalars,
		const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const COORD_TYPE pointX[DIM3],
		const COORD_TYPE plane_normal[DIM3],
		NUM_TYPE & num_singular_vals,
		EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE isovert_coords[DIM3],
    bool & flag_coord_on_plane)
{
  COORD_TYPE A[DIM3*DIM3];
  COORD_TYPE B[DIM3];
  compute_A_B(num_vert, err_tolerance, isovalue, vert_scalars, vert_coords,
              vert_grads, A, B);

	compute_sharp_point_on_plane_lindstrom_3x3
    (A, B, err_tolerance, pointX, plane_normal, num_singular_vals, 
     singular_vals, isovert_coords, flag_coord_on_plane);
}


// Calculate the sharp iso vertex using SVD,
// and the lindstrom approach
// this is called from svd_compute_sharp_vertex_for_cube in sharpiso_feature.cxx
/// @param pointX Compute vertex closest to pointX.
void svd_calculate_sharpiso_vertex_using_lindstrom
(		const bool useLindstrom2,
		const COORD_TYPE * vert_coords, const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars, const NUM_TYPE num_vert,
		const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3], 
		NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
    COORD_TYPE isoVertcoords[DIM3])
{

	if (num_vert == 0)
		return;
	// Initialize variables
	for (int i = 0; i < DIM3; i++) {
		singular_vals[i] = 0.0; // set the default singular values to zero.
		isoVertcoords[i] = 0.5; // set the default value to cube center.
	}

	/// Find point x (3 coordinates) such that:
	/// (g_i) cdot (x - p_i) + s_i = isovalue
	// singular values
	MatrixXf singular_values;

	//Compute A where A is g_i's
	//call the version which normalizes the gradients
	MatrixXf A(num_vert, DIM3);
	compute_A_normalize(vert_grads, num_vert, A);
	//Compute B where B is isovalue - s_i + g_i*p_i;
	RowVectorXf B(num_vert);
	compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert,
			isovalue, B);

	RowVectorXf eigen_pX(DIM3);
	eigen_pX << pointX[0], pointX[1], pointX[2];
	// check which one to use lindstrom2 or lindstrom
	if (useLindstrom2 == false)
		// Compute the cube vertex
		compute_cube_vertex
		(A, B, singular_values, err_tolerance,
				num_singular_vals, eigen_pX, isoVertcoords);
	else
		compute_cube_vertex_lind2
		(A, B, singular_values, err_tolerance,
				num_singular_vals, eigen_pX, isoVertcoords);

	//set up singular values. convert from eigen data type to floating type.
	for (int i = 0; i < num_singular_vals; i++) {
		singular_vals[i] = singular_values(i);
	}

}

// Calculate the sharp iso vertex using SVD and lindstrom approach.
// Input is isosurface-edge intersections.
/// @param pointX Compute vertex closest to pointX.
void svd_calculate_sharpiso_vertex_using_lindstrom(
		const bool useLindstrom2,
		const COORD_TYPE * edgeI_coords, 
    const GRADIENT_COORD_TYPE * edgeI_normal_coord,
		const NUM_TYPE num_intersections,
		const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
		const COORD_TYPE pointX[DIM3], 
		NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
    COORD_TYPE * isoVertcoords)
{
  std::vector<SCALAR_TYPE> edgeI_scalar(num_intersections);

  for (NUM_TYPE i = 0; i < num_intersections; i++) 
    { edgeI_scalar[i] = isovalue; }

  svd_calculate_sharpiso_vertex_using_lindstrom
    (useLindstrom2, edgeI_coords, edgeI_normal_coord,
     &(edgeI_scalar[0]), num_intersections, isovalue,
     err_tolerance, pointX, 
     num_singular_vals, singular_vals, isoVertcoords);
}

// calculate the sharp iso vertex using SVD, but normalize the
// gradients
void svd_calculate_sharpiso_vertex_unit_normals
(const COORD_TYPE * vert_coords,
 const GRADIENT_COORD_TYPE * vert_grads,
 const SCALAR_TYPE * vert_scalars, const NUM_TYPE num_vert,
 const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
 NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
 COORD_TYPE isoVertcoords[DIM3], GRADIENT_COORD_TYPE edge_direction[DIM3]) 
{
	if (num_vert == 0)
		return;


	// Initialize variables
	for (int i = 0; i < DIM3; i++) {
		singular_vals[i] = 0.0; // set the default singular values to zero.
		isoVertcoords[i] = 0.5; // set the default value to cube center.
	}

	/// Find point x (3 coordinates) such that:
	/// (g_i) cdot (x - p_i) + s_i = isovalue
	// singular values
	MatrixXf singular_values;

	//Compute A where A is g_i's
	//call the version which normalizes the gradients
	MatrixXf A(num_vert, DIM3);
	compute_A_normalize(vert_grads, num_vert, A);


	//Compute B where B is isovalue - s_i + g_i*p_i;
	RowVectorXf B(num_vert);
	compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert,
			isovalue, B);
	//Compute A inverse using svd
	//MatrixXf inA = compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals);
	MatrixXf inA(DIM3, num_vert);

	compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals, inA);
	//set up singular values. convert from eigen data type to floating type.
	for (int i = 0; i < num_singular_vals; i++) {
		singular_vals[i] = singular_values(i);
	}

	//Compute X as Ainverse times B

	RowVectorXf X(DIM3);
	compute_X(inA, B, X);

	//copy X to isoVertcoords
	for (int i = 0; i < DIM3; i++) {
		isoVertcoords[i] = X(i);
	}

	//if num of singular values is 2 then it must return a direction.
	if (num_singular_vals == 2) {

    RowVector3f edge_dir;
    compute_edge_direction(inA, A, edge_dir);

    for (int d = 0; d < DIM3; d++)
      { edge_direction[d] = edge_dir[d]; }

		normalize(edge_direction, edge_direction);
	}
}

// SVD calculate the sharp vertex using only top 2 singular values
void svd_calculate_sharpiso_vertex_2_svals(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars, const NUM_TYPE num_vert,
		const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords, GRADIENT_COORD_TYPE * edge_direction) {

	if (num_vert == 0)
		return;

	// Initialize variables
	for (int i = 0; i < DIM3; i++) {
		singular_vals[i] = 0.0; // set the default singular values to zero.
		isoVertcoords[i] = 0.5; // set the default value to cube center.
	}

	/// Find point x (3 coordinates) such that:
	/// (g_i) cdot (x - p_i) + s_i = isovalue
	// singular values
	MatrixXf singular_values;

	//Compute A where A is g_i's
	MatrixXf A(num_vert, DIM3);
	compute_A(vert_grads, num_vert, A);
	//Compute B where B is isovalue - s_i + g_i*p_i;
	RowVectorXf B(num_vert);
	compute_B(vert_coords, vert_grads, vert_scalars, num_vert, isovalue, B);

	//Compute A inverse using svd
	MatrixXf inverseA(DIM3, num_vert);
	compute_A_inverse_top_2(A, err_tolerance, singular_values,
			num_singular_vals, inverseA);

	//set up singular values. convert from eigen data type to floating type.
	for (int i = 0; i < num_singular_vals; i++) {
		singular_vals[i] = singular_values(i);
	}

	//Compute X as Ainverse times B

	RowVectorXf X(DIM3);
	compute_X(inverseA, B, X);

	//copy X to isoVertcoords
	for (int i = 0; i < 3; i++) {
		isoVertcoords[i] = X(i);
	}

	//if num of singular values is 2 then it must return a direction.
	if (num_singular_vals == 2) {

		RowVectorXf dir;
		MatrixXf I(3, 3);
		I << 1, 0, 0, 0, 1, 0, 0, 0, 1;

		//Obtaining all solutions of a linear system
		//x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

		//Calculate w
		//RowVectorXf w = calculate_w (inverseA, A, I);
		RowVectorXf w(DIM3);
		calculate_w(inverseA, A, I, w);
		//Calculate [I - A_pseudo_inv A ]w
		dir = (I - inverseA * A) * w.transpose();

		edge_direction[0] = dir[0];
		edge_direction[1] = dir[1];
		edge_direction[2] = dir[2];
		normalize(edge_direction, edge_direction);
	}
}

// normalized version
// SVD calculate the sharp vertex using only top 2 singular values
void svd_calculate_sharpiso_vertex_2_svals_unit_normals
(const COORD_TYPE * vert_coords, const GRADIENT_COORD_TYPE * vert_grads,
 const SCALAR_TYPE * vert_scalars, const NUM_TYPE num_vert,
 const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
 NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
 COORD_TYPE * isoVertcoords, GRADIENT_COORD_TYPE * edge_direction) 
{

	if (num_vert == 0)
		return;

	// Initialize variables
	for (int i = 0; i < DIM3; i++) {
		singular_vals[i] = 0.0; // set the default singular values to zero.
		isoVertcoords[i] = 0.5; // set the default value to cube center.
	}

	/// Find point x (3 coordinates) such that:
	/// (g_i) cdot (x - p_i) + s_i = isovalue
	// singular values
	MatrixXf singular_values;

	//Compute A where A is g_i's
	MatrixXf A(num_vert, DIM3);
	compute_A_normalize(vert_grads, num_vert, A);
	//Compute B where B is isovalue - s_i + g_i*p_i;
	RowVectorXf B(num_vert);
	compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert,
			isovalue, B);
	//Compute A inverse using svd
	MatrixXf inverseA(DIM3, num_vert);
	compute_A_inverse_top_2(A, err_tolerance, singular_values,
			num_singular_vals, inverseA);

	//set up singular values. convert from eigen data type to floating type.
	for (int i = 0; i < num_singular_vals; i++) {
		singular_vals[i] = singular_values(i);
	}

	//Compute X as Ainverse times B
	RowVectorXf X(DIM3);
	compute_X(inverseA, B, X);

	//copy X to isoVertcoords
	for (int i = 0; i < 3; i++) {
		isoVertcoords[i] = X(i);
	}

	//if num of singular values is 2 then it must return a direction.
	if (num_singular_vals == 2) {

		RowVectorXf dir;
		MatrixXf I(3, 3);
		I << 1, 0, 0, 0, 1, 0, 0, 0, 1;

		//Obtaining all solutions of a linear system
		//x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

		//Calculate w
		//RowVectorXf w = calculate_w (inverseA, A, I);
		RowVectorXf w(DIM3);
		calculate_w(inverseA, A, I, w);
		//Calculate [I - A_pseudo_inv A ]w
		dir = (I - inverseA * A) * w.transpose();

		edge_direction[0] = dir[0];
		edge_direction[1] = dir[1];
		edge_direction[2] = dir[2];
		normalize(edge_direction, edge_direction);
	}
}

//FUNCTION Compute A, where A is g_i's
void compute_A(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert,
		MatrixXf & A) {
	for (int i = 0; i < num_vert; i++) {
		for (int j = 0; j < DIM3; j++) {
			int k = DIM3 * i + j;
			A(i, j) = vert_grads[k];
		}
	}
}


//Compute magnitude of a vector of dimension 'dim'

void calculate_magnitude(const GRADIENT_COORD_TYPE *vec, const int dim,
		GRADIENT_COORD_TYPE & result) {
	GRADIENT_COORD_TYPE mag = 0.0;
	for (int i = 0; i < dim; i++) {
		mag = mag + vec[i] * vec[i];
	}
	result = sqrt(mag);
}
//FUNCTION Compute A, where A is g_i's
// normalize the gradients
void compute_A_normalize(const GRADIENT_COORD_TYPE * vert_grads,
		const int num_vert, MatrixXf & A) {
	GRADIENT_COORD_TYPE vert_grads_normalized[DIM3] = { 0.0 };
	GRADIENT_COORD_TYPE magnitude = 0.0;
	for (int i = 0; i < num_vert; i++) {
		//compute magnitude of the gradients
		calculate_magnitude(&(vert_grads[DIM3 * i]), DIM3, magnitude);
		//normalize(&(vert_grads[DIM*i]), vert_grads_normalized);
		for (int j = 0; j < DIM3; j++) {
			int k = DIM3 * i + j;
			if (magnitude > TOLERANCE)
				A(i, j) = vert_grads[k] / magnitude;
			else
				A(i, j) = 0.0;
		}
	}
}
;

// FUNCTION to compute the dot product
void compute_dot_pdt(const GRADIENT_COORD_TYPE * A, const COORD_TYPE * B,
		const int num, float &result) {
	float sum(0.0);
	for (int i = 0; i < num; i++) {
		sum += A[i] * B[i];
	}
	result = sum;
}

// FUNCTION to compute B where B is given as isovalue - c_i + g_i*p_i
void compute_B(const COORD_TYPE *vert_cooords,
		const GRADIENT_COORD_TYPE *vert_grads, const SCALAR_TYPE *vert_scalars,
		const int num_vert, const SCALAR_TYPE isovalue, RowVectorXf &B) {
	for (int i = 0; i < num_vert; i++) {
		float dtpdt(0.0);
		compute_dot_pdt(&(vert_grads[DIM3 * i]), &(vert_cooords[DIM3 * i]),
				DIM3, dtpdt);
		B(i) = isovalue - vert_scalars[i] + dtpdt;
	}

}

// FUNCTION to compute B where B is given as isovalue - c_i + g_i*p_i
// Normalize the gradients

void compute_B_normalize
(const COORD_TYPE *vert_coords,
		const GRADIENT_COORD_TYPE *vert_grads, const SCALAR_TYPE *vert_scalars,
		const int num_vert, const SCALAR_TYPE isovalue, RowVectorXf & B) {
	GRADIENT_COORD_TYPE magnitude;
	for (int i = 0; i < num_vert; i++) {
		//compute magnitude of the gradients
		calculate_magnitude(&(vert_grads[DIM3 * i]), DIM3, magnitude);
		float dtpdt;
		compute_dot_pdt(vert_grads + i*DIM3, vert_coords + i*DIM3, DIM3, dtpdt);
		B(i) = isovalue - vert_scalars[i] + dtpdt;

		if (magnitude > TOLERANCE)
			B(i) = B(i) / magnitude;
		else
			B(i) = 0.0;
	}
}

// Compute sigma using only large singular values
// Precondition: sigma is a square kxk matrix where k=singular_values.rows().
void compute_sigma
(const MatrixXf & singular_values,
		const EIGENVALUE_TYPE error_tolerance,
		const int max_num_singular,
		MatrixXf & sigma_plus,
		MatrixXf & sigma,
		int & num_large_singular)
{
	int num_sval = singular_values.rows();

	// Initialize
	num_large_singular = 0;
	sigma.setZero(num_sval, num_sval);
	sigma_plus.setZero(num_sval, num_sval);

	EIGENVALUE_TYPE max_sval = singular_values.maxCoeff();
	EIGENVALUE_TYPE scaled_error_tolerance = error_tolerance * max_sval;

	if (scaled_error_tolerance <= 0) {
		return;
	}

	int max_num_sval = std::min(num_sval, max_num_singular);

	for (int i = 0; i < max_num_sval; i++) {
		if (singular_values(i) > scaled_error_tolerance) {
			num_large_singular++;
			sigma_plus(i, i) = 1.0 / singular_values(i);
			sigma(i,i) = singular_values(i);
		}
	}

}

// FUNCTION compute the pseudo inverse of A
// helper function to compute_A_inverse.
// it calculates the sigma interms of the given tolerance
// debug change this to error type.
void compute_A_pseudoinverse(const MatrixXf &A, MatrixXf & singular_values,
		const float err_tolerance, int & num_large_singular,
		MatrixXf &pseudoinverseA) {
	// Compute the singular values for the matrix
	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
	singular_values = svd.singularValues();
	int num_sval = singular_values.rows();

	// Compute sigma using only large singular values
	MatrixXf sigma(num_sval, num_sval);
	MatrixXf sigma_plus(num_sval, num_sval);
	compute_sigma(singular_values, err_tolerance, num_sval, sigma_plus, sigma,
			num_large_singular);

	pseudoinverseA = svd.matrixV() * sigma_plus.transpose()
															* svd.matrixU().transpose();
}

// Compute_point for edge based dual contouring
// applying the formula used in  lindstrom
// need not be the centroid , may even be the cube_center
void compute_cube_vertex
(const MatrixXf &A, const RowVectorXf &b,
		MatrixXf &singular_values, const float err_tolerance,
		int & num_large_sval, const RowVectorXf &centroid, float * sharp_point)
{
	// Compute the singular values for the matrix
	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);

	singular_values = svd.singularValues();
	int num_sval = singular_values.rows();

	// Compute sigma using only large singular values
	MatrixXf sigma(num_sval, num_sval);
	MatrixXf sigma_plus(num_sval, num_sval);
	compute_sigma
	(singular_values, err_tolerance, num_sval, sigma_plus, sigma,
			num_large_sval);

	MatrixXf point = centroid.transpose() +
			svd.matrixV() * sigma_plus * svd.matrixU().transpose() *
			(b.transpose() - A * centroid.transpose());

	for (int i = 0; i < 3; i++) { sharp_point[i] = point(i); }
}

// Compute_point for edge based dual contouring
// applying the formula used in  lindstrom2
// sharp point = centroid + V*Sigma * U^t *(b - A'*centroid)
// where A' = U*Sigma *V^t
void compute_cube_vertex_lind2
(const MatrixXf &A, const RowVectorXf &b,
		MatrixXf &singular_values,
		const float err_tolerance,
		int & num_large_sval,
		const RowVectorXf &centroid,
		float * sharp_point)
{
	// Compute the singular values for the matrix
	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);

	singular_values = svd.singularValues();
	int num_sval = singular_values.rows();

	// Compute sigma using only large singular values
	MatrixXf sigma(num_sval, num_sval);
	MatrixXf sigma_plus(num_sval, num_sval);
	compute_sigma(singular_values, err_tolerance, num_sval, sigma_plus, sigma,
			num_large_sval);
	MatrixXf A2 = svd.matrixU()*sigma*svd.matrixV().transpose();
	//
	MatrixXf point = centroid.transpose() +
			svd.matrixV() * sigma_plus * svd.matrixU().transpose() *
			(b.transpose() - A2 * centroid.transpose());
	for (int i = 0; i < 3; i++) { sharp_point[i] = point(i); }
}

///
// FUNCTION compute the pseudo inverse of A using the TOP 2 singular values.
// helper function to compute_A_inverse.
// it calculates the sigma interms of the given tolerance
void compute_A_pseudoinverse_top_2(const MatrixXf &A,
		MatrixXf &singular_values, const float err_tolerance,
		int & num_large_singular, MatrixXf & inA) {
	// Compute the singular values for the matrix
	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
	singular_values = svd.singularValues();
	int num_sval = singular_values.rows();
	// Compute sigma using only large singular values.
	// Use at most two singular values.
	MatrixXf sigma_plus(num_sval, num_sval);
	MatrixXf sigma(num_sval, num_sval);
	compute_sigma(singular_values, err_tolerance, 2, sigma_plus, sigma, num_large_singular);

	inA = svd.matrixV() * sigma_plus.transpose() * svd.matrixU().transpose();
}

// FUNCTION compute the inverse of A,
// accepts as input the matrix A , it calculates the singular values and the number of singular
// values above the user set tolerance

void compute_A_inverse(const MatrixXf &A, const EIGENVALUE_TYPE err_tolerance,
		MatrixXf &singularValues, NUM_TYPE & num_singular_vals,
		MatrixXf &pseudoInverseA) {
	compute_A_pseudoinverse(A, singularValues, err_tolerance,
			num_singular_vals, pseudoInverseA);
}

// FUNCTION compute the inverse of A,
// accepts as input the matrix A , it calculates the singular values and the number of singular
// values above the user set tolerance

void compute_A_inverse_top_2(const MatrixXf &A,
		const EIGENVALUE_TYPE err_tolerance, MatrixXf &singularValues,
		NUM_TYPE & num_singular_vals, MatrixXf& pseudoInverseA) {
	compute_A_pseudoinverse_top_2(A, singularValues, err_tolerance,
			num_singular_vals, pseudoInverseA);
}

// Calculate the magnitude
SCALAR_TYPE calculate_mag(const RowVectorXf &res) {
	SCALAR_TYPE sum = 0.0;
	for (int i = 0; i < 3; i++) {
		sum += res(i) * res(i);
	}
	return sqrt(sum);
}

// FUNCTION compute w
void calculate_w
(const MatrixXf & inA, const MatrixXf & A, const MatrixXf & I,
 RowVectorXf & w) 
{
	vector<RowVectorXf> e;
	RowVectorXf e1(3);
	e1 << 1, 0, 0;
	e.push_back(e1);
	e1 << 0, 1, 0;
	e.push_back(e1);
	e1 << 0, 0, 1;
	e.push_back(e1);

	RowVectorXf res;
	NUM_TYPE index;
	SCALAR_TYPE mag = 0.0;
	SCALAR_TYPE max = -1.0;
	for (int i = 0; i < e.size(); i++) {
		res = (I - inA * A) * e[i].transpose();

		mag = calculate_mag(res);
		if (mag > max) {
			max = mag;
			index = i;
		}
	}
	w = e[index];
}

////FUNCTION DEFINITION

//// Inputs:
//// grid_vertex_coords
//// gris_vertex_scalars
//// grid_vertex_gradients
//// Number of grid vertex
//// Isovalue
//// EigenValue Tolerance

///This function calculates the
//singular values and the isovert predicted by the singular value.
//if the number of singular values is 2 , it also calculates the direction of the sharp edge.
//It does NOT compute the final vertex, which is computed under
//'svd_compute_sharp_vertex_neighborhood_S'
///

void svd_calculate_sharpiso_vertex(const COORD_TYPE * vert_coords,
		const GRADIENT_COORD_TYPE * vert_grads,
		const SCALAR_TYPE * vert_scalars, const NUM_TYPE num_vert,
		const SCALAR_TYPE isovalue, const EIGENVALUE_TYPE err_tolerance,
		NUM_TYPE & num_singular_vals, EIGENVALUE_TYPE singular_vals[DIM3],
		COORD_TYPE * isoVertcoords, GRADIENT_COORD_TYPE * edge_direction) {

	if (num_vert == 0)
		return;

	// Initialize variables
	for (int i = 0; i < DIM3; i++) {
		singular_vals[i] = 0.0; // set the default singular values to zero.
		isoVertcoords[i] = 0.5; // set the default value to cube center.
	}

	//Find point x (3 coordinates) such that:
	//(g_i) cdot (x - p_i) + s_i = isovalue
	//singular values
	MatrixXf singular_values;

	//Compute A where A is g_i's
	MatrixXf A(num_vert, DIM3);
	compute_A(vert_grads, num_vert, A);

	//Compute B where B is isovalue - s_i + g_i*p_i;
	RowVectorXf B(num_vert);
	compute_B(vert_coords, vert_grads, vert_scalars, num_vert, isovalue, B);

	//Compute A inverse using svd
	MatrixXf inA(DIM3, num_vert);
	compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals, inA);


	//set up singular values. convert from eigen data type to floating type.
	for (int i = 0; i < num_singular_vals; i++) {
		singular_vals[i] = singular_values(i);
	}

	//Compute X as Ainverse times B
	RowVectorXf X(DIM3);
	compute_X(inA, B, X);

	//copy X to isoVertcoords
	for (int i = 0; i < 3; i++) {
		isoVertcoords[i] = X(i);
	}

	//if num of singular values is 2 then it must return a direction.
	if (num_singular_vals == 2) {

		RowVectorXf dir;
		MatrixXf I(3, 3);
		I << 1, 0, 0, 0, 1, 0, 0, 0, 1;

		//Obtaining all solutions of a linear system
		//x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

		//Calculate w
		//RowVectorXf w = calculate_w (inA, A, I);
		RowVectorXf w(DIM3);
		calculate_w(inA, A, I, w);
		//Calculate [I - A_pseudo_inv A ]w
		dir = (I - inA * A) * w.transpose();

		edge_direction[0] = dir[0];
		edge_direction[1] = dir[1];
		edge_direction[2] = dir[2];
		normalize(edge_direction, edge_direction);
	}

}


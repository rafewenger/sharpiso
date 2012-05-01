/*
*  sharpiso_svd.cxx
*  SHARPISO
*
*  Created by arindam bhattacharya on 11/28/11.
*  Copyright 2011 Ohio State University. All rights reserved.
*
*/

#include<iostream>
#include<vector>
#include <cmath>

#include "sharpiso_svd.h"

using namespace std;
using namespace SHARPISO;

// FUNCTION PROTOTYPES

// Compute A where A is g_i's
void compute_A
(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert, MatrixXf &);

// Compute A where A is g_i's
void compute_A_normalize
(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert, MatrixXf &);

// Compute B where B is isovalue - s_i;
void compute_B
(const COORD_TYPE *vert_cooords,
 const GRADIENT_COORD_TYPE *vert_grads,
 const  SCALAR_TYPE *vert_scalars, const int num_vert,
 const SCALAR_TYPE isovalue,
 RowVectorXf&);

// Compute B where B is isovalue - s_i;
void compute_B_normalize
(const COORD_TYPE *vert_cooords,
 const GRADIENT_COORD_TYPE *vert_grads,
 const  SCALAR_TYPE *vert_scalars, const int num_vert,
 const SCALAR_TYPE isovalue,
 RowVectorXf&);

// Compute A inverse using svd
void  compute_A_inverse_top_2
(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance,
 MatrixXf &singularValues, NUM_TYPE & num_singular_vals,
 MatrixXf& );

// Normalize a vector

void normalize(const GRADIENT_COORD_TYPE *intial,GRADIENT_COORD_TYPE  *normalized)
{
  double sum(0.0),mag(0.0);
  for (int i=0; i<DIM3; i++) {
    sum = sum + intial[i]*intial[i];
  }
  if (sum>0) {
    mag = sqrt(sum);
    for (int j=0; j<3; j++) {
      normalized[j] = intial[j] / mag;
    }
  }
}

// calculate the sharp iso vertex using SVD, 
// and the lindstrom approach
// this is called from svd_compute_sharp_vertex_for_cube in sharpiso_feature.cxx
void svd_calculate_sharpiso_vertex_using_lindstrom
(const COORD_TYPE * vert_coords,
 const GRADIENT_COORD_TYPE * vert_grads,
 const SCALAR_TYPE * vert_scalars,
 const NUM_TYPE  num_vert,
 const SCALAR_TYPE isovalue,
 const EIGENVALUE_TYPE err_tolerance,
 NUM_TYPE & num_singular_vals,
 EIGENVALUE_TYPE singular_vals[DIM3],
 COORD_TYPE * cubecenter,
 COORD_TYPE * isoVertcoords)
{
  if (num_vert == 0)
    return;
  // Initialize variables
  for (int i=0; i<DIM3; i++)
  {
    singular_vals[i] = 0.0;   // set the default singular values to zero.
    isoVertcoords[i] = 0.5;   // set the default value to cube center.
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
  compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert, isovalue, B);
  RowVectorXf eigen_cubecenter(DIM3);
  eigen_cubecenter<<  cubecenter[0],cubecenter[1],cubecenter[2];
  
  // compute the cube vertex
  compute_cube_vertex
    ( A, B, singular_values, err_tolerance, num_singular_vals, eigen_cubecenter, isoVertcoords); 
  
  //set up singular values. convert from eigen data type to floating type.
  for (int i=0; i<num_singular_vals; i++)
  { singular_vals[i]  = singular_values(i); }

}


// calculate the sharp iso vertex using SVD, but normalize the
// gradients
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
 GRADIENT_COORD_TYPE * ray_direction)
{
  if (num_vert == 0)
    return;

  // Initialize variables
  for (int i=0; i<DIM3; i++)
  {
    singular_vals[i] = 0.0;   // set the default singular values to zero.
    isoVertcoords[i] = 0.5;   // set the default value to cube center.
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
  compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert, isovalue, B);
  //Compute A inverse using svd
  //MatrixXf inA = compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals);
  MatrixXf inA(DIM3, num_vert);



  compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals, inA);
  //set up singular values. convert from eigen data type to floating type.
  for (int i=0; i<num_singular_vals; i++)
  { singular_vals[i]  = singular_values(i); }

  //Compute X as Ainverse times B

  RowVectorXf X(DIM3);
  compute_X(inA, B, X);

  //copy X to isoVertcoords
  for (int i=0; i<DIM3; i++)
  { isoVertcoords[i] = X(i); }

  //if num of singular values is 2 then it must return a direction.
  if(num_singular_vals == 2){

    RowVectorXf dir;
    MatrixXf I(3,3);
    I<< 1,0,0, 0,1,0, 0,0,1;

    //Obtaining all solutions of a linear system
    //x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

    //Calculate w
    //RowVectorXf w = calculate_w (inA, A, I);
    RowVectorXf w (DIM3);
    calculate_w (inA, A, I, w);

    //Calculate [I - A_pseudo_inv A ]w
    dir = ( I - inA * A ) * w.transpose();

    ray_direction[0] = dir[0];
    ray_direction[1] = dir[1];
    ray_direction[2] = dir[2];
    normalize(ray_direction, ray_direction);
  }
}

// SVD calculate the sharp vertex using only top 2 singular values
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
 GRADIENT_COORD_TYPE * ray_direction)
{

  if (num_vert == 0)
    return;

  // Initialize variables
  for (int i=0; i<DIM3; i++)
  {
    singular_vals[i] = 0.0;   // set the default singular values to zero.
    isoVertcoords[i] = 0.5;   // set the default value to cube center.
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
  MatrixXf inverseA( DIM3, num_vert);
  compute_A_inverse_top_2(A, err_tolerance, singular_values, num_singular_vals, inverseA);

  //set up singular values. convert from eigen data type to floating type.
  for (int i=0; i<num_singular_vals; i++)
  { singular_vals[i]  = singular_values(i); }

  //Compute X as Ainverse times B

  RowVectorXf X(DIM3);
  compute_X(inverseA, B, X);

  //copy X to isoVertcoords
  for (int i=0; i<3; i++)
  { isoVertcoords[i] = X(i); }

  //if num of singular values is 2 then it must return a direction.
  if(num_singular_vals == 2){

    RowVectorXf dir;
    MatrixXf I(3,3);
    I<< 1,0,0, 0,1,0, 0,0,1;

    //Obtaining all solutions of a linear system
    //x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

    //Calculate w
    //RowVectorXf w = calculate_w (inverseA, A, I);
    RowVectorXf w (DIM3);
    calculate_w (inverseA, A, I, w);
    //Calculate [I - A_pseudo_inv A ]w
    dir = ( I - inverseA * A ) * w.transpose();

    ray_direction[0] = dir[0];
    ray_direction[1] = dir[1];
    ray_direction[2] = dir[2];
    normalize(ray_direction, ray_direction);
  }
}


// normalized version
// SVD calculate the sharp vertex using only top 2 singular values
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
 GRADIENT_COORD_TYPE * ray_direction)
{

  if (num_vert == 0)
    return;

  // Initialize variables
  for (int i=0; i<DIM3; i++)
  {
    singular_vals[i] = 0.0;   // set the default singular values to zero.
    isoVertcoords[i] = 0.5;   // set the default value to cube center.
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
  compute_B_normalize(vert_coords, vert_grads, vert_scalars, num_vert, isovalue, B);
  //Compute A inverse using svd
  MatrixXf inverseA( DIM3, num_vert);
  compute_A_inverse_top_2(A, err_tolerance, singular_values, num_singular_vals, inverseA);

  //set up singular values. convert from eigen data type to floating type.
  for (int i=0; i<num_singular_vals; i++)
  { singular_vals[i]  = singular_values(i); }

  //Compute X as Ainverse times B
  RowVectorXf X(DIM3);
  compute_X(inverseA, B, X);

  //copy X to isoVertcoords
  for (int i=0; i<3; i++)
  { isoVertcoords[i] = X(i); }

  //if num of singular values is 2 then it must return a direction.
  if(num_singular_vals == 2){

    RowVectorXf dir;
    MatrixXf I(3,3);
    I<< 1,0,0, 0,1,0, 0,0,1;

    //Obtaining all solutions of a linear system
    //x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

    //Calculate w
    //RowVectorXf w = calculate_w (inverseA, A, I);
    RowVectorXf w (DIM3);
    calculate_w (inverseA, A, I, w);
    //Calculate [I - A_pseudo_inv A ]w
    dir = ( I - inverseA * A ) * w.transpose();

    ray_direction[0] = dir[0];
    ray_direction[1] = dir[1];
    ray_direction[2] = dir[2];
    normalize(ray_direction, ray_direction);
  }
}






//FUNCTION Compute A, where A is g_i's
void compute_A(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert, MatrixXf & A)
{
  for (int i=0; i<num_vert; i++) {
    for (int j=0; j<DIM3; j++) {
      int k = DIM3*i + j;
      A(i,j) = vert_grads[k];
    }
  }
};

//Compute magnitude of a vector of dimension 'dim'

void calculate_magnitude
(const GRADIENT_COORD_TYPE *vec, const int dim, GRADIENT_COORD_TYPE & result)
{
  GRADIENT_COORD_TYPE mag = 0.0;
  for (int i=0; i<dim; i++) {
    mag = mag + vec[i]*vec[i];
  }
  result = sqrt(mag);
}
//FUNCTION Compute A, where A is g_i's
// normalize the gradients
void compute_A_normalize
(const GRADIENT_COORD_TYPE * vert_grads,
 const int num_vert,
 MatrixXf & A)
{
  GRADIENT_COORD_TYPE vert_grads_normalized[DIM3]={0.0};
  GRADIENT_COORD_TYPE magnitude = 0.0;
  for (int i=0; i<num_vert; i++)
  {
    //compute magnitude of the gradients
    calculate_magnitude (&(vert_grads[DIM3*i]), DIM3, magnitude);
    //normalize(&(vert_grads[DIM*i]), vert_grads_normalized);
    for (int j=0; j<DIM3; j++)
    {
      int k = DIM3*i + j;
      if (magnitude > TOLERANCE)
        A(i,j) = vert_grads[k] / magnitude;
      else
        A(i,j) = 0.0;
    }
  }
};


// FUNCTION to compute the dot product
void compute_dot_pdt
(const GRADIENT_COORD_TYPE * A, const COORD_TYPE * B , const int num, float &result )
{
  float sum(0.0);
  for (int i=0; i<num; i++) {
    sum+=A[i]*B[i];
  }
  result=sum;
}



// FUNCTION to compute B where B is given as isovalue - c_i + g_i*p_i
void compute_B
(const COORD_TYPE *vert_cooords, const GRADIENT_COORD_TYPE *vert_grads,
 const  SCALAR_TYPE *vert_scalars, const int num_vert, const SCALAR_TYPE isovalue,
 RowVectorXf &B)
{
  for (int i=0; i<num_vert; i++) {
    float dtpdt(0.0);
    compute_dot_pdt(&(vert_grads[DIM3*i]), &(vert_cooords[DIM3*i]), DIM3, dtpdt );
    B(i) = isovalue - vert_scalars[i] + dtpdt;
  }

}

// FUNCTION to compute B where B is given as isovalue - c_i + g_i*p_i
// Normalize the gradients

void compute_B_normalize
(const COORD_TYPE *vert_cooords, const GRADIENT_COORD_TYPE *vert_grads,
 const  SCALAR_TYPE *vert_scalars, const int num_vert, const SCALAR_TYPE isovalue,
 RowVectorXf & B)
{ 
  GRADIENT_COORD_TYPE vert_grads_normalized[DIM3]={0.0};
  GRADIENT_COORD_TYPE magnitude=0.0;
  for (int i=0; i<num_vert; i++)
  {
    //compute magnitude of the gradients
    calculate_magnitude (&(vert_grads[DIM3*i]), DIM3, magnitude);
    float dtpdt(0.0);
    compute_dot_pdt(&(vert_grads[DIM3*i]), &(vert_cooords[DIM3*i]), DIM3, dtpdt );
    B(i) = isovalue - vert_scalars[i] + dtpdt;

    if (magnitude > TOLERANCE)
      B(i) = B(i) / magnitude;
    else
      B(i) = 0.0;
  }
}

// Compute sigma using only large singular values
// Precondition: sigma is a square kxk matrix where k=singular_values.rows().
void compute_sigma(const MatrixXf & singular_values,
                   const EIGENVALUE_TYPE error_tolerance,
                   const int max_num_singular,
                   MatrixXf & sigma,
                   int & num_large_singular)
{
  int num_sval = singular_values.rows();

  // Initialize
  num_large_singular = 0;
  sigma.setZero(num_sval, num_sval);

  EIGENVALUE_TYPE max_sval = singular_values.maxCoeff();
  EIGENVALUE_TYPE scaled_error_tolerance = error_tolerance * max_sval;

  if (scaled_error_tolerance <= 0) { return; }

  int max_num_sval = std::min(num_sval, max_num_singular);

  for (int i = 0; i < max_num_sval; i++) {
    if (singular_values(i) > scaled_error_tolerance) {
      num_large_singular++;
      sigma(i,i) = 1.0/singular_values(i);
    }
  }

}

// FUNCTION compute the pseudo inverse of A
// helper function to compute_A_inverse.
// it calculates the sigma interms of the given tolerance
// debug change this to error type.
void compute_A_pseudoinverse
(const MatrixXf A,
 MatrixXf & singular_values,
 const float  err_tolerance,
 int & num_large_singular,
 MatrixXf &pseudoinverseA)
{
  // Compute the singular values for the matrix
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  singular_values  = svd.singularValues();
  int num_sval = singular_values.rows();

  // Compute sigma using only large singular values
  MatrixXf sigma(num_sval, num_sval);
  compute_sigma(singular_values, err_tolerance, num_sval,
                sigma, num_large_singular);

  pseudoinverseA = svd.matrixV()*sigma.transpose()*svd.matrixU().transpose();
}

// compute_point for edge based dual contouring
// applying the formula used in  lindstrom
// need not be the centroid , may even be the cube_center
void compute_cube_vertex
(const MatrixXf A,
 const RowVectorXf b,
 MatrixXf &singular_values,
 const float  err_tolerance,
 int & num_large_sval,
 const RowVectorXf centroid,
 float * sharp_point)
{
  // Compute the singular values for the matrix
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  singular_values  = svd.singularValues();
  int num_sval = singular_values.rows();

  // Compute sigma using only large singular values
  MatrixXf sigma(num_sval, num_sval);
  compute_sigma(singular_values, err_tolerance, num_sval,
                sigma, num_large_sval);

  MatrixXf point = centroid.transpose() + svd.matrixV()*sigma*
    svd.matrixU().transpose()*(b.transpose() - A*centroid.transpose());

  for (int i=0;i<3;i++)
    { sharp_point[i]=point(i); }
}

// FUNCTION compute the pseudo inverse of A using the TOP 2 singular values.
// helper function to compute_A_inverse.
// it calculates the sigma interms of the given tolerance
void  compute_A_pseudoinverse_top_2
(const MatrixXf A,
 MatrixXf &singular_values,
 const float  err_tolerance,
 int & num_large_singular,
 MatrixXf & inA)
{
  // Compute the singular values for the matrix
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  singular_values  = svd.singularValues();
  int num_sval = singular_values.rows();

  // Compute sigma using only large singular values.
  // Use at most two singular values.
  MatrixXf sigma(num_sval, num_sval);
  compute_sigma(singular_values, err_tolerance, 2,
                sigma, num_large_singular);

  inA =  svd.matrixV()*sigma.transpose()*svd.matrixU().transpose();
}

// FUNCTION compute the inverse of A,
// accepts as input the matrix A , it calculates the singular values and the number of singular
// values above the user set tolerance

void compute_A_inverse(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance,
                       MatrixXf &singularValues, NUM_TYPE & num_singular_vals, MatrixXf &pseudoInverseA )
{
  compute_A_pseudoinverse(A, singularValues, err_tolerance,
    num_singular_vals, pseudoInverseA);
}


// FUNCTION compute the inverse of A,
// accepts as input the matrix A , it calculates the singular values and the number of singular
// values above the user set tolerance

void compute_A_inverse_top_2
(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance,
 MatrixXf &singularValues, NUM_TYPE & num_singular_vals,
 MatrixXf& pseudoInverseA )
{
  // in this case we are actually computing the pseudo inverse
  compute_A_pseudoinverse_top_2(A, singularValues, err_tolerance,
    num_singular_vals, pseudoInverseA);
}


// FUNCTION compute x which computes the position x based on x - Ainverse times B
void compute_X(const MatrixXf Inv_A, RowVectorXf B, RowVectorXf &X)
{
  // compute the vector X
  X = Inv_A*B.transpose();
}

// Calculate the magnitude
SCALAR_TYPE calculate_mag(const RowVectorXf &res)
{
  SCALAR_TYPE sum=0.0;
  for (int i=0; i<3; i++) {
    sum += res(i)*res(i);
  }
  return sqrt(sum);
}

// FUNCTION compute w
void calculate_w
(const MatrixXf & inA,
 const MatrixXf & A,
 const MatrixXf &I,
 RowVectorXf& w)
{
  vector<RowVectorXf> e;
  RowVectorXf e1(3);
  e1<<1,0,0;
  e.push_back(e1);
  e1<<0,1,0;
  e.push_back(e1);
  e1<<0,0,1;
  e.push_back(e1);

  RowVectorXf res;
  NUM_TYPE index;
  SCALAR_TYPE mag = 0.0;
  SCALAR_TYPE max = -1.0;
  for (int i=0; i<e.size(); i++) {
    res = (I - inA*A)*e[i].transpose();

    mag = calculate_mag(res);
    if (mag>max) {
      max=mag;
      index=i;
    }
  }
  w=e[index];
}


////FUNCTION DEFINITION

//// Inputs:
//// grid_vertex_coords
//// gris_vertex_scalars
//// grid_vertex_gradients
//// Number of grid vertex
//// Isovalue
//// EigenValue Tolerance

///* this function calculates the
//singular values and the isovert predicted by the singular value.
//if the number of singular values is 2 , it also calculates the direction of the ray.
//It does NOT compute the final vertex, which is computed under
//'svd_compute_sharp_vertex_neighborhood_S'
//*/

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
 GRADIENT_COORD_TYPE * ray_direction)
{
  if (num_vert == 0)
    return;

  // Initialize variables
  for (int i=0; i<DIM3; i++)
  {
    singular_vals[i] = 0.0;   // set the default singular values to zero.
    isoVertcoords[i] = 0.5;   // set the default value to cube center.
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
  MatrixXf inA(DIM3,num_vert);
  compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals, inA);

  //set up singular values. convert from eigen data type to floating type.
  for (int i=0; i<num_singular_vals; i++)
  { singular_vals[i]  = singular_values(i); }

  //Compute X as Ainverse times B
  RowVectorXf X(DIM3);
  compute_X(inA, B, X);

  //copy X to isoVertcoords
  for (int i=0; i<3; i++)
  { isoVertcoords[i] = X(i); }

  //if num of singular values is 2 then it must return a direction.
  if(num_singular_vals == 2){

    RowVectorXf dir;
    MatrixXf I(3,3);
    I<< 1,0,0, 0,1,0, 0,0,1;

    //Obtaining all solutions of a linear system
    //x = A_pseudo_inv b + [I - A_pseudo_inv A ]w

    //Calculate w
    //RowVectorXf w = calculate_w (inA, A, I);
    RowVectorXf w(DIM3);
    calculate_w(inA, A, I, w);
    //Calculate [I - A_pseudo_inv A ]w
    dir = ( I - inA * A ) * w.transpose();

    ray_direction[0] = dir[0];
    ray_direction[1] = dir[1];
    ray_direction[2] = dir[2];
    normalize(ray_direction, ray_direction);
  }

}


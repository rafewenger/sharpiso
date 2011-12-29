/*
 *  sharpiso_svd.cxx
 *  SHARPISO
 *
 *  Created by arindam bhattacharya on 11/28/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */

#include "sharpiso_svd.h"

#include<iostream>
using namespace std;
const int DIM(3);
  //FUNCTION PROTOTYPES

  //Compute A where A is g_i's
MatrixXf compute_A(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert);
  //Compute B where B is isovalue - s_i;
RowVectorXf compute_B(const COORD_TYPE *vert_cooords, const GRADIENT_COORD_TYPE *vert_grads,const  SCALAR_TYPE *vert_scalars, const int num_vert, 
                      const SCALAR_TYPE isovalue);
  //Compute A inverse using svd
MatrixXf compute_A_inverse(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance, 
                           MatrixXf &singularValues, NUM_TYPE & num_singular_vals );
  //Compute X as Ainverse times B
RowVectorXf compute_X(const MatrixXf Inv_A, RowVectorXf B);

  //FUNNCTION DEFINITION 
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
  //sharp iso eigen as the include file.
  //eigen vaue type and eigen value coord type !


void svd_calculate_sharpiso_vertex

(const COORD_TYPE * vert_coords, 
 const GRADIENT_COORD_TYPE * vert_grads,
 const SCALAR_TYPE * vert_scalars,
 const NUM_TYPE  num_vert,
 const SCALAR_TYPE isovalue,
 const EIGENVALUE_TYPE err_tolerance,
 NUM_TYPE & num_singular_vals,
 EIGENVALUE_TYPE * singular_vals,
 COORD_TYPE * isoVertcoords){
  

    /// Find point x (3 coordinates) such that:
    /// (g_i) cdot (x - p_i) + s_i = isovalue
    // singular values
  
  MatrixXf singular_values;

    //Compute A where A is g_i's
  MatrixXf A = compute_A(vert_grads, num_vert);
  
    //Compute B where B is isovalue - s_i + g_i*p_i;
  RowVectorXf B =  compute_B(vert_coords, vert_grads, vert_scalars, num_vert, isovalue);
 
    //Compute A inverse using svd
  MatrixXf inA = compute_A_inverse(A, err_tolerance, singular_values, num_singular_vals);

    //set up singular values. convert from egein data type to normal type.
  for (int i=0; i<num_singular_vals; i++)
    singular_vals[i]  = singular_values(i);
  
    //Compute X as Ainverse times B
  RowVectorXf X = compute_X(inA, B);


    //copy X to isoVertcoords
  if (num_singular_vals == 3 ) {
    for (int i=0; i<3; i++) {
      isoVertcoords[i] = X(i);
    }
  }
}



  //FUNCTION Compute A, where A is g_i's

MatrixXf compute_A(const GRADIENT_COORD_TYPE * vert_grads, const int num_vert){
    //declare the A matrix;
  MatrixXf A(num_vert, DIM);
  
  for (int i=0; i<num_vert; i++) {
    for (int j=0; j<DIM; j++) {
      int k = DIM*i + j;
      A(i,j) = vert_grads[k];
    }
  }
  return A;
};


  //FUNCTION to compute the dot product 
GRADIENT_COORD_TYPE compute_dot_pdt(const GRADIENT_COORD_TYPE * A, const COORD_TYPE * B , const int num)
{
  double sum(0.0);
  for (int i=0; i<num; i++) {
    sum+=A[i]*B[i];
  }
  return sum;
}

  //FUNCTION to compute B where B is given as isovalue - c_i + g_i*p_i
RowVectorXf compute_B(const COORD_TYPE *vert_cooords, const GRADIENT_COORD_TYPE *vert_grads,const  SCALAR_TYPE *vert_scalars, const int num_vert, 
const SCALAR_TYPE isovalue)
{
  RowVectorXf B(num_vert);
  for (int i=0; i<num_vert; i++) {
    
    B(i) = isovalue - vert_scalars[i] + compute_dot_pdt(&(vert_grads[DIM*i]), &(vert_cooords[DIM*i]), DIM ); 
  }
  return B;
}
  //FUNCTION compute the pseudo inverse of A 
  //helper function to compute_A_inverse.
  //it calcualates the sigma interms of the given tolerance  
MatrixXf compute_A_pseudoinverse(const MatrixXf A, MatrixXf &singular_values,
                                 const float  err_tolerance, int &num_singular_vals)
{
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
    //compute the singular values for the matrix
  singular_values  = svd.singularValues();
  
    //Compute the sigma interms of the tolerance
  MatrixXf sigma(3,3);
  sigma << 0,0,0,0,0,0,0,0,0;
  
  for (int i=0; i<DIM; i++) {
    if (singular_values(i)>err_tolerance) {
        //increment the number of singular values of A
      num_singular_vals++;
        //sigma(i,i) is updated to 1 / the singulat value,
        //only if it is more than the error tolerance
      sigma(i,i) = 1.0/singular_values(i);
    }
  }
  return svd.matrixV()*sigma.transpose()*svd.matrixU().transpose();
}

  //FUNCTION compute the inverse of A,
/* accepts as input the matrix A , it calculates the singular values and the number of singular
 * values above the user set tolerance
 */
MatrixXf compute_A_inverse(const MatrixXf A, const EIGENVALUE_TYPE  err_tolerance, 
                           MatrixXf &singularValues, NUM_TYPE & num_singular_vals )
{
  MatrixXf pseudoInverseA = compute_A_pseudoinverse(A, singularValues, err_tolerance, num_singular_vals); 
  return pseudoInverseA;
}


  //FUNCTION compute x which computes the position x based on x - Ainverse times B
RowVectorXf compute_X(const MatrixXf Inv_A, RowVectorXf B)
{
    //compute the vector X
   RowVectorXf x = Inv_A*B.transpose();
  return x;
  
}

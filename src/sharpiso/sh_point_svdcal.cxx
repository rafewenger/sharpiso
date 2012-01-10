/*
 *  sh_point_svdcal.cxx
 *  SHARPPOINT
 *
 *  Created by arindam bhattacharya on 11/6/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */

#include "sh_point_svdcal.h"
#include "sh_point_datastruct.h"
#include "sharpiso_findIntersect.h"
#include <vector>

using namespace sh_cube;
/*
 Normalize a vector
 */
 
void normalize3D(GRADIENT_COORD_TYPE intial[], GRADIENT_COORD_TYPE normalized[])
{
  double sum(0.0),mag(0.0);
	for (int i=0; i<3; i++) {
    sum = sum + intial[i]*intial[i];
  }
  if (sum>0) {
    mag = sqrt(sum);
    for (int j=0; j<3; j++) {
      normalized[j] = intial[j] / mag;
    }
  }
}


/*
 Check if the point is inside the cube.
 */
bool inCube(COORD_TYPE *pt){
  if (((pt[0]>= 0.0) && (pt[0]<= 1.0))
      &&((pt[1]>= 0.0) && (pt[1]<= 1.0))
      &&((pt[2]>= 0.0) && (pt[2]<= 1.0))) 
    {
    return true;
    }
  else {
    return false;  
  }
  
};

/*
 clamp to 0 /1 based on a threshold 't'.
 */
void clamp (COORD_TYPE *coord, const SCALAR_TYPE t)
{
  for (int i=0; i<3; i++) {
    if(coord[i]< 0.0)
      {
      if (coord[i] >= -t) {
        coord[i] = 0.0;
      }
      }
    if (coord[i] > 1.0) {
      if (coord[i] <= 1.0+t) {
        coord[i] = 1.0;
      }
    } 
  }
}

/*
 set sharp point to cube center
 */ 
void setCubeCenter(COORD_TYPE *shpoint)
{
  for (int i=0; i<3; i++) {
    shpoint[i]=0.5;
  }
};

/*
 set sharp point to cube centroid.
 */

void setCubeCentroid(const CUBE &cb, NUM_TYPE index[], COORD_TYPE *shpoint)
{
  for( int k=0; k<3; k++){
    double temp_k = 0.0;    
    for(int i=0; i<12; i++){
      if (cb.edges[i].is_intersect) {
        temp_k=temp_k+cb.edges[i].pt_intersect.pos[k];
      }
      }
      //Divide by the number of intersects to find the centroid.
    temp_k = (temp_k / double(cb.ne_intersect));
    shpoint[k] = temp_k;
    }  
}

  //calculate the sum for the entire row
inline void cal_sum(const RowVectorXf &res1, vector<float> &sum)
{
  float tempSum = res1(0) + res1(1) +res1(2);
  sum.push_back(tempSum); // sum squared
}

  //calcualte  magnitude for a single vector
float CalMag( RowVectorXf &res){
  float sum = 0;
  
  for (int i=0; i<3; i++){
    sum = sum +  res(i)*res(i);
  }
  return sqrt(sum);
  
};

/*
 Function to Calculate the w
 */
RowVectorXf Cal_w(MatrixXf &pseudoInverse,const MatrixXf  &m,const MatrixXf  &I)
{
  vector<RowVectorXf> e;
  RowVectorXf e1(3);
  e1 << 1,0,0;
  e.push_back(e1);
  
  e1 << 0,1,0;
  e.push_back(e1);
  
  e1 << 0,0,1;
  e.push_back(e1);
  
  RowVectorXf res;
  int index;
  float Mag = 0;
  float Max = -1;

  for (int i=0; i<e.size(); i++){
    res = (I - pseudoInverse*m)*e[i].transpose();
    Mag = CalMag(res);

    if (Mag > Max)
      {
      Max = Mag;
      index = i;
      }
  }
  return e[index];
}
/*
 Compute pseudo inverse of sigma
 */
void computePseudoSigmaInv(MatrixXf &PseudoInvSigma, MatrixXf &singular_vals,
                           int& num_svals ,const double & err )
{
  for (int i=0; i<singular_vals.rows(); i++) {
    if (singular_vals(i) > err) {
      num_svals++;
      PseudoInvSigma(i,i) = 1.0/singular_vals(i);
    }
  }
};


/*
 compute the pseudo inverse of m
 */
MatrixXf computePseudoinv
( const MatrixXf &m, const double err, int & num_svals,float eigenvalues[DIM3])
{
  JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);

  MatrixXf singular_vals;
  singular_vals = svd.singularValues();
  const int num_singular_values = singular_vals.rows();

  for (int i=0; i <DIM3;i++) 
    { eigenvalues[i] = 0; }

  for (int i=0; i < singular_vals.rows(); i++)
    { eigenvalues[i] = singular_vals(i); }

    // compute the pseudo inverse.
  MatrixXf PseudoInvSigma =
    MatrixXf::Zero(num_singular_values,num_singular_values);

  computePseudoSigmaInv(PseudoInvSigma, singular_vals, num_svals , err );

  MatrixXf pseudoInverse;
  pseudoInverse = svd.matrixV()*PseudoInvSigma.transpose()*svd.matrixU().transpose();
  return pseudoInverse;
};

/*
Find the sharp point in-cube.
*/
void sh_cube::findPoint
(const CUBE &cb,
const SCALAR_TYPE err,
float eigenvalues[DIM3],
int &num_large_eigenvalues,
SVD_INFO &svd_debug_info,
COORD_TYPE *shpoint)
{ 
 //set m 
  int ind[cb.ne_intersect]; //ind (index) keeep tracks of the edges which are intersected

  int i=0;
  MatrixXf m(cb.ne_intersect, cb.dim);
  for (int k=0; k<cb.num_edges; k++) {
    if (cb.edges[k].is_intersect) {
      for (int j=0; j<3; j++) {
        m(i,j) = cb.edges[k].pt_intersect.grads[j];
      }
      ind[i]=k;
      i++; 
    }
  }

  //set b
  RowVectorXf b(cb.ne_intersect);
  for (int i=0; i<cb.ne_intersect; i++) {
    b(i)=m(i,0)*cb.edges[ind[i]].pt_intersect.pos[0]+
    m(i,1)*cb.edges[ind[i]].pt_intersect.pos[1]+
    m(i,2)*cb.edges[ind[i]].pt_intersect.pos[2];
  }  
    //compute pseudo inverse of m. and also return the number of sing vals.
  int num_svals = 0;

  MatrixXf pseudoInv_m = computePseudoinv(m, err, num_svals, eigenvalues);
  num_large_eigenvalues = num_svals;

  // x = m inv *b, mx=b
  RowVectorXf x;
  x = pseudoInv_m*b.transpose();

  if (num_svals == 3 ) {
      //set shpoint to x
    shpoint[0] = x(0); 
    shpoint[1] = x(1);
    shpoint[2] = x(2);
      //clamp the shpoint
    clamp (shpoint, clamp_threshold);
    
    bool isInsideCube = inCube(shpoint);
    if (!isInsideCube) {
        //check centroid.
      setCubeCentroid(cb, ind, shpoint);
	// using centroid 
	svd_debug_info.location = CENTROID;
	
      bool isInsideCube2 = inCube(shpoint);
      if (!isInsideCube2) {
      svd_debug_info.location = CUBE_CENTER; //set using centre
        setCubeCenter(shpoint);
      }
    }
  }
  else if(num_svals == 2){
    RowVectorXf tempdir;
    MatrixXf I(3,3);

    I<< 1,0,0, 0,1,0, 0,0,1;
    RowVectorXf w = Cal_w(pseudoInv_m, m, I );
    tempdir = (I-pseudoInv_m*m)*w.transpose();
    GRADIENT_COORD_TYPE dir[DIM3];
    dir[0]=tempdir(0);
    dir[1]=tempdir(1);
    dir[2]=tempdir(2);
    normalize3D(dir, dir);
   
      //given x and ray direction , check if this intersects the identity cube.
    SCALAR_TYPE intersect[3]={0.0};
    SCALAR_TYPE p[3]={0.0};
    p[0]=x(0);
    p[1]=x(1);
    p[2]=x(2);
    bool isIntersect = calculate_point_intersect(p, dir, intersect);
    
    svd_debug_info.ray_direction[0] = dir[0];
    svd_debug_info.ray_direction[1] = dir[1];
    svd_debug_info.ray_direction[2] = dir[2];
    svd_debug_info.ray_initial_point[0] = p[0];
    svd_debug_info.ray_initial_point[1] = p[1];
    svd_debug_info.ray_initial_point[2] = p[2];
    
    if (isIntersect) {
    svd_debug_info.ray_intersect_cube = true;
    /* BIGGER CLAMP.  TEMPORARILY DISABLED.
    //intersect is on a bigger cube so we clamp it.
        clamp(intersect, clamp_threshold);
        //test if the intersect is within the cube.
      bool  isInCube = inCube(intersect);
      if (isInCube) {
        for (int i=0; i<3; i++) {
          shpoint[i] = intersect[i];
        }
      else {
      svd_debug_info.location = CENTROID ; 
        setCubeCentroid(cb, ind, shpoint);
      }
      */
       clamp(intersect, clamp_threshold);
        for (int i=0; i<3; i++) {
          shpoint[i] = intersect[i];
      }
    }
    else {
     svd_debug_info.location = CENTROID;
      setCubeCentroid(cb, ind, shpoint);
      /* BIGGER CLAMP.  TEMPORARILY DISABLED.
      bool isInsideCube2 = inCube(shpoint);
      if (!isInsideCube2) {
        setCubeCenter(shpoint);
      }
      */
    }
  }
  else {
     // One or less singular values 
    setCubeCentroid(cb, ind, shpoint);
     svd_debug_info.location = CENTROID;
    /* BIGGER CLAMP.  TEMPORARILY DISABLED.
    bool isInsideCube2 = inCube(shpoint);
    if (!isInsideCube2) {
      setCubeCenter(shpoint);
    }
    */
  }
}

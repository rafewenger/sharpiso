//FindIntersect.cxx
//

#include <iostream>
#include <cmath>
#include <limits>

#include "sharpiso_findIntersect.h"
using namespace std;

const int THRESHOLD_CLAMP = 0.0001; //
int findmax(const float dir[]);


/*
 FindIntersect  new , 
 accepts as inputs  a point p[], a direction dir[].It returns a boolean which is TRUE 
 if the Ray intersects the cube. It returns FALSE if the ray does not intersect the cube.
 If the bool is true , intersect[] returns the MID-point of intersection of the ray and the cueb.
 */
 bool calculate_point_intersect
(const SCALAR_TYPE *p, const SCALAR_TYPE *dir, SCALAR_TYPE *intersect)
{

  int ind  = findmax(dir);
    //find t0 and t1
  SCALAR_TYPE t0 = ( -1.0*p[ind])/(dir[ind]);
  SCALAR_TYPE t1 = (1.0 - p[ind])/(dir[ind]);
  SCALAR_TYPE min_coord_j;
  SCALAR_TYPE max_coord_j;
  
  for (int j=0; j<3; j++)
    {
    if(j != ind)
      {
      min_coord_j = p[j] + t0*dir[j];
      max_coord_j = p[j] + t1*dir[j];
      
        //swap
      if (min_coord_j > max_coord_j)
        {
        swap(t0, t1);
        swap(min_coord_j, max_coord_j);
        }
      if (min_coord_j - 1.0 < THRESHOLD_CLAMP ) {
        min_coord_j = 1.0;
      }
      if (max_coord_j-0.0 > -THRESHOLD_CLAMP ) {
        max_coord_j=0.0;
      }
     
      if(min_coord_j > 1.0)
        { 
          return false ;
        }
      else if(max_coord_j < 0.0)
        {
        return false;
        }
      else
        {
        
        if ((min_coord_j <= 0.0) && (max_coord_j > 0.0 ))
          {
          t0  = (-p[j])/dir[j];
          }
        if ((min_coord_j < 1.0) &&(max_coord_j >= 1.0 ))
          {
          t1  = (1.0 - p[j])/dir[j];
          }
        }
      }
    }
  
  SCALAR_TYPE endPt0[3];
  SCALAR_TYPE endPt1[3];
  for (int i=0; i<3; i++)
    {
    SCALAR_TYPE temp = p[i] + t0*dir[i];
    endPt0[i] = temp;
    }
    //t1
  for (int i=0; i<3; i++)
    {
    SCALAR_TYPE temp = p[i] + t1*dir[i];
    endPt1[i] = temp;
    }
  /* /// used for debug purpose 
  cout <<" pt 1 x:["<< endPt0[0] <<"] y:["<< endPt0[1] <<"] z:["<< endPt0[2] <<"]" 
  <<" t0: "<<t0<<endl;
  cout <<" pt 2 x:["<< endPt1[0] <<"] y:["<< endPt1[1] <<"] z:["<< endPt1[2] <<"]"
  <<" t1: "<<t1<<endl;
  */
    //Find the point of intersection.
  for (int i=0; i<3; i++) {
    intersect[i] = (endPt0[i] + endPt1[i])/2.0;
  }
  return true;
  
};
/*
The complex form of the cube , line intersection which  intersects a larger cube
it takes an extra parameter t which is the threshold.
*/
/*
bool calculate_point_intersect_cmplx(const double *p, const double *dir, 
double *intersect, const double t)
{
   //debug
    //cout<<"The point is : ["<<p[0]<<"] ["<<p[1]<<"] ["<<p[2]<<"]"<<endl;
    //cout<<"The direction is:["<<dir[0]<<"] ["<<dir[1]<<"] ["<<dir[2]<<"]"<<endl;
    //find index for max
    cout<<"complex function "<<endl;
  int ind  = findmax(dir);
    //find t0 and t1
  SCALAR_TYPE t0 = (-t -1.0*p[ind])/(dir[ind]);
  SCALAR_TYPE t1 = (1.0 + t - p[ind])/(dir[ind]);
  SCALAR_TYPE min_coord_j;
  SCALAR_TYPE max_coord_j;
  
  cout <<"index is : "<<ind<<endl;
  cout <<" t0 "<<t0<<endl;
  cout <<" t1 "<<t1<<endl;
  
  for (int j=0; j<3;j++)
    {
    if(j != ind)
      {
      cout <<" j :"<<j<<endl;
      cout <<" pj ("<<p[j]<<") t0 *dir ("<<t0*dir[j]<<endl;
      min_coord_j = p[j] + t0*dir[j];
      max_coord_j = p[j] + t1*dir[j];
      cout <<" pj ("<<p[j]<<") t1 *dir ("<<t1*dir[j]<<endl;
      cout <<"min_coord_j: "<<min_coord_j<<endl;
      cout <<"max_coord_j: "<<max_coord_j<<endl;
      
        //swap
      if (min_coord_j > max_coord_j)
        {
        swap(t0, t1);
        swap(min_coord_j, max_coord_j);
        }
        
     // if (min_coord_j-1.0 < THRESHOLD_CLAMP ) {
      //  min_coord_j=1.0;
     // }
     // if (max_coord_j-0.0 > -THRESHOLD_CLAMP ) {
     //   max_coord_j=0.0;
     // }
     
      cout<<"after swap"<<endl;
      cout <<"min_coord_j: "<<min_coord_j<<endl;
      cout <<"max_coord_j: "<<max_coord_j<<endl;
      if(min_coord_j > 1.0 + t)
        { 
          cout <<"no intersection  min coord is bigger than 1 ["<<(min_coord_j - 1.0)<<"]" <<endl;
          return false ;
        }
      else if(max_coord_j < -t )
        {
        cout <<"no_intersection max coord is less than 0 ["<<(max_coord_j - 0.0)<<"] "<<endl;
        return false;
        }
      else
        {
        
        if ((min_coord_j <= 0.0 - t) && (max_coord_j > 0.0 - t ))
          {
          cout <<"((min_coord_j < 0.0) && (max_coord_j > 0.0 )) "<<endl;
          t0  = (-t -p[j])/dir[j];
          }
        if ((min_coord_j < 1.0 + t) &&(max_coord_j >= 1.0 - t ))
          {
          cout <<"((min_coord_j < 1.0) &&(max_coord_j > 1.0 ))"<<endl;
          t1  = (1.0 + t - p[j])/dir[j];
          }
        }
      }
    }
  
  SCALAR_TYPE endPt0[3];
  SCALAR_TYPE endPt1[3];
  for (int i=0; i<3; i++)
    {
    SCALAR_TYPE temp = p[i] + t0*dir[i];
    endPt0[i] = temp;
    }
    //t1
  for (int i=0; i<3; i++)
    {
    SCALAR_TYPE temp = p[i] + t1*dir[i];
    endPt1[i] = temp;
    }
  
  cout <<" pt 1 x:["<< endPt0[0] <<"] y:["<< endPt0[1] <<"] z:["<< endPt0[2] <<"]" 
  <<" t0: "<<t0<<endl;
  cout <<" pt 2 x:["<< endPt1[0] <<"] y:["<< endPt1[1] <<"] z:["<< endPt1[2] <<"]"
  <<" t1: "<<t1<<endl;
  
    //Find the point of intersection.
  for (int i=0; i<3; i++) {
    intersect[i] = (endPt0[i] + endPt1[i])/2.0;
  }
  
  cout << " intersect "<<intersect[0]<<" "<<intersect[1]<<" "<<intersect[2]<<endl;
  return true;
  
};
*/

/////
int findmax(const SCALAR_TYPE dir[])
{
  int ind=0;
  SCALAR_TYPE max = abs(dir[ind]);
  for (int j=0;j<3;j++)
    {
    if (abs(dir[j]) >= max)
      {
      max = abs(dir[j]);
      ind = j;
      }
    }
  return ind;
}


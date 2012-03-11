//FindIntersect.cxx
//

#include <iostream>
#include <cmath>
#include <limits>

#include "sharpiso_findIntersect.h"
#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_linear_alg.txx"

using namespace std;
using namespace SHARPISO;


// HELPER FUNCTIONS:
// Findmax direction
int findmax(const float dir[]);



// FindIntersect  ,
// accepts as inputs  a point p[], a direction dir[].It returns a boolean which is TRUE
// if the Ray intersects the cube. It returns FALSE if the ray does not intersect the cube.
// If the bool is true , intersect[] returns the MID-point of intersection of the ray and the cueb.

bool SHARPISO::calculate_point_intersect
(const SCALAR_TYPE *p,
 const SCALAR_TYPE *dir,
 SCALAR_TYPE *intersect)
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
    //Find the point of intersection.
    for (int i=0; i<3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i])/2.0;
    }
    return true;
    
};


// Caculate point intesect between the ray and the cube.
// This version takes in the coordinate of the cube index and translates,
// to find the intersection with the unit cube and translate it back.

bool SHARPISO::calculate_point_intersect
(const COORD_TYPE cube_coord[],
 const SCALAR_TYPE *original_pt,
 const SCALAR_TYPE *dir,
 SCALAR_TYPE *intersect)
{
	SCALAR_TYPE p[DIM3] ={0.0};
	for (int i=0; i < DIM3; i++){
		p[i] = original_pt[i] - cube_coord[i];
	}
    
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
    
    SCALAR_TYPE endPt0[DIM3];
    SCALAR_TYPE endPt1[DIM3];
    for (int i=0; i<DIM3; i++)
    {
        SCALAR_TYPE temp = p[i] + t0*dir[i];
        endPt0[i] = temp;
    }
    //t1
    for (int i=0; i<DIM3; i++)
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
    for (int i=0; i<DIM3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i])/2.0;
        intersect[i] += cube_coord[i];
    }
    return true;
    
};




// Calculate the intersection with a larger cube.
// Takes an extra parameter of how big the cube is 

bool SHARPISO::calculate_point_intersect_complex
(const COORD_TYPE cube_coord[],
 const SCALAR_TYPE *original_pt,
 const SCALAR_TYPE *dir,
 const float th,
 SCALAR_TYPE *intersect)
{ 
	SCALAR_TYPE p[DIM3] ={0.0};
    
	for (int i = 0; i < DIM3; i++){
		p[i] = original_pt[i] - cube_coord[i];
	}
    
    int ind  = findmax(dir);    
    //find t0 and t1
    SCALAR_TYPE t0 = (-th -p[ind])/(dir[ind]);
    SCALAR_TYPE t1 = (1.0 + th - p[ind])/(dir[ind]);
    SCALAR_TYPE min_coord_j;
    SCALAR_TYPE max_coord_j;
    
    for (int j = 0; j < DIM3; j++)
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
            
            if (( min_coord_j - (1.0 + th)) > EPSILON)  
            {
                return false ;
            }
            else if (((0.0 - th) - max_coord_j) > EPSILON)
            {
                return false;
            }
            else
            {
                if ((min_coord_j <= (0.0 - th))
                    && (max_coord_j > (0.0 - th)))
                {  t0 = (-th - p[j])/dir[j]; }
                
                
                if ((min_coord_j < (1.0 + th)) &&
                    (max_coord_j >= (1.0 + th)))
                { t1  = (1.0 + th - p[j])/dir[j]; }
            }
        }
    }
    
    SCALAR_TYPE endPt0[3];
    SCALAR_TYPE endPt1[3];
    for (int i=0; i<3; i++)
    {
        endPt0[i] = p[i] + t0*dir[i];
        endPt1[i] = p[i] + t1*dir[i];
    }
    /*
     /// used for debug purpose
     cout <<" pt 1 ["<< endPt0[0] <<" "<< endPt0[1] <<" "<< endPt0[2] <<"]"
     <<" t0: "<<t0<<endl;
     cout <<" pt 2 ["<< endPt1[0] <<" "<< endPt1[1] <<"  "<< endPt1[2] <<"]"
     <<" t1: "<<t1<<endl;
     */
    
    //Find the point of intersection.
    for (int i = 0; i < DIM3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i] ) / 2.0;
        intersect[i] = intersect[i] + cube_coord[i];
    }
    return true;    
};


//  Given a point and a ray direction and the cord of the cube index 
//  Find a point on the ray which is at the closest distance to the cube-center
//  Note calculate the cube center from the cube coord (index) by adding 0.5 0.5 0.5
void SHARPISO::compute_closest_point_to_cube_center
(const COORD_TYPE cube_coord[],
 const COORD_TYPE coord[],
 const COORD_TYPE ray_direction[],
 COORD_TYPE closest_point[DIM3])
{
  COORD_TYPE ray_direction_normalized[DIM3]={0.0};
  bool is_ray_direc_zero = false;
  // Normalize the ray_direction
  normalize_vector_3D
  (ray_direction, ray_direction_normalized, is_ray_direc_zero, 0.0001);
 
  // find the cube center
  COORD_TYPE center_offset[DIM3] = {0.5, 0.5, 0.5};
  COORD_TYPE cube_center[DIM3];
  IJK::add_coord(DIM3, center_offset, cube_coord, cube_center);
  
  
  // the vector perpendicular is given by 
  // c =a -(a dot b^)* b^
  COORD_TYPE a[DIM3]={0.0};  
    //set a and b is given  by the ray direction
  IJK::subtract_coord(DIM3, cube_center, coord, a);

  COORD_TYPE c[DIM3]={0.0};
  // dot product
  float dotpdt=0.0;
  IJK::compute_inner_product(DIM3, a, ray_direction_normalized, dotpdt);
  
  // compute c and invert it at same time 
  for (int d=0; d<DIM3; d++) 
    c[d]=-(a[d]-dotpdt*ray_direction_normalized[d]);
  
  
  //debug
 /*
  cout <<" cube_coord  "<<cube_coord[0]<<","<<cube_coord[1]<<","<<cube_coord[2]<<endl;
  cout <<" cube_center "<<cube_center[0]<<","<<cube_center[1]<<","<<cube_center[2]<<endl;
  cout <<" coord       "<<coord[0]<<","<<coord[1]<<","<<coord[2]<<endl;
  cout <<" a           "<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
  cout <<" dot pdt     "<<dotpdt<<endl;
  cout <<" c           "<<c[0]<<","<<c[1]<<","<<c[2]<<endl;
  cout <<" ray dir     "<<ray_direction_normalized[0]<<","<<ray_direction_normalized[1]<<","
  <<ray_direction_normalized[2]<<endl;
  */
  IJK::add_coord(DIM3, cube_center, c, closest_point);
}



// HELPER FUNCTIONS:
// find the max direction
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


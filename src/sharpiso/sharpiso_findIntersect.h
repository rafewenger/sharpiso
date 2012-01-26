//FindIntersect.h
//Find Intersect of the ray and the cube given a point and a direction.


#ifndef _SHARPISO_FINDINTERSECT_
#define _SHARPISO_FINDINTERSECT_

#include<iostream>
#include <Eigen/Dense>
#include "sharpiso_types.h"
#include "sharpiso_eigen.h"


using namespace Eigen;
//using namespace SHARPISO;

namespace SHARPISO {
    
    const int THRESHOLD_CLAMP = 0.001; // clamp threshold 
    
    // FindIntersect , 
    // accepts as inputs  a point p[], a direction dir[].It returns a boolean which is TRUE 
    // if the Ray intersects the cube. It returns FALSE if the ray does not intersect the cube.
    // If the bool is true , intersect returns the MID-point of intersection of the ray and the cueb.
    
    
    bool calculate_point_intersect
    (const SCALAR_TYPE * point,
     const SCALAR_TYPE *dir,
     SCALAR_TYPE *intersect);
    
    // Separate version which translate back and forth to find the intersection with the unit cube.
    
    bool calculate_point_intersect
    (const COORD_TYPE cube_coord[],
     const SCALAR_TYPE *p,
     const SCALAR_TYPE *dir,
     SCALAR_TYPE *intersect);
    
    // Calculate the intersection with a larger cube.
    
    bool calculate_point_intersect_complex
    (const COORD_TYPE cube_coord[],
     const SCALAR_TYPE *original_pt,
     const SCALAR_TYPE *dir,
     const float th,  
     SCALAR_TYPE *intersect);
    
};
#endif

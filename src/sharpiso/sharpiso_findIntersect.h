// FindIntersect.h
// Find intersection of a ray and a cube.


#ifndef _SHARPISO_FINDINTERSECT_
#define _SHARPISO_FINDINTERSECT_

#include<iostream>
#include <Eigen/Dense>
#include "sharpiso_types.h"
#include "sharpiso_eigen.h"

namespace SHARPISO {
    
    const float THRESHOLD_CLAMP = 0.001; // clamp threshold 
    const float EPSILON = 0.0001;
    
    /// Compute the midpoint of the intersection of the ray and the cube.
    /// Return true if the ray intersects the cube
    bool calculate_point_intersect
    (const SCALAR_TYPE * point, const SCALAR_TYPE *dir, SCALAR_TYPE *intersect);
    
    /// Separate version which translate back and forth 
    ///   to find the intersection with the unit cube.
    bool calculate_point_intersect
    (const COORD_TYPE cube_coord[], const SCALAR_TYPE *p,
     const SCALAR_TYPE *dir, SCALAR_TYPE *intersect);
    
    // Calculate the intersection with an enlarged cube.
    bool calculate_point_intersect_complex
    (const COORD_TYPE cube_coord[], const SCALAR_TYPE *original_pt,
     const SCALAR_TYPE *dir, const float th,  SCALAR_TYPE *intersect);
    
};
#endif

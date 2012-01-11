//FindIntersect.h
//Find Intersect of the ray and the cube given a point and a direction.
/*
# include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
#include "sharpiso_types.h"
*/
#include<iostream>
#include <Eigen/Dense>
#include "sharpiso_types.h"
#include "sharpiso_eigen.h"

using namespace std;
using namespace Eigen;
using namespace SHARPISO;
/*
FindIntersect , 
accepts as inputs  a point p[], a direction dir[].It returns a boolean which is TRUE 
if the Ray intersects the cube. It returns FALSE if the ray does not intersect the cube.
If the bool is true , intersect returns the MID-point of intersection of the ray and the cueb.
*/

bool calculate_point_intersect
(const SCALAR_TYPE * point, const SCALAR_TYPE *dir, SCALAR_TYPE *intersect);


//Separate version which translate back and forth to find the intersection with the unit cube.

bool calculate_point_intersect
(const COORD_TYPE cube_coord[], const SCALAR_TYPE *p,
const SCALAR_TYPE *dir, SCALAR_TYPE *intersect);

//debug
/*
//Function calculates the intersection between a larger cube and the ray, it takes  a extra parameter
//which decides how big the larger cube is.

bool calculate_point_intersect_cmplx(const double *p, const double *dir, 
double *intersect, const double t)
*/
/*bool calculate_point_intersect_cmplx(const double *p, const double *dir, 
double *intersect, const double t);*/
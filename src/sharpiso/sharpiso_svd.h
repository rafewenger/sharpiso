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
 EIGENVALUE_TYPE * singular_vals,
 COORD_TYPE * isoVertcoords);
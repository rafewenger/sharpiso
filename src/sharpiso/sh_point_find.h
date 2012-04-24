/*
 *  sh_point_find.h
 *  SHARPPOINT
 *
 *  Created by arindam bhattacharya on 11/6/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */
//new code adding to the sharpinfo 


#include <iostream>
#include "sh_point_datastruct.h"
#include "sharpiso_types.h"
#include "sharpiso_feature.h"

using namespace sh_cube;
using namespace SHARPISO;
/*
 For each cube this function returns a double[3] 'point'
 which the dual-vertex for this cube
 input :
 gradients, 
 scalar_values,
 isovalue, 
 EIGEN_VALUE_CUTOFF for svd calculation 
 point(which is also its out put)
 */

bool shFindPoint
(const GRADIENT_COORD_TYPE gradients[],
 const  SCALAR_TYPE  scalar_vals[],
 const  SCALAR_TYPE  isovalue, const  bool use_cmplx_interp, 
 const  SCALAR_TYPE EIGEN_VALUE_CUTOFF, float eigenvalues[DIM3],
 int &num_large_eigenvalues, 
 SVD_INFO &svd_debug_info, 
 COORD_TYPE *shpoint);




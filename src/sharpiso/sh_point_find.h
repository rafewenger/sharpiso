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

/// Get gradients at edge intersection points determined by edge endpoints.
/// Use sharp formula for computing gradient at intersection.
void get_edgeI_sharp_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & point_coord,
 std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 std::vector<SCALAR_TYPE> & scalar,
 NUM_TYPE & num_gradients);


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

void setup_edgeIntercepts
(CUBE &cb,
 const double isovalue,
 const bool use_cmplx_interp);




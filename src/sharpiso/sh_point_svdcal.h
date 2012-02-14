/*
 *  sh_point_svdcal.h
 *  SHARPPOINT
 *
 *  Created by arindam bhattacharya on 11/6/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */

#ifndef _SHARPISO_SVD_CAL_
#define _SHARPISO_SVD_CAL_

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "sh_point_datastruct.h"
#include "sharpiso_types.h"
#include "sharpiso_feature.h"

using namespace Eigen;
using namespace SHARPISO;


// Accepts a cube as an input and calculates the singular values , the num of large singular values
// and the sharp point in the cube.

namespace sh_cube
{
    ///Find the sharp point in-cube.
    void findPoint
    (const CUBE &cb,
     const SCALAR_TYPE EIGEN_VALUE_CUTOFF,
     float eigenvalues[DIM3],
     int &num_large_eigenvalues,
     SVD_INFO &svd_debug_info,
     COORD_TYPE *shpoint);
}

#endif
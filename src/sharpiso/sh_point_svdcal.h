/*
 *  sh_point_svdcal.h
 *  SHARPPOINT
 *
 *  Created by arindam bhattacharya on 11/6/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "sh_point_datastruct.h"
#include "sharpiso_types.h"

using namespace std;
using namespace sh_cube;
using namespace Eigen;
using namespace SHARPISO;


// Accepts a cube as an input and calculates the singular values , the num of large singular values
// and the sharp point in the cube.

void findPoint(const CUBE &cb, const SCALAR_TYPE err, float eigenvalues[DIM3],
int &num_large_eigenvalues, COORD_TYPE *shpoint);
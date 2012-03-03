#ifndef _ANISOGRADINFO_
#define _ANISOGRADINFO_



#include "ijkscalar_grid.txx"
#include "sharpiso_grids.h"

using namespace SHARPISO;

class ANISOINFO_TYPE{
    
public:
    
    int iv1;
    int flag_aniso;
    bool flag_k;
    bool flag_normals;
    int iter;
    bool flag_m;
    GRADIENT_COORD_TYPE *normals;
    
    GRADIENT_COORD_TYPE mX[DIM3];
    GRADIENT_COORD_TYPE mX_prev_vert_X[DIM3];
    
    GRADIENT_COORD_TYPE mY[DIM3];
    GRADIENT_COORD_TYPE mY_prev_vert_Y[DIM3];
    
    GRADIENT_COORD_TYPE mZ[DIM3];
    GRADIENT_COORD_TYPE mZ_prev_vert_Z[DIM3];
    
    
    bool print_c;
    GRADIENT_COORD_TYPE c[DIM3];
    GRADIENT_COORD_TYPE gradientS_d[DIM3];
    int dirc;
    GRADIENT_COORD_TYPE   fwd_diff_d_normals[DIM3];
    GRADIENT_COORD_TYPE fwd_diff_d;
    bool print_gradientH_d_normals;
    int gradH_d_normals_direc;
    GRADIENT_COORD_TYPE  gradientH_d_Normals[9];
    
    SCALAR_TYPE K[DIM3];
    SCALAR_TYPE gK[DIM3];
};


void compute_curvature 
 (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
  const GRADIENT_GRID & gradient_grid,
  const float mu,
  const float lambda,
  const VERTEX_INDEX iv,
  ANISOINFO_TYPE &aniso_info);

void print_info
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const ANISOINFO_TYPE aniso_info);

#endif
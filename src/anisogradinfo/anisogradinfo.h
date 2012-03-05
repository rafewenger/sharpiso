#ifndef _ANISOGRADINFO_
#define _ANISOGRADINFO_



#include "ijkscalar_grid.txx"
#include "sharpiso_grids.h"

using namespace SHARPISO;

// **************************************************
// ANISOINFO_TYPE
// **************************************************

class ANISOINFO_TYPE {

 protected:
  void Init();

 public:
    
  int iv1;
  int flag_aniso;
  bool flag_k;
  bool flag_m;
  bool flag_print_scalar;
  bool flag_print_normal;
  bool flag_print_neighbors;
  bool flag_print_c;
  bool flag_print_w;
  bool flag_print_gradS;
  bool flag_print_gradN;
  bool flag_print_cdiffN;
  bool flag_print_fdiffN;
  int half_edge_direction;

  int iter;
  GRADIENT_COORD_TYPE *normals;

  GRADIENT_COORD_TYPE mX[DIM3];
  GRADIENT_COORD_TYPE mX_prev_vert_X[DIM3];
    
  GRADIENT_COORD_TYPE mY[DIM3];
  GRADIENT_COORD_TYPE mY_prev_vert_Y[DIM3];
    
  GRADIENT_COORD_TYPE mZ[DIM3];
  GRADIENT_COORD_TYPE mZ_prev_vert_Z[DIM3];
    
  GRADIENT_COORD_TYPE c[DIM3];
  GRADIENT_COORD_TYPE w[DIM3];
  GRADIENT_COORD_TYPE w_orth[DIM3];

  GRADIENT_COORD_TYPE fwd_diff_d_normals[DIM3];
  GRADIENT_COORD_TYPE fwd_diff_d;

  GRADIENT_COORD_TYPE  gradientH_d_Normals[DIM3*DIM3];
  GRADIENT_COORD_TYPE  gradientS_d[DIM3];

  GRADIENT_COORD_TYPE cdiffN[DIM3*DIM3];

  SCALAR_TYPE K[DIM3];
  SCALAR_TYPE gK[DIM3];

 public:
  ANISOINFO_TYPE() { Init(); }
};


void compute_curvature 
 (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
  const GRADIENT_GRID & gradient_grid,
  const float mu,
  const float lambda,
  const VERTEX_INDEX iv,
  ANISOINFO_TYPE &aniso_info);

void print_aniso_info
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const ANISOINFO_TYPE aniso_info);

void print_gradient_info
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv);

#endif


#include <iostream>
#include <iomanip>

#include "anisogradinfo.h"
#include "anisograd.h"
#include "anisograd_operators.h"

using namespace std;

// **************************************************
// COMPUTE
// **************************************************

void compute_curvature 
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const float mu,
 const float lambda,
 const VERTEX_INDEX iv1,
 ANISOINFO_TYPE  &aniso_info)
{
  const int hdir = aniso_info.half_edge_direction;
  int icube = iv1;
  aniso_info.iv1 = iv1;

  // Compute the central difference.

  // Compute M d for direction 'd' for  vertex iv1
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 0, aniso_info.mX);
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 1, aniso_info.mY);
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 2, aniso_info.mZ);

  // compute the curvature k for a vertex iv1
  compute_curvature_iv 
    (scalar_grid, gradient_grid, iv1, aniso_info.K);
    
  for (int d=0; d<DIM3; d++)
    compute_g_x(mu, aniso_info.K[d],aniso_info.flag_aniso, aniso_info.gK[d]);

  compute_c_d
    ( scalar_grid, gradient_grid, iv1, hdir, aniso_info.c);
    
  compute_forward_difference_d_normals
    (gradient_grid, iv1, hdir, aniso_info.fwd_diff_d_normals);
    
  compute_forward_difference_d
    (scalar_grid, iv1, hdir, aniso_info.fwd_diff_d);

  compute_gradH_d_scalar_grid
    (scalar_grid, iv1, hdir, aniso_info.gradientS_d);
    
  // Compute operator gradH for the direction dir
  // for the Normal field
  compute_gradH_d_normals 
    (gradient_grid,iv1, hdir, aniso_info.gradientH_d_Normals);

  // Compute gradient of normals using central difference
  compute_gradient_normals(gradient_grid, iv1, aniso_info.cdiffN);
}


void compute_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv,
 GRADIENT_COORD_TYPE cdiff[DIM3],
 GRADIENT_COORD_TYPE fdiff[DIM3])
{
  for (int d = 0; d < DIM3; d++) {
    compute_central_difference_d(scalar_grid, iv, d, cdiff[d]);
    compute_forward_difference_d(scalar_grid, iv, d, fdiff[d]);
  }
}

// **************************************************
// PRINT
// **************************************************

template <typename CTYPE>
void aniso_print_coord(std::ostream & out, const CTYPE coord[DIM3])
{
  const int OWIDTH = 10;

  out << "(";
  for (int i = 0; i < DIM3; i++) {
    out << setw(10);
    out << coord[i];
    if (i+1 < DIM3) { out << ","; };
  }
  out << ")";
}

void print_aniso_info
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const ANISOINFO_TYPE aniso_info)
{
  const int hdir = aniso_info.half_edge_direction;

  if (aniso_info.flag_k){
    cout <<"  K: ";
    aniso_print_coord(cout, aniso_info.K);

    cout <<" Gk: ";
    aniso_print_coord(cout, aniso_info.gK);
    cout << endl;
  }

  if (aniso_info.flag_m) {
    cout << "  m.x: ";
    aniso_print_coord(cout, aniso_info.mX);
    cout << endl;
    cout << "  m.y: ";
    aniso_print_coord(cout, aniso_info.mY);
    cout << endl;
    cout << "  m.z: ";
    aniso_print_coord(cout, aniso_info.mZ);
    cout << endl;
  }
    
  if (aniso_info.flag_print_c) {
    cout <<"  c (direction " << hdir << "):  ";
    aniso_print_coord(cout, aniso_info.c);
    cout << endl;
  }

  if (aniso_info.flag_print_gradS) {

    cout << "  gradientS normals (direction " 
         << hdir << "): ";
    aniso_print_coord(cout, aniso_info.gradientS_d);
    cout << endl;
  }

  if (aniso_info.flag_print_cdiffN) {
    cout << "  gradient normals (central difference):" << endl;
    cout << "    index 0: ";
    aniso_print_coord(cout, aniso_info.cdiffN);
    cout << endl;
    cout << "    index 1: ";
    aniso_print_coord(cout, aniso_info.cdiffN+DIM3);
    cout << endl;
    cout << "    index 2: ";
    aniso_print_coord(cout, aniso_info.cdiffN+2*DIM3);
    cout << endl;
  }

  if (aniso_info.flag_print_fdiffN) {
    cout << "  gradient normals (forward difference, direction "
         << hdir << " ): ";
    aniso_print_coord(cout, aniso_info.fwd_diff_d_normals);
    cout << endl;
  }

  if (aniso_info.flag_print_gradN) {

    cout << "  gradientH_d normals (direction " 
         << hdir << "): " << endl;
    cout << "    index 0: ";
    aniso_print_coord(cout, aniso_info.gradientH_d_Normals);
    cout << endl;
    cout << "    index 1: ";
    aniso_print_coord(cout, aniso_info.gradientH_d_Normals+DIM3);
    cout << endl;
    cout << "    index 2: ";
    aniso_print_coord(cout, aniso_info.gradientH_d_Normals+2*DIM3);
    cout << endl;
  }

  if (aniso_info.flag_normals){
    cout << "  Vertex normal: ";
    aniso_print_coord(cout, aniso_info.normals);
    cout << endl;
        
    for (int i=0; i<DIM3; i++) {
      VERTEX_INDEX prev = scalar_grid.PrevVertex(aniso_info.iv1, i);
      VERTEX_INDEX next = scalar_grid.NextVertex(aniso_info.iv1, i);
      const GRADIENT_COORD_TYPE *gr_prev;
      const GRADIENT_COORD_TYPE *gr_next;
      gr_prev=gradient_grid.VectorPtrConst(prev);
      cout << "    Direction " << i << ". ";
      cout <<"  Prev vert: ";
      aniso_print_coord(cout, gr_prev);
      gr_next=gradient_grid.VectorPtrConst(next);
      cout <<"  Next vert: ";
      aniso_print_coord(cout, gr_next);
      cout << endl;
    }
    
  }
    
}

void print_gradient_info
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX iv)
{
  GRADIENT_COORD_TYPE cdiff[DIM3];
  GRADIENT_COORD_TYPE normalized_cdiff[DIM3];
  GRADIENT_COORD_TYPE fdiff[DIM3];
  GRADIENT_COORD_TYPE normalized_fdiff[DIM3];

  compute_gradients(scalar_grid, iv, cdiff, fdiff);
  std::copy(cdiff, cdiff+DIM3, normalized_cdiff);
  std::copy(fdiff, fdiff+DIM3, normalized_fdiff);
  normalize(normalized_cdiff, DIM3);
  normalize(normalized_fdiff, DIM3);

  cout << "Gradients computed from scalar grid:" << endl;
  cout << "  central difference: ";
  aniso_print_coord(cout, cdiff);
  cout << " normalized: ";
  aniso_print_coord(cout, normalized_cdiff);
  cout << endl;
  cout << "  forward difference: ";
  aniso_print_coord(cout, fdiff);
  cout << " normalized: ";
  aniso_print_coord(cout, normalized_fdiff);
  cout << endl;
}


// **************************************************
// ANISOINFO_TYPE MEMBER FUNCTIONS
// **************************************************

void ANISOINFO_TYPE::Init()
{
  flag_k = false;
  flag_normals = false;
  flag_m = false;
  flag_print_c = false;
  flag_print_gradS = false;
  flag_print_gradN = false;
  flag_print_cdiffN = false;
  flag_print_fdiffN = false;
  half_edge_direction = 0;
}

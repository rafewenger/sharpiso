
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
  int icube = iv1;
  aniso_info.iv1 = iv1;

  // Compute M d for direction 'd' for  vertex iv1
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 0, aniso_info.mX);
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 1, aniso_info.mY);
  compute_m_d( scalar_grid, gradient_grid, icube, iv1, 2, aniso_info.mZ);
    
  // compute the curvature k for a vertex iv1
  compute_curvature_iv 
    (scalar_grid, gradient_grid, iv1, aniso_info.K);
    
  for (int d=0; d<DIM3; d++)
    compute_g_x(mu, aniso_info.K[d],aniso_info.flag_aniso, aniso_info.gK[d]);
    
  compute_gradH_d_scalar_grid
    (scalar_grid, iv1, aniso_info.dirc, aniso_info.gradientS_d);
    
  compute_c_d
    ( scalar_grid, gradient_grid, iv1, aniso_info.dirc, aniso_info.c);
    
  compute_forward_difference_d_normals
    ( gradient_grid, iv1, aniso_info.dirc, aniso_info.fwd_diff_d_normals);
    
  compute_forward_difference_d
    (scalar_grid, iv1, aniso_info.dirc, aniso_info.fwd_diff_d);
    
  // Compute operator gradH for the direction 'd'
  // for the Normal field
  compute_gradH_d_normals 
    (gradient_grid,iv1, aniso_info.gradH_d_normals_direc, 
     aniso_info.gradientH_d_Normals);
    
};


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

void print_info
(
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const ANISOINFO_TYPE aniso_info)
{
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
    
  if (aniso_info.print_c) {
    cout <<"  c (direction" << aniso_info.dirc << "):  ";
    aniso_print_coord(cout, aniso_info.c);
    cout << endl;
  }

  if (aniso_info.print_gradientH_d_normals) {

    cout << "  gradientH_d normals (direction " 
         << aniso_info.gradH_d_normals_direc << "): " << endl;
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
    
  if(aniso_info.flag_normals){
    aniso_print_coord(cout, aniso_info.normals);
    cout << endl;
        
    for (int i=0; i<DIM3; i++) {
      VERTEX_INDEX prev = scalar_grid.PrevVertex(aniso_info.iv1, i);
      VERTEX_INDEX next = scalar_grid.NextVertex(aniso_info.iv1, i);
      const GRADIENT_COORD_TYPE *gr_prev;
      const GRADIENT_COORD_TYPE *gr_next;
      gr_prev=gradient_grid.VectorPtrConst(prev);
      cout <<"prev vert in direction "<<i<<" ";
      cout <<"("<<gr_prev[0]<<","<<gr_prev[1]<<","<<gr_prev[2]<<")"<<endl;
      gr_next=gradient_grid.VectorPtrConst(next);
      cout <<"next vert in direction "<<i<<" ";
      cout <<"("<<gr_next[0]<<","<<gr_next[1]<<","<<gr_next[2]<<")"<<endl;
    }
    
  }
    
}

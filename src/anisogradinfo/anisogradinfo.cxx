#include "anisogradinfo.h"

#include "anisograd.h"
#include "anisograd_operators.h"
#include <iomanip>


using namespace std;
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
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 0, aniso_info.mX);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 1, aniso_info.mY);
    
    compute_m_d
    ( scalar_grid, gradient_grid, icube, iv1, 2, aniso_info.mZ);
    
    // compute the curvature k for a vertex iv1
    compute_curvature_iv 
     (scalar_grid, gradient_grid, iv1, aniso_info.K);
    
    for (int d=0; d<DIM3; d++)
        compute_g_x(mu, aniso_info.K[d],aniso_info.flag_aniso, aniso_info.gK[d]);
    
};



void debug_print(const ANISOINFO_TYPE aniso_info)
{
    cout <<"K vertex( "<< aniso_info.iv1<<") iter "<<aniso_info.iter<<" is ";
    
    cout <<" ("<<setw(10)<<aniso_info.K[0]<<","<<setw(10)<<aniso_info.K[1]<<","<<setw(10)<<aniso_info.K[2]<<")";
    
    cout <<" Gk ";
    
    cout <<" ("<<setw(10)<<aniso_info.gK[0]<<","<<setw(10)<<aniso_info.gK[1]<<","<<setw(10)<<aniso_info.gK[2]<<")"<<endl;
    if (aniso_info.flag_m) {
        cout <<" m "<<endl;
        cout <<"x("<<aniso_info.mX[0]<<" "<<aniso_info.mX[1]<<" "<<aniso_info.mX[2]<<")"<<endl;
        cout <<"y("<<aniso_info.mY[0]<<" "<<aniso_info.mY[1]<<" "<<aniso_info.mY[2]<<")"<<endl;
        cout <<"z("<<aniso_info.mZ[0]<<" "<<aniso_info.mZ[1]<<" "<<aniso_info.mZ[2]<<")"<<endl;
        
    }
}
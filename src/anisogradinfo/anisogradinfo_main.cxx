/// \file anisogradinfo_main.cxx
/// compute gradient vectors from scalar data
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
 
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <iostream>
#include <vector>



#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"

#include "anisograd.h"
#include "sharpiso_types.h"
#include "anisogradinfo.h"

using namespace IJK;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
bool report_time_flag = false;
bool flag_gzip = false;
bool flag_cdiff = false;
bool flag_iso = false;
float large_magnitude = 2.0;

bool flag_norm_after = false;
bool flag_norm_before = false;
bool flag_m = false;
bool flag_mprev = false;
bool flag_c = false;
bool flag_w = false;
bool flag_print_curvature = false;
bool flag_icube = false;
float mu(0.1);
float lambda(1.0);
int num_iter = 10;
VERTEX_INDEX icube = 0;

bool debug = true;

using namespace std;

// local subroutines
void memory_exhaustion();
void usage_error();
void parse_command_line(int argc, char **argv);

ANISOINFO_TYPE aniso_info;

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);
    
  vector<GRADIENT_COORD_TYPE> mag_list;
    
  IJK::ERROR error;
    
  try {
        
    std::set_new_handler(memory_exhaustion);
        
    parse_command_line(argc, argv);
        
    SHARPISO_SCALAR_GRID full_scalar_grid;
    GRID_NRRD_IN<int,int> nrrd_in;
    NRRD_DATA<int,int> nrrd_header;
     
        
    nrrd_in.ReadScalarGrid
      (scalar_filename, full_scalar_grid, nrrd_header, error);
        
    if (nrrd_in.ReadFailed()) { throw error; }
        
    GRADIENT_GRID gradient_grid;
        
    // compute the central gradients first 
    compute_gradient_central_difference(full_scalar_grid, icube, gradient_grid);
            
    normalize_and_store_gradient_magnitudes(full_scalar_grid, gradient_grid, mag_list);
            
    if (flag_iso) {

      // Calculate the anisotropic diff of the gradients.
      const int dimension = full_scalar_grid.Dimension();
      gradient_grid.SetSize(full_scalar_grid, dimension);
      for (int k=0; k<num_iter; k++) {
        aniso_info.iter=k;
        if(flag_icube) {
          aniso_info.flag_aniso = 0;
          compute_curvature
            (full_scalar_grid, gradient_grid, mu, lambda, icube, aniso_info);
          print_info(full_scalar_grid, gradient_grid, aniso_info);
        }
                    
        anisotropic_diff_iter_k
          (full_scalar_grid, mu, lambda, k, 0, icube, dimension, gradient_grid);
      }            
    }
    else {
      cout << "Anisostropic gradients called."<<endl;
      const int dimension = full_scalar_grid.Dimension();
      gradient_grid.SetSize(full_scalar_grid, dimension);
      for (int k=0; k<num_iter; k++) {
        aniso_info.iter=k;
        cout <<"iteration "<<k<<endl;
        if(flag_icube) {
          aniso_info.flag_aniso = 1;
          aniso_info.normals = gradient_grid.VectorPtr(icube);
          compute_curvature
            (full_scalar_grid, gradient_grid, mu, lambda, icube, aniso_info);
          print_info(full_scalar_grid, gradient_grid,aniso_info);
        }
        anisotropic_diff_iter_k
          (full_scalar_grid, mu, lambda, k, 1, icube, dimension, gradient_grid);
      }
    }
        
    if (report_time_flag) {
            
      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);
            
      cout << endl;
      cout << "Elapsed time = " << total_elapsed_time << " seconds." 
           << endl;
    };
  }
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };
    
}


// **************************************************
// MISC ROUTINES
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
    
  while (iarg < argc && argv[iarg][0] == '-') {
    if (string(argv[iarg]) == "-time") 
      { report_time_flag = true;   }
    else if (string(argv[iarg]) == "-gzip")
      { flag_gzip = true; }
    else if (string(argv[iarg]) == "-cdiff")
      { flag_cdiff = true; }
    else if (string(argv[iarg]) == "-iso")
      { flag_iso = true; }
    else if (string(argv[iarg]) == "-k")
      { aniso_info.flag_k = true; }
    else if (string(argv[iarg]) == "-m")
      { aniso_info.flag_m = true; }
    else if (string(argv[iarg]) == "-n")
      { aniso_info.flag_normals = true; }
    else if (string(argv[iarg]) == "-c")
      {
        aniso_info.print_c = true;
        iarg++;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%d", &aniso_info.dirc);
      }
    else if (string(argv[iarg]) == "-gradH_d_Normals")
      {
        aniso_info.print_gradientH_d_normals= true;
        iarg++;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%d", &aniso_info.gradH_d_normals_direc);
      }
    else if (string(argv[iarg]) == "-mu")
      {
        iarg++;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%f", &mu);
      }
    else if (string(argv[iarg]) == "-lambda")
      {
        iarg++;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%f", &lambda);
      }
    else if (string(argv[iarg]) == "-num_iter")
      {
        iarg++;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%d", &num_iter);
      }
    else if (string(argv[iarg]) == "-icube")
      {
        iarg++;
        flag_icube = true;
        if (iarg >= argc) { usage_error(); };
        sscanf(argv[iarg], "%d", &icube);
      }
    else 
      { usage_error(); }
    iarg++;
  }
    
  if (iarg+1 != argc) { usage_error(); };
    
  scalar_filename = argv[iarg];
    
}

void usage_msg()
{
  cerr <<"Usage: anisogradinfo [options]  {scalar nrrd file}"<<endl;
  cerr <<"                 [-gzip] [-time]"<<endl;
  cerr <<"                 [-icube]    cube  index "<<endl; 
  cerr <<"                 [-cdiff]    central difference"<<endl;
  cerr <<"                 [-iso]      isotropic diffusion "<<endl;
  cerr <<"                 [-mu]       extent of anisotropic diffusion"<< endl; 
  cerr <<"                 [-lambda]   extent of diffusion in each iteration " <<endl; 
  cerr <<"                 [-num_iter] number of iterations "<<endl;
  cerr <<"                  -k prints the curvature"<<endl;
  cerr <<"                  -n prints normals"<<endl;
  cerr <<"                  -m prints the m  values for the vertex"<<endl;
  cerr <<"                  -c <d>  d is the direction{0 1 2} prints the c for that direction "<<endl;
  cerr <<"                  -gradH_d_Normals <d> d is the direction "<<endl;
  cerr << endl;
}

void usage_error()
{
  usage_msg();
  exit(100);
}

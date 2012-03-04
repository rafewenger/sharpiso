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
#include <string>
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
bool flag_iso = false;
float large_magnitude = 2.0;

bool flag_norm_after = false;
bool flag_norm_before = false;
bool flag_m = false;
bool flag_mprev = false;
bool flag_c = false;
bool flag_w = false;
bool flag_print_curvature = false;
bool flag_vertex = false;
float mu(0.1);
float lambda(1.0);
int num_iter = 10;
VERTEX_INDEX vertex_index(0);

bool debug = true;

using namespace std;

// local subroutines
void memory_exhaustion();
void usage_error(), help();
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
    compute_gradient_central_difference(full_scalar_grid, vertex_index, gradient_grid);
            
    normalize_and_store_gradient_magnitudes(full_scalar_grid, gradient_grid, mag_list);
            
    if (flag_iso) {

      // Calculate the anisotropic diff of the gradients.
      const int dimension = full_scalar_grid.Dimension();
      gradient_grid.SetSize(full_scalar_grid, dimension);
      for (int k=0; k<num_iter; k++) {
        aniso_info.iter=k;
        if(flag_vertex) {
          aniso_info.flag_aniso = 0;
          compute_curvature
            (full_scalar_grid, gradient_grid, mu, lambda, vertex_index, aniso_info);
          print_aniso_info(full_scalar_grid, gradient_grid, aniso_info);
        }
                    
        anisotropic_diff_iter_k
          (full_scalar_grid, mu, lambda, k, 0, vertex_index, dimension, gradient_grid);
      }            
    }
    else {
      print_gradient_info(full_scalar_grid, vertex_index);
      cout << endl;

      cout << "Anisostropic gradients called."<<endl;

      const int dimension = full_scalar_grid.Dimension();
      gradient_grid.SetSize(full_scalar_grid, dimension);
      for (int k=0; k<num_iter; k++) {
        aniso_info.iter=k;
        cout <<"iteration "<<k<<endl;
        if(flag_vertex) {
          aniso_info.flag_aniso = 1;
          aniso_info.normals = gradient_grid.VectorPtr(vertex_index);
          compute_curvature
            (full_scalar_grid, gradient_grid, mu, lambda, vertex_index, aniso_info);
          print_aniso_info(full_scalar_grid, gradient_grid,aniso_info);
        }
        anisotropic_diff_iter_k
          (full_scalar_grid, mu, lambda, k, 1, vertex_index, dimension, gradient_grid);
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

    std::string s = string(argv[iarg]);

    if (s == "-time") 
      { report_time_flag = true;   }
    else if (s == "-iso")
      { flag_iso = true; }
    else if (s == "-k")
      { aniso_info.flag_k = true; }
    else if (s == "-m")
      { aniso_info.flag_m = true; }
    else if (s == "-n")
      { aniso_info.flag_normals = true; }
    else if (s == "-hdir") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &aniso_info.half_edge_direction);
    }
    else if (s == "-c") 
      { aniso_info.flag_print_c = true; }
    else if (s == "-gradS") 
      { aniso_info.flag_print_gradS = true; }
    else if (s == "-gradN") 
      { aniso_info.flag_print_gradN = true; }
    else if (s == "-cdiffN")
      { aniso_info.flag_print_cdiffN = true; }
    else if (s == "-fdiffN")
      { aniso_info.flag_print_fdiffN = true; }
    else if (s == "-mu") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%f", &mu);
    }
    else if (s == "-lambda") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%f", &lambda);
    }
    else if (s == "-num_iter") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &num_iter);
    }
    else if (s == "-vertex") {
      iarg++;
      flag_vertex = true;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &vertex_index);
    }
    else if (s == "-help") {
      help();
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
  cerr <<"Usage: anisogradinfo [OPTIONS]  {scalar nrrd file}"<<endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -time | -vertex <iv> | -iso | -mu |-lambda | -num_iter"
       << endl;
  cerr << "  -k | -n | -m | -hdir <dir> | -c | -gradS | -gradN | -cdiffN | -fdiffN" 
       << endl;
  cerr << "  -help" << endl;
}

void help()
{
  cerr << "Usage: anisogradinfo [OPTIONS]  {scalar nrrd file}"<<endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -time: Output running time." << endl;
  cerr << "  -vertex <iv>: Set vertex index to <iv>." << endl;
  cerr << "  -iso: Apply isotropic diffusion." << endl;
  cerr << "  -mu: Exponential anisotropic diffusion parameter." << endl;
  cerr << "  -lambda: Diffusion parameter." << endl;
  cerr << "  -num_iter: Number of iterations." << endl;
  cerr << "  -k: Print curvature." << endl;
  cerr << "  -n: Print normal vectors." << endl;
  cerr << "  -m: Prints m vectors." << endl;
  cerr << "  -hdir <d>: Offset point by half edge in direction d." << endl;
  cerr << "  -c: Print c vector." << endl;
  cerr << "  -gradS: Print gradients of scalar function."<<endl;
  cerr << "  -gradN: Print gradients of normals." << endl;
  cerr << "  -cdiffN: Print gradient of normals (central difference)." 
       << endl;
  cerr << "  -fdiffN: Print gradient of normals (forward difference)." 
       << endl;
  cerr << endl;

  exit(15);
}

void usage_error()
{
  usage_msg();
  exit(100);
}

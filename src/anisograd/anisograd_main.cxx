/// \file anisograd_main.cxx
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
#include "isodual3D_datastruct.h"

using namespace IJK;
using namespace ISODUAL3D;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
bool report_time_flag = false;
bool flag_gzip = false;
bool flag_cdiff = false;
bool flag_iso = false;
float large_magnitude = 2.0;

float mu(0.1);
float lambda(0.1);
int num_iter = 15;
VERTEX_INDEX icube = 0;

bool debug = true;

using namespace std;

// local subroutines
void memory_exhaustion();
void usage_error();
void parse_command_line(int argc, char **argv);

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
    time_t start_time;
    time(&start_time);
    
    vector<GRADIENT_TYPE> mag_list;
    
    IJK::ERROR error;
    
    try {
        
        std::set_new_handler(memory_exhaustion);
        
        parse_command_line(argc, argv);
        
        ISODUAL_SCALAR_GRID full_scalar_grid;
        GRID_NRRD_IN<int,int> nrrd_in;
        NRRD_DATA<int,int> nrrd_header;
        
        nrrd_in.ReadScalarGrid
        (scalar_filename, full_scalar_grid, nrrd_header, error);
        
        if (nrrd_in.ReadFailed()) { throw error; }
        
        GRADIENT_GRID gradient_grid;
        if(debug){
            cout <<" lambda " << lambda <<endl;
            cout <<" mu " << mu <<endl;
            cout <<" num_iter " << num_iter << endl;
            cout <<" icube "<<icube<<endl;
        }
        if (flag_cdiff) {
            compute_gradient_central_difference(full_scalar_grid, icube, gradient_grid);
        }
        else
        {
            // compute the central gradients first 
            compute_gradient_central_difference(full_scalar_grid, icube, gradient_grid);
            
            // Normalize the gradients and also 
            // store the magnitudes so that they can be later added back
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
                GRADIENT_TYPE  * N = gradient_grid.VectorPtr(iv);
                GRADIENT_TYPE   mag = 0.0;
                vector_magnitude (N, DIM3, mag);
                
                if (mag > 0.0001)
                {
                    mag_list.push_back(mag);
                    normalize (N, DIM3);
                    gradient_grid.Set(iv, N);
                }
                else 
                {
                    mag_list.push_back(0.0);
                }
            }
            
            
            
            if (flag_iso) {
                // Calculate the anisotropic diff of the gradients.
                cout << "isotropic diffusion called "<<endl;
                
                anisotropic_diff (full_scalar_grid,  mu, lambda, num_iter, 0, icube, gradient_grid);
            }
            else
            {
                // Calculate the anisotropic diff of the gradients.
                anisotropic_diff (full_scalar_grid,  mu, lambda, num_iter, 1, icube, gradient_grid);
            }
            
            
            // debug
            // reset the gradients to be normalized
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
                GRID_COORD_TYPE coord[DIM3];
                
                GRADIENT_TYPE  * N = gradient_grid.VectorPtr(iv);
                GRADIENT_TYPE   mag = 0.0;
                vector_magnitude (N, DIM3, mag);
                
                if (mag > 0.0001)
                {
                    normalize (N, DIM3);
                    gradient_grid.Set(iv, N);
                }
            }
            
            
            //reset the magnitudes of the gradients
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
                GRADIENT_TYPE  * N = gradient_grid.VectorPtr(iv);
                for (int i=0; i<DIM3; i++) {
                    N[i] = N[i]*mag_list[iv];
                }
                gradient_grid.Set(iv, N);
            }    
        }
        if (flag_gzip) {
            write_vector_grid_nrrd_gzip(gradient_filename, gradient_grid);
        }
        else {
            write_vector_grid_nrrd(gradient_filename, gradient_grid);
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
            if (iarg >= argc) { usage_error(); };
            sscanf(argv[iarg], "%d", &icube);
        }
        else 
        { usage_error(); }
        iarg++;
    }
    
    if (iarg+2 != argc) { usage_error(); };
    
    scalar_filename = argv[iarg];
    gradient_filename = argv[iarg+1];
}

void usage_msg()
{
    cerr <<"Usage: anisograd [options]  {scalar nrrd file} {gradient nrrd file}"<<endl;
    cerr <<"                 [-gzip] [-time]"<<endl;
    cerr <<"                 [-icube]    cube  index "<<endl; 
    cerr <<"                 [-cdiff]    central difference"<<endl;
    cerr <<"                 [-iso]      isotropic diffusion "<<endl;
    cerr <<"                 [-mu]       extent of anisotropic diffusion"<< endl; 
    cerr <<"                 [-lambda]   extent of diffusion in each iteration " <<endl; 
    cerr <<"                 [-num_iter] number of iterations "<<endl;
    cerr << endl;
}

void usage_error()
{
    usage_msg();
    exit(100);
}

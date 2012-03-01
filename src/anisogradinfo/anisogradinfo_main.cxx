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
bool flag_k =false;
bool flag_c = false;
bool flag_w = false;

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
        
        // debug
        //ISODUAL_SCALAR_GRID full_scalar_grid;
        SHARPISO_SCALAR_GRID full_scalar_grid;
        GRID_NRRD_IN<int,int> nrrd_in;
        NRRD_DATA<int,int> nrrd_header;
        
        nrrd_in.ReadScalarGrid
        (scalar_filename, full_scalar_grid, nrrd_header, error);
        
        if (nrrd_in.ReadFailed()) { throw error; }
        
        GRADIENT_GRID gradient_grid;
       
        if(false)
        {
            cout <<" lambda " << lambda <<endl;
            cout <<" mu " << mu <<endl;
            cout <<" num_iter " << num_iter << endl;
        }
        if (flag_cdiff) 
        {
            cout <<"num of iteration " << num_iter <<endl;
            compute_gradient_central_difference
            (full_scalar_grid, icube, gradient_grid);
        }
        else
        {
            cout <<" lambda (" << lambda<<")" <<endl;
            cout <<" mu " << mu <<")"<<endl;
            cout <<" num of iteration " << num_iter << endl;
            cout <<" Computing central gradients ..."<<endl;
            // compute the central gradients first 
            compute_gradient_central_difference(full_scalar_grid, icube, gradient_grid);
            
            // Normalize the gradients and also 
            // store the magnitudes so that they can be later added back
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
                GRADIENT_COORD_TYPE  * N = gradient_grid.VectorPtr(iv);
                GRADIENT_COORD_TYPE   mag = 0.0;
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
                cout << "Isotropic diffusion called "<<endl;
                
                anisotropic_diff
                (full_scalar_grid,  mu, lambda, num_iter, 0, icube, gradient_grid);
            }
            else
            {
                cout << "Anisostropic gradients called."<<endl;
                // Calculate the anisotropic diff of the gradients.
                anisotropic_diff
                (full_scalar_grid,  mu, lambda, num_iter, 1, icube, gradient_grid);
            }
            
            // reset the gradients to be normalized
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
                GRADIENT_COORD_TYPE  * N = gradient_grid.VectorPtr(iv);
                GRADIENT_COORD_TYPE   mag = 0.0;
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
                GRADIENT_COORD_TYPE  * N = gradient_grid.VectorPtr(iv);
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
        else if (string(argv[iarg]) == "-norm_before")
        { flag_norm_before = true; }
        else if (string(argv[iarg]) == "-norm_after")
        { flag_norm_after = true; }
        else if (string(argv[iarg]) == "-m")
        { flag_m = true; }
        else if (string(argv[iarg]) == "-mprev")
        { flag_mprev = true; }
        else if (string(argv[iarg]) == "-k")
        { flag_k = true; }
        else if (string(argv[iarg]) == "-c")
        { flag_c = true; }
        else if (string(argv[iarg]) == "-w")
        { flag_w = true; }
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
    cerr <<"Usage: anisogradinfo [options]  {scalar nrrd file} {gradient nrrd file}"<<endl;
    cerr <<"                 [-gzip] [-time]"<<endl;
    cerr <<"                 [-icube]    cube  index "<<endl; 
    cerr <<"                 [-cdiff]    central difference"<<endl;
    cerr <<"                 [-iso]      isotropic diffusion "<<endl;
    cerr <<"                 [-mu]       extent of anisotropic diffusion"<< endl; 
    cerr <<"                 [-lambda]   extent of diffusion in each iteration " <<endl; 
    cerr <<"                 [-num_iter] number of iterations "<<endl;
    cerr <<"                 -norm_before"<<endl;
    cerr <<"                 -norm_after"<<endl;
    cerr <<"                 -m"<<endl;
    cerr <<"                 -mprev"<<endl;
    cerr <<"                 -k"<<endl;
    cerr <<"                 -c"<<endl;
    cerr <<"                 -w"<<endl;
    cerr << endl;
}

void usage_error()
{
    usage_msg();
    exit(100);
}

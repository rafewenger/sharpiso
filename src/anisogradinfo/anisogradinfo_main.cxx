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
        
        ANISO_INFO_TYPE aniso_info;
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
                cout << "Isotropic diffusion called "<<endl;
                
                anisotropic_diff_debug 
                (full_scalar_grid,  mu, lambda, num_iter, 0, icube, gradient_grid, aniso_info);
            }
            else
            {
                cout << "Anisostropic gradients called."<<endl;
                // Calculate the anisotropic diff of the gradients.
                anisotropic_diff_debug
                (full_scalar_grid,  mu, lambda, num_iter, 1, icube, gradient_grid, aniso_info);
            }
            
            // reset the gradients to be normalized
            for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
            {
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

/*******************/
//DEBUG
void debug_print(ANISO_INFO_TYPE & aniso_info)
{
    cout.precision(5);
    cout << " vertex " << aniso_info.iv;
    cout << " num iteration " << aniso_info.num_iter<<endl;
    
    if(flag_norm_before){
    cout << " norm before  (";
    
    cout <<aniso_info.Normals[0]<<" " << aniso_info.Normals[1]<<" " <<aniso_info.Normals[2]<<") "<<endl;
    }
    if(flag_norm_after)
    {
        cout << " norm after  (";
        cout <<aniso_info.Normals2[0] <<" "<< aniso_info.Normals2[1]<<" " <<aniso_info.Normals2[2]<<") "<<endl;

    }
    if(flag_m)
    {
        cout << " mX (";
        cout <<aniso_info.mX[0] <<" "<< aniso_info.mX[1]<<" " <<aniso_info.mX[2]<<") "<<endl;
        cout << " mY (";
        cout <<aniso_info.mY[0] <<" "<< aniso_info.mY[1]<<" " <<aniso_info.mY[2]<<") "<<endl;
        cout << " mZ (";
        cout <<aniso_info.mZ[0] <<" "<< aniso_info.mZ[1]<<" " <<aniso_info.mZ[2]<<") "<<endl;
    }
    if(flag_mprev)
    {
        cout << " mX prev vert "<<aniso_info.prev_vert[0]<<" (";
        cout <<aniso_info.mX_prev_vert_X[0] <<" "<< aniso_info.mX_prev_vert_X[1]<<" "<<aniso_info.mX_prev_vert_X[2]<<") "<<endl;
        cout << " mY prev vert "<<aniso_info.prev_vert[1]<<" (";
        cout <<aniso_info.mY_prev_vert_Y[0] <<" "<< aniso_info.mY_prev_vert_Y[1]<<" " <<aniso_info.mY_prev_vert_Y[2]<<") "<<endl;
        cout << " mZ prev vert "<<aniso_info.prev_vert[2]<<" (";
        cout <<aniso_info.mZ_prev_vert_Z[0] <<" "<< aniso_info.mZ_prev_vert_Z[1]<<" " <<aniso_info.mZ_prev_vert_Z[2]<<") "<<endl;
    }
    if(flag_k)
    {
        cout << " k  (";
        cout <<aniso_info.K[0] <<" "<< aniso_info.K[1]<<" " <<aniso_info.K[2]<<") "<<endl;

        cout << " k prev {"<< aniso_info.prev_vert[0]<<" "<< aniso_info.prev_vert[1]<<" "<<aniso_info.prev_vert[2]<<"} (";
        cout <<aniso_info.K[0] <<" "<< aniso_info.K[1]<<" " <<aniso_info.K[2]<<") "<<endl;

    }
    if (flag_c) {
        cout <<" c \n";
        for (int i=0; i<DIM3; i++) {
            cout <<aniso_info.c[DIM3*i + 0]<<","<<aniso_info.c[DIM3*i + 1]<<","<<aniso_info.c[DIM3*i + 2]<<"\n";
        }
        cout <<"\n";
        cout <<" c prev\n";
        for (int i=0; i<DIM3; i++) {
            cout <<aniso_info.c_prev[DIM3*i + 0]<<","<<aniso_info.c_prev[DIM3*i + 1]<<","<<aniso_info.c_prev[DIM3*i + 2]<<"\n";
        }
        cout <<"\n";
    }
    
    if(flag_w)
    {
        cout <<" w \n ("<<aniso_info.w[0]<<" "<<aniso_info.w[1]<<" "<<aniso_info.w[2]<<")"<<endl;
       
        cout <<" w' \n ("<<aniso_info.w_dash[0]<<" "<<aniso_info.w_dash[1]<<" "<<aniso_info.w_dash[2]<<")"<<endl;
    }
}
/*
 class ANISO_INFO_TYPE{
 public :
 VERTEX_INDEX iv;
 int num_iter;
 GRADIENT_TYPE Normals[DIM3];
 GRADIENT_TYPE Normals2[DIM3];
 GRADIENT_TYPE mX[DIM3];
 GRADIENT_TYPE mX_prev_vert_X[DIM3];
 
 GRADIENT_TYPE mY[DIM3];
 GRADIENT_TYPE mY_prev_vert_Y[DIM3];
 
 GRADIENT_TYPE mZ[DIM3];
 GRADIENT_TYPE mZ_prev_vert_Z[DIM3];
 VERTEX_INDEX prev_vert[DIM3];
 
 SCALAR_TYPE K[DIM3];
 SCALAR_TYPE gK[DIM3];
 
 SCALAR_TYPE Kprev[DIM3];
 SCALAR_TYPE gKprev[DIM3];
 
 GRADIENT_TYPE   gradHNd[DIM9*DIM3];
 GRADIENT_TYPE   gradHNd_prev[DIM9*DIM3];
 GRADIENT_TYPE   c[DIM3*DIM3];
 GRADIENT_TYPE   c_prev[DIM3*DIM3];
 GRADIENT_TYPE   w[DIM3];
 GRADIENT_TYPE   wN ;
 GRADIENT_TYPE   w_dash[DIM3];
 void mycopy(const GRADIENT_TYPE temp[], GRADIENT_TYPE temp1[],int  n)
 {
 for(int i=0;i<n;i++)
 {
 temp1[i] = temp[i];
 }
 };
 
 };
 */

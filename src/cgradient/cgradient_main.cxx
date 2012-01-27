/// \file cgradient_main.cxx
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

#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"

#include "cgradient.h"
#include "isodual3D_datastruct.h"

using namespace IJK;
using namespace ISODUAL3D;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
bool report_time_flag = false;
bool flag_gzip = false;

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
    compute_gradient_central_difference(full_scalar_grid, gradient_grid);

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
  cerr << "Usage: cgradient [-gzip] [-time] {scalar nrrd file} {gradient nrrd file}"
       << endl;
}

void usage_error()
{
  usage_msg();
  exit(100);
}

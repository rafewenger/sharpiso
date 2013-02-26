/// \file grad2hermite.cxx
/// Convert gradients to hermite data
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2013 Rephael Wenger

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


#include <fstream>
#include <iostream>

#include "ijkIO.txx"
#include "ijkgrid_nrrd.txx"
#include "ijkstring.txx"

#include "ijkNrrd.h"
#include "ijkgrid_macros.h"

#include "sharpiso_types.h"
#include "sharpiso_grids.h"
#include "sharpiso_get_gradients.h"

using namespace IJK;
using namespace SHARPISO;
using namespace std;

// global variables
char * scalar_filename(NULL);
char * gradient_filename(NULL);
char * normals_filename(NULL);
SCALAR_TYPE isovalue;
GRADIENT_COORD_TYPE max_small_magnitude(0.0001);

// Computation routine.
void compute_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & edgeI_coord,
 std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord);

// I/O routine.
void write_off_file
(const char * filename, 
 const std::vector<COORD_TYPE> & coord,
 const std::vector<GRADIENT_COORD_TYPE> & normal_coord);

// local subroutines
void memory_exhaustion();
void parse_command_line(int argc, char **argv);
void usage_error();
string remove_nrrd_suffix(const string & filename);
void construct_gradient_filename
(const char * scalar_filename, std::string & gradient_filename);

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  SHARPISO_SCALAR_GRID scalar_grid;
  GRADIENT_GRID gradient_grid;
  std::vector<COORD_TYPE> edgeI_coord;
  std::vector<GRADIENT_COORD_TYPE> edgeI_normal_coord;
  IJK::ERROR error;

  try {
    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    GRID_NRRD_IN<int,int> nrrd_in;
    NRRD_DATA<int,int> nrrd_header;

    nrrd_in.ReadScalarGrid
      (scalar_filename, scalar_grid, nrrd_header, error);
    if (nrrd_in.ReadFailed()) { throw error; }

    if (gradient_filename == NULL) {
      std::string gradient_filename;
      construct_gradient_filename(scalar_filename, gradient_filename);
      nrrd_in.ReadVectorGrid
        (gradient_filename.c_str(), gradient_grid, error);
    }
    else {
      nrrd_in.ReadVectorGrid
        (gradient_filename, gradient_grid, error);
    }
    if (nrrd_in.ReadFailed()) { throw error; }

    if (!gradient_grid.CompareSize(scalar_grid)) {
      error.AddMessage("Input error. Grid mismatch.");
      error.AddMessage
        ("  Dimension or axis sizes of gradient grid and scalar grid do not match.");
        throw error;
    }

    compute_edgeI(scalar_grid, gradient_grid, isovalue, 
                  edgeI_coord, edgeI_normal_coord);

    write_off_file(normals_filename, edgeI_coord, edgeI_normal_coord);
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
// COMPUTE GRID EDGE-ISOSURFACE INTERSECTIONS
// **************************************************

void compute_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 std::vector<COORD_TYPE> & edgeI_coord,
 std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord)
{
  COORD_TYPE p[DIM3];
  GRADIENT_COORD_TYPE normal[DIM3];

  IJK_FOR_EACH_GRID_EDGE(iv0, edge_dir, scalar_grid, VERTEX_INDEX) {

    VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, edge_dir);
    if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

      NUM_TYPE num_coord = edgeI_coord.size();
      edgeI_coord.resize(num_coord+DIM3);
      edgeI_normal_coord.resize(num_coord+DIM3);

      compute_isosurface_grid_edge_intersection
        (scalar_grid, gradient_grid, isovalue,
         iv0, iv1, edge_dir, max_small_magnitude, 
         &(edgeI_coord.front())+num_coord,
         &(edgeI_normal_coord.front())+num_coord);
    }
  }
}


// **************************************************
// I/O ROUTINES
// **************************************************

void write_off_file
(const char * filename, 
 const std::vector<COORD_TYPE> & coord,
 const std::vector<GRADIENT_COORD_TYPE> & normal_coord)
{
  std::string outfilename;
  IJK::PROCEDURE_ERROR error("write_off_file");

  if (filename != NULL) {
    outfilename = filename;
  }
  else {
    outfilename = remove_nrrd_suffix(scalar_filename);
    outfilename += ".normals.off";
  }

  ofstream out(outfilename.c_str(), ios::out);
  if (!out.good()) {
    error.AddMessage
      ("Error.  Unable to open file ", outfilename, ".");
    throw error;
  };

  // Dummy array containing no simplex vertices.
  std::vector<VERTEX_INDEX> simplex_vert;

  ijkoutNormalsOFF(out, DIM3, NUM_VERT_PER_TRI,
                   coord, normal_coord, simplex_vert);

  out.close();

  cout << "Wrote output to file: " << outfilename << endl;
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

  if (argc < (iarg+2) || (iarg+4) < argc) 
    { usage_error(); };

  if (!string2val(argv[iarg], isovalue)) {
    cerr << "Input error.  Illegal isovalue: " << isovalue << "." << endl;
    usage_error();
  }
  iarg++;

  scalar_filename = argv[iarg];
  iarg++;

  if (iarg < argc) {
    gradient_filename = argv[iarg];
    iarg++;
  }

  if (iarg < argc)
    { normals_filename = argv[iarg]; }
}

void usage_msg()
{
  cerr << "Usage: grad2hermite {isovalue} {scalar nrrd file} [gradient nrrd file] [normals off filename]"
       << endl;
}

void usage_error()
{
  usage_msg();
  exit(100);
}

#ifdef _WIN32
const char PATH_DELIMITER = '\\';
#else
const char PATH_DELIMITER = '/';
#endif

string remove_nrrd_suffix(const string & filename)
{
  string prefix, suffix;

  // create output filename
  string fname = filename;

#ifndef _WIN32
  // remove path from file name
  split_string(fname, PATH_DELIMITER, prefix, suffix);
  if (suffix != "") { fname = suffix; }
#endif

  // construct output filename
  split_string(fname, '.', prefix, suffix);
  if (suffix == "nrrd" || suffix == "nhdr") { fname = prefix; }
  else { fname = filename; }

  return(fname);
}

// Construct gradient filename from scalar filename.
void construct_gradient_filename
(const char * scalar_filename, std::string & gradient_filename)
{
	std::string prefix;
	std::string suffix;

	split_string(scalar_filename, '.', prefix, suffix);

	gradient_filename = prefix + ".grad." + suffix;
}

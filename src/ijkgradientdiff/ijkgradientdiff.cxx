/// \file ijkgradientdiff.cxx
/// Find difference between two gradient vector fields
/// Version v0.1.3

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011-2015 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include "ijkcoord.txx"
#include "ijkgrid_nrrd.txx"
#include "ijkprint.txx"
#include "ijkscalar_grid.txx"
#include "ijkstring.txx"
#include "ijkvector_grid.txx"

using namespace std;
using namespace IJK;

// types
typedef float SCALAR_TYPE;
typedef int AXIS_SIZE_TYPE;
typedef int LENGTH_TYPE;
typedef int VERTEX_INDEX;
typedef int GRID_COORD_TYPE;
typedef float COORD_TYPE;
typedef float GRADIENT_COORD_TYPE;
typedef float ANGLE_TYPE;
typedef double WEIGHT_TYPE;
typedef int NUM_TYPE;

typedef enum { IGNORE_SMALL_MAG_ZERO, IGNORE_SMALL_MAG_ONE,
               IGNORE_SMALL_MAG_EITHER, IGNORE_SMALL_MAG_BOTH }
  IGNORE_SMALL_MAG;

typedef IJK::GRID_PLUS<NUM_TYPE, AXIS_SIZE_TYPE, VERTEX_INDEX, NUM_TYPE> 
INFO_GRID;
typedef IJK::SCALAR_GRID<INFO_GRID, SCALAR_TYPE> SGRID;
typedef IJK::VECTOR_GRID<INFO_GRID, LENGTH_TYPE, GRADIENT_COORD_TYPE> 
GRADIENT_GRID;

// global variables
char * gradient0_filename(NULL);
char * gradient1_filename(NULL);
char * scalar_filename(NULL);
ANGLE_TYPE min_angle(0);
GRADIENT_COORD_TYPE min_magnitude(0);
VERTEX_INDEX vertex_index(0);
bool flag_diff_all(true);
bool flag_report_each_diff(true);
bool flag_format_numeric(false);
IGNORE_SMALL_MAG ignore_small_mag(IGNORE_SMALL_MAG_BOTH);


// output routines
void output_gradient_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1);
void output_gradient_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1,
 const VERTEX_INDEX vertex_index);


// misc routines
void split_string(const string & s, const char c,
                  string & prefix, string & suffix);
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  IJK::PROCEDURE_ERROR error("ijkgradientdiff");

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in_gradient0;
    GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in_gradient1;
    GRADIENT_GRID gradient0_grid;
    GRADIENT_GRID gradient1_grid;

    nrrd_in_gradient0.ReadVectorGrid
      (gradient0_filename, gradient0_grid, error);
    if (nrrd_in_gradient0.ReadFailed()) { throw error; }

    nrrd_in_gradient1.ReadVectorGrid
      (gradient1_filename, gradient1_grid, error);
    if (nrrd_in_gradient1.ReadFailed()) { throw error; }

    GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in_scalar;
    SGRID scalar_grid;


    if (flag_diff_all) {
      output_gradient_diff(gradient0_grid, gradient1_grid);
    }
    else {
      output_gradient_diff(gradient0_grid, gradient1_grid, vertex_index);
    }
  }
  catch(IJK::ERROR error) {
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
// Compute routines
// **************************************************

template <typename DTYPE, typename CTYPE0, typename CTYPE1,
          typename CTYPE2>
void compute_cos_angle
(const DTYPE dimension, const CTYPE0 v0[], const CTYPE1 v1[],
 CTYPE2 & cos_angle)
{
  CTYPE2 magnitude0, magnitude1;
  IJK::compute_magnitude(dimension, v0, magnitude0);
  IJK::compute_magnitude(dimension, v1, magnitude1);

  cos_angle = 1;
  if (magnitude0 == 0 || magnitude1 == 0) { return; }

  cos_angle = 0;
  CTYPE2 scale = magnitude0*magnitude1;
  for (DTYPE d = 0; d < dimension; d++) {
    cos_angle += (v0[d]*v1[d]/scale);
  }
}


// **************************************************
// Output Information
// **************************************************

// Output gradients and angle between gradients.
template <typename DTYPE>
void output_gradients
(const DTYPE dimension, const VERTEX_INDEX iv,
 const GRADIENT_COORD_TYPE * v0, const GRADIENT_COORD_TYPE * v1,
 const ANGLE_TYPE cos_angle)
{
  cout << "Vertex " << iv;
  cout << " grad0: ";
  IJK::print_list(cout, v0, dimension);
  cout << " grad1: ";
  IJK::print_list(cout, v1, dimension);
  cout << " angle: "<< std::acos(cos_angle)*180.0/M_PI;
  cout << endl;
}

// Output gradients
template <typename DTYPE>
void output_gradients
(const DTYPE dimension, const VERTEX_INDEX iv,
 const GRADIENT_COORD_TYPE * v0, const GRADIENT_COORD_TYPE * v1)
{
  cout << "Vertex " << iv;
  cout << " grad0: ";
  IJK::print_list(cout, v0, dimension);
  cout << " grad1: ";
  IJK::print_list(cout, v1, dimension);
  cout << endl;
}

// Output coordinates of gradients and angle between gradients.
// Only numeric values
template <typename DTYPE>
void output_gradients_numeric
(const DTYPE dimension, const VERTEX_INDEX iv,
 const GRADIENT_COORD_TYPE * v0, const GRADIENT_COORD_TYPE * v1,
 const ANGLE_TYPE cos_angle)
{
  cout << iv << " ";
  IJK::print_list_values(cout, v0, dimension, ' ');
  IJK::print_list_values(cout, v1, dimension, ' ');
  cout << " " << std::acos(cos_angle)*180.0/M_PI;
  cout << endl;
}

bool output_gradient_angular_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1,
 const GRADIENT_COORD_TYPE cos_angle_threshold, const VERTEX_INDEX iv,
 const GRADIENT_COORD_TYPE magnitude0, const GRADIENT_COORD_TYPE magnitude1,
 GRADIENT_COORD_TYPE & cos_angle);

void output_gradient_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1)
{
  const NUM_TYPE dimension = gradient0.Dimension();
  const GRADIENT_COORD_TYPE cos_angle_threshold 
    = std::cos(min_angle*M_PI/180.0);
  GRADIENT_COORD_TYPE cos_angle;
  IJK::PROCEDURE_ERROR error("error");
  
  if (!gradient0.Check
      (gradient1, "gradient0", "gradient1", error)) 
    { throw error; }

  NUM_TYPE icount = 0;
  GRADIENT_COORD_TYPE min_cos_angle = 1.0;
  VERTEX_INDEX max_angle_iv = 0;

  if (ignore_small_mag == IGNORE_SMALL_MAG_BOTH) {

    for (VERTEX_INDEX iv = 0; iv < gradient0.NumVertices(); iv++) {
      GRADIENT_COORD_TYPE magnitude0 = gradient0.ComputeMagnitude(iv);
      GRADIENT_COORD_TYPE magnitude1 = gradient1.ComputeMagnitude(iv);

      if (magnitude0 >= min_magnitude || magnitude1 >= min_magnitude) {

        if (output_gradient_angular_diff
            (gradient0, gradient1, cos_angle_threshold, iv, 
             magnitude0, magnitude1, cos_angle)) {
          icount++; 
          if (cos_angle < min_cos_angle) {
            min_cos_angle = cos_angle; 
            max_angle_iv = iv;
          }
        }
      }
    }
  }
  else if (ignore_small_mag == IGNORE_SMALL_MAG_EITHER) {

    for (VERTEX_INDEX iv = 0; iv < gradient0.NumVertices(); iv++) {
      GRADIENT_COORD_TYPE magnitude0 = gradient0.ComputeMagnitude(iv);
      GRADIENT_COORD_TYPE magnitude1 = gradient1.ComputeMagnitude(iv);

      if (magnitude0 >= min_magnitude && magnitude1 >= min_magnitude) {

        if (output_gradient_angular_diff
            (gradient0, gradient1, cos_angle_threshold, iv, 
             magnitude0, magnitude1, cos_angle)) {
          icount++; 
          if (cos_angle < min_cos_angle) { 
            min_cos_angle = cos_angle; 
            max_angle_iv = iv;
          }
        }
      }
    }
  }
  else if (ignore_small_mag == IGNORE_SMALL_MAG_ZERO) {

    for (VERTEX_INDEX iv = 0; iv < gradient0.NumVertices(); iv++) {
      GRADIENT_COORD_TYPE magnitude0 = gradient0.ComputeMagnitude(iv);
      GRADIENT_COORD_TYPE magnitude1 = gradient1.ComputeMagnitude(iv);

      if (magnitude0 >= min_magnitude) {

        if (output_gradient_angular_diff
            (gradient0, gradient1, cos_angle_threshold, iv, 
             magnitude0, magnitude1, cos_angle)) {
          icount++; 
          if (cos_angle < min_cos_angle) { 
            min_cos_angle = cos_angle; 
            max_angle_iv = iv;
          }
        }
      }
    }
  }
  else if (ignore_small_mag == IGNORE_SMALL_MAG_ONE) {

    for (VERTEX_INDEX iv = 0; iv < gradient0.NumVertices(); iv++) {
      GRADIENT_COORD_TYPE magnitude0 = gradient0.ComputeMagnitude(iv);
      GRADIENT_COORD_TYPE magnitude1 = gradient1.ComputeMagnitude(iv);

      if (magnitude1 >= min_magnitude) {

        if (output_gradient_angular_diff
            (gradient0, gradient1, cos_angle_threshold, iv, 
             magnitude0, magnitude1, cos_angle)) {
          icount++;
          if (cos_angle < min_cos_angle) {
            min_cos_angle = cos_angle; 
            max_angle_iv = iv;
          }
        }
      }
    }
  }
  else {
    error.AddMessage("Programming error.  Illegal value for ignore_small_mag.");
    throw error;
  }

  if (flag_format_numeric) {
    if (!flag_report_each_diff) {
      cout << icount << endl;
    }
  }
  else {
    cout << "Total: " << icount
         << " differences between gradients." << endl;
    cout << "Max difference (vertex " << max_angle_iv << "): " 
         << std::acos(min_cos_angle)*180.0/M_PI << endl;
  }
}

// Output gradients with large angular difference.
// Return true if gradients have large angular difference.
bool output_gradient_angular_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1,
 const GRADIENT_COORD_TYPE cos_angle_threshold, const VERTEX_INDEX iv,
 const GRADIENT_COORD_TYPE magnitude0, const GRADIENT_COORD_TYPE magnitude1,
 GRADIENT_COORD_TYPE & cos_angle)
{
  const NUM_TYPE dimension = gradient0.Dimension();
  const GRADIENT_COORD_TYPE * v0 = gradient0.VectorPtrConst(iv);
  const GRADIENT_COORD_TYPE * v1 = gradient1.VectorPtrConst(iv);

  compute_cos_angle(dimension, v0, v1, cos_angle);

  if (cos_angle < cos_angle_threshold) {

    if (flag_report_each_diff) {
      if (flag_format_numeric) 
        { output_gradients_numeric(dimension, iv, v0, v1, cos_angle); }
      else 
        { output_gradients(dimension, iv, v0, v1, cos_angle); }
    }

    return (true);
  }
  else if (magnitude0 == 0 || magnitude1 == 0) {
    // Special handling for 0 magnitudes.
    if (magnitude0 != magnitude1) {
      if (flag_report_each_diff) {
        if (flag_format_numeric) {
          // Output 0 as default angle.
          const float default_cos_angle = 1.0;
          output_gradients_numeric
            (dimension, iv, v0, v1, default_cos_angle);
        }
        else
          { output_gradients(dimension, iv, v0, v1); }
      }

      return(true);
    }
  }

  return(false);
}

void output_gradient_diff
(const GRADIENT_GRID & gradient0, const GRADIENT_GRID & gradient1,
 const VERTEX_INDEX vertex_index)
{
  const NUM_TYPE dimension = gradient0.Dimension();
  const GRADIENT_COORD_TYPE * v0 = gradient0.VectorPtrConst(vertex_index);
  const GRADIENT_COORD_TYPE * v1 = gradient1.VectorPtrConst(vertex_index);

  ANGLE_TYPE cos_angle;
  compute_cos_angle(dimension, v0, v1, cos_angle);
  output_gradients(dimension, vertex_index, v0, v1, cos_angle);
}

// **************************************************
// Misc Routines
// **************************************************

float get_float(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  float x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

int get_int(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  int x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];

    if (s == "-minmag") {
      min_magnitude = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-angle") {
      min_angle = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-ignore_small_mag") {
      iarg++;
      if (iarg >= argc) {
        cerr << "Usage error.  Missing parameter for -ignore_small_mag."
             << endl;
        usage_error();
      }
      string s2 = argv[iarg];
      if (s2 == "0" || s2 == "zero") 
        { ignore_small_mag = IGNORE_SMALL_MAG_ZERO; }
      else if (s2 == "1" || s2 == "one")
        { ignore_small_mag = IGNORE_SMALL_MAG_ONE; }
      else if (s2 == "either")
        { ignore_small_mag = IGNORE_SMALL_MAG_EITHER; }
      else if (s2 == "both")
        { ignore_small_mag = IGNORE_SMALL_MAG_BOTH; }
      else {
        cerr << "Usage error.  Illegal parameter for -ignore_small_mag."
             << endl;
        usage_error();
      }
    }
    else if (s == "-vertex") {
      vertex_index = get_int(iarg, argc, argv);
      iarg++;
      flag_diff_all = false;
    }
    else if (s == "-total") {
      flag_report_each_diff = false;
    }
    else if (s == "-fnum") {
      flag_format_numeric = true;
    }
    else if (s == "-help") 
      { help(); }
    else {
      cerr << "Usage error. Illegal parameter: " << argv[iarg] << endl;
      cerr << endl;
      usage_error();
    };
    
    iarg++;
  };

  if (iarg+2 != argc) { usage_error(); }

  gradient0_filename = argv[iarg];
  gradient1_filename = argv[iarg+1];
}

void usage_error()
{
  cerr << "Usage: ijkgradientdiff [OPTIONS]"
       << " <gradient0 filename> <gradient1 filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -minmag <mag> | -angle <angle>" << endl;
  cerr << "  -ignore_small_mag {0|1|either|both}" << endl;
  cerr << "  -vertex <iv> | -total | -fnum | -help" << endl;
  exit(10);
}

void help()
{
  cerr << "Usage: ijkgradientdiff [OPTIONS]"
       << " <gradient0 filename> <gradient1 filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -minmag <mag>:  Ignore gradients with magnitude under <mag>."
       << endl;
  cerr << "  -angle <angle>: Report angles greater than <angle>." << endl;
  cerr << "  -ignore_small_mag 0:  Ignore if gradient0 is at most <mag>."
       << endl;
  cerr << "  -ignore_small_mag 1:  Ignore if gradient1 is at most <mag>."
       << endl;
  cerr << "  -ignore_small_mag either:" << endl
       << "     Ignore if either gradient0 or gradient1 is at most <mag>."
       << endl;
  cerr << "  -ignore_small_mag both:" << endl
       << "     Ignore if both gradient0 and gradient1 is at most <mag>."
       << endl;
  cerr << "  -vertex <iv>:   Print gradient differences for vertex <iv>." 
       << endl;
  cerr << "  -total:         Print only total number of gradient differences."
       << endl;
  cerr << "  -fnum:          Numeric format.  Print only numeric values."
       << endl;
  cerr << "  -help:          Print this help message." << endl;

  exit(10);
}

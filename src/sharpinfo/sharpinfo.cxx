/// \file sharpinfo.cxx
/// Output information about sharp isosurface vertices and edges.
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

#include "ijkcoord.txx"
#include "ijkgrid_macros.h"
#include "ijkgrid_nrrd.txx"

#include "sharpiso_feature.h"
#include "sharpiso_eigen.h"
#include "sharpiso_types.h"
#include "sharpiso_scalar.txx"

using namespace SHARPISO;

using namespace std;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
SCALAR_TYPE isovalue(0);
GRADIENT_COORD_TYPE max_small_mag(0.0);
EIGENVALUE_TYPE max_small_eigenvalue(0.1);
VERTEX_INDEX cube_index(0);
AXIS_SIZE_TYPE subgrid_axis_size(3);
bool flag_list_gradients(false);
bool flag_list_subgrid(false);
std::vector<COORD_TYPE> location;
bool flag_location_set(false);
bool flag_isovalue_set(false);
bool flag_use_neighboring_facets(false);
bool flag_use_selective_neighbors(true);
bool flag_list_eigen(false);

bool flag_svd_gradients(true); //default
bool flag_svd_edges_simple(false);
bool flag_svd_edges_cmplx(false);
COORD_TYPE cube_offset(0);

// check routines
bool check_gradient_grid
(const GRADIENT_GRID & gradient_grid, IJK::ERROR & error);
bool check_input_grids
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid, IJK::ERROR & error);

// output routines
void output_cube_coordinates
(std::ostream & output,
 const SHARPISO_SCALAR_GRID & scalar_grid,
 const VERTEX_INDEX icube);
void output_gradients
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_points);
void output_svd_results
(std::ostream & output, const COORD_TYPE sharp_coord[DIM3],
 const EIGENVALUE_TYPE eigenvalues[DIM3], const NUM_TYPE num_large_eigenvalues,
 const EIGENVALUE_TYPE eigenvalue_tolerance, SVD_INFO & svd_info);
void output_subgrid_results
(std::ostream & output, const COORD_TYPE sharp_coord[DIM3],
 const SCALAR_TYPE scalar_stdev, const SCALAR_TYPE max_abs_scalar_error);
void output_cube_subgrid_scalar_errors
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar, const NUM_TYPE num_points,
 const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
 const NUM_TYPE subgrid_axis_size);
void output_gradient_based_scalars
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar, const NUM_TYPE num_points,
 const std::vector<COORD_TYPE> & location);
void output_cube_eigenvalues
(std::ostream & output,
 const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE max_zero_mag,
 const EIGENVALUE_TYPE eigenvalue_tolerance);



// local subroutines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  GRID_COORD_TYPE cube_coord[DIM3];
  IJK::ERROR error;
  // accumulate information from svd_calls
  SVD_INFO svd_info;
  svd_info.ray_intersect_cube = false;
  svd_info.location = NONE;

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    SHARPISO_SCALAR_GRID scalar_grid;
    IJK::GRID_NRRD_IN<NUM_TYPE,AXIS_SIZE_TYPE> nrrd_in_scalar;
    nrrd_in_scalar.ReadScalarGrid
      (scalar_filename, scalar_grid,  error);

    GRADIENT_GRID gradient_grid;
    IJK::GRID_NRRD_IN<NUM_TYPE,AXIS_SIZE_TYPE> nrrd_in_gradient;

    nrrd_in_gradient.ReadVectorGrid
      (gradient_filename, gradient_grid, error);

    if (!check_gradient_grid(gradient_grid, error))
      { throw error; };

    if (!check_input_grids(scalar_grid, gradient_grid, error))
      { throw error; }

    NUM_TYPE num_gradients = 0;
    std::vector<COORD_TYPE> point_coord;
    std::vector<GRADIENT_COORD_TYPE> gradient_coord;
    std::vector<SCALAR_TYPE> scalar;

    if (flag_use_neighboring_facets) {

      if (flag_use_selective_neighbors) {
        OFFSET_CUBE_111 cube_111(cube_offset);
        get_selected_cube_neighbor_gradients
          (scalar_grid, gradient_grid, cube_index, max_small_mag, isovalue,
           point_coord, gradient_coord, scalar, num_gradients, cube_111);
      }
      else {
        get_large_cube_neighbor_gradients
          (scalar_grid, gradient_grid, cube_index, max_small_mag,
           point_coord, gradient_coord, scalar, num_gradients);
      }
    }
    else {
      get_large_cube_gradients
        (scalar_grid, gradient_grid, cube_index, max_small_mag,
         point_coord, gradient_coord, scalar, num_gradients);
    }

    if (flag_list_gradients) {
      output_gradients
        (cout, point_coord, gradient_coord, scalar, num_gradients);
      cout << endl;
    }

    if (flag_isovalue_set) {

      COORD_TYPE sharp_coord[DIM3];
      EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
      NUM_TYPE num_large_eigenvalues(0);

      if (flag_svd_edges_simple){

        svd_compute_sharp_vertex_in_cube_edge_based_simple
          (scalar_grid, gradient_grid, cube_index, isovalue,
           max_small_mag, max_small_eigenvalue, sharp_coord, eigenvalues,
           num_large_eigenvalues, svd_info);
      }
      else if (flag_svd_edges_cmplx){

        svd_compute_sharp_vertex_in_cube_edge_based_cmplx
          (scalar_grid, gradient_grid, cube_index, isovalue,
           max_small_mag, max_small_eigenvalue, sharp_coord, eigenvalues,
           num_large_eigenvalues, svd_info);
      }
      else if(flag_svd_gradients) {

        if (flag_use_neighboring_facets &&
            flag_use_selective_neighbors) {
          OFFSET_CUBE_111 cube_111(cube_offset);

          // Use only cube gradients
          svd_compute_sharp_vertex_neighborhood
            (scalar_grid, gradient_grid, cube_index, isovalue,
             max_small_mag, max_small_eigenvalue, sharp_coord, eigenvalues,
             num_large_eigenvalues, svd_info, cube_111);
        }
        else {
          // Use only cube gradients
          svd_compute_sharp_vertex_in_cube
            (scalar_grid, gradient_grid, cube_index, isovalue,
             max_small_mag, max_small_eigenvalue, sharp_coord, eigenvalues,
             num_large_eigenvalues, svd_info);
        }

      }

      output_svd_results
        (cout, sharp_coord, eigenvalues, num_large_eigenvalues,
         max_small_eigenvalue, svd_info);
      cout << endl;

      if (flag_list_subgrid) {
        scalar_grid.ComputeCoord(cube_index, cube_coord);
        output_cube_subgrid_scalar_errors
          (cout, point_coord, gradient_coord, scalar, num_gradients,
           cube_coord, isovalue, subgrid_axis_size);
        cout << endl;
      }

      SCALAR_TYPE scalar_stdev;
      SCALAR_TYPE max_abs_scalar_error;
      if (flag_use_neighboring_facets) {
        if (flag_use_selective_neighbors) {
          OFFSET_CUBE_111 cube_111(cube_offset);
          subgrid_compute_sharp_vertex_neighborhood
            (scalar_grid, gradient_grid, cube_index, isovalue,
             max_small_mag, subgrid_axis_size,
             sharp_coord, scalar_stdev, max_abs_scalar_error,
             cube_111);
        }
        else {
          scalar_grid.ComputeCoord(cube_index, cube_coord);
          subgrid_calculate_iso_vertex_in_cube
            (point_coord, gradient_coord, scalar,
             num_gradients, cube_coord, isovalue, subgrid_axis_size,
             sharp_coord, scalar_stdev, max_abs_scalar_error);
        }
      }
      else {
        subgrid_compute_sharp_vertex_in_cube
          (scalar_grid, gradient_grid, cube_index, isovalue,
           max_small_mag, subgrid_axis_size,
           sharp_coord, scalar_stdev, max_abs_scalar_error);
      }

      output_subgrid_results
        (cout, sharp_coord, scalar_stdev, max_abs_scalar_error);
      cout << endl;
    }

    if (flag_location_set) {
      output_gradient_based_scalars
        (cout, point_coord, gradient_coord, scalar, num_gradients,
         location);
    }

    if (flag_list_eigen) {
      output_cube_eigenvalues
        (cout, scalar_grid, gradient_grid, isovalue,
         max_small_mag, max_small_eigenvalue);
    }
  }
  catch (IJK::ERROR error) {
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
// OUTPUT ROUTINES
// **************************************************

void output_cube_coordinates
(std::ostream & output,
 const SHARPISO_SCALAR_GRID & scalar_grid,
 const VERTEX_INDEX icube)
{
  GRID_COORD_TYPE coord[DIM3];

  scalar_grid.ComputeCoord(icube, coord);

  output << "Cube " << icube << ".";
  output << "  Coordinates: ";
  IJK::ijkgrid_output_coord(output, DIM3, coord);
  output << ".";
}

void output_gradients
(std::ostream & output,
 const COORD_TYPE * point_coord,
 const GRADIENT_COORD_TYPE * gradient_coord,
 const SCALAR_TYPE * scalar,
 const NUM_TYPE num_points)
{
  using namespace std;

  for (NUM_TYPE i = 0; i < num_points; i++) {
    output << "Point " << i << " ";
    IJK::ijkgrid_output_coord(output, DIM3, point_coord + i*DIM3);
    output << ":  Scalar " << scalar[i];
    output << " Gradient ";
    IJK::ijkgrid_output_coord(output, DIM3, gradient_coord + i*DIM3);
    GRADIENT_COORD_TYPE magnitude;
    IJK::compute_magnitude(DIM3, gradient_coord + i*DIM3, magnitude);
    output << "  Magnitude " << magnitude;
    output << endl;
  }
}

void output_gradients
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_points)
{
  output_gradients(output, &(point_coord[0]), &(gradient_coord[0]),
                   &(scalar[0]), num_points);

}

void output_svd_results
(std::ostream & output,
 const COORD_TYPE sharp_coord[DIM3],
 const EIGENVALUE_TYPE eigenvalues[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const EIGENVALUE_TYPE eigenvalue_tolerance,
 SVD_INFO & svd_info)
{
  output << "Sharp coordinates: ";
  IJK::ijkgrid_output_coord(output, DIM3, sharp_coord);
  output << endl;
  output << "Eigenvalues: ";
  IJK::ijkgrid_output_coord(output, DIM3, eigenvalues);
  output << endl;
  output << "Number of large eigenvalues (>= " << eigenvalue_tolerance
         << "): "
         << num_large_eigenvalues << endl;

  if(num_large_eigenvalues == 2) {
    output << "ray direction " << svd_info.ray_direction[0] << " "
           << svd_info.ray_direction[1] << " "
           << svd_info.ray_direction[2] << endl;
    output << "ray initial point "
           << svd_info.ray_initial_point[0] << " "
           << svd_info.ray_initial_point[1] << " "
           << svd_info.ray_initial_point[2] << endl;
    output << "ray intersected cube? "
           << svd_info.ray_intersect_cube << endl;
  }
  output << "used none[0]/centroid[1]/cube-center[2]: "
         << svd_info.location << endl;
}

void output_subgrid_results
(std::ostream & output,
 const COORD_TYPE sharp_coord[DIM3],
 const SCALAR_TYPE scalar_stdev,
 const SCALAR_TYPE max_abs_scalar_error)
{
  output << "Sharp coordinates: ";
  IJK::ijkgrid_output_coord(output, DIM3, sharp_coord);
  output << endl;
  output << "Scalar standard deviation: " << scalar_stdev << endl;
  output << "Max absolute scalar error: " << max_abs_scalar_error << endl;
}


void output_cube_subgrid_scalar_errors
(std::ostream & output,
 const COORD_TYPE * point_coord, const GRADIENT_COORD_TYPE * gradient_coord,
 const SCALAR_TYPE * scalar, const NUM_TYPE num_points,
 const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
 const NUM_TYPE subgrid_axis_size)
{
  COORD_TYPE coord[DIM3];
  COORD_TYPE center_coord[DIM3];
  IJK::PROCEDURE_ERROR error("output_cube_subgrid_scalar_errors");

  if (subgrid_axis_size < 1) {
    error.AddMessage
      ("Programming error. Subgrid axis size must be at least 1.");
    error.AddMessage("  Subgrid axis size = ", subgrid_axis_size, ".");
    throw error;
  }

  // Compute center coordinate
  for (NUM_TYPE d = 0; d < DIM3; d++)
    { center_coord[d] = cube_coord[d] + 0.5; }

  const COORD_TYPE h = 1.0/(subgrid_axis_size+1);

  for (NUM_TYPE ix = 0; ix < subgrid_axis_size; ix++) {
    coord[0] = cube_coord[0] + (ix+1)*h;
    for (NUM_TYPE iy = 0; iy < subgrid_axis_size; iy++) {
      coord[1] = cube_coord[1] + (iy+1)*h;
      for (NUM_TYPE iz = 0; iz < subgrid_axis_size; iz++) {

        coord[2] = cube_coord[2] + (iz+1)*h;
        SCALAR_TYPE s, stdev_squared, max_abs_error;

        compute_gradient_based_scalar_diff
          (coord, isovalue, point_coord, gradient_coord, scalar, num_points,
           stdev_squared, max_abs_error);

        output << "Coord: ";
        IJK::ijkgrid_output_coord(cerr, DIM3, coord);
        output << "  stdev: " << std::sqrt(stdev_squared);
        output << "  max abs error: " << max_abs_error;
        COORD_TYPE dist2center_squared;
        IJK::compute_distance_squared
          (DIM3, coord, center_coord, dist2center_squared);
        cerr << "  dist to center: " << std::sqrt(dist2center_squared);
        cerr << endl;
      }
    }
  }

}

void output_cube_subgrid_scalar_errors
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar, const NUM_TYPE num_points,
 const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
 const NUM_TYPE subgrid_axis_size)
{
  output_cube_subgrid_scalar_errors
    (output, &(point_coord[0]), &(gradient_coord[0]),
     &(scalar[0]), num_points, cube_coord, isovalue, subgrid_axis_size);
}

void output_gradient_based_scalars
(std::ostream & output,
 const COORD_TYPE * point_coord,
 const GRADIENT_COORD_TYPE * gradient_coord,
 const SCALAR_TYPE * scalar, const NUM_TYPE num_points,
 const COORD_TYPE * location)
{
  output << "Location: ";
  IJK::ijkgrid_output_coord(output, DIM3, location);
  output << endl;

  for (int i = 0; i < num_points; i++) {
    output << "  Point ";
    IJK::ijkgrid_output_coord(output, DIM3, point_coord+i*DIM3);
    output << ".  Scalar " << scalar[i];
    output << ".  Gradient ";
    IJK::ijkgrid_output_coord(output, DIM3, gradient_coord+i*DIM3);
    output << ".";

    SCALAR_TYPE s =
      compute_gradient_based_scalar
      (location, point_coord+i*DIM3, gradient_coord+i*DIM3, scalar[i]);

    output << "  Location scalar: " << s << "." << endl;
  }
}

void output_gradient_based_scalars
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar, const NUM_TYPE num_points,
 const std::vector<COORD_TYPE> & location)
{
  output_gradient_based_scalars
    (output, &(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
     num_points, &(location[0]));
}

void output_cube_eigenvalues
(std::ostream & output,
 const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE max_zero_mag,
 const EIGENVALUE_TYPE eigenvalue_tolerance)
{
  COORD_TYPE sharp_coord[DIM3];
  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
  NUM_TYPE num_large_eigenvalues(0);
  SVD_INFO svd_info;

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {

    svd_compute_sharp_vertex_in_cube
      (scalar_grid, gradient_grid, icube, isovalue,
       max_zero_mag, eigenvalue_tolerance, sharp_coord,
       eigenvalues, num_large_eigenvalues,
       svd_info);

    output_cube_coordinates(output, scalar_grid, icube);

    output << " Num eigen: " << num_large_eigenvalues
           << " Eigen: ";
    IJK::ijkgrid_output_coord(output, DIM3, eigenvalues);
    output << endl;
  }
}

// **************************************************
// CHECK ROUTINES
// **************************************************

bool check_input_grids
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid, IJK::ERROR & error)
{
  if (!check_gradient_grid(gradient_grid, error))
    { return(false); }

  IJK::ERROR size_error;
  if (!gradient_grid.Check
      (scalar_grid, "Gradient grid", "Scalar grid", size_error)) {
    error.AddMessage
      ("Scalar grid (file ", scalar_filename,
       ") and gradient grid (file ", gradient_filename, ") do not match.");
    for (int i = 0; i < size_error.NumMessages(); i++)
      { error.AddMessage(size_error.Message(i)); }

    return(false);
  }

  return(true);
}

bool check_gradient_grid
(const GRADIENT_GRID & gradient_grid, IJK::ERROR & error)
{
  if (gradient_grid.VectorLength() != gradient_grid.Dimension()) {
    error.AddMessage
      ("Error in gradient grid (file ", gradient_filename, ").");
    error.AddMessage
      ("  Vector length of gradient grid should equal volume dimension.");
    error.AddMessage
      ("  Volume dimension = ", gradient_grid.Dimension(), ".");
    error.AddMessage
      ("  Vector length = ", gradient_grid.VectorLength(), ".");

    return(false);
  }

  return(true);
}

// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_error()
{
  cerr << "Usage: sharpinfo [OPTIONS] <scalar filename> <gradient filename>"
       << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -isovalue <isovalue> | -cube <cube_index>" << endl;
  cerr << "  -neighbor <cube_index> | -coord \"point coord\"" << endl;
  cerr << "  -svd_grad | -svd_edge_simple | -svd_edge_cmplx"<<endl;
  cerr << "  -listg | -list_subgrid" << endl;
  exit(10);
}

template <typename CTYPE>
void get_coord(const char * s, vector<CTYPE> & coord)
{
  istringstream coord_string;

  coord.clear();

  string s2 = s;
  // remove trailing blanks from s2
  size_t pos = 0;
  for (size_t i = 0; i < s2.length(); i++) {
    if (!isspace(s2[i])) { pos = i+1; }
  }
  if (pos < s2.length()) { s2.erase(pos); };

  coord_string.str(s2);
  while (coord_string.good()) {
    CTYPE c;
    coord_string >> c;
    coord.push_back(c);
  }

  if (coord_string.fail() && !coord_string.eof()) {
    cerr << "Error reading coordinates: "
         << "\"" << s << "\"" << endl;
    cerr << "  Non-numeric character in coordinate string." << endl;
    exit(600);
  }

}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  bool flag_none = false;
  while (iarg < argc && argv[iarg][0] == '-') {

    std::string s = argv[iarg];
    if (s == "-cube") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &cube_index);
    }
    else if (s == "-isovalue") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%f", &isovalue);
      flag_isovalue_set = true;
    }
    else if (s == "-coord") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      get_coord(argv[iarg], location);
      flag_location_set = true;
    }
    else if (s == "-listg") {
      flag_list_gradients = true;
    }
    else if (s == "-list_subgrid") {
      flag_list_subgrid = true;
    }
    else if (s == "-list_eigen") {
      flag_list_eigen = true;
    }
    else if(s == "-svd_grad")
      {
        flag_svd_gradients = true;
      }
    else if(s == "-svd_edge_simple")
      {
        flag_svd_edges_simple = true;
      }
    else if(s == "-svd_edge_cmplx")
      {
        flag_svd_edges_cmplx = true;
      }
    else if (s == "-neighbor") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &cube_index);
      flag_use_neighboring_facets = true;
      flag_use_selective_neighbors = true;
    }
    else if (s == "-neighborL") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%d", &cube_index);
      flag_use_neighboring_facets = true;
      flag_use_selective_neighbors = false;
    }
    else if (s == "-offset") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      sscanf(argv[iarg], "%f", &cube_offset);
    }
    else if (s == "-help") {
      help();
    }
    else {
      cerr << "Option error. Unknown option: " << argv[iarg] << endl;
      cerr << endl;
      usage_error();
    }

    iarg++;
  }

  if (iarg >= argc) {
    cerr << "Error. Missing scalar and gradient file names." << endl;
    usage_error();
  }

  if (iarg + 1 >= argc) {
    cerr << "Error. Missing gradient file name." << endl;
    usage_error();
  }

  if (iarg + 2 < argc) {
    cerr << "Error. Command line has more than two input file names." << endl;
    usage_error();
  }

  scalar_filename = argv[iarg];
  gradient_filename = argv[iarg+1];

  if (!flag_isovalue_set && !flag_location_set && !flag_list_gradients) {
    cerr << "Error.  Option  -isovalue or -coord or -listg must be specified."
         << endl;
    usage_error();
    exit(15);
  }

  if (!flag_isovalue_set && flag_list_subgrid) {
    cerr << "Error.  Option -list_subgrid cannot be used without -isovalue."
         << endl;
    exit(15);
  }

  if (cube_offset <= -1 || cube_offset > 1) {
    cerr << "Error in option -cube_offset." << endl;
    cerr << "  Cube offset must be greater than -1 and at most 1." << endl;
    exit(15);
  }
}

void help()
{
  cerr << "Usage: sharpinfo [OPTIONS] <scalar filename> <gradient filename>"
       << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -isovalue <isovalue>:  Compute isosurface vertex for given <isovalue>." << endl;
  cerr << "  -cube <cube_index>:  Compute isosurface vertex for cube <cube_index>." << endl;
  cerr << "      Default is cube 0." << endl;
  cerr << "  -neighbor <cube_index>:  Use gradients from cube and" << endl
       << "             neighbors of cube <cube_index>." << endl;
  cerr << "  -svd_grad: Compute using svd directly on gradients."<<endl;
  cerr << "  -svd_edge_simple: Interpolate intersection points/normals and"
       << endl
       << "                   apply svd." << endl;
  cerr << "  -svd_edge_cmplx:  Compute edge-isosurface intersection points/normals"
       << endl
       << "                   using gradient assignment and apply svd." << endl;
  cerr << "  -coord \"point_coord\":  Compute scalar values at coordinate point_coord." << endl;
  cerr << "  -listg: List gradients." << endl;
  cerr << "  -list_subgrid:  List all scalar values at vertices of subgrid." << endl;
  exit(15);
}

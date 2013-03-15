/// \file sharpinfo.cxx
/// Output information about sharp isosurface vertices and edges.
/// Version 0.1.0

/*
 IJK: Isosurface Jeneration Kode
 Copyright (C) 2011-2013 Rephael Wenger

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


#include <iomanip>
#include <iostream>

#include "ijkcoord.txx"
#include "ijkgrid_macros.h"
#include "ijkgrid_nrrd.txx"
#include "ijkinterpolate.txx"
#include "ijkprint.txx"
#include "ijkstring.txx"

#include "sharpiso_feature.h"
#include "sharpiso_intersect.h"
#include "sharpiso_eigen.h"
#include "sharpiso_types.h"
#include "sharpiso_scalar.txx"
#include "sharpiso_findIntersect.h"
#include "sh_point_find.h"

using namespace SHARPISO;
using namespace IJK;

using namespace std;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
SCALAR_TYPE isovalue(0);
VERTEX_INDEX cube_index(0);
VERTEX_INDEX facet_v0(0);
int facet_orth_dir(0);
VERTEX_INDEX vertex_index(0);
std::vector<GRID_COORD_TYPE> cube_coord;
std::vector<GRID_COORD_TYPE> vertex_coord;
AXIS_SIZE_TYPE subgrid_axis_size(3);
bool flag_list_gradients(false);
bool flag_list_subgrid(false);
bool flag_list_edgeI(false);
std::vector<COORD_TYPE> location;
bool flag_location_set(false);
bool flag_isovalue_set(false);
bool flag_centroid(false);
SHARP_ISOVERT_PARAM sharpiso_param;
bool flag_list_eigen(false);
bool flag_svd_gradients(true); //default
bool flag_edge_intersect_interpolate(false);
bool flag_edge_intersect_sharp(false);
bool flag_dist2vert(false);
bool flag_cube_set(false);
bool flag_vertex_set(false);
bool flag_use_lindstrom(false);
bool flag_use_lindstrom2(false);
bool flag_use_lindstrom_fast(false);
bool flag_edge_intersection(false);
bool flag_subgrid(false);
bool flag_output_param(true);
bool flag_facet_set(false);

// compute isosurface vertices
void compute_iso_vertex_using_svd
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_nonzero_eigenvalues,
 SVD_INFO & svd_info);
void compute_iso_vertex_using_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GET_GRADIENTS_PARAM & get_gradients_param,
 const NUM_TYPE subgrid_axis_size,
 COORD_TYPE sharp_coord[DIM3],
 SCALAR_TYPE & scalar_stdev,
 SCALAR_TYPE & max_abs_scalar_error);


// check routines
bool check_gradient_grid
(const GRADIENT_GRID & gradient_grid, IJK::ERROR & error);
bool check_input_grids
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid, IJK::ERROR & error);
bool check_cube_coord
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const std::vector<GRID_COORD_TYPE> & cube_coord, IJK::ERROR & error);
bool check_vertex_coord
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const std::vector<GRID_COORD_TYPE> & vertex_coord, IJK::ERROR & error);
bool check_facet
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const VERTEX_INDEX facet_v0, const NUM_TYPE facet_orth_dir,
 IJK::ERROR & error);

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
 const NUM_TYPE num_points,
 const SCALAR_TYPE isovalue,
 const COORD_TYPE cube_center[DIM3]);
void output_svd_results
(std::ostream & output, 
 const SHARPISO_GRID & grid,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE sharp_coord[DIM3],
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
 const SHARP_ISOVERT_PARAM & sharpiso_param);
void output_edge_intersections
(std::ostream & output,
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE zero_tolerance);
void output_facet_svd_results
(std::ostream & output,
 const GRID_COORD_TYPE facet_coord[DIM3],
 const NUM_TYPE facet_orth_dir,
 const NUM_TYPE sharp_vertex_location,
 const COORD_TYPE sharp_coord[DIM3],
 const EIGENVALUE_TYPE eigenvalues[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const EIGENVALUE_TYPE eigenvalue_tolerance);

void output_centroid_results
(std::ostream & output, const COORD_TYPE coord[DIM3]);
void output_dist2vert
(std::ostream & output,
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1);
void report_sharp_param(const SHARP_ISOVERT_PARAM & sharp_param);

// local subroutines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  IJK::ERROR error;
  // accumulate information from svd_calls
  SVD_INFO svd_info;
  svd_info.location = LOC_NONE;

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    SHARPISO_SCALAR_GRID scalar_grid;
    IJK::GRID_NRRD_IN<NUM_TYPE,AXIS_SIZE_TYPE> nrrd_in_scalar;
    nrrd_in_scalar.ReadScalarGrid
      (scalar_filename, scalar_grid,  error);
    if (nrrd_in_scalar.ReadFailed()) { throw(error); }

    GRADIENT_GRID gradient_grid;
    IJK::GRID_NRRD_IN<NUM_TYPE,AXIS_SIZE_TYPE> nrrd_in_gradient;

    nrrd_in_gradient.ReadVectorGrid
      (gradient_filename, gradient_grid, error);
    if (nrrd_in_gradient.ReadFailed()) { throw(error); }

    if (!check_gradient_grid(gradient_grid, error))
      { throw error; };

    if (!check_input_grids(scalar_grid, gradient_grid, error))
      { throw error; }

    if (cube_coord.size() > 0) {
      if (!check_cube_coord(scalar_grid, cube_coord, error))
        { throw error; }

      cube_index = scalar_grid.ComputeVertexIndex(cube_coord);
      flag_cube_set = true;
    }

    if (vertex_coord.size() > 0) {
      if (!check_vertex_coord(scalar_grid, vertex_coord, error))
        { throw error; }

      vertex_index = scalar_grid.ComputeVertexIndex(vertex_coord);
    }

    if (flag_facet_set) {
      if (!check_facet(scalar_grid, facet_v0, facet_orth_dir, error))
        { throw error; }
    }

    NUM_TYPE num_gradients = 0;
    std::vector<COORD_TYPE> point_coord;
    std::vector<GRADIENT_COORD_TYPE> gradient_coord;
    std::vector<SCALAR_TYPE> scalar;

    if (flag_list_gradients || flag_isovalue_set || flag_location_set) {
      GRADIENT_COORD_TYPE max_small_mag =
        sharpiso_param.max_small_magnitude;
      GRADIENT_COORD_TYPE max_small_eigenvalue =
        sharpiso_param.max_small_eigenvalue;

      OFFSET_CUBE_111 cube_111
        (sharpiso_param.grad_selection_cube_offset);

      if (flag_edge_intersect_sharp) {
        get_edgeI_sharp_gradients
          (scalar_grid, gradient_grid, cube_index, isovalue,
           point_coord, gradient_coord, scalar, num_gradients);
      }
      else {
        get_gradients
          (scalar_grid, gradient_grid, cube_index, isovalue,
           sharpiso_param, cube_111, sharpiso_param.flag_sort_gradients,
           point_coord, gradient_coord, scalar, num_gradients);
      }

      if (flag_list_gradients) {
        COORD_TYPE cube_center[DIM3];
  
        scalar_grid.ComputeCoord(cube_index, cube_center);

        for (int d = 0; d < DIM3; d++)
          { cube_center[d] += 0.5; }

        output_gradients
          (cout, point_coord, gradient_coord, scalar, num_gradients, 
           isovalue, cube_center);
        cout << endl;

        if (flag_facet_set) {
          std::vector<COORD_TYPE> point_coord2;
          std::vector<GRADIENT_COORD_TYPE> gradient_coord2;
          std::vector<SCALAR_TYPE> scalar2;
          int num_gradients2;

          get_two_cube_gradients
            (scalar_grid, gradient_grid, facet_v0, facet_orth_dir, 
             isovalue, sharpiso_param,
             point_coord2, gradient_coord2, scalar2, num_gradients2);

          cout << "Gradients around facet:" << endl;
          output_gradients
          (cout, point_coord2, gradient_coord2, scalar2, num_gradients2, 
           isovalue, cube_center);
          cout << endl;
        }

      }

      if (flag_isovalue_set) {

        COORD_TYPE sharp_coord[DIM3];
        EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
        NUM_TYPE num_large_eigenvalues(0);

        output_cube_coordinates(cout, scalar_grid, cube_index);
        cout << endl << endl;

        if (flag_centroid) {
          COORD_TYPE coord[DIM3];
          compute_edgeI_centroid
            (scalar_grid, gradient_grid, isovalue, cube_index, 
             sharpiso_param.use_sharp_edgeI, coord);

          output_centroid_results(cout, coord);
          cout << endl;
        }
        else {

          if (flag_facet_set) {
            GRID_COORD_TYPE facet_coord[DIM3];
            GRADIENT_COORD_TYPE line_direction[DIM3];
            NUM_TYPE sharp_vertex_location;

            svd_compute_sharp_vertex_near_facet
              (scalar_grid, gradient_grid, facet_v0, facet_orth_dir,
               isovalue, sharpiso_param, 
               sharp_vertex_location, sharp_coord, line_direction, 
               eigenvalues, num_large_eigenvalues);

            scalar_grid.ComputeCoord(facet_v0, facet_coord);

            output_facet_svd_results
              (cout, facet_coord, facet_orth_dir, sharp_vertex_location,
               sharp_coord, eigenvalues, num_large_eigenvalues,
               max_small_eigenvalue);
            cout << endl;
          }

          if (flag_cube_set) {
            GRID_COORD_TYPE cube_coord[DIM3];

            compute_iso_vertex_using_svd
              (scalar_grid, gradient_grid, cube_index, isovalue,
               sharp_coord, eigenvalues,
               num_large_eigenvalues, svd_info);

            scalar_grid.ComputeCoord(cube_index, cube_coord);

            output_svd_results
              (cout, scalar_grid, cube_coord, sharp_coord, 
               eigenvalues, num_large_eigenvalues,
               max_small_eigenvalue, svd_info);
            cout << endl;
          }

          if (flag_output_param) {
            report_sharp_param(sharpiso_param);
            cout << endl;
          }
        }

        if (flag_list_edgeI) {
          output_edge_intersections
            (cout, scalar_grid, gradient_grid, cube_index, isovalue,
             sharpiso_param.zero_tolerance);
          cout << endl;
        }

        if (flag_list_subgrid) {
          GRID_COORD_TYPE cube_coord[DIM3];

          scalar_grid.ComputeCoord(cube_index, cube_coord);
          output_cube_subgrid_scalar_errors
            (cout, point_coord, gradient_coord, scalar, num_gradients,
             cube_coord, isovalue, subgrid_axis_size);
          cout << endl;
        }

        if (flag_subgrid || flag_list_subgrid) {

          SCALAR_TYPE scalar_stdev;
          SCALAR_TYPE max_abs_scalar_error;
          compute_iso_vertex_using_subgrid
            (scalar_grid, gradient_grid, cube_index, isovalue,
             sharpiso_param, subgrid_axis_size,
             sharp_coord, scalar_stdev, max_abs_scalar_error);

          output_subgrid_results
            (cout, sharp_coord, scalar_stdev, max_abs_scalar_error);
          cout << endl;
        }
      }


      if (flag_location_set) {
        output_gradient_based_scalars
          (cout, point_coord, gradient_coord, scalar, num_gradients,
           location);
      }

    }

    if (flag_list_eigen) {

      output_cube_eigenvalues
        (cout, scalar_grid, gradient_grid, isovalue, sharpiso_param);

    }

    if (flag_dist2vert) {
      output_dist2vert(cout, scalar_grid, gradient_grid, vertex_index);
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
// COMPUTE ISOSURFACE VERTEX
// **************************************************

void compute_iso_vertex_using_svd
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 COORD_TYPE sharp_coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
  SIGNED_COORD_TYPE grad_selection_cube_offset =
    sharpiso_param.grad_selection_cube_offset;
  GRADIENT_COORD_TYPE max_small_mag =
    sharpiso_param.max_small_magnitude;
  GRADIENT_COORD_TYPE max_small_eigenvalue =
    sharpiso_param.max_small_eigenvalue;


  if (flag_edge_intersect_interpolate) {
    svd_compute_sharp_vertex_edgeI_interpolate_gradients
      (scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
       sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
  }
  else if (flag_edge_intersect_sharp) {
    svd_compute_sharp_vertex_edgeI_sharp_gradient
      (scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
       sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
  }
  else {

    OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

    if (sharpiso_param.use_lindstrom) {
      svd_compute_sharp_vertex_for_cube_lindstrom
        (scalar_grid, gradient_grid, cube_index, isovalue,
         sharpiso_param, cube_111,
         sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
    }
    else {
      svd_compute_sharp_vertex_for_cube_lc_intersection
        (scalar_grid, gradient_grid, cube_index, isovalue,
         sharpiso_param, cube_111,
         sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
    }
  }

}


void compute_iso_vertex_using_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GET_GRADIENTS_PARAM & get_gradients_param,
 const NUM_TYPE subgrid_axis_size,
 COORD_TYPE sharp_coord[DIM3],
 SCALAR_TYPE & scalar_stdev,
 SCALAR_TYPE & max_abs_scalar_error)
{
  SIGNED_COORD_TYPE grad_selection_cube_offset =
    sharpiso_param.grad_selection_cube_offset;

  OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

  subgrid_compute_sharp_vertex_in_cube
    (scalar_grid, gradient_grid, cube_index, isovalue,
     get_gradients_param, cube_111, subgrid_axis_size,
     sharp_coord, scalar_stdev, max_abs_scalar_error);
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
 const NUM_TYPE num_points,
 const SCALAR_TYPE isovalue,
 const COORD_TYPE cube_center[DIM3])
{
  using namespace std;

  for (NUM_TYPE i = 0; i < num_points; i++) {
    GRADIENT_COORD_TYPE magnitude;
    IJK::compute_magnitude(DIM3, gradient_coord + i*DIM3, magnitude);

    COORD_TYPE distance;
    compute_signed_distance_to_gfield_plane
      (gradient_coord+i*DIM3, point_coord+i*DIM3, scalar[i], cube_center,
       isovalue, distance);

    output << "Point " << setw(2)  << i << " ";
    IJK::ijkgrid_output_coord(output, DIM3, point_coord + i*DIM3);
    output << ":  Scalar " << scalar[i];
    output << " Grad ";
    IJK::ijkgrid_output_coord(output, DIM3, gradient_coord + i*DIM3);
    output << "  Mag " << magnitude;
    output << "  Dist " << distance;
    output << endl;
  }
}

void output_gradients
(std::ostream & output,
 const std::vector<COORD_TYPE> & point_coord,
 const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
 const std::vector<SCALAR_TYPE> & scalar,
 const NUM_TYPE num_points,
 const SCALAR_TYPE isovalue,
 const COORD_TYPE cube_center[DIM3])
{
  output_gradients(output, &(point_coord[0]), &(gradient_coord[0]),
                   &(scalar[0]), num_points, isovalue, cube_center);
 
}

void output_svd_results
(std::ostream & output,
 const SHARPISO_GRID & grid,
 const GRID_COORD_TYPE cube_coord[DIM3],
 const COORD_TYPE sharp_coord[DIM3],
 const EIGENVALUE_TYPE eigenvalues[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const EIGENVALUE_TYPE eigenvalue_tolerance,
 SVD_INFO & svd_info)
{
  COORD_TYPE closest_point[DIM3];

  bool use_only_cube_gradients = sharpiso_param.use_only_cube_gradients;
  bool use_selected_gradients = sharpiso_param.use_selected_gradients;

  
  output << "SVD: Sharp coordinates ";
  if (svd_info.location == LOC_SVD)
    { output << "(svd): "; }
  else if (svd_info.location == CENTROID)
    { output << "(centroid): "; }
  else if (svd_info.location == CUBE_CENTER)
    { output << "(cube center): "; }
  else
    { output << "(unknown): "; }

  IJK::ijkgrid_output_coord(output, DIM3, sharp_coord);
  output << endl;

  output << "Central point for distance: ";
  IJK::ijkgrid_output_coord(output, DIM3, svd_info.central_point);
  output << endl;

  output << "Eigenvalues: ";
  IJK::ijkgrid_output_coord(output, DIM3, eigenvalues);
  output << endl;

  if (eigenvalues[0] > eigenvalue_tolerance) {
    EIGENVALUE_TYPE normalized_eigenvalues[DIM3];
    multiply_coord_3D(1/eigenvalues[0], eigenvalues, normalized_eigenvalues);
    output << "Normalized eigenvalues: ";
    IJK::ijkgrid_output_coord(output, DIM3, normalized_eigenvalues);
    output << endl;
  }

  output << "Number of large eigenvalues (>= " << eigenvalue_tolerance
         << "): "
         << num_large_eigenvalues << endl;

  if(num_large_eigenvalues == 2 && !flag_use_lindstrom) {
    output << "Ray: ";
    print_coord3D(output, svd_info.ray_initial_point);
    output << " + t";
    print_coord3D(output, svd_info.ray_direction);
    output << endl;

    compute_closest_point_on_line
      (svd_info.central_point, 
       svd_info.ray_initial_point, svd_info.ray_direction,
       sharpiso_param.zero_tolerance, closest_point);
    output <<"  Closest (L2) point on ray to central point: ";
    print_coord3D(output, closest_point);
    output << endl;

    compute_closest_point_on_line_linf
      (svd_info.central_point, 
       svd_info.ray_initial_point, svd_info.ray_direction,
       sharpiso_param.zero_tolerance, closest_point);
    output <<"  Closest (Linf) point on ray to central point: ";
    print_coord3D(output, closest_point);
    output << endl;

    compute_closest_point_on_line
      (svd_info.central_point, 
       svd_info.ray_initial_point, svd_info.ray_direction,
       sharpiso_param.zero_tolerance, closest_point);
  }

  if (svd_info.flag_conflict) {
    const VERTEX_INDEX cube2 = svd_info.cube_containing_coord;
    COORD_TYPE cube2_coord[DIM3];

    grid.ComputeCoord(cube2, cube2_coord);
    output << "Conflict with cube: ";
    print_coord3D(output, cube2_coord);
    output << " (Index: " << cube2 << ")" << endl;
  }

}

void output_facet_svd_results
(std::ostream & output,
 const GRID_COORD_TYPE facet_coord[DIM3],
 const NUM_TYPE facet_orth_dir,
 const NUM_TYPE sharp_vertex_location,
 const COORD_TYPE sharp_coord[DIM3],
 const EIGENVALUE_TYPE eigenvalues[DIM3],
 const NUM_TYPE num_large_eigenvalues,
 const EIGENVALUE_TYPE eigenvalue_tolerance)
{
  COORD_TYPE closest_point[DIM3];

  bool use_only_cube_gradients = sharpiso_param.use_only_cube_gradients;
  bool use_selected_gradients = sharpiso_param.use_selected_gradients;


  output << "Facet: ";
  IJK::ijkgrid_output_coord(output, DIM3, facet_coord);
  output << " Orth dir: " << facet_orth_dir;
  output << endl;

  if (sharp_vertex_location < 0) 
    { output << "Sharp feature left/below facet." << endl; }
  else if (sharp_vertex_location > 0) 
    { output << "Sharp feature right/above facet." << endl; }
  else
    { output << "Sharp feature intersects facet." << endl; }

  output << "SVD: Sharp coordinates ";

  IJK::ijkgrid_output_coord(output, DIM3, sharp_coord);
  output << endl;
  output << "Eigenvalues: ";
  IJK::ijkgrid_output_coord(output, DIM3, eigenvalues);
  output << endl;

  if (eigenvalues[0] > eigenvalue_tolerance) {
    EIGENVALUE_TYPE normalized_eigenvalues[DIM3];
    multiply_coord_3D(1/eigenvalues[0], eigenvalues, normalized_eigenvalues);
    output << "Normalized eigenvalues: ";
    IJK::ijkgrid_output_coord(output, DIM3, normalized_eigenvalues);
    output << endl;
  }

  output << "Number of large eigenvalues (>= " << eigenvalue_tolerance
         << "): "
         << num_large_eigenvalues << endl;
}

void output_centroid_results
(std::ostream & output, const COORD_TYPE coord[DIM3])
{
  output << "Centroid: coordinates ";
  print_coord3D(output, coord);
  output << endl;
}

void output_subgrid_results
(std::ostream & output,
 const COORD_TYPE sharp_coord[DIM3],
 const SCALAR_TYPE scalar_stdev,
 const SCALAR_TYPE max_abs_scalar_error)
{
    output << "SUBGRID: Sharp coordinates: ";
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
    COORD_TYPE vdiff[DIM3];
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
                COORD_TYPE distance;
                IJK::compute_distance(DIM3, coord, center_coord, distance);
                cerr << "  dist to center: " << distance << endl;
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

void output_edge_intersections
(std::ostream & output,
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const GRADIENT_COORD_TYPE zero_tolerance)
{
  const int dimension = scalar_grid.Dimension();
  COORD_TYPE vcoord0[DIM3];
  COORD_TYPE vcoord1[DIM3];
  COORD_TYPE pcoord[DIM3];
  GRADIENT_COORD_TYPE grad[DIM3];
  SCALAR_TYPE s, s_split, t_split;
  bool flag_no_split;

  for (int d = 0; d < dimension; d++) {
    for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
      VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
      VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);

      if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        scalar_grid.ComputeCoord(iv0, vcoord0);
        scalar_grid.ComputeCoord(iv1, vcoord1);

        output << "Edge: ";
        IJK::ijkgrid_output_coord(output, DIM3, vcoord0);
        IJK::ijkgrid_output_coord(output, DIM3, vcoord1);

        if (sharpiso_param.use_sharp_edgeI) {

          compute_isosurface_grid_edge_intersection
            (scalar_grid, gradient_grid, isovalue,
             iv0, iv1, d, pcoord);

          output << " Intersection point: ";
          IJK::ijkgrid_output_coord(output, DIM3, pcoord);

          cout << endl << "     ";

          compute_gradient_change_on_edge
            (scalar_grid, gradient_grid, iv0, iv1, d,
             zero_tolerance, flag_no_split, t_split, s_split);

          if (flag_no_split) {
            output << "  No split." << endl;
          }
          else {
            output << "  split scalar: " << s_split;
            output << "  split t: " << t_split << endl;
          }
        }
        else {
          IJK::linear_interpolate_coord
            (DIM3, s0, vcoord0, s1, vcoord1, isovalue, pcoord);

          output << " Intersection point: ";
          IJK::ijkgrid_output_coord(output, DIM3, pcoord);

          output << endl;
        }
      }
    }
  }
}


/// Output distance from vert to plane defined by gradient field at neighbors.
void output_dist2vert
(std::ostream & output, const GRADIENT_COORD_TYPE gfield_gradient[DIM3],
 const GRID_COORD_TYPE gfield_point[DIM3],
 const SCALAR_TYPE gfield_point_scalar,
 const GRID_COORD_TYPE point[DIM3],
 const SCALAR_TYPE plane_scalar)
{
  GRADIENT_COORD_TYPE distance;
  GRADIENT_COORD_TYPE gradient_magnitude;

  IJK::compute_magnitude(DIM3, gfield_gradient, gradient_magnitude);

  // Skip zero magnitude gradients.
  if (gradient_magnitude <= 0) { return; }

  compute_signed_distance_to_gfield_plane
    (gfield_gradient, gfield_point, gfield_point_scalar, point, plane_scalar,
     distance);

  output << "  Point ";
  IJK::ijkgrid_output_coord(output, DIM3, gfield_point);
  output << ".  Scalar " << gfield_point_scalar;
  output << ".  Gradient ";
  IJK::ijkgrid_output_coord(output, DIM3, gfield_gradient);
  output << ".";

  output << "  Distance plane to vertex: " << distance << "." << endl;
}

/// Output distance from vert to plane defined by gradient field at neighbors.
void output_dist2vert
(std::ostream & output,
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const VERTEX_INDEX iv1)
{
  const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);
  GRID_COORD_TYPE coord0[DIM3], coord1[DIM3], coord2[DIM3];

  scalar_grid.ComputeCoord(iv1, coord1);

  output << "Vertex " << iv1 << " coordinates: ";
  IJK::ijkgrid_output_coord(output, DIM3, coord1);
  output << ".  Scalar " << s1 << "." << endl;

  for (int d = 0; d < DIM3; d++) {

    if (coord1[d] > 0) {
      VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
      scalar_grid.ComputeCoord(iv0, coord0);

      output_dist2vert(output, gradient_grid.VectorPtrConst(iv0),
                       coord0, scalar_grid.Scalar(iv0),
                       coord1, s1);
    }

    if (coord1[d] + 1 < scalar_grid.AxisSize(d)) {
      VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
      scalar_grid.ComputeCoord(iv2, coord2);

      output_dist2vert(output, gradient_grid.VectorPtrConst(iv2),
                       coord2, scalar_grid.Scalar(iv2),
                       coord1, s1);
    }

  }
}


void output_cube_eigenvalues
(std::ostream & output,
 const SHARPISO_SCALAR_GRID & scalar_grid,
 const GRADIENT_GRID & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param)
{
  COORD_TYPE sharp_coord[DIM3];
  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
  NUM_TYPE num_large_eigenvalues(0);
  SVD_INFO svd_info;

  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    sharpiso_param.grad_selection_cube_offset;

  OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {

    if (sharpiso_param.use_lindstrom) {
      svd_compute_sharp_vertex_for_cube_lindstrom
        (scalar_grid, gradient_grid, icube, isovalue,
         sharpiso_param, cube_111,
         sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
    }
    else {
      svd_compute_sharp_vertex_for_cube_lc_intersection
        (scalar_grid, gradient_grid, icube, isovalue,
         sharpiso_param, cube_111,
         sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);
    }

    output_cube_coordinates(output, scalar_grid, icube);

    output << " Num eigen: " << num_large_eigenvalues
           << " Eigen: ";
    IJK::ijkgrid_output_coord(output, DIM3, eigenvalues);
    output << endl;
  }
}

void report_sharp_param(const SHARP_ISOVERT_PARAM & sharp_param)
{
  cout << "Gradient selection cube offset: "
       << sharp_param.grad_selection_cube_offset << endl;
  cout << "Maximum small eigenvalue: "
       << sharp_param.max_small_eigenvalue << endl;
  cout << "Max (Linf) distance from cube to isosurface vertex: "
       << sharp_param.max_dist << endl;
  cout << "Max small gradient magnitude: "
       << sharp_param.max_small_magnitude << endl;
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

bool check_cube_coord
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const std::vector<GRID_COORD_TYPE> & cube_coord, IJK::ERROR & error)
{
    if (scalar_grid.Dimension() != cube_coord.size()) {
        error.AddMessage("Error: Option -cc followed by incorrect number of coordinates.");
        error.AddMessage("  Option -cc followed by ", cube_coord.size(), " coordinates.");
        error.AddMessage("  Option -cc should be followed by ", scalar_grid.Dimension(),
                         " coordinates.");
        return(false);
    }

    for (int d = 0; d < scalar_grid.Dimension(); d++) {
        if (cube_coord[d] < 0) {
            error.AddMessage("Error: Illegal negative coordinate following -cc.");
            return(false);
        }

        if (cube_coord[d]+1 >= scalar_grid.AxisSize(d)) {
            error.AddMessage("Error: Coordinate following -cc is not a valid cube coordinate.");
            return(false);
        }
    }

    return(true);
}

bool check_vertex_coord
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const std::vector<GRID_COORD_TYPE> & vertex_coord, IJK::ERROR & error)
{
    if (scalar_grid.Dimension() != vertex_coord.size()) {
        error.AddMessage("Error: Option -vc followed by incorrect number of coordinates.");
        error.AddMessage("  Option -vc followed by ", vertex_coord.size(), " coordinates.");
        error.AddMessage("  Option -vc should be followed by ", scalar_grid.Dimension(),
                         " coordinates.");
        return(false);
    }

    for (int d = 0; d < scalar_grid.Dimension(); d++) {
        if (vertex_coord[d] < 0) {
            error.AddMessage("Error: Illegal negative coordinate following -vc.");
            return(false);
        }

        if (vertex_coord[d] >= scalar_grid.AxisSize(d)) {
            error.AddMessage("Error: Coordinate following -vc is not a valid vertex coordinate.");
            return(false);
        }
    }

    return(true);
}

bool check_facet
(const SHARPISO_SCALAR_GRID & scalar_grid,
 const VERTEX_INDEX facet_v0, const NUM_TYPE facet_orth_dir,
 IJK::ERROR & error)
{
  GRID_COORD_TYPE coord[DIM3];

  if (facet_v0 < 0) {
    error.AddMessage("Illegal facet vertex following -facet.");
    error.AddMessage("Facet vertex index must be non-negative.");
    return(false);
  }

  if (facet_v0+1 >= scalar_grid.NumVertices()) {
    error.AddMessage("Illegal facet vertex following -facet.");
    error.AddMessage("Facet vertex index must be less than ",
                     scalar_grid.NumVertices(), ".");
    return(false);
  }

  scalar_grid.ComputeCoord(facet_v0, coord);

  if (coord[facet_orth_dir] <= 0 ||
      coord[facet_orth_dir]+1 >= scalar_grid.AxisSize(facet_orth_dir)) {

    error.AddMessage
      ("Illegal facet specification.  Facet not in grid interior.");
    return(false);
  }
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
    cerr << "  -isovalue <isovalue> | -cube <cube_index> | -cc \"cube coordinates\""
    << endl;
    cerr << "  [-centroid | -gradC | -gradN | -gradCS | -gradNS |" << endl;
    cerr << "   -gradIE | -gradIES | -gradIEDir | -gradNIE | -gradNIES |" 
         << endl;
    cerr << "   -gradCD | -gradCDdup | -gradES | -gradEC ]" << endl;
    cerr << "  -subgrid | -lindstrom | -rayI" << endl;
    cerr << "  [-allow_conflict | -clamp_conflict | -centroid_conflict]" 
         << endl;
    cerr << "  -clamp_far | [-recompute_eigen2 | -no_recompute_eigen2]" << endl;
    cerr << "  [-sharp_edgeI | -interpolate_edgeI]" << endl;
    cerr << "  [-dist2center | -dist2centroid]" << endl;
    cerr << "  [-removeg | -no_removeg]"
         << "  [-centroid_eigen1 | -no_centroid_eigen1] [-no_Linf]" << endl;
    cerr << "  [-no_round | -round <n>]" << endl;
    cerr << "  -coord \"point coord\"" << endl;
    cerr << "  -facet <vertex_index> <orth_dir>" << endl;
    cerr << "  -dist2vert | -vertex <vertex_index>"
         << " | -vc \"vertex coordinates\"" << endl;
    cerr << "  -max_eigen <value>" << endl;
    cerr << "  -gradS_offset <value> | -max_dist <value> | max_mag <value>" 
         << endl;
    cerr << "  -listg | -list_subgrid | -list_edgeI" << endl;
    cerr << "  -help | -no_output_param" << endl;
    exit(10);
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

int get_two_int(const int iarg, const int argc, char **argv,
                int & x0, int & x1)
{
  const int num_arg = 2;

  if (iarg+num_arg >= argc) {
    cerr << "Usage error. Missing arguments for option "
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!IJK::string2val(argv[iarg+1], x0)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  if (!IJK::string2val(argv[iarg+2], x1)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg+2] << endl;
    cerr << "Note: " << argv[iarg] << " takes two integer arguments."
         << endl;
    exit(51);
  }
}

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

void parse_command_line(int argc, char **argv)
{
  bool use_only_cube_gradients(true);
  bool use_selected_gradients(true);
  bool use_intersected_edge_endpoint_gradients(false);
  bool use_gradients_determining_edge_intersections(false);
  bool select_based_on_grad_dir(false);

  if (argc == 1) { usage_error(); }

  int iarg = 1;
  bool flag_none = false;
  while (iarg < argc && argv[iarg][0] == '-') {

    std::string s = argv[iarg];

    if (s == "-cube") {
      cube_index = get_int(iarg, argc, argv);
      iarg++;
      flag_cube_set = true;
    }
    else if (s == "-facet") {
      get_two_int(iarg, argc, argv, facet_v0, facet_orth_dir);
      flag_facet_set = true;
      iarg = iarg+2;
    }
    else if (s == "-dist2vert") {
      flag_dist2vert = true;
    }
    else if (s == "-vertex") {
      vertex_index = get_int(iarg, argc, argv);
      iarg++;
      flag_vertex_set = true;
    }
    else if (s == "-cc") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      IJK::string2vector(argv[iarg], cube_coord);
    }
    else if (s == "-vc") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      IJK::string2vector(argv[iarg], vertex_coord);
    }
    else if (s == "-isovalue") {
      isovalue = get_float(iarg, argc, argv);
      iarg++;
      flag_isovalue_set = true;
    }
    else if (s == "-coord") {
      iarg++;
      if (iarg >= argc) { usage_error(); };
      IJK::string2vector(argv[iarg], location);
      flag_location_set = true;
    }
    else if (s == "-listg") {
      flag_list_gradients = true;
    }
    else if (s == "-list_subgrid") {
      flag_list_subgrid = true;
    }
    else if (s == "-list_edgeI") {
      flag_list_edgeI = true;
    }
    else if (s == "-list_eigen") {
      flag_list_eigen = true;
    }
    else if (s == "-sortg") {
      sharpiso_param.flag_sort_gradients = true;
    }
    else if (s == "-gradES") {
      flag_edge_intersect_interpolate = true;
      flag_use_lindstrom = true;    // *** CURRENTLY ONLY LINDSTROM
      use_gradients_determining_edge_intersections = true;
      use_selected_gradients = false;
      flag_edge_intersection = true;
    }
    else if (s == "-gradEC") {
      flag_edge_intersect_sharp = true;
      flag_use_lindstrom = true;    // *** CURRENTLY ONLY LINDSTROM
      use_selected_gradients = false;
      use_gradients_determining_edge_intersections = true;
      flag_edge_intersection = true;
      sharpiso_param.use_sharp_edgeI = true;
    }
    else if (s == "-centroid") {
      flag_centroid = true;
    }
    else if (s == "-gradC") {
      use_only_cube_gradients = true;
      use_selected_gradients = false;
    }
    else if (s == "-gradN") {
      use_only_cube_gradients = false;
      use_selected_gradients = false;
      use_intersected_edge_endpoint_gradients = false;
    }
    else if (s == "-gradCS") {
      use_only_cube_gradients = true;
      use_selected_gradients = true;
      use_intersected_edge_endpoint_gradients = false;
    }
    else if (s == "-gradNS") {
      use_only_cube_gradients = false;
      use_selected_gradients = true;
      use_intersected_edge_endpoint_gradients = false;
    }
    else if (s == "-gradIE") {
      use_only_cube_gradients = true;
      use_selected_gradients = false;
      use_intersected_edge_endpoint_gradients = true;
    }
    else if (s == "-gradIES") {
      use_only_cube_gradients = true;
      use_selected_gradients = true;
      use_intersected_edge_endpoint_gradients = true;
    }
    else if (s == "-gradIEDir") {
      use_only_cube_gradients = true;
      use_selected_gradients = true;
      select_based_on_grad_dir = true;
      use_intersected_edge_endpoint_gradients = true;
    }
    else if (s == "-gradNIE") {
      use_only_cube_gradients = false;
      use_selected_gradients = false;
      use_intersected_edge_endpoint_gradients = true;
    }
    else if (s == "-gradNIES") {
      use_only_cube_gradients = false;
      use_selected_gradients = true;
      use_intersected_edge_endpoint_gradients = true;
    }
    else if (s == "-gradCD") {
      use_only_cube_gradients = true;
      use_selected_gradients = false;
      use_intersected_edge_endpoint_gradients = false;
      use_gradients_determining_edge_intersections = true;
      flag_edge_intersection = false;
    }
    else if (s == "-gradCDdup") {
      use_only_cube_gradients = true;
      use_selected_gradients = false;
      use_intersected_edge_endpoint_gradients = false;
      use_gradients_determining_edge_intersections = true;
      flag_edge_intersection = false;
      sharpiso_param.allow_duplicates = true;
    }
    else if (s == "-subgrid") {
      flag_subgrid = true;
    }
    else if (s == "-sharp_edgeI") {
      sharpiso_param.use_sharp_edgeI = true;
    }
    else if (s == "-gradS_offset") {
      sharpiso_param.grad_selection_cube_offset = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-lindstrom") {
      flag_use_lindstrom = true;
    }
    else if (s == "-lindstrom2") {
      flag_use_lindstrom = true;
      flag_use_lindstrom2 = true;
    }
    else if (s == "-lindstrom_fast"){
    	flag_use_lindstrom = true;
      flag_use_lindstrom_fast = true;
    }
    else if (s == "-rayI") {
      flag_use_lindstrom = false;
    }
    else if (s == "-max_dist") {
      sharpiso_param.max_dist = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-max_mag") {
      sharpiso_param.max_small_magnitude = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-no_round") {
      sharpiso_param.flag_round = false;
    }
    else if (s == "-round") {
      sharpiso_param.flag_round = true;
      sharpiso_param.round_denominator = get_int(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-allow_conflict") {
      sharpiso_param.flag_allow_conflict = true;
    }
    else if (s == "-clamp_conflict") {
      sharpiso_param.flag_clamp_conflict = true;
    }
    else if (s == "-centroid_conflict") {
      sharpiso_param.flag_clamp_conflict = false;
    }
    else if (s == "-clamp_far") {
      sharpiso_param.flag_clamp_far = true;
    }
    else if (s == "-dist2centroid") {
      sharpiso_param.flag_dist2centroid = true;
    }
    else if (s == "-dist2center") {
      sharpiso_param.flag_dist2centroid = false;
    }
    else if (s == "-no_Linf") {
      sharpiso_param.use_Linf_dist = false;
    }
    else if (s == "-max_eigen") {
      sharpiso_param.max_small_eigenvalue = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-no_output_param") {
      flag_output_param = false;
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

  if (!flag_isovalue_set && !flag_location_set && !flag_list_gradients &&
      !flag_dist2vert) {
    cerr << "Error.  Option  -isovalue or -coord or -listg or -dist2vert must be specified."
         << endl;
    usage_error();
    exit(15);
  }

  if (!flag_isovalue_set && flag_list_subgrid) {
    cerr << "Error.  Option -list_subgrid cannot be used without -isovalue."
         << endl;
    exit(15);
  }

  if (!flag_isovalue_set && flag_edge_intersect_interpolate) {
    cerr << "Error.  Option -gradES cannot be used without -isovalue."
         << endl;
    exit(15);
  }

  if (!flag_isovalue_set && flag_edge_intersect_sharp) {
    cerr << "Error.  Option -gradEC cannot be used without -isovalue."
         << endl;
    exit(15);
  }

  if (!flag_isovalue_set && flag_list_gradients && use_selected_gradients) {
    cerr << "Error. Option -isovalue required when listing selected gradients."
         << endl;
    exit(15);
  }

  if (flag_dist2vert && (!flag_vertex_set && vertex_coord.size() == 0)) {
    cerr << "Error. Options -vertex or -vc required when using -dist2vert."
         << endl;
    exit(15);
  }

  if (!flag_isovalue_set && flag_location_set && use_selected_gradients) {
    cerr << "Error. Option -isovalue required when using selected gradients."
         << endl;
    exit(15);
  }

  if (sharpiso_param.grad_selection_cube_offset <= -1) {
    cerr << "Error in option -gradS_offset." << endl;
    cerr << "  Grad selection cube offset must be at least -1." << endl;
    exit(15);
  }

  if (sharpiso_param.grad_selection_cube_offset > 1) {
    cerr << "Error in option -gradS_offset." << endl;
    cerr << "  Grad selection cube offset must be at most 1." << endl;
    exit(15);
  }

  sharpiso_param.use_only_cube_gradients = use_only_cube_gradients;
  sharpiso_param.use_selected_gradients = use_selected_gradients;
  sharpiso_param.use_intersected_edge_endpoint_gradients = 
    use_intersected_edge_endpoint_gradients;
  sharpiso_param.use_gradients_determining_edge_intersections =
    use_gradients_determining_edge_intersections;
  sharpiso_param.use_lindstrom = flag_use_lindstrom;
  sharpiso_param.use_lindstrom2 = flag_use_lindstrom2;
  sharpiso_param.use_lindstrom_fast = flag_use_lindstrom_fast;
  sharpiso_param.select_based_on_grad_dir = select_based_on_grad_dir;
}

void help()
{
  cerr << "Usage: sharpinfo [OPTIONS] <scalar filename> <gradient filename>"
       << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -isovalue <isovalue>:  Compute isosurface vertex for given <isovalue>." << endl;
  cerr << "  -cube <cube_index>:  Compute isosurface vertex for cube <cube_index>." << endl;
  cerr << "           Default is cube 0." << endl;
  cerr << "  -cc \"cube coordinates\":  Compute isosurface vertex for cube"
       << endl
       << "           at given coordinates." << endl;
  cerr << "  -centroid:  Return centroid of isosurface-edge intersections."
       << endl;
  cerr << "  -gradC:  Use only cube gradients." << endl;
  cerr << "  -gradN:  Use gradients from cube and neighboring cubes." << endl;
  cerr << "  -gradCS: Use selected cube gradients." << endl;
  cerr << "           Isosurfaces from selected gradients must intersect the cube." << endl;
  cerr << "  -gradNS: Use selected gradients from cube and neighboring cubes."
       << endl;
  cerr << "           Isosurfaces from selected gradients must intersect the cube." << endl;
  cerr << "  -gradIE: Use gradients at endpoints of intersected cube edges."
       << endl;
  cerr << "  -gradIES: Use gradients at endpoints of intersected cube edges."
       << endl;
  cerr << "           Isosurfaces from selected gradients must intersect the cube." << endl;
  cerr << "  -gradIEDir: Use gradients at endpoints of intersected cube edges."
       << endl;
  cerr << "           Select gradients based on direction along edge."
       << endl;
  cerr << "  -gradCD: Use gradients determining intersections of isosurface"
       << endl
       << "           and cube edges." << endl;
  cerr << "  -gradCDdup: Use gradients determining intersections of isosurface"
       << endl
       << "           and cube edges." << endl;
  cerr << "           Allow duplicates." << endl;
  cerr << "  -gradNIE: Use gradients at endpoints of intersected cube edges"
       << endl
       << "              and intersected cube edges in neighboring cubes."
       << endl;
  cerr << "  -gradNIES: Use gradients at endpoints of intersected cube edges"
       << endl
       << "              and intersected cube edges in neighboring cubes."
       << endl;
  cerr << "           Isosurfaces from selected gradients must intersect the cube." << endl;
  cerr << "  -neighbor <cube_index>:  Use gradients from cube and" << endl
       << "           neighbors of cube <cube_index>." << endl;
  cerr << "  -gradES:    Interpolate intersection points/normals and"
       << endl
       << "           apply svd." << endl;
  cerr << "  -gradEC:    Compute edge-isosurface intersection points/normals"
       << endl
       << "                using gradient assignment and apply svd." << endl;
  cerr << "  -subgrid:   Output isosurface vertex based on subgrid." << endl;
  cerr << "  -lindstrom: Use Lindstrom's equation for calculating point on sharp feature." << endl;
  cerr << "  -sharp_edgeI: Use sharp formula for calculating"
       << endl
       << "                intersections of isosurface and grid edges." << endl;
  cerr << "  -interpolate_edgeI: Use linear interpolation for calculating"
       << endl
       << "                intersections of isosurface and grid edges." << endl;
  cout << "  -dist2center:  Use distance to center in lindstrom." << endl;
  cout << "  -dist2centroid:  Use distance to centroid of isourface-edge"
       << "                   intersections in lindstrom." << endl;
  cerr << "  -edgeI:     Output edge intersection information." << endl;
  cerr << "  -rayI:      Use intersection of ray and cubes for calculating point on sharp feature." << endl;
  cerr << "  -coord \"point_coord\":  Compute scalar values at coordinate point_coord." << endl;
  cerr << "  -facet <vertex_index> <orth_dir>: Set grid facet." << endl;
  cerr << "     Lower/left vertex of facet is <vertex_index>." << endl;
  cerr << "     Facet is orthgonal to direction <orth_dir>." << endl;
  cerr << "  -cube_offset2: Cube offset for intersection calculations."
       << endl;
  cerr << "                 Initial value set to 0.3." << endl;
  cerr << "  -max_eigen <V>: Normalized eigenvalues below V are set to zero." 
       << endl;
  cerr << "  -max_dist <V>:  Maximum distance from cube to isosurface vertex."
       << endl;
  cerr << "  -max_mag {max}:  Set maximum small gradient magnitude to max." 
       << endl;
  cerr << "           Gradients with magnitude below max are ignored." << endl;
  cerr << "  -listg: List gradients." << endl;
  cerr << "  -sortg: Sort gradients by distance of isoplane to cube center."
       << endl;
  cerr << "  -list_subgrid:  List all scalar values at vertices of subgrid." << endl;
  cerr << "  -list_edgeI: List intersections of isosurface and grid edges."
       << endl;
  cerr << "  -dist2vert: Distance from vertex to planes defined by gradients" 
       << endl
       << "              at neighboring vertices." << endl;
  cerr << "  -vertex <vertex_index>:  Set vertex." << endl;
  cerr << "  -vc \"vertex coordinates\":  Vertex coordinates."
       << endl;
  cout << "  -no_round:  Don't round coordinates." << endl;
  cout << "  -round <n>: Round coordinates to nearest 1/n." << endl;
  cout << "              Suggest using n=16,32,64,... or 2^k for some k."
       << endl;
  cout << "  -no_output_param:  Do not ouput parameters for sharpiso calculations." << endl;
  exit(15);
}

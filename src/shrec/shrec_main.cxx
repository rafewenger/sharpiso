/// \file shrec_main.cxx
/// generate isosurface using dual contouring algorithm
/// Version 0.1.0

/*
  Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

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


#include <iostream>

#include "shrec.h"
#include "shrecIO.h"

#include "ijkmesh.txx"
#include "ijkmesh_cpp11.txx"
#include "ijkmesh_geom.txx"

using namespace IJK;
using namespace SHREC;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_isosurface
(const INPUT_INFO & input_info, const SHREC_DATA & shrec_data,
 SHREC_TIME & shrec_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  SHREC_TIME shrec_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};
  INPUT_INFO input_info;
  IJK::ERROR error;
  bool flag_gradient(false);
  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv, input_info);

    SHARPISO_SCALAR_GRID full_scalar_grid;
    NRRD_INFO nrrd_info;
    read_nrrd_file
      (input_info.scalar_filename, full_scalar_grid,  nrrd_info, io_time);

    GRADIENT_GRID full_gradient_grid;
    NRRD_INFO nrrd_gradient_info;
    std::vector<COORD_TYPE> edgeI_coord;
    std::vector<GRADIENT_COORD_TYPE> edgeI_normal_coord;

    if (input_info.GradientsRequired()) {

      string gradient_filename;

      if (input_info.gradient_filename == NULL) {
        construct_gradient_filename
          (input_info.scalar_filename, gradient_filename);
      }
      else {
        gradient_filename = string(input_info.gradient_filename);
      }

      read_nrrd_file(gradient_filename.c_str(), full_gradient_grid,
                     nrrd_gradient_info);
      flag_gradient = true;

      if (!full_gradient_grid.CompareSize(full_scalar_grid)) {
        error.AddMessage("Input error. Grid mismatch.");
        error.AddMessage
          ("  Dimension or axis sizes of gradient grid and scalar grid do not match.");
        throw error;
      }
    }
    else if (input_info.NormalsRequired()) {

      if (input_info.normal_filename == NULL) {
        error.AddMessage("Programming error.  Missing normal filename.");
        throw error;
      }

      read_off_file
        (input_info.normal_filename, edgeI_coord, edgeI_normal_coord);
    }

    if (!check_input(input_info, full_scalar_grid, error))
      { throw(error); };

    // copy nrrd_info into input_info
    set_input_info(nrrd_info, input_info);

    // set DUAL datastructures and flags
    SHREC_DATA shrec_data;
    shrec_data.grad_selection_cube_offset = 0.1;

    if (flag_gradient) {
      shrec_data.SetGrids
        (full_scalar_grid, full_gradient_grid,
         input_info.flag_subsample, input_info.subsample_resolution,
         input_info.flag_supersample, input_info.supersample_resolution);
    }
    else
    {
      shrec_data.SetScalarGrid
        (full_scalar_grid, 
         input_info.flag_subsample, input_info.subsample_resolution,
         input_info.flag_supersample, input_info.supersample_resolution);

      if (input_info.NormalsRequired()) {
        shrec_data.SetEdgeI(edgeI_coord, edgeI_normal_coord);
      }

    }
    // Note: shrec_data.SetScalarGrid or shrec_data.SetGrids
    //       must be called before set_shrec_data.
    set_shrec_data(input_info, shrec_data, shrec_time);

    if (input_info.flag_output_param) 
      { report_shrec_param(shrec_data); }
 

    report_num_cubes(full_scalar_grid, input_info, shrec_data);
    construct_isosurface(input_info, shrec_data, shrec_time, io_time);

    if (input_info.report_time_flag) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(input_info, io_time, shrec_time, total_elapsed_time);
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
/**
* Construct Isosurface 
* Rescaling using spacing now done when ISO coordinated  are computed. 
*/
void construct_isosurface
(const INPUT_INFO & input_info, const SHREC_DATA & shrec_data,
 SHREC_TIME & shrec_time, IO_TIME & io_time)
{
  const int dimension = shrec_data.ScalarGrid().Dimension();
  const int num_cube_vertices = compute_num_cube_vertices(dimension);
  const int num_cubes = shrec_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;
  for (unsigned int i = 0; i < input_info.isovalue.size(); i++) {

    ISOVERT isovert;

    const SCALAR_TYPE isovalue = input_info.isovalue[i];

    DUAL_ISOSURFACE dual_isosurface;
    SHREC_INFO shrec_info(dimension);
    shrec_info.grid.num_cubes = num_cubes;

    dual_contouring
      (shrec_data, isovalue, dual_isosurface, isovert, shrec_info);
    shrec_time.Add(shrec_info.time);

    OUTPUT_INFO output_info;
    set_output_info(input_info, i, output_info);

    int grow_factor = 1;
    int shrink_factor = 1;
    if (input_info.flag_subsample)
      { grow_factor = input_info.subsample_resolution; }
    if (input_info.flag_supersample)
      { shrink_factor = input_info.supersample_resolution; }

    if (shrec_data.flag_convert_quad_to_tri) {

      VERTEX_INDEX_ARRAY quad_vert(dual_isosurface.quad_vert);
      VERTEX_INDEX_ARRAY quad_vert2;
      DUAL_ISOSURFACE isosurface_tri_mesh;
      isosurface_tri_mesh.vertex_coord = dual_isosurface.vertex_coord;
      isosurface_tri_mesh.tri_vert = dual_isosurface.tri_vert;

      IJK::reorder_quad_vertices(quad_vert);

      triangulate_quad_sharing_multiple_edges
        (quad_vert, isosurface_tri_mesh.tri_vert, quad_vert2);

      if (shrec_data.quad_tri_method == SPLIT_MAX_ANGLE) {

        // *** CREATE create_dual_tri IN shrec.cxx ***
        triangulate_quad_split_max_angle
          (DIM3, isosurface_tri_mesh.vertex_coord, quad_vert2,
           shrec_data.max_small_magnitude, isosurface_tri_mesh.tri_vert);
      }
      else {
        triangulate_quad(quad_vert2, isosurface_tri_mesh.tri_vert);
      }
	   
      output_dual_isosurface
        (output_info, shrec_data, isosurface_tri_mesh, isovert,
         shrec_info, io_time);
    }
    else {
      output_dual_isosurface
        (output_info, shrec_data, dual_isosurface, isovert,
         shrec_info, io_time);
    }
  }

}

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


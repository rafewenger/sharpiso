/// \file isodual3D.cxx
/// generate isosurface using dual contouring algorithm
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

#include "isodual3DIO.h"
#include "isodual3D.h"

using namespace IJK;
using namespace ISODUAL3D;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_isosurface
(const IO_INFO & io_info, const ISODUAL_DATA & isodual_data,
 ISODUAL_TIME & isodual_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  ISODUAL_TIME isodual_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;
  bool flag_gradient(false);

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv, io_info);

    ISODUAL_SCALAR_GRID full_scalar_grid;
    NRRD_INFO nrrd_info;
    read_nrrd_file
      (io_info.scalar_filename, full_scalar_grid,  nrrd_info, io_time);

    GRADIENT_GRID full_gradient_grid;
    NRRD_INFO nrrd_gradient_info;

    if (io_info.GradientsRequired()) {

      string gradient_filename;

      if (io_info.gradient_filename == NULL) {
        construct_gradient_filename
          (io_info.scalar_filename, gradient_filename);
      }
      else {
        gradient_filename = string(io_info.gradient_filename);
      }

      read_nrrd_file(gradient_filename.c_str(), full_gradient_grid,
                     nrrd_gradient_info);
      flag_gradient = true;

      if (!full_gradient_grid.CompareSize(full_scalar_grid)) {
        error.AddMessage("Input error. Grid mismatch.");
        error.AddMessage("  Dimension or axis sizes of gradient grid and scalar grid do not match.");
        throw error;
      }
    }

    if (!check_input(io_info, full_scalar_grid, error))
      { throw(error); };

    // copy nrrd_info into io_info
    set_io_info(nrrd_info, io_info);

    // set DUAL datastructures and flags
    ISODUAL_DATA isodual_data;
    isodual_data.grad_selection_cube_offset = 0.1;
    isodual_data.ray_intersection_cube_offset = 0.3;

    if (flag_gradient) {
      isodual_data.SetGrids
        (full_scalar_grid, full_gradient_grid,
         io_info.flag_subsample, io_info.subsample_resolution,
         io_info.flag_supersample, io_info.supersample_resolution);
    }
    else
    {
      isodual_data.SetScalarGrid
        (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
         io_info.flag_supersample, io_info.supersample_resolution);
    }
    // Note: isodual_data.SetScalarGrid or isodual_data.SetGrids
    //       must be called before set_mesh_data.
    set_isodual_data(io_info, isodual_data, isodual_time);
    report_num_cubes(full_scalar_grid, io_info, isodual_data);

    construct_isosurface(io_info, isodual_data, isodual_time, io_time);

    if (io_info.report_time_flag) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, isodual_time, total_elapsed_time);
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

void construct_isosurface
(const IO_INFO & io_info, const ISODUAL_DATA & isodual_data,
 ISODUAL_TIME & isodual_time, IO_TIME & io_time)
{
  const int dimension = isodual_data.ScalarGrid().Dimension();
  const int num_cube_vertices = compute_num_cube_vertices(dimension);
  const int num_cubes = isodual_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;
  for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue = io_info.isovalue[i];

    DUAL_ISOSURFACE dual_isosurface(num_cube_vertices);
    ISODUAL_INFO isodual_info(dimension);
    isodual_info.grid.num_cubes = num_cubes;

    dual_contouring(isodual_data, isovalue, dual_isosurface, isodual_info);
    isodual_time.Add(isodual_info.time);

    OUTPUT_INFO output_info;
    set_output_info(io_info, i, output_info);

    VERTEX_INDEX num_poly = dual_isosurface.NumIsoPoly();

    int grow_factor = 1;
    int shrink_factor = 1;
    if (io_info.flag_subsample)
      { grow_factor = io_info.subsample_resolution; }
    if (io_info.flag_supersample)
      { shrink_factor = io_info.supersample_resolution; }

    rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
                         dual_isosurface.vertex_coord);

    output_dual_isosurface
      (output_info, isodual_data, dual_isosurface, isodual_info, io_time);
  }
}

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


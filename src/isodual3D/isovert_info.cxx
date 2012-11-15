/// \file isovert_info.cxx
/// Print isovert information
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012 Rephael Wenger

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
#include <iomanip>

#include "isodual3D.h"
#include "isodual3DIO.h"
#include "isodual3D_decimate.h"
#include "isodual3D_extract.h"
#include "isodual3D_isovert.h"

#include "ijkcoord.txx"
#include "ijkmesh.txx"
#include "ijkmesh_geom.txx"
#include "ijkprint.txx"

using namespace IJK;
using namespace ISODUAL3D;

using namespace std;

// local subroutines
void memory_exhaustion();
void print_isovert_info
(const INPUT_INFO & input_info, const ISODUAL_DATA & isodual_data);
void compute_eigenvalues
(const SCALAR_TYPE isovalue, 
 const ISODUAL_DATA & isodual_data, const GRID_CUBE & gcube,
 EIGENVALUE_TYPE eigenvalues[DIM3], 
 EIGENVALUE_TYPE normalized_eigenvalues[DIM3],
 NUM_TYPE num_large_eigenvalues);
void out_gcube(std::ostream & out, const SCALAR_TYPE isovalue,
               const ISODUAL_DATA & isodual_data, const GRID_CUBE & gcube);
void out_gcube_type(std::ostream & out, const GRID_CUBE_FLAG flag);

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  ISODUAL_TIME isodual_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};
  INPUT_INFO input_info;
  IJK::ERROR error;
  bool flag_gradient(false);
  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv, input_info);

    ISODUAL_SCALAR_GRID full_scalar_grid;
    NRRD_INFO nrrd_info;
    read_nrrd_file
      (input_info.scalar_filename, full_scalar_grid,  nrrd_info, io_time);

    GRADIENT_GRID full_gradient_grid;
    NRRD_INFO nrrd_gradient_info;

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

    if (!check_input(input_info, full_scalar_grid, error))
      { throw(error); };

    // copy nrrd_info into input_info
    set_input_info(nrrd_info, input_info);

    // set DUAL datastructures and flags
    ISODUAL_DATA isodual_data;
    isodual_data.grad_selection_cube_offset = 0.1;

    if (flag_gradient) {
      isodual_data.SetGrids
        (full_scalar_grid, full_gradient_grid,
         input_info.flag_subsample, input_info.subsample_resolution,
         input_info.flag_supersample, input_info.supersample_resolution);
    }
    else
    {
      isodual_data.SetScalarGrid
        (full_scalar_grid, input_info.flag_subsample, input_info.subsample_resolution,
         input_info.flag_supersample, input_info.supersample_resolution);
    }
    // Note: isodual_data.SetScalarGrid or isodual_data.SetGrids
    //       must be called before set_mesh_data.
    set_isodual_data(input_info, isodual_data, isodual_time);
    if (input_info.flag_output_param) 
      { report_isodual_param(isodual_data); }

    report_num_cubes(full_scalar_grid, input_info, isodual_data);

    print_isovert_info(input_info, isodual_data);
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

void print_isovert_info
(const INPUT_INFO & input_info, const ISODUAL_DATA & isodual_data)
{
  const int dimension = isodual_data.ScalarGrid().Dimension();
  const int num_cube_vertices = compute_num_cube_vertices(dimension);
  const int num_cubes = isodual_data.ScalarGrid().ComputeNumCubes();
  GRID_COORD_TYPE cube_coord[DIM3];
  IJK::BOX<COORD_TYPE> box(DIM3);

  box.SetAllMinCoord(0);
  box.SetMaxCoord(isodual_data.ScalarGrid().AxisSize());
  if (input_info.minc.size() == DIM3) 
    { box.SetMinCoord(vector2pointer(input_info.minc)); };
  if (input_info.maxc.size() == DIM3) 
    { box.SetMaxCoord(vector2pointer(input_info.maxc)); };

  for (unsigned int i = 0; i < input_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue = input_info.isovalue[i];

    ISOVERT isovert;
    DUAL_ISOSURFACE dual_isosurface;
    ISODUAL_INFO isodual_info;

    isovert.linf_dist_threshold = isodual_data.linf_dist_thresh_merge_sharp;
    compute_dual_isovert
      (isodual_data.ScalarGrid(), isodual_data.GradientGrid(),
       isovalue, isodual_data, isovert);

    extract_dual_isopoly
      (isodual_data.ScalarGrid(), isovalue, isovert, 
       dual_isosurface, isodual_info);

    const NUM_TYPE num_gcube = isovert.gcube_list.size();
    std::vector<VERTEX_INDEX> gcube_map(num_gcube);
    decimate_dual_isopoly(isovert, dual_isosurface, gcube_map);

    for (int i = 0; i < isovert.gcube_list.size(); i++) {
      VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
      isodual_data.ScalarGrid().ComputeCoord(cube_index, cube_coord);

      if (box.Contains(cube_coord)) {
        out_gcube(cout, isovalue, isodual_data, isovert.gcube_list[i]); 

        VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
        ISO_VERTEX_INDEX isov = isovert.sharp_ind_grid.Scalar(cube_index);
        if (gcube_map[isov] != isov) {
          cout << "      Mapped to isovert for cube: "
               << isovert.gcube_list[gcube_map[isov]].cube_index << endl;
        }
      }

     }
  }
}

void out_gcube
(std::ostream & out, const SCALAR_TYPE isovalue, 
 const ISODUAL_DATA & isodual_data, const GRID_CUBE & gcube)
{
  const VERTEX_INDEX cube_index = gcube.cube_index;
  GRID_COORD_TYPE cube_coord[DIM3];
  COORD_TYPE isovert_coord[DIM3];
  EIGENVALUE_TYPE eigenvalues[DIM3];
  EIGENVALUE_TYPE normalized_eigenvalues[DIM3];
  NUM_TYPE num_large_eigenvalues;

  isodual_data.ScalarGrid().ComputeCoord(cube_index, cube_coord);

  compute_eigenvalues
    (isovalue, isodual_data, gcube, eigenvalues, normalized_eigenvalues,
     num_large_eigenvalues);

  out << "Cube: " << setw(6) << cube_index << " ";
  print_coord3D(out, cube_coord);
  out << ". Type ";
  out_gcube_type(out, gcube.flag);
  out << ".  Isovert: ";
  IJK::print_coord3D(out, gcube.isovert_coord);
  out << ". Linf dist: " << gcube.linf_dist << ".";
  out << endl;
  out << "  ";
  out << "    # eigen: " << int(gcube.num_eigen);
  out << ".  Eigenvalues: ";
  print_coord3D(out, eigenvalues);
  out << ".  Normalized eigenvalues: ";
  print_coord3D(out, normalized_eigenvalues);
  out << "." << endl;
}

void out_gcube_type(std::ostream & out, const GRID_CUBE_FLAG flag)
{
  out << setw(9);
  switch(flag) {
  case (AVAILABLE_GCUBE):
    out << "Available";
    break;

  case (SELECTED_GCUBE):
    out << "Selected";
    break;

  case (COVERED_GCUBE):
    out << "Covered";
    break;

  case (UNAVAILABLE_GCUBE):
    out << "Unavailable";
    break;

  case (SMOOTH_GCUBE):
    out << "Smooth";
    break;

  default:
    out << "Unknown";
    break;
  }
}

void compute_eigenvalues
(const SCALAR_TYPE isovalue, 
 const ISODUAL_DATA & isodual_data, const GRID_CUBE & gcube,
 EIGENVALUE_TYPE eigenvalues[DIM3], 
 EIGENVALUE_TYPE normalized_eigenvalues[DIM3],
 NUM_TYPE num_large_eigenvalues)
{
  const VERTEX_INDEX cube_index = gcube.cube_index;
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
			isodual_data.grad_selection_cube_offset;
  COORD_TYPE isovert_coord[DIM3];
	OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);
	SVD_INFO svd_info;

  svd_compute_sharp_vertex_for_cube
    (isodual_data.ScalarGrid(), isodual_data.GradientGrid(), cube_index, 
     isovalue, isodual_data, cube_111, isovert_coord,
     eigenvalues, num_large_eigenvalues, svd_info);

  if (num_large_eigenvalues != gcube.num_eigen) {
    cerr << "Warning: Number of large eigenvalues mismatch for cube "
         << cube_index << endl;
    cerr << "  gcube.num_eigen: " << gcube.num_eigen << ".";
    cerr << "  svd_compute...() num_large_eigen: "
         << num_large_eigenvalues << endl;
  }

  if (!is_coord_equal(DIM3, isovert_coord, gcube.isovert_coord)) {
    cerr << "Warning: Number of large eigenvalues mismatch for cube "
         << cube_index << endl;
    cerr << "  gcube.isovert_coord[] = ";
    print_coord3D(cerr, gcube.isovert_coord);
    cerr << ".  svd_compute...() isovert_coord[] = ";
    print_coord3D(cerr, isovert_coord);
    cerr << "." << endl;
  }

  if (eigenvalues[0] > isodual_data.max_small_eigenvalue) {
    multiply_coord_3D(1.0/eigenvalues[0], eigenvalues, normalized_eigenvalues);
  }
  else {
    set_coord_3D(0, normalized_eigenvalues);
  }
}


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


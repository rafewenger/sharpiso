/// \file isodual3DIO.cxx
/// IO routines for isodual3D

/*
 IJK: Isosurface Jeneration Kode
 Copyright (C) 2011,2012 Rephael Wenger

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

#include <assert.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "isodual3DIO.h"
#include "ijkIO.txx"
#include "ijkmesh.txx"
#include "ijkstring.txx"

using namespace IJK;
using namespace ISODUAL3D;

using namespace std;

// **************************************************
// PARSE COMMAND LINE
// **************************************************

// local namespace
namespace {

  typedef enum
    {SUBSAMPLE_PARAM, SUPERSAMPLE_PARAM,
     GRADIENT_PARAM, POSITION_PARAM, TRIMESH_PARAM, UNIFORM_TRIMESH_PARAM,
     MAX_EIGEN_PARAM, MAX_DIST_PARAM, GRAD_S_OFFSET_PARAM, RAY_I_OFFSET_PARAM,
     HELP_PARAM, OFF_PARAM, IV_PARAM, OUTPUT_PARAM_PARAM,
     OUTPUT_FILENAME_PARAM, STDOUT_PARAM,
     NOWRITE_PARAM, SILENT_PARAM, TIME_PARAM, UNKNOWN_PARAM} PARAMETER;
  const char * parameter_string[] =
    {"-subsample", "-supersample",
     "-gradient", "-position", "-trimesh", "-uniform_trimesh",
     "-max_eigen", "-max_dist", "-gradS_offset", "-rayI_offset",
     "-help", "-off", "-iv", "-out_param",
     "-o", "-stdout",
     "-nowrite", "-s", "-time", "-unknown"};

  PARAMETER get_parameter_token(char * s)
  // convert string s into parameter token
  {
    for (int i = 0; i < int(UNKNOWN_PARAM); i++)
      if (strcmp(parameter_string[i], s) == 0)
        return(PARAMETER(i));
    return(UNKNOWN_PARAM);
  }

  INTERPOLATION_TYPE get_interpolation_type(char * s)
  // convert string s into parameter token
  {
    INTERPOLATION_TYPE type = LINEAR_INTERPOLATION;

    if (strcmp(s, "linear") == 0)
      { type = LINEAR_INTERPOLATION; }
    else if (strcmp(s, "multilinear") == 0)
      { type = MULTILINEAR_INTERPOLATION; }
    else {
      cerr << "Error in input parameter -interpolate.  Illegal interpolation type: "
           << s << "." << endl;
      exit(1030);
    }

    return(type);
  }

  // Set vertex position method and flags.
  void set_vertex_position_method(const char * s, IO_INFO & io_info)
  {
    const string str = s;

    // set default method
    io_info.vertex_position_method = CENTROID_EDGE_ISO;
    io_info.use_only_cube_gradients = false;
    io_info.use_selected_gradients = true;
    io_info.use_intersected_edge_endpoint_gradients = false;

    if (str == "cube_center") {
      io_info.vertex_position_method = CUBECENTER;
    }
    else if (str == "centroid") {
      io_info.vertex_position_method = CENTROID_EDGE_ISO;
    }
    else if (str == "gradC"){
      io_info.vertex_position_method = GRADIENT_POSITIONING;
      io_info.use_only_cube_gradients = true;
      io_info.use_selected_gradients = false;
    }
    else if (str == "gradN"){
      io_info.vertex_position_method = GRADIENT_POSITIONING;
      io_info.use_only_cube_gradients = false;
      io_info.use_selected_gradients = false;
    }
    else if (str == "gradCS"){
      io_info.vertex_position_method = GRADIENT_POSITIONING;
      io_info.use_only_cube_gradients = true;
      io_info.use_selected_gradients = true;
    }
    else if (str == "gradNS"){
      io_info.vertex_position_method = GRADIENT_POSITIONING;
      io_info.use_only_cube_gradients = false;
      io_info.use_selected_gradients = true;
    }
    else if (str == "gradIE"){
      io_info.vertex_position_method = GRADIENT_POSITIONING;
      io_info.use_only_cube_gradients = true;
      io_info.use_selected_gradients = false;
      io_info.use_intersected_edge_endpoint_gradients = true;
    }
    else if (str == "gradES"){
      io_info.vertex_position_method = EDGE_SIMPLE;
    }
    else if (str == "gradEC"){
      io_info.vertex_position_method = EDGE_COMPLEX;
    }
    else {
      cerr << "Error in input parameter -position.  Illegal position method: "
           << str << "." << endl;
      exit(1030);
    }
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

}

void ISODUAL3D::parse_command_line(int argc, char **argv, IO_INFO & io_info)
// parse command line
// control parameters, followed by one or more isovalues,
// followed by input file name
{
  if (argc == 1) { usage_error(); };

  // Set default
  io_info.ray_intersection_cube_offset = 0.3;

  int iarg = 1;
  bool is_vertex_position_method_set = false;
  while (iarg < argc && argv[iarg][0] == '-') {
    PARAMETER param = get_parameter_token(argv[iarg]);
    if (param == UNKNOWN_PARAM) break;

    switch(param) {

    case SUBSAMPLE_PARAM:
      io_info.subsample_resolution = get_int(iarg, argc, argv);
      iarg++;
      io_info.flag_subsample = true;
      break;

    case SUPERSAMPLE_PARAM:
      io_info.supersample_resolution = get_int(iarg, argc, argv);
      iarg++;
      io_info.flag_supersample = true;
      break;

    case GRADIENT_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.gradient_filename = argv[iarg];
      break;

    case POSITION_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      set_vertex_position_method(argv[iarg], io_info);

      is_vertex_position_method_set = true;
      break;

    case TRIMESH_PARAM:
      io_info.flag_convert_quad_to_tri = true;
      io_info.quad_tri_method = SPLIT_MAX_ANGLE;
      break;

    case UNIFORM_TRIMESH_PARAM:
      io_info.flag_convert_quad_to_tri = true;
      io_info.quad_tri_method = UNIFORM_TRI;
      break;

    case MAX_EIGEN_PARAM:
      io_info.max_small_eigenvalue = get_float(iarg, argc, argv);
      iarg++;
      break;

    case MAX_DIST_PARAM:
      io_info.max_dist = get_float(iarg, argc, argv);
      iarg++;
      break;

    case GRAD_S_OFFSET_PARAM:
      io_info.grad_selection_cube_offset = get_float(iarg, argc, argv);
      iarg++;
      break;

    case RAY_I_OFFSET_PARAM:
      io_info.ray_intersection_cube_offset = get_float(iarg, argc, argv);
      iarg++;
      break;

    case OFF_PARAM:
      io_info.output_format = OFF;
      break;

    case IV_PARAM:
      io_info.output_format = IV;
      break;

    case OUTPUT_PARAM_PARAM:
      io_info.flag_output_param = true;
      break;

    case OUTPUT_FILENAME_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.output_filename = argv[iarg];
      break;

    case STDOUT_PARAM:
      io_info.use_stdout = true;
      break;

    case NOWRITE_PARAM:
      io_info.nowrite_flag = true;
      break;

    case SILENT_PARAM:
      io_info.flag_silent = true;
      break;

    case TIME_PARAM:
      io_info.report_time_flag = true;
      break;

    case HELP_PARAM:
      help();
      break;
    };

    iarg++;
  };

  // remaining parameters should be list of isovalues followed
  // by input file name

  // check for more parameter tokens
  for (int j = iarg; j < argc; j++) {
    if (get_parameter_token(argv[j]) != UNKNOWN_PARAM) {
      // argv[iarg] is not an isovalue
      cerr << "Error. Illegal parameter: " << argv[iarg] << endl;
      usage_error();
    }
  }

  if (iarg+2 > argc) {
    cerr << "Error.  Missing input isovalue or input file name." << endl;
    usage_error();
  };

  // store isovalues
  for (int j = iarg; j+1 < argc; j++) {
    io_info.isovalue_string.push_back(argv[j]);
    SCALAR_TYPE value;

    istringstream input_string(argv[j]);
    input_string >> value;

    if (input_string.fail()) {
      cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue."
           << endl;
      usage_error();
    };

    io_info.isovalue.push_back(value);
  }

  io_info.scalar_filename = argv[argc-1];

  if (!is_vertex_position_method_set && io_info.gradient_filename == NULL) {
    io_info.vertex_position_method = GRADIENT_POSITIONING;
  }

  if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
    cerr << "Error.  Subsample resolution must be an integer greater than 1."
         << endl;
    exit(230);
  };

  if (io_info.output_filename != NULL && io_info.use_stdout) {
    cerr << "Error.  Can't use both -o and -stdout parameters."
         << endl;
    exit(230);
  };

  if (io_info.flag_subsample && io_info.flag_supersample) {
    cerr << "Error.  Can't use both -subsample and -supersample parameters."
         << endl;
    exit(555);
  }
}

// Check input information/flags.
bool ISODUAL3D::check_input
(const IO_INFO & io_info,
 const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 IJK::ERROR & error)
{
  // Construct isosurface
  if (io_info.isovalue.size() > 1 && io_info.use_stdout) {
    error.AddMessage
      ("Error.  Cannot use stdout for more than one isovalue.");
    return(false);
  }

  if (io_info.isovalue.size() > 1 && io_info.output_filename != NULL) {
    error.AddMessage
      ("Error.  Cannot specify output file for more than one isovalue.");
    return(false);
  }

  return(true);
}

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

void ISODUAL3D::read_nrrd_file
(const char * input_filename, ISODUAL_SCALAR_GRID & scalar_grid,
 NRRD_INFO & nrrd_info)
{
  int dimension = 0;
  Nrrd *nin;

    // get scalar field data from nrrd file
  nin = nrrdNew();
  if (nrrdLoad(nin, input_filename, NULL)) {
    char *err = biffGetDone(NRRD);
    cerr << "Error reading: " << input_filename << endl;
    cerr << "  Error: " << err << endl;
    exit(35);
  };
  dimension = nin->dim;

  if (dimension < 1) {
    cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
    exit(20);
  };

  IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);
  size_t size[NRRD_DIM_MAX];
  double grid_spacing[NRRD_DIM_MAX];

  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size);
  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSpacing, grid_spacing);

  nrrd_info.grid_spacing.clear();
  for (int d = 0; d < dimension; d++) {
    axis_size[d] = size[d];
    nrrd_info.grid_spacing.push_back(grid_spacing[d]);
  }

  scalar_grid.SetSize(dimension, axis_size.PtrConst());
  nrrd2scalar(nin, scalar_grid.ScalarPtr());

  nrrdNuke(nin);

  nrrd_info.dimension = dimension;
}

void ISODUAL3D::read_nrrd_file
(const char * input_filename, GRADIENT_GRID & gradient_grid,
 NRRD_INFO & nrrd_info)
{
  int nrrd_dimension = 0;
  int dimension = 0;
  Nrrd *nin;

    // get scalar field data from nrrd file
  nin = nrrdNew();
  if (nrrdLoad(nin, input_filename, NULL)) {
    char *err = biffGetDone(NRRD);
    cerr << "Error reading: " << input_filename << endl;
    cerr << "  Error: " << err << endl;
    exit(35);
  };
  nrrd_dimension = nin->dim;

  if (nrrd_dimension < 2) {
    cerr << "Illegal nrrd dimension for gradient nrrd file." << endl;
    cerr << "  Dimension in nrrd file must be at least 2." << endl;
    exit(20);
  };

  dimension = nrrd_dimension-1;

  IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);
  size_t size[NRRD_DIM_MAX];
  double grid_spacing[NRRD_DIM_MAX];

  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size);
  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSpacing, grid_spacing);

  if (size[0] != dimension) {
    cerr << "Illegal gradient nrrd file." << endl;
    cerr << "  Gradient nrrd file should have size[0] (gradient vector length)"
    << endl;
    cerr << "     equal to (nrrd_dimension-1) (volume dimension)."
    << endl;
    cerr << "  size[0] = " << size[0]
    << ".  (nrrd_dimension-1) = " << nrrd_dimension-1 << "." << endl;
    exit(25);
  }

  nrrd_info.grid_spacing.clear();
  for (int d = 0; d < dimension; d++) {
    axis_size[d] = size[d+1];
    nrrd_info.grid_spacing.push_back(grid_spacing[d+1]);
  }

  int gradient_length = dimension;
  gradient_grid.SetSize(dimension, axis_size.PtrConst(), gradient_length);
  nrrd2scalar(nin, gradient_grid.VectorPtr());

  nrrdNuke(nin);

  nrrd_info.dimension = dimension;
}

void ISODUAL3D::read_nrrd_file
(const char * input_filename, ISODUAL_SCALAR_GRID & scalar_grid,
 NRRD_INFO & nrrd_info, IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  read_nrrd_file(input_filename, scalar_grid, nrrd_info);
  io_time.read_nrrd_time = wall_time.getElapsed();
}

  // **************************************************
  // PATH_DELIMITER
  // **************************************************

namespace {

#ifdef _WIN32
  const char PATH_DELIMITER = '\\';
#else
  const char PATH_DELIMITER = '/';
#endif
}


  // **************************************************
  // OUTPUT ISOSURFACE
  // **************************************************

void ISODUAL3D::output_dual_isosurface
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const DUAL_ISOSURFACE & dual_isosurface,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  if (output_info.flag_output_tri_mesh) {
    output_tri_isosurface
    (output_info, isodual_data,
     dual_isosurface.vertex_coord, dual_isosurface.tri_vert,
     isodual_info, io_time);
  }
  else {
    output_quad_isosurface
    (output_info, isodual_data,
     dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
     isodual_info, io_time);
  }
}

/// Output isosurface of quadrilaterals.
void ISODUAL3D::output_quad_isosurface
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const vector<COORD_TYPE> & vertex_coord,
 const vector<VERTEX_INDEX> & quad_vert,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, isodual_data,
                    vertex_coord, quad_vert, isodual_info);
  }

  if (!output_info.nowrite_flag) {
    write_dual_quad_mesh(output_info, vertex_coord, quad_vert, io_time);
  }
}

/// Output isosurface of triangles.
void ISODUAL3D::output_tri_isosurface
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const vector<COORD_TYPE> & vertex_coord,
 const vector<VERTEX_INDEX> & tri_vert,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, isodual_data,
                    vertex_coord, tri_vert, isodual_info);
  }

  if (!output_info.nowrite_flag) {
    write_dual_tri_mesh(output_info, vertex_coord, tri_vert, io_time);
  }
}

void ISODUAL3D::output_dual_isosurface_color
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const DUAL_ISOSURFACE & dual_isosurface,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  output_dual_isosurface_color
  (output_info, isodual_data,
   dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
   front_color, back_color, isodual_info, io_time);
}

void ISODUAL3D::output_dual_isosurface_color
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, isodual_data,
                    vertex_coord, slist, isodual_info);
  }

  if (!output_info.nowrite_flag) {
    write_dual_mesh_color
    (output_info, vertex_coord, slist, front_color, back_color, io_time);
  }
}

void ISODUAL3D::output_dual_isosurface_color_alternating
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const DUAL_ISOSURFACE & dual_isosurface,
 const ISODUAL_INFO & isodual_info, IO_TIME & io_time)
{
  const VERTEX_INDEX num_poly = dual_isosurface.NumIsoPoly();

  IJK::ARRAY<COLOR_TYPE> front_color(4*num_poly);
  IJK::ARRAY<COLOR_TYPE> back_color(4*num_poly);
  set_color_alternating
  (isodual_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
   front_color.Ptr());
  set_color_alternating
  (isodual_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
   back_color.Ptr());

  output_dual_isosurface_color
  (output_info, isodual_data, dual_isosurface.vertex_coord,
   dual_isosurface.isopoly_vert,
   front_color.PtrConst(), back_color.PtrConst(),
   isodual_info, io_time);
}

  // **************************************************
  // WRITE_DUAL_MESH
  // **************************************************

void ISODUAL3D::write_dual_quad_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_isopoly;
  const bool use_stdout = output_info.use_stdout;
  ofstream output_file;
  ERROR error_mcube("write_dual_mesh");

    // Output vertices in counter-clockwise order around quadrilateral.
  const bool flag_reorder_quad_vertices = true;

  string ofilename = output_info.output_filename;

  switch (output_info.output_format) {

    case OFF:
      if (!use_stdout) {
        output_file.open(ofilename.c_str(), ios::out);
        if (dimension == 3) {
          ijkoutQuadOFF(output_file, dimension, vertex_coord, plist,
                        flag_reorder_quad_vertices);
        }
        else {
          ijkoutOFF(output_file, dimension, numv_per_simplex,
                    vertex_coord, plist);
        }
        output_file.close();
      }
      else {
        if (dimension == 3) {
          ijkoutQuadOFF(dimension, vertex_coord, plist,
                        flag_reorder_quad_vertices);
        }
        else {
          ijkoutOFF(dimension, numv_per_simplex, vertex_coord, plist);
        }
      };
      break;

    case IV:
      if (dimension == 3) {
        if (!use_stdout) {
          output_file.open(ofilename.c_str(), ios::out);
          ijkoutIV(output_file, dimension, vertex_coord, plist);
          output_file.close();
        }
        else {
          ijkoutIV(dimension, vertex_coord, plist);
        }
      }
      else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
      break;

    default:
      throw error_mcube("Illegal output format.");
      break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}

void ISODUAL3D::write_dual_quad_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_quad_mesh(output_info, vertex_coord, plist);

  io_time.write_time += wall_time.getElapsed();
}

void ISODUAL3D::write_dual_mesh_color
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_isopoly;
  const bool use_stdout = output_info.use_stdout;

  ofstream output_file;
  ERROR error_mcube("write_dual_mesh_color");

  string ofilename = output_info.output_filename;

  switch (output_info.output_format) {

    case OFF:
      if (!use_stdout) {
        output_file.open(ofilename.c_str(), ios::out);
        ijkoutColorFacesOFF(output_file, dimension, numv_per_simplex,
                            vertex_coord, plist, front_color, back_color);
        output_file.close();
      }
      else {
        ijkoutColorFacesOFF(std::cout, dimension, numv_per_simplex,
                            vertex_coord, plist, front_color, back_color);
      };
      break;

    case IV:
      if (dimension == 3) {
        if (!use_stdout) {
          output_file.open(ofilename.c_str(), ios::out);
          ijkoutIV(output_file, dimension, vertex_coord, plist);
          output_file.close();
        }
        else {
          ijkoutOFF(dimension, vertex_coord, plist);
        }
      }
      else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
      break;

    default:
      throw error_mcube("Illegal output format.");
      break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}

void ISODUAL3D::write_dual_mesh_color
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_mesh_color(output_info, vertex_coord, plist,
                        front_color, back_color);

  io_time.write_time += wall_time.getElapsed();
}

/// Convert dual quadrilateral meshes to triangles and write mesh.
/// @param output_info Output information.
/// @param vertex_coord List of vertex coordinates.
/// @param quad_vert[] List of quadrilateral vertices.
///        quad_vert[4*i+k] is k'th quadrilateral vertices.
///        Quadrilateral vertices are listed in order:
///            Lower-Left, Lower-Right, Upper-Left, Upper-Right
void ISODUAL3D::write_dual_tri_mesh
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & tri_vert)
{
  const int NUMV_PER_QUAD = 4;
  const int NUMV_PER_TRI = 3;
  const int dimension = output_info.dimension;
  /* OBSOLETE
  const int numv_per_poly = output_info.num_vertices_per_isopoly;
  */
  const bool use_stdout = output_info.use_stdout;
  ofstream output_file;
  /* OBSOLETE
  std::vector<VERTEX_INDEX> tri_vert;
  */
  PROCEDURE_ERROR error("write_dual_tri_mesh");

  if (dimension != 3) {
    error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
    error.AddMessage("   Routine only allowed for dimension 3.");
    throw error;
  }

  /* OBSOLETE
  if (numv_per_poly != NUMV_PER_QUAD) {
    error.AddMessage
      ("Programming error.  Illegal number of vertices per quad: ",
       numv_per_poly, ".");
    error.AddMessage("  Should be ", NUMV_PER_QUAD, " vertices per quad.");
    throw error;
  }
  */

  /* OBSOLETE
  // Convert quadrilaterals to triangles
  convert_quad_to_tri(quad_vert, tri_vert);
  */

  string ofilename = output_info.output_filename;

  switch (output_info.output_format) {

    case OFF:
      if (!use_stdout) {
        output_file.open(ofilename.c_str(), ios::out);
        ijkoutOFF(output_file, dimension, NUMV_PER_TRI,
                  vertex_coord, tri_vert);
        output_file.close();
      }
      else {
        ijkoutOFF(dimension, NUMV_PER_TRI, vertex_coord, tri_vert);
      };
      break;

    case IV:
      ijkoutIV(dimension, vertex_coord, tri_vert);
      break;

    default:
      throw error("Illegal output format.");
      break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}

void ISODUAL3D::write_dual_tri_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord,
 const vector<VERTEX_INDEX> & quad_vert,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_tri_mesh(output_info, vertex_coord, quad_vert);

  io_time.write_time += wall_time.getElapsed();
}

// **************************************************
// RESCALE ROUTINES
// **************************************************

namespace {

  void grow_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
  for (unsigned int i = 0; i < vertex_coord.size(); i++) {
    vertex_coord[i] = scale * vertex_coord[i];
  };
  }

  void shrink_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
  for (unsigned int i = 0; i < vertex_coord.size(); i++) {
    vertex_coord[i] = vertex_coord[i]/scale;
  };
  }

  bool unit_spacing(const std::vector<COORD_TYPE> & spacing)
    // return true if spacing not defined or spacing along all axes equals 1.0
  {
  for (unsigned int d = 0; d < spacing.size(); d++) {
    if (!AIR_EXISTS(spacing[d])) { return(true); }
    else if (spacing[d] != 1.0) { return(false); };
  }

  return(true);
  }

  void rescale_coord(const std::vector<COORD_TYPE> & grid_spacing,
                     std::vector<COORD_TYPE> & vertex_coord)
  {
  const int dimension = grid_spacing.size();

  if (unit_spacing(grid_spacing)) { return; }

  const VERTEX_INDEX numv = vertex_coord.size()/dimension;
  for (int iv = 0; iv < numv; iv++) {
    for (int d = 0; d < dimension; d++) {
      vertex_coord[iv*dimension+d] *= grid_spacing[d];
    }
  };
  }

}

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
void ISODUAL3D::rescale_vertex_coord
(const OUTPUT_INFO & output_info, vector<COORD_TYPE> & vertex_coord)
{
  const int grow_factor = output_info.grow_factor;
  const int shrink_factor = output_info.shrink_factor;
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (grow_factor <= 0) {
    error.AddMessage("Illegal grow factor ", grow_factor, ".");
    error.AddMessage("  Grow factor must be a positive integer");
  }

  if (shrink_factor <= 0) {
    error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
    error.AddMessage("  Shrink factor must be a positive integer");
  }

  if (output_info.dimension != output_info.grid_spacing.size()) {
    error.AddMessage("Size of grid spacing array does not equal volume dimension.");
    error.AddMessage("  Grid spacing array has ",
                     output_info.grid_spacing.size(), " elements.");
    error.AddMessage("  Volume dimension = ", output_info.dimension, ".");
  }

  if (output_info.grow_factor != 1)
    { grow_coord(output_info.grow_factor, vertex_coord); };

  if (output_info.shrink_factor != 1)
    { shrink_coord(output_info.shrink_factor, vertex_coord); };

  rescale_coord(output_info.grid_spacing, vertex_coord);
}

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
void ISODUAL3D::rescale_vertex_coord
(const int grow_factor, const int shrink_factor,
 const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord)
{
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (grow_factor <= 0) {
    error.AddMessage("Illegal grow factor ", grow_factor, ".");
    error.AddMessage("  Grow factor must be a positive integer");
  }

  if (shrink_factor <= 0) {
    error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
    error.AddMessage("  Shrink factor must be a positive integer");
  }

  if (vertex_coord.size() == 0) { return; };

  if (grid_spacing.size() < 1) {
    error.AddMessage("Illegal size ", grid_spacing.size(),
                     " of array grid spacing.");
    error.AddMessage("Size must equal vertex dimension.");
    throw error;
  }

  if (grow_factor != 1)
    { grow_coord(grow_factor, vertex_coord); };

  if (shrink_factor != 1)
    { shrink_coord(shrink_factor, vertex_coord); };

  rescale_coord(grid_spacing, vertex_coord);
}

  // **************************************************
  // REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

void ISODUAL3D::report_num_cubes
(const ISODUAL_GRID & full_scalar_grid, const IO_INFO & io_info,
 const ISODUAL_DATA & isodual_data)
{
  const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
  const int num_cubes_in_isodual_data =
  isodual_data.ScalarGrid().ComputeNumCubes();

  if (!io_info.use_stdout && !io_info.flag_silent) {

    if (io_info.flag_subsample) {
        // subsampled grid
      cout << num_grid_cubes << " grid cubes.  "
      << num_cubes_in_isodual_data << " subsampled grid cubes." << endl;
    }
    else if (io_info.flag_supersample) {
        // supersample grid
      cout << num_grid_cubes << " grid cubes.  "
      << num_cubes_in_isodual_data << " supersampled grid cubes." << endl;
    }
    else {
        // use full_scalar_grid
      cout << num_grid_cubes << " grid cubes." << endl;
    }
  }

}

void ISODUAL3D::report_isodual_param(const ISODUAL_PARAM & isodual_param)
{
  cout << "Gradient selection cube offset: "
       << isodual_param.grad_selection_cube_offset << endl;
  cout << "Ray intersection cube offset: "
       << isodual_param.ray_intersection_cube_offset << endl;
  cout << "Maximum small eigenvalue: "
       << isodual_param.max_small_eigenvalue << endl;
  cout << "Max (Linf) distance from cube to isosurface vertex: "
       << isodual_param.max_dist << endl;
  cout << "Max small gradient magnitude: "
       << isodual_param.max_small_magnitude << endl;
  cout << endl;
}


void ISODUAL3D::report_iso_info
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const vector<COORD_TYPE> & vertex_coord,
 const vector<VERTEX_INDEX> & plist,
 const ISODUAL_INFO & isodual_info)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_isopoly;

  const char * indent4 = "    ";
  string grid_element_name = "cubes";
  if (dimension == 2) { grid_element_name = "squares"; };

  VERTEX_INDEX numv = (vertex_coord.size())/dimension;
  VERTEX_INDEX num_poly = (plist.size())/numv_per_simplex;
  VERTEX_INDEX num_grid_cubes = isodual_info.grid.num_cubes;
  VERTEX_INDEX num_non_empty_cubes = isodual_info.scalar.num_non_empty_cubes;

  float percent = 0.0;
  if (num_grid_cubes > 0)
    { percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
  int ipercent = int(100*percent);
  cout << "  Isovalue " << output_info.isovalue[0] << ".  "
       << numv << " isosurface vertices.  ";
  if (output_info.flag_output_tri_mesh) {
    cout << 2*num_poly << " isosurface triangles." << endl;
  }
  else {
    cout << num_poly << " isosurface quadrilaterals." << endl;
  }
}

  // **************************************************
  // REPORT TIMING INFORMATION
  // **************************************************

void ISODUAL3D::report_isodual_time
(const IO_INFO & io_info, const ISODUAL_TIME & isodual_time,
 const char * mesh_type_string)
{
  cout << "CPU time to run Marching Cubes: "
  << isodual_time.total << " seconds." << endl;
  cout << "    Time to extract " << mesh_type_string << " triangles: "
  << isodual_time.extract << " seconds." << endl;
  cout << "    Time to merge identical "
  << mesh_type_string << " vertices: "
  << isodual_time.merge << " seconds." << endl;
  cout << "    Time to position "
  << mesh_type_string << " vertices: "
  << isodual_time.position << " seconds." << endl;
}


void ISODUAL3D::report_time
(const IO_INFO & io_info, const IO_TIME & io_time,
 const ISODUAL_TIME & isodual_time, const double total_elapsed_time)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * INTERVAL_VOLUME_STRING = "interval volume";
  const char * mesh_type_string = NULL;

  mesh_type_string = ISOSURFACE_STRING;

  cout << "Time to read file " << io_info.scalar_filename << ": "
  << io_time.read_nrrd_time << " seconds." << endl;

  cout << "Time to read " << mesh_type_string << " lookup tables: "
  << io_time.read_table_time << " seconds." << endl;

  report_isodual_time(io_info, isodual_time, mesh_type_string);
  if (!io_info.nowrite_flag) {
    cout << "Time to write "
    << mesh_type_string << ": "
    << io_time.write_time << " seconds." << endl;
  };
  cout << "Total elapsed time: " << total_elapsed_time
  << " seconds." << endl;
}

  // **************************************************
  // USAGE/HELP MESSAGES
  // **************************************************

  // local namespace
namespace {

  void usage_msg(std::ostream & out)
  {
  out << "Usage: isodual3D [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
  }


  void options_msg()
  {
  cerr << "OPTIONS:" << endl;
  cerr << "  [-subsample S] [-supersample S]" << endl;
  cerr << "  [-position {centroid|cube_center|gradC|gradN|gradCS|gradNS|gradIE|gradES|gradEC}]" << endl;
  cerr << "  [-gradient {gradient_nrrd_filename}]" << endl;
  cerr << "  [-max_eigen {max}]" << endl;
  cerr << "  [-max_dist {D}] [-gradS_offset {offset}] [-rayI_offset {offset}]" 
       << endl;
  cerr << "  [-off|-iv] [-o {output_filename}] [-stdout]"
  << endl;
  cerr << "  [-help] [-s] [-out_param] [-nowrite] [-time]" << endl;
  }

}

void ISODUAL3D::usage_error()
{
  usage_msg(cerr);
  options_msg();
  exit(10);
}

void ISODUAL3D::help()
{
  usage_msg(cout);
  cout << endl;
  cout << "isodual3D - Marching cubes isosurface generation algorithm." << endl;
  cout << endl;
  cout << "OPTIONS:" << endl;

  cout << "  -subsample S: Subsample grid at every S vertices." << endl;
  cout << "                S must be an integer greater than 1." << endl;
  cout << "  -position {method}: Isosurface vertex position method." << endl;
  cout << "  -position centroid: Position isosurface vertices at centroid of"
       << endl;
  cout << "                      intersection of grid edges and isosurface."
       << endl;
  cout << "  -position cube_center: Position isosurface vertices at cube centers." << endl;
  cout << "  -position gradC: Position isosurface vertices using cube vertex gradients" << endl;
  cout << "                   and singular value decomposition (svd)." << endl;
  cout << "  -position gradN: Position isosurface vertices using svd on"
       << endl;
  cout << "                   vertex gradients of cube and cube neighbors."
       << endl;
  cout << "  -position gradCS: Position isosurface vertices using selected"
       << endl;
  cout << "                    cube vertex gradients and svd." << endl;
  cout << "  -position gradNS: Position isosurface vertices using svd"
       << endl;
  cout << "       on selected vertex gradients of cube and cube neighbors."
       << endl;
  cout << "  -position gradIE: Position isosurface vertices using svd"
       << endl;
  cout << "       on gradients at endpoints of intersected cube edges."
       << endl;
  cout << "  -position gradES: Position using isosurface vertices using svd"
       << endl;
  cout << "       on gradients at isosurface-edge intersections."
       << endl;
  cout << "       Use simple interpolation to compute isosurface-edge intersections." << endl;
  cout << "  -position gradEC: Position isosurface vertices using svd"
       << endl;
  cout << "       on computed gradients at isosurface-edge intersections."
       << endl;
  cout << "       Use complex choice to compute isosurface-edge intersections." << endl;
  cout << "  -gradient {gradient_nrrd_filename}: Read gradients from gradient nrrd file." << endl;
  cerr << "  -max_eigen {max}: Set maximum small eigenvalue to max."
       << endl;
  cerr << "  -max_dist {D}:    Set max Linf distance from cube to isosurface vertex." 
       << endl;
  cerr << "  -gradS_offset {offset}: Set cube offset for gradient selection to offset."
       << endl;
  cerr << "  -rayI_offset {offset}:  Set cube offset for ray intersection to offset."
       << endl;
  cout << "  -trimesh:   Output triangle mesh." << endl;
  cout << "  -off: Output in geomview OFF format. (Default.)" << endl;
  cout << "  -iv: Output in OpenInventor .iv format." << endl;
  cout << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
  cout << "  -stdout: Write isosurface to standard output." << endl;
  cout << "  -nowrite: Don't write isosurface." << endl;
  cout << "  -time: Output running time." << endl;
  cout << "  -s: Silent mode." << endl;
  cout << "  -help: Print this help message." << endl;
  cout << "  -out_param: Print isodual3D parameters." << endl;
  exit(20);
}


  // **************************************************
  // CLASS IO_INFO
  // **************************************************

  /// IO information
void ISODUAL3D::IO_INFO::Init()
{
  isovalue.clear();
  isovalue_string.clear();
  scalar_filename = NULL;
  gradient_filename = NULL;
  output_filename = NULL;
  isotable_directory = "";
  output_format = OFF;
  report_time_flag = false;
  use_stdout = false;
  nowrite_flag = false;
  flag_silent = false;
  flag_subsample = false;
  subsample_resolution = 2;
  flag_supersample = false;
  supersample_resolution = 2;
  flag_color_alternating = false;  // color simplices in alternating cubes
  region_length = 1;
  max_small_eigenvalue = 0.1;
  flag_output_param = false;
}

  // **************************************************
  // class OUTPUT_INFO
  // **************************************************

void ISODUAL3D::OUTPUT_INFO::Init()
{
  output_filename = "";
  dimension = 3;
  num_vertices_per_isopoly = 4;
  isovalue[0] = 0;
  isovalue[1] = 0;
  nowrite_flag = false;
  use_stdout = false;
  flag_silent = false;
  flag_output_tri_mesh = false;
  output_format = OFF;
  grow_factor = 1;
  shrink_factor = 1;
  grid_spacing.resize(3,1);
}

namespace {

  /// split string at last occurrence of character c into prefix and suffix
  void split_string(const string & s, const char c,
                    string & prefix, string & suffix)
  {
  string::size_type i = s.rfind(c);
  if (i == string::npos) {
    prefix = s;
    suffix = "";
  }
  else {
    if (i > 0) { prefix = s.substr(0,i); }
    else { prefix = ""; };

    if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
    else { suffix = ""; };
  }
  }

  string construct_output_filename
  (const IO_INFO & io_info, const int i)
  {
  string prefix, suffix;

    // create output filename
  string fname = string(io_info.scalar_filename);

#ifndef _WIN32
    // remove path from file name
  split_string(fname, PATH_DELIMITER, prefix, suffix);
  if (suffix != "") { fname = suffix; }
#endif

  string ofilename;

    // construct output filename
  split_string(fname, '.', prefix, suffix);
  if (suffix == "nrrd" || suffix == "nhdr") { ofilename = prefix; }
  else { ofilename = string(io_info.scalar_filename); }

  ofilename += string(".") + string("isov=") + io_info.isovalue_string[i];

  switch (io_info.output_format) {
    case OFF:
      ofilename += ".off";
      break;

    case IV:
      ofilename += ".iv";
      break;
  }

  return(ofilename);
  }

}

  // **************************************************
  // SET ROUTINES
  // **************************************************

void ISODUAL3D::set_isodual_data
(const IO_INFO & io_info, ISODUAL_DATA & isodual_data, ISODUAL_TIME & isodual_time)
{
  PROCEDURE_ERROR error("set_isodual_data");

  if (!isodual_data.IsScalarGridSet()) {
    error.AddMessage
      ("Programming error. Scalar field must be set before set_isodual_data is called.");
    throw error;
  }

  // Set data structures in isodual_data
  isodual_data.SetVertexPositionMethod(io_info.vertex_position_method);
  isodual_data.SetUseSelectedGradients(io_info.use_selected_gradients);
  isodual_data.SetUseOnlyCubeGradients(io_info.use_only_cube_gradients);
  isodual_data.SetUseIntersectedEdgeEndpointGradients
    (io_info.use_intersected_edge_endpoint_gradients);
  isodual_data.max_small_eigenvalue = io_info.max_small_eigenvalue;
  isodual_data.max_dist = io_info.max_dist;
  isodual_data.grad_selection_cube_offset =
    io_info.grad_selection_cube_offset;
  isodual_data.ray_intersection_cube_offset =
    io_info.ray_intersection_cube_offset;
  isodual_data.flag_convert_quad_to_tri = io_info.flag_convert_quad_to_tri;
  isodual_data.quad_tri_method = io_info.quad_tri_method;
}

void ISODUAL3D::set_io_info
(const NRRD_INFO & nrrd_info, IO_INFO & io_info)
{
  io_info.grid_spacing.clear();
  for (int d = 0; d < nrrd_info.dimension; d++) {
    io_info.grid_spacing.push_back(nrrd_info.grid_spacing[d]);
  }
}

void ISODUAL3D::set_output_info
(const IO_INFO & io_info,
 const int i, OUTPUT_INFO & output_info)
{
  output_info.nowrite_flag = io_info.nowrite_flag;
  output_info.use_stdout = io_info.use_stdout;
  output_info.flag_silent = io_info.flag_silent;
  /* OBSOLETE
  output_info.flag_output_tri_mesh = io_info.flag_output_tri_mesh;
  */
  output_info.flag_output_tri_mesh = io_info.flag_convert_quad_to_tri;

  output_info.grow_factor = 1;
  if (io_info.flag_subsample)
    { output_info.grow_factor = io_info.subsample_resolution; }

  output_info.shrink_factor = 1;
  if (io_info.flag_supersample)
    { output_info.shrink_factor = io_info.supersample_resolution; }

  output_info.grid_spacing.clear();
  output_info.grid_spacing.resize(io_info.grid_spacing.size());
  for (unsigned int j = 0; j < io_info.grid_spacing.size(); j++)
    { output_info.grid_spacing[j] = io_info.grid_spacing[j]; }

  output_info.output_format = io_info.output_format;
  output_info.isovalue[0] = io_info.isovalue[i];
  if (i+1 < int(io_info.isovalue.size()))
    { output_info.isovalue[1] = io_info.isovalue[i+1]; };

  if (io_info.output_filename != NULL) {
    output_info.output_filename = string(io_info.output_filename);
  }
  else {
    output_info.output_filename =
    construct_output_filename(io_info, i);
  }
}

void ISODUAL3D::set_color_alternating
(const ISODUAL_GRID & grid, const vector<VERTEX_INDEX> & cube_list,
 COLOR_TYPE * color)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<GRID_COORD_TYPE> coord(dimension);

  const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
  const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
  const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

  VERTEX_INDEX icube = 0;
  int parity = 0;
  COLOR_TYPE * color_ptr = color;
  for (unsigned int i = 0; i < cube_list.size(); i++) {
    int new_cube = cube_list[i];
    if (icube != new_cube) {
      icube = new_cube;
      grid.ComputeCoord(icube, coord.Ptr());
      int sum = 0;
      for (int d = 0; d < dimension; d++)
        { sum += coord[d]; }
      parity = sum%2;
    }

    if (parity == 0)
      { std::copy(red, red+3, color_ptr); }
    else
      { std::copy(blue, blue+3, color_ptr); }

      // set opacity
    color_ptr[3] = 1.0;
    color_ptr += 4;
  }

}

// **************************************************
// NRRD INFORMATION
// **************************************************

/// Construct gradient filename from scalar filename.
void ISODUAL3D::construct_gradient_filename
(const char * scalar_filename, std::string & gradient_filename)
{
  std::string prefix;
  std::string suffix;

  split_string(scalar_filename, '.', prefix, suffix);

  gradient_filename = prefix + ".grad." + suffix;
}

// **************************************************
// NRRD INFORMATION
// **************************************************

NRRD_INFO::NRRD_INFO()
{
  Clear();
}

NRRD_INFO::~NRRD_INFO()
{
  Clear();
}

void NRRD_INFO::Clear()
{
  dimension = 0;
}

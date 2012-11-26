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
#include "sharpiso_get_gradients.h"
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

typedef enum {
  SUBSAMPLE_PARAM,
  GRADIENT_PARAM, NORMAL_PARAM, POSITION_PARAM, POS_PARAM, 
  TRIMESH_PARAM, UNIFORM_TRIMESH_PARAM,
  MAX_EIGEN_PARAM, MAX_DIST_PARAM, GRAD_S_OFFSET_PARAM, 
  MAX_MAG_PARAM, SNAP_DIST_PARAM, 
  SHARP_EDGEI_PARAM, INTERPOLATE_EDGEI_PARAM,
	REPOSITION_PARAM, NO_REPOSITION_PARAM, SEPDIST_PARAM,
	ALLOW_CONFLICT_PARAM,
	CLAMP_CONFLICT_PARAM, CENTROID_CONFLICT_PARAM,
  MERGE_CONFLICT_PARAM, MERGE_SHARP_PARAM, MERGE_SHARP_LINF_THRES_PARAM,
	CLAMP_FAR_PARAM, CENTROID_FAR_PARAM,
	RECOMPUTE_EIGEN2_PARAM, NO_RECOMPUTE_EIGEN2_PARAM,
	REMOVEG_PARAM, NO_REMOVEG_PARAM,
	RESELECT_GRAD_PARAM, NO_RESELECT_GRAD_PARAM,
  DIST2CENTER_PARAM, DIST2CENTROID_PARAM,
	CENTROID_EIGEN1_PARAM, NO_CENTROID_EIGEN1_PARAM,
	LINF_PARAM, NO_LINF_PARAM,
	USE_LINDSTROM_PARAM,
	USE_LINDSTROM2_PARAM,
	USE_LINDSTROM_FAST,
	SINGLE_ISOV_PARAM, MULTI_ISOV_PARAM,
	SEP_NEG_PARAM, SEP_POS_PARAM, RESOLVE_AMBIG_PARAM,
	ROUND_PARAM, NO_ROUND_PARAM,
  MINC_PARAM, MAXC_PARAM,
	HELP_PARAM, OFF_PARAM, IV_PARAM, OUTPUT_PARAM_PARAM,
	OUTPUT_FILENAME_PARAM, STDOUT_PARAM,
	NOWRITE_PARAM, OUTPUT_INFO_PARAM, SILENT_PARAM,
	TIME_PARAM, UNKNOWN_PARAM} PARAMETER;
	const char * parameter_string[] =
	{ "-subsample",
    "-gradient", "-normal", "-position", "-pos", "-trimesh", "-uniform_trimesh",
    "-max_eigen", "-max_dist", "-gradS_offset", "-max_mag", "-snap_dist",
    "-sharp_edgeI", "-interpolate_edgeI",
    "-reposition", "-no_reposition", "-sepdist",
    "-allow_conflict", "-clamp_conflict", "-centroid_conflict", 
    "-merge_conflict", "-merge_sharp","-merge_linf_th",
    "-clamp_far", "-centroid_far",
    "-recompute_eigen2", "-no_recompute_eigen2",
    "-removeg", "-no_removeg",
    "-reselectg", "-no_reselectg",
    "-dist2center", "-dist2centroid",
    "-centroid_eigen1", "-no_centroid_eigen1",
    "-Linf", "-no_Linf",
    "-lindstrom",	"-lindstrom2","-lindstrom_fast",
    "-single_isov", "-multi_isov",
    "-sep_neg", "-sep_pos", "-resolve_ambig",
    "-round", "-no_round",
    "-minc", "-maxc",
    "-help", "-off", "-iv", "-out_param",
    "-o", "-stdout",
    "-nowrite", "-info", "-s", "-time", "-unknown"};

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
	void set_vertex_position_method(const char * s, INPUT_INFO & input_info)
	{
		const string str = s;

		if (str == "cube_center") {
			input_info.vertex_position_method = CUBECENTER;
		}
		else if (str == "centroid") {
			input_info.vertex_position_method = CENTROID_EDGE_ISO;
		}
		else if (str == "edgeIinterp" || str == "gradES"){
			input_info.vertex_position_method = EDGEI_INTERPOLATE;
		}
		else if (str == "edgeIgrad" || str == "gradEC"){
			input_info.vertex_position_method = EDGEI_GRADIENT;
		}
    else {
      input_info.vertex_position_method = GRADIENT_POSITIONING;

      GRAD_SELECTION_METHOD grad_selection_method
        = get_grad_selection_method(str);

      if (grad_selection_method == UNKNOWN_GRAD_SELECTION_METHOD) {
        cerr << "Error in input parameter -position.  Illegal position method: "
             << str << "." << endl;
        exit(1030);
      }

      input_info.SetGradSelectionMethod(grad_selection_method);
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

  /// Get string and convert to list of arguments.
  /// Does not modify iarg.
  template <typename ETYPE>
  void get_multiple_arguments
  (const int iarg, const int argc, char **argv, vector<ETYPE> & v)
  {
    if (iarg+1 >= argc) { 
      cerr << "Usage error. Missing argument for option " 
           << argv[iarg] << " and missing file name." << endl;
      usage_error();
    }

    if (!IJK::string2vector(argv[iarg+1], v)) {
      cerr << "Error in argument for option: " << argv[iarg] << endl;
      cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
      exit(50);
    }
  }

}

void ISODUAL3D::parse_command_line
(int argc, char **argv, INPUT_INFO & input_info)
// parse command line
// control parameters, followed by one or more isovalues,
// followed by input file name
{
	if (argc == 1) { usage_error(); };

	int iarg = 1;
	bool is_vertex_position_method_set = false;
  bool is_use_sharp_edgeI_set = false;
  bool is_conflict_set = false;
	while (iarg < argc && argv[iarg][0] == '-') {
		PARAMETER param = get_parameter_token(argv[iarg]);
		if (param == UNKNOWN_PARAM) break;

		switch(param) {

		case SUBSAMPLE_PARAM:
			input_info.subsample_resolution = get_int(iarg, argc, argv);
			iarg++;
			input_info.flag_subsample = true;
			break;

		case GRADIENT_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			input_info.gradient_filename = argv[iarg];
			break;

		case NORMAL_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			input_info.normal_filename = argv[iarg];
      input_info.vertex_position_method = EDGEI_INPUT_DATA;

			is_vertex_position_method_set = true;
			break;

		case POSITION_PARAM:
		case POS_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			set_vertex_position_method(argv[iarg], input_info);

			is_vertex_position_method_set = true;
			break;

		case TRIMESH_PARAM:
			input_info.flag_convert_quad_to_tri = true;
			input_info.quad_tri_method = SPLIT_MAX_ANGLE;
			break;

		case UNIFORM_TRIMESH_PARAM:
			input_info.flag_convert_quad_to_tri = true;
			input_info.quad_tri_method = UNIFORM_TRI;
			break;

		case MAX_EIGEN_PARAM:
			input_info.max_small_eigenvalue = get_float(iarg, argc, argv);
			iarg++;
			break;

		case MAX_DIST_PARAM:
			input_info.max_dist = get_float(iarg, argc, argv);
			iarg++;
			break;

		case MAX_MAG_PARAM:
			input_info.max_small_magnitude = get_float(iarg, argc, argv);
			iarg++;
			break;

		case SNAP_DIST_PARAM:
			input_info.snap_dist = get_float(iarg, argc, argv);
			iarg++;
			break;

    case SHARP_EDGEI_PARAM:
      input_info.use_sharp_edgeI = true;
      is_use_sharp_edgeI_set = true;
      break;

    case INTERPOLATE_EDGEI_PARAM:
      input_info.use_sharp_edgeI = false;
      is_use_sharp_edgeI_set = true;
      break;

		case GRAD_S_OFFSET_PARAM:
			input_info.grad_selection_cube_offset = get_float(iarg, argc, argv);
			iarg++;
			break;

		case SEPDIST_PARAM:
			input_info.separation_distance = get_float(iarg, argc, argv);
			iarg++;
			break;

		case REPOSITION_PARAM:
			input_info.flag_reposition = true;
			break;

		case NO_REPOSITION_PARAM:
			input_info.flag_reposition = false;
			break;

		case LINF_PARAM:
			input_info.use_Linf_dist = true;
			break;

		case NO_LINF_PARAM:
			input_info.use_Linf_dist = false;
			break;

		case USE_LINDSTROM_PARAM:
			input_info.use_lindstrom =true;
			break;

		case USE_LINDSTROM2_PARAM:
			input_info.use_lindstrom = true;
			input_info.use_lindstrom2 = true;
			break;

		case USE_LINDSTROM_FAST:
					input_info.use_lindstrom = true;
					input_info.use_lindstrom_fast = true;
					break;

		case SINGLE_ISOV_PARAM:
			input_info.allow_multiple_iso_vertices = false;
			break;

		case MULTI_ISOV_PARAM:
			input_info.allow_multiple_iso_vertices = true;
			break;

		case SEP_NEG_PARAM:
			input_info.flag_separate_neg = true;
      input_info.flag_resolve_ambiguous_facets = false;
			input_info.allow_multiple_iso_vertices = true;
			break;

		case SEP_POS_PARAM:
			input_info.flag_separate_neg = false;
      input_info.flag_resolve_ambiguous_facets = false;
			input_info.allow_multiple_iso_vertices = true;
			break;

		case RESOLVE_AMBIG_PARAM:
			input_info.flag_resolve_ambiguous_facets = true;
			input_info.allow_multiple_iso_vertices = true;
			break;

		case ALLOW_CONFLICT_PARAM:
			input_info.flag_allow_conflict = true;
      is_conflict_set = true;
			break;

		case CLAMP_CONFLICT_PARAM:
			input_info.flag_clamp_conflict = true;
      is_conflict_set = true;
			break;

		case CENTROID_CONFLICT_PARAM:
			input_info.flag_clamp_conflict = false;
      is_conflict_set = true;
			break;

		case MERGE_CONFLICT_PARAM:
			input_info.flag_merge_conflict = true;
      break;

		case MERGE_SHARP_PARAM:
			input_info.linf_dist_thresh_merge_sharp = 0.8;
			input_info.flag_merge_sharp = true;
      break;
      //
		case MERGE_SHARP_LINF_THRES_PARAM:
					input_info.linf_dist_thresh_merge_sharp = get_float(iarg, argc, argv);
					iarg++;
					input_info.flag_merge_sharp = true;
		      break;
      //

		case CLAMP_FAR_PARAM:
			input_info.flag_clamp_far = true;
			break;

		case CENTROID_FAR_PARAM:
			input_info.flag_clamp_far = false;
			break;

		case RECOMPUTE_EIGEN2_PARAM:
			input_info.flag_recompute_eigen2 = true;
			break;

		case NO_RECOMPUTE_EIGEN2_PARAM:
			input_info.flag_recompute_eigen2 = false;
			break;

		case REMOVEG_PARAM:
			input_info.flag_remove_gradients = true;
			break;

		case NO_REMOVEG_PARAM:
			input_info.flag_remove_gradients = false;
			break;

		case RESELECT_GRAD_PARAM:
			input_info.flag_reselect_gradients = true;
			break;

		case NO_RESELECT_GRAD_PARAM:
			input_info.flag_reselect_gradients = false;
			break;

    case DIST2CENTER_PARAM:
      input_info.flag_dist2centroid = false;
      break;

    case DIST2CENTROID_PARAM:
      input_info.flag_dist2centroid = true;
      break;

		case CENTROID_EIGEN1_PARAM:
			input_info.flag_centroid_eigen1 = true;
			break;

		case NO_CENTROID_EIGEN1_PARAM:
			input_info.flag_centroid_eigen1 = false;
			break;

		case ROUND_PARAM:
			input_info.flag_round = true;
			input_info.round_denominator = get_int(iarg, argc, argv);
			iarg++;
			break;

		case NO_ROUND_PARAM:
			input_info.flag_round = false;
			break;

    case MINC_PARAM:
      get_multiple_arguments(iarg, argc, argv, input_info.minc);
      iarg++;
      break;

    case MAXC_PARAM:
      get_multiple_arguments(iarg, argc, argv, input_info.maxc);
      iarg++;
      break;

		case OFF_PARAM:
			input_info.output_format = OFF;
			break;

		case IV_PARAM:
			input_info.output_format = IV;
			break;

		case OUTPUT_PARAM_PARAM:
			input_info.flag_output_param = true;
			break;

		case OUTPUT_FILENAME_PARAM:
			iarg++;
			if (iarg >= argc) usage_error();
			input_info.output_filename = argv[iarg];
			break;

		case STDOUT_PARAM:
			input_info.use_stdout = true;
			break;

		case NOWRITE_PARAM:
			input_info.nowrite_flag = true;
			break;

		case OUTPUT_INFO_PARAM:
			input_info.flag_output_alg_info = true;
			break;

		case SILENT_PARAM:
			input_info.flag_silent = true;
			break;

		case TIME_PARAM:
			input_info.report_time_flag = true;
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
		input_info.isovalue_string.push_back(argv[j]);
		SCALAR_TYPE value;

		istringstream input_string(argv[j]);
		input_string >> value;

		if (input_string.fail()) {
			cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue."
					<< endl;
			usage_error();
		};

		input_info.isovalue.push_back(value);
	}

	input_info.scalar_filename = argv[argc-1];

	if (!is_vertex_position_method_set && input_info.gradient_filename == NULL) {
		input_info.vertex_position_method = GRADIENT_POSITIONING;
    input_info.SetGradSelectionMethod(GRAD_NS);
	}

  if (!is_conflict_set && input_info.flag_merge_sharp) {
    // Set merge_sharp defaults.
    input_info.flag_allow_conflict = true;
    input_info.flag_clamp_conflict = false;
  }

  if (is_use_sharp_edgeI_set) {
    if (input_info.use_sharp_edgeI) {
      if (input_info.vertex_position_method == EDGEI_INTERPOLATE) {
        cerr << "Error.  Cannot use -interpolate_edgeI with -position gradES."
             << endl;
        exit(230);
      }
    }
    else {
      if (input_info.vertex_position_method == EDGEI_GRADIENT) {
        cerr << "Error.  Cannot use -interpolate_edgeI with -position gradEC."
             << endl;
        exit(230);
      }
    }
  }

  if (input_info.vertex_position_method == EDGEI_GRADIENT) {
    input_info.use_sharp_edgeI = true;
  }

  if (input_info.vertex_position_method == EDGEI_INTERPOLATE) {
    input_info.use_sharp_edgeI = false;
  }

	if (input_info.flag_subsample && input_info.subsample_resolution <= 1) {
		cerr << "Error.  Subsample resolution must be an integer greater than 1."
				<< endl;
		exit(230);
	};

	if (input_info.output_filename != NULL && input_info.use_stdout) {
		cerr << "Error.  Can't use both -o and -stdout parameters."
				<< endl;
		exit(230);
	};

	if (input_info.flag_subsample && input_info.flag_supersample) {
		cerr << "Error.  Can't use both -subsample and -supersample parameters."
				<< endl;
		exit(555);
	}

	if (input_info.flag_round) {
		if (input_info.round_denominator < 1) {
			cerr << "Error.  Illegal -round <n> parameter. Integer <n> must be positive." << endl;
			exit(560);
		}
	}
}

// Check input information/flags.
bool ISODUAL3D::check_input
(const INPUT_INFO & input_info,
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		IJK::ERROR & error)
{
	// Construct isosurface
	if (input_info.isovalue.size() > 1 && input_info.use_stdout) {
		error.AddMessage
		("Error.  Cannot use stdout for more than one isovalue.");
		return(false);
	}

	if (input_info.isovalue.size() > 1 && input_info.output_filename != NULL) {
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
// READ OFF FILE
// **************************************************

void ISODUAL3D::read_off_file
(const char * input_filename,
 std::vector<COORD_TYPE> & coord, std::vector<GRADIENT_COORD_TYPE> & normal, 
 std::vector<int> & simplex_vert)
{
  int dimension, mesh_dimension;
  IJK::PROCEDURE_ERROR error("read_off_file");

  ifstream in(input_filename, ios::in);
  if (!in.good()) {
    error.AddMessage("Error.  Unable to open file ", input_filename, ".");
    throw error;
  };

  ijkinOFF(in, dimension, mesh_dimension, coord, normal, simplex_vert);

  if (dimension != DIM3) {
    error.AddMessage("Error.  Vertices in OFF file have dimension ",
                     dimension, ".");
    error.AddMessage("  Dimension should be ", DIM3, ".");
    throw error;
  }

  in.close();
}

// Read off file.  Ignore simplex vert.
void ISODUAL3D::read_off_file
(const char * input_filename,
 std::vector<COORD_TYPE> & coord, std::vector<GRADIENT_COORD_TYPE> & normal)
{
  std::vector<int> simplex_vert;

  read_off_file(input_filename, coord, normal, simplex_vert);
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
	if (!output_info.use_stdout && !output_info.flag_silent) {
		report_iso_info3D
      (output_info, isodual_data, dual_isosurface, isodual_info);
	}

	if (!output_info.nowrite_flag) {
    bool flag_reorder_quad_vertices = true;
		write_dual_mesh3D
      (output_info, dual_isosurface, flag_reorder_quad_vertices, io_time);
	}
}


// **************************************************
// WRITE_DUAL_MESH
// **************************************************

// Write dual mesh.
void ISODUAL3D::write_dual_mesh3D
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord, 
 const std::vector<VERTEX_INDEX> & tri_vert,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const bool flag_reorder_quad_vertices)
{
  const string output_filename = output_info.output_filename;
	ofstream output_file;
  IJK::PROCEDURE_ERROR error("write_dual_mesh");

  if (output_filename == "") {
    error.AddMessage("Programming error. Missing output filename.");
    throw error;
  }

  output_file.open(output_filename.c_str(), ios::out);
  if (!output_file.good()) {
    cerr << "Unable to open output file " << output_filename << "." << endl;
    exit(65);
  };

  if (flag_reorder_quad_vertices) {
    std::vector<VERTEX_INDEX> quad_vert2(quad_vert);
    IJK::reorder_quad_vertices(quad_vert2);
    ijkoutOFF(output_file, DIM3, vertex_coord, 
              tri_vert, NUM_VERT_PER_TRI, quad_vert2, NUM_VERT_PER_QUAD);
  }
  else {
    ijkoutOFF(output_file, DIM3, vertex_coord, 
              tri_vert, NUM_VERT_PER_TRI, quad_vert, NUM_VERT_PER_QUAD);
  }

  output_file.close();

	if (!output_info.flag_silent)
		cout << "Wrote output to file: " << output_filename << endl;
}

// Write dual mesh.
// Time write.
void ISODUAL3D::write_dual_mesh3D
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord, 
 const std::vector<VERTEX_INDEX> & tri_vert,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const bool flag_reorder_quad_vertices,
 IO_TIME & io_time)
{
	ELAPSED_TIME wall_time;

	write_dual_mesh3D(output_info, vertex_coord, tri_vert, quad_vert, 
                    flag_reorder_quad_vertices);

	io_time.write_time += wall_time.getElapsed();
}

// Write dual mesh.
// Version with DUAL_ISOSURFACE parameter.
void ISODUAL3D::write_dual_mesh3D
(const OUTPUT_INFO & output_info,
 const DUAL_ISOSURFACE & dual_isosurface,
 const bool flag_reorder_quad_vertices,
 IO_TIME & io_time)
{
  write_dual_mesh3D
    (output_info, dual_isosurface.vertex_coord,
     dual_isosurface.tri_vert, dual_isosurface.quad_vert,
     flag_reorder_quad_vertices, io_time);
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
(const ISODUAL_GRID & full_scalar_grid, const INPUT_INFO & input_info,
		const ISODUAL_DATA & isodual_data)
{
	const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
	const int num_cubes_in_isodual_data =
			isodual_data.ScalarGrid().ComputeNumCubes();

	if (!input_info.use_stdout && !input_info.flag_silent) {

		if (input_info.flag_subsample) {
			// subsampled grid
			cout << num_grid_cubes << " grid cubes.  "
					<< num_cubes_in_isodual_data << " subsampled grid cubes." << endl;
		}
		else if (input_info.flag_supersample) {
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
  const VERTEX_POSITION_METHOD vpos_method = 
    isodual_param.vertex_position_method;

	cout << "Gradient selection cube offset: "
			<< isodual_param.grad_selection_cube_offset << endl;
	cout << "Maximum small eigenvalue: "
			<< isodual_param.max_small_eigenvalue << endl;
	cout << "Max (Linf) distance from cube to isosurface vertex: "
			<< isodual_param.max_dist << endl;
	cout << "Max small gradient magnitude: "
			<< isodual_param.max_small_magnitude << endl;

  cout << "Vertex positioning: ";

  switch (vpos_method) {

  case CUBECENTER:
    cout << "cube_center" << endl;
    break;

  case CENTROID_EDGE_ISO:
    cout << "centroid" << endl;
    break;

  case EDGEI_INTERPOLATE:
    cout << "edgeIinterp" << endl;
    break;

  case EDGEI_GRADIENT:
    cout << "edgeIgrad" << endl;
    break;

  case GRADIENT_POSITIONING:
    {
      string s;

      GRAD_SELECTION_METHOD grad_selection_method =
        isodual_param.GradSelectionMethod();
      get_grad_selection_string(grad_selection_method, s);

      cout << s << endl;
    }
    break;

  default:
    cout << "Unkown" << endl;
    break;
  }

  if (isodual_param.use_sharp_edgeI) {
    cout << "Using sharp formula to calculate intersections of isosurface and grid edges." << endl;
  }
  else {
    cout << "Using interpolation to calculate intersections of isosurface and grid edges." << endl;
  }

  if (isodual_param.use_lindstrom) {
    cout << "Using Lindstrom formula." << endl;
  }

  if (isodual_param.flag_dist2centroid) {
    cout << "Using distance to centroid." << endl;
  }
  else {
    cout << "Using distance to center." << endl;
  }

  if (isodual_param.flag_allow_conflict) {
    cout << "Allow isosurface vertex conflicts." << endl;
  }
  else {
    if (isodual_param.flag_clamp_conflict) 
      { cout << "Resolve conflict by clamping coordinates." << endl; }
    else
      { cout << "Resolve conflict by reverting to centroid." << endl; }
  }

  if (isodual_param.allow_multiple_iso_vertices) {
    cout << "Allow multiple isosurface vertices per grid cube." << endl;

    if (isodual_param.flag_resolve_ambiguous_facets) {
      cout << "Resolve ambiguous facets." << endl;
    }
    else if (isodual_param.flag_separate_neg) {
      cout << "Separate negative vertices." << endl;
    }
    else {
      cout << "Separate positive vertices." << endl;
    }
  }
  else {
    cout << "One isosurface vertex per grid cube." << endl;
  }

	cout << endl;
}


void ISODUAL3D::report_iso_info3D
(const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
 const DUAL_ISOSURFACE & dual_isosurface,
 const ISODUAL_INFO & isodual_info)
{
	const char * indent4 = "    ";
	string grid_element_name = "cubes";

	VERTEX_INDEX numv = (dual_isosurface.vertex_coord.size())/DIM3;
	VERTEX_INDEX num_tri = 
    (dual_isosurface.tri_vert.size())/NUM_VERT_PER_TRI;
	VERTEX_INDEX num_quad = 
    (dual_isosurface.quad_vert.size())/NUM_VERT_PER_QUAD;
	VERTEX_INDEX num_grid_cubes = isodual_info.grid.num_cubes;
	VERTEX_INDEX num_non_empty_cubes = isodual_info.scalar.num_non_empty_cubes;

	float percent = 0.0;
	if (num_grid_cubes > 0)
	{ percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
	int ipercent = int(100*percent);
	cout << "  Isovalue " << output_info.isovalue << ".  "
			<< numv << " isosurface vertices.  ";
  if (num_tri > 0) {
		cout << num_tri << " isosurface triangles." << endl;
	}
  if (num_quad > 0) {
		cout << num_quad << " isosurface quadrilaterals." << endl;
  }
  if (num_tri+num_quad == 0) {
		cout << "No isosurface polygons." << endl;
  }


	if (output_info.flag_output_alg_info) {
    VERTEX_POSITION_METHOD vpos_method = output_info.VertexPositionMethod();

    if (vpos_method == GRADIENT_POSITIONING ||
        vpos_method == EDGEI_INTERPOLATE ||
        vpos_method == EDGEI_GRADIENT) {

      cout << endl;
      if (!output_info.flag_merge_sharp) {
        cout << "  Number of conflicts: "
             << isodual_info.sharpiso.num_conflicts << endl;
      }
      cout << "  Number of sharp corners: "
           << isodual_info.sharpiso.num_sharp_corners << endl;
      cout << "  Number of isosurface vertices on sharp edges: "
           << isodual_info.sharpiso.num_sharp_edges << endl;
      cout << "  Number of smooth isosurface vertices: "
           << isodual_info.sharpiso.num_smooth_vertices << endl;
      if (!output_info.flag_merge_sharp && output_info.use_Linf_dist) {
        cout << "  Number of vertices at min Linf distance to cube center: "
             << isodual_info.sharpiso.num_Linf_iso_vertex_locations << endl;
      }

      if (output_info.flag_merge_conflict) {
        cout << "  Number of edge collapses: " 
             << isodual_info.sharpiso.num_edge_collapses << endl;
      }

      if (output_info.flag_reposition) {
        cout << "  Number of repositioned isosurface vertices: "
             << isodual_info.sharpiso.num_repositioned_vertices << endl;
      }

      if (output_info.flag_merge_sharp) {
        cout << "  Number of merged isosurface vertices: "
             << isodual_info.sharpiso.num_merged_iso_vertices << endl;
      }

      if (output_info.allow_multiple_iso_vertices) {

        if (output_info.flag_resolve_ambiguous_facets) {
          cout << "  Number of non-ambiguous cubes: "
               << isodual_info.sharpiso.num_cube_not_ambiguous << endl;
          cout << "  Number of separate positive cubes: "
               << isodual_info.sharpiso.num_cube_separate_pos << endl;
          cout << "  Number of separate negative cubes: "
               << isodual_info.sharpiso.num_cube_separate_neg << endl;
          cout << "  Number of unresolved ambiguous cubes: "
               << isodual_info.sharpiso.num_cube_unresolved_ambiguity << endl;
        }

        cout << "  Number of cubes with single isov: "
             << isodual_info.sharpiso.num_cube_single_isov << endl;
        cout << "  Number of cubes with multi isov: "
             << isodual_info.sharpiso.num_cube_multi_isov << endl;
      }

      cout << endl;
		}

    if (vpos_method == CENTROID_EDGE_ISO) {
      if (output_info.allow_multiple_iso_vertices) {
        cout << endl;
        cout << "  Number of cubes with single isov: "
             << isodual_info.sharpiso.num_cube_single_isov << endl;
        cout << "  Number of cubes with multi isov: "
             << isodual_info.sharpiso.num_cube_multi_isov << endl;
        cout << endl;
      }
    }
	}
}

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

void ISODUAL3D::report_isodual_time
(const INPUT_INFO & input_info, const ISODUAL_TIME & isodual_time,
		const char * mesh_type_string)
{
	cout << "CPU time to run Marching Cubes: "
			<< isodual_time.total << " seconds." << endl;

  if (input_info.flag_merge_sharp) {

    cout << "    Time to position "
         << mesh_type_string << " vertices: "
         << isodual_time.position << " seconds." << endl;

    cout << "    Time to extract " << mesh_type_string << " triangles: "
         << isodual_time.extract << " seconds." << endl;

    cout << "    Time to merge sharp "
         << mesh_type_string << " vertices: "
         << isodual_time.merge_sharp << " seconds." << endl;
  }
  else {

    cout << "    Time to extract " << mesh_type_string << " triangles: "
         << isodual_time.extract << " seconds." << endl;

    cout << "    Time to merge identical "
         << mesh_type_string << " vertices: "
         << isodual_time.merge_identical << " seconds." << endl;

    cout << "    Time to position "
         << mesh_type_string << " vertices: "
         << isodual_time.position << " seconds." << endl;
  }

}


void ISODUAL3D::report_time
(const INPUT_INFO & input_info, const IO_TIME & io_time,
		const ISODUAL_TIME & isodual_time, const double total_elapsed_time)
{
	const char * ISOSURFACE_STRING = "isosurface";
	const char * INTERVAL_VOLUME_STRING = "interval volume";
	const char * mesh_type_string = NULL;

	mesh_type_string = ISOSURFACE_STRING;

	cout << "Time to read file " << input_info.scalar_filename << ": "
			<< io_time.read_nrrd_time << " seconds." << endl;

	cout << "Time to read " << mesh_type_string << " lookup tables: "
			<< io_time.read_table_time << " seconds." << endl;

	report_isodual_time(input_info, isodual_time, mesh_type_string);
	if (!input_info.nowrite_flag) {
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
	cerr << "  [-subsample S]" << endl;
	cerr << "  [-position {centroid|cube_center|gradC|gradN|gradCS|gradNS|"
			<< endl
			<< "              gradIE|gradIES|gradIEDir|gradCD|gradNIE|gradNIES|"
			<< endl
			<< "              gradES|gradEC}]" << endl;
	cerr << "  [-gradient {gradient_nrrd_filename}]" << endl;
	cerr << "  [-single_isov | -multi_isov | -sep_pos | -sep_neg | -resolve_ambig]"
			<< endl;
	cerr << "  [-max_eigen {max}]" << endl;
	cerr << "  [-max_dist {D}] [-gradS_offset {offset}] [-max_mag {M}] [-snap_dist {D}]" << endl;
  cerr << "  [-sharp_edgeI | -interpolate_edgeI]" << endl;
	cerr << "  [-reposition | -no_reposition] [-sepdist {dist}]" << endl;
	cerr << "  [-lindstrom]" << endl;
	cerr << "  [-allow_conflict |-clamp_conflict | -centroid_conflict] [-merge_conflict]" << endl;
  cerr << "  [-merge_sharp]"<<"[-merge_linf_th]" << endl;
	cerr << "  [-clamp_far] [-centroid_far]" << endl;
	cerr << "  [-recompute_eigen2 | -no_recompute_eigen2]" << endl;
	cerr << "  [-Linf | -no_Linf]" << endl;
	cerr << "  [-removeg | -no_removeg] [-reselectg | -no_reselectg]" << endl;
  cerr << "  [-dist2center | -dist2centroid]" << endl;
	cerr << "  [-centroid_eigen1 | -no_centroid_eigen1]" << endl;
	cerr << "  [-no_round | -round <n>]" << endl;
	cerr << "  [-off|-iv] [-o {output_filename}] [-stdout]"
			<< endl;
	cerr << "  [-help] [-s] [-out_param] [-info] [-nowrite] [-time]" << endl;
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
	cout << "  -position gradIES: Position isosurface vertices using svd"
			<< endl;
	cout << "       on selected gradients at endpoints of intersected cube edges."
			<< endl;
	cout << "  -position gradIEDir: Position isosurface vertices using svd"
			<< endl;
	cout << "       on selected gradients at endpoints of intersected cube edges."
			<< endl;
	cerr << "           Select gradients based on direction along edge."
			<< endl;
	cout << "  -position gradCD: Position isosurface vertices using svd"
			<< endl;
	cout << "       using gradients determining intersections of cube edge"
			<< endl
			<< "       and isosurface." << endl;
	cout << "  -position gradNIE: Position isosurface vertices using svd"
			<< endl;
	cout << "       on gradients at endpoints of intersected cube" << endl;
	cerr << "       and cube neighbor edges." << endl;
	cout << "  -position gradNIES: Position isosurface vertices using svd"
			<< endl;
	cout << "       on selected gradients at endpoints of intersected cube"
			<< endl;
	cerr << "       and cube neighbor edges." << endl;
	cout << "  -position gradES: Position using isosurface vertices using svd"
			<< endl;
	cout << "       on gradients at isosurface-edge intersections."
			<< endl;
	cout << "       Use simple interpolation to compute isosurface-edge intersections." << endl;
	cout << "  -position gradEC: Position isosurface vertices using svd"
			<< endl;
	cout << "       on computed gradients at isosurface-edge intersections."
			<< endl;
	cout << "       Use endpoint gradients to compute isosurface-edge intersections." << endl;
	cout << "  -gradient {gradient_nrrd_filename}: Read gradients from gradient nrrd file." << endl;
	cout << "  -single_isov: Each intersected cube generates a single isosurface vertex." << endl;
	cout << "  -multi_isov:  An intersected cube may generate multiple isosurface vertices."  << endl;
	cout << "  -sep_pos:     Use dual isosurface table separating positive "
			<< "                grid vertices." << endl;
	cout << "  -sep_neg:     Use dual isosurface table separating negative "
			<< "                grid vertices." << endl;
	cout << "  -resolve_ambig:  Selectively resolve ambiguities." << endl
			<< "       Note: Not all ambiguities may be resolved." << endl
			<< "             Unresolved ambiguities may create non-manifold regions."
			<< endl;
	cerr << "  -max_eigen {max}: Set maximum small eigenvalue to max."
			<< endl;
	cerr << "  -max_dist {D}:    Set max Linf distance from cube to isosurface vertex."
			<< endl;
	cerr << "  -max_mag {max}:  Set maximum small gradient magnitude to max."
			<< endl;
	cerr << "           Gradients with magnitude below max are ignored." << endl;
	cerr << "  -snap_dist {D}:  Snap points within distance D to cube."
       << endl;
	cerr << "  -gradS_offset {offset}: Set cube offset for gradient selection to offset."
			<< endl;
	cout << "  -lindstrom:   Use Lindstrom's equation to compute sharp point."
			<< endl;
	cout << "  -allow_conflict:  Allow more than one isosurface vertex in a cube."
			<< endl;
	cout << "  -clamp_conflict:  Settle conflicts by clamping to cube."
			<< " (Default.)"  << endl;
	cout << "  -centroid_conflict:  Settle conflicts by using centroid."
			<< endl;
	cout << "  -merge_conflict:  Settle conflicting isosurface vertices."
       << endl;
	cout << "  -clamp_far: Clamp isosurface vertices at distance greater"
			<< " than max_dist." << endl;
	cout << "  -centroid_far: Revert to centroid when an isosurface vertex is"
			<< endl
			<< "                 at distance greater than max_dist." << endl;
  cout << "  -sharp_edgeI:  Use sharp formula for computing intersections" << endl
       << "                 of isosurface and grid edges." << endl;
  cout << "  -interpolate_edgeI:  Interpolate intersections of isosurface and grid edges."
       << endl;
	cout << "  -recompute_eigen2:  Recompute with only 2 eigenvalues to settle conflicts." << endl;
	cout << "  -no_recompute_eigen2:  Don't recompute with only 2 eigenvalues."
			<< endl;
  cout << "  -dist2center:  Use distance to center in lindstrom." << endl;
  cout << "  -dist2centroid:  Use distance to centroid of isourface-edge"
       << "                   intersections in lindstrom." << endl;
	cout << "  -Linf:     Use Linf metric to resolve conflicts." << endl;
	cout << "  -no_Linf:  Don't use Linf metric to resolve conflicts." << endl;
	cout << "  -removeg:  Remove gradients to resolve conflicts." << endl;
	cout << "  -no_removeg: Don't remove gradients to resolve conflicts." << endl;
	cout << "  -centroid_eigen1: Use centroid with one large eigenvalue." << endl;
	cout << "  -no_centroid_eigen1:  Don't use centroid." << endl;
	cout << "  -no_round:  Don't round coordinates." << endl;
	cout << "  -round <n>: Round coordinates to nearest 1/n." << endl;
	cout << "              Suggest using n=16,32,64,... or 2^k for some k."
			<< endl;
	cout << "  -trimesh:   Output triangle mesh." << endl;
	cout << "  -off: Output in geomview OFF format. (Default.)" << endl;
	cout << "  -iv: Output in OpenInventor .iv format." << endl;
	cout << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
	cout << "  -stdout: Write isosurface to standard output." << endl;
	cout << "  -nowrite: Don't write isosurface." << endl;
	cout << "  -time: Output running time." << endl;
	cout << "  -s: Silent mode." << endl;
	cout << "  -out_param: Print isodual3D parameters." << endl;
	cout << "  -info: Print algorithm information." << endl;
	cout << "  -help: Print this help message." << endl;
	exit(20);
}


// **************************************************
// CLASS IO_INFO
// **************************************************

/// IO information
void ISODUAL3D::IO_INFO::Init()
{
	dimension = 3;
	scalar_filename = NULL;
	gradient_filename = NULL;
	output_filename = NULL;
	output_format = OFF;
	report_time_flag = false;
	use_stdout = false;
	nowrite_flag = false;
	flag_output_alg_info = false;
	flag_silent = false;
	flag_subsample = false;
	subsample_resolution = 2;
	flag_supersample = false;
	supersample_resolution = 2;
	flag_color_alternating = false;  // color simplices in alternating cubes
	region_length = 1;
	max_small_eigenvalue = 0.1;
	flag_output_param = false;

	grid_spacing.resize(3,1);
}

void ISODUAL3D::IO_INFO::Set(const IO_INFO & io_info)
{
	*this = io_info;
}

// **************************************************
// CLASS INPUT_INFO
// **************************************************

/// Clear input information
void ISODUAL3D::INPUT_INFO::Init()
{
	Clear();
}

/// Clear input information
void ISODUAL3D::INPUT_INFO::Clear()
{
	isovalue.clear();
	isovalue_string.clear();
	isotable_directory = "";
}

// **************************************************
// class OUTPUT_INFO
// **************************************************

void ISODUAL3D::OUTPUT_INFO::Init()
{
	output_filename = "";
	isovalue = 0;
	grow_factor = 1;
	shrink_factor = 1;

	SetOutputTriMesh(false);
}

namespace {

string construct_output_filename
(const INPUT_INFO & input_info, const int i)
{
	string prefix, suffix;

	// create output filename
	string fname = string(input_info.scalar_filename);

#ifndef _WIN32
	// remove path from file name
	split_string(fname, PATH_DELIMITER, prefix, suffix);
	if (suffix != "") { fname = suffix; }
#endif

	string ofilename;

	// construct output filename
	split_string(fname, '.', prefix, suffix);
	if (suffix == "nrrd" || suffix == "nhdr") { ofilename = prefix; }
	else { ofilename = string(input_info.scalar_filename); }

	ofilename += string(".") + string("isov=") + input_info.isovalue_string[i];

	switch (input_info.output_format) {
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
(const INPUT_INFO & input_info, ISODUAL_DATA & isodual_data, ISODUAL_TIME & isodual_time)
{
	PROCEDURE_ERROR error("set_isodual_data");

	if (!isodual_data.IsScalarGridSet()) {
		error.AddMessage
		("Programming error. Scalar field must be set before set_isodual_data is called.");
		throw error;
	}

	// Set data structures in isodual_data
	isodual_data.Set(input_info);
}

void ISODUAL3D::set_input_info
(const NRRD_INFO & nrrd_info, INPUT_INFO & input_info)
{
	input_info.grid_spacing.clear();
	for (int d = 0; d < nrrd_info.dimension; d++) {
		input_info.grid_spacing.push_back(nrrd_info.grid_spacing[d]);
	}
}

void ISODUAL3D::set_output_info
(const INPUT_INFO & input_info,
		const int i, OUTPUT_INFO & output_info)
{
	// Set data structures in isodual_data
	output_info.IO_INFO::Set(input_info);

	output_info.SetOutputTriMesh(input_info.flag_convert_quad_to_tri);

	output_info.grow_factor = 1;
	if (input_info.flag_subsample)
	{ output_info.grow_factor = input_info.subsample_resolution; }

	output_info.shrink_factor = 1;
	if (input_info.flag_supersample)
	{ output_info.shrink_factor = input_info.supersample_resolution; }

	output_info.isovalue = input_info.isovalue[i];

	if (input_info.output_filename != NULL) {
		output_info.output_filename = string(input_info.output_filename);
	}
	else {
		output_info.output_filename =
				construct_output_filename(input_info, i);
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

// **************************************************
// Class OUTPUT_INFO member functions
// **************************************************

void ISODUAL3D::OUTPUT_INFO::SetOutputTriMesh(const bool flag)
{
	flag_output_tri_mesh = flag;

	if (flag_output_tri_mesh)
	{ num_vertices_per_isopoly = 3; }
	else
	{ num_vertices_per_isopoly = 4; }
}

bool ISODUAL3D::OUTPUT_INFO::OutputTriMesh() const
{
	return(flag_output_tri_mesh);
}

NUM_TYPE ISODUAL3D::OUTPUT_INFO::NumVerticesPerIsopoly() const
{
	return(num_vertices_per_isopoly);
}


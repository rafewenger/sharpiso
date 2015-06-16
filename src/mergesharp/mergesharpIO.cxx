/// \file mergesharpIO.cxx
/// IO routines for mergesharp

/*
  Copyright (C) 2011-2015 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3 of the License, or any later version.

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

#include "ijkcoord.txx"
#include "ijkgrid_nrrd.txx"
#include "ijkIO.txx"
#include "ijkmesh.txx"
#include "ijkprint.txx"
#include "ijkstring.txx"

#include "sharpiso_array.txx"
#include "mergesharpIO.h"
#include "sharpiso_get_gradients.h"

using namespace IJK;
using namespace MERGESHARP;

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
    GRAD2HERMITE_PARAM, GRAD2HERMITE_INTERPOLATE_PARAM,
    MAX_EIGEN_PARAM, MAX_DIST_PARAM, 
    GRAD_S_OFFSET_PARAM, MIN_GRAD_S_OFFSET_PARAM,
    MAX_MAG_PARAM, SNAP_DIST_PARAM, MAX_GRAD_DIST_PARAM,
    MIN_TRIANGLE_ANGLE_PARAM, MIN_NORMAL_ANGLE_PARAM,
    SHARP_EDGEI_PARAM, INTERPOLATE_EDGEI_PARAM,
    ALLOW_CONFLICT_PARAM,
    CLAMP_CONFLICT_PARAM, CENTROID_CONFLICT_PARAM,
    MERGE_SHARP_PARAM, NO_MERGE_SHARP_PARAM, 
    MERGE_SHARP_LINF_THRES_PARAM,
    CLAMP_FAR_PARAM, CENTROID_FAR_PARAM,
    RECOMPUTE_ISOVERT, NO_RECOMPUTE_ISOVERT,
    RECOMPUTE_USING_ADJACENT, NO_RECOMPUTE_USING_ADJACENT,
    CHECK_TRIANGLE_ANGLE, NO_CHECK_TRIANGLE_ANGLE,
    DIST2CENTER_PARAM, DIST2CENTROID_PARAM,
    LINF_PARAM, NO_LINF_PARAM,

    // DEPRECATED
    USE_LINDSTROM_PARAM,
    USE_LINDSTROM2_PARAM,
    USE_LINDSTROM_FAST,
    NO_LINDSTROM_PARAM,

    SINGLE_ISOV_PARAM, MULTI_ISOV_PARAM,
    SPLIT_NON_MANIFOLD_PARAM, SELECT_SPLIT_PARAM,
    SEP_NEG_PARAM, SEP_POS_PARAM, RESOLVE_AMBIG_PARAM,
    CHECK_DISK_PARAM, NO_CHECK_DISK_PARAM,
    MANIFOLD_PARAM,
    ROUND_PARAM, NO_ROUND_PARAM,
    KEEPV_PARAM,
    MINC_PARAM, MAXC_PARAM,
    MAP_EXTENDED_PARAM, NO_MAP_EXTENDED_PARAM,
    COLLAPSE_TRIANGLES_PARAM, NO_COLLAPSE_TRIANGLES_PARAM,
    SELECT_MOD3_PARAM, SELECT_MOD6_PARAM,
    HELP_PARAM, OFF_PARAM, IV_PARAM, 
    OUTPUT_FILENAME_PARAM, STDOUT_PARAM, NOWRITE_PARAM, 
    OUTPUT_PARAM_PARAM, OUTPUT_INFO_PARAM, 
    OUTPUT_SELECTED_PARAM, OUTPUT_SHARP_PARAM, OUTPUT_ACTIVE_PARAM,
    OUTPUT_MAP_TO_SELF_PARAM, 
    OUTPUT_MAP_TO_PARAM, OUTPUT_NEIGHBORS_PARAM,
    OUTPUT_ISOVERT_PARAM,
    WRITE_ISOV_INFO_PARAM, SILENT_PARAM, TIME_PARAM, 
    UNKNOWN_PARAM} PARAMETER;
  const char * parameter_string[] =
    { "-subsample",
      "-gradient", "-normal", "-position", "-pos", 
      "-trimesh", "-uniform_trimesh",
      "-grad2hermite", "-grad2hermiteI",
      "-max_eigen", "-max_dist", "-gradS_offset", "-min_gradS_offset", 
      "-max_mag", "-snap_dist", "-max_grad_dist",
      "-min_triangle_angle", "-min_normal_angle",
      "-sharp_edgeI", "-interpolate_edgeI",
      "-allow_conflict", "-clamp_conflict", "-centroid_conflict", 
      "-merge_sharp","-no_merge_sharp", "-merge_linf_th",
      "-clamp_far", "-centroid_far",
      "-recompute_isovert", "-no_recompute_isovert",
      "-recompute_using_adjacent", "-no_recompute_using_adjacent",
      "-check_triangle_angle", "-no_check_triangle_angle",
      "-dist2center", "-dist2centroid",
      "-Linf", "-no_Linf",

      // DEPRECATED
      "-lindstrom", "-lindstrom2","-lindstrom_fast", "-no_lindstrom",

      "-single_isov", "-multi_isov", "-split_non_manifold", "-select_split",
      "-sep_neg", "-sep_pos", "-resolve_ambig", 
      "-check_disk", "-no_check_disk",
      "-manifold",
      "-round", "-no_round",
      "-keepv",
      "-minc", "-maxc",
      "-map_extended", "-no_map_extended",
      "-collapse_triangles", "-no_collapse_triangles",
      "-select_mod3", "-select_mod6",
      "-help", "-off", "-iv", 
      "-o", "-stdout", "-nowrite", 
      "-out_param", "-info", "-out_selected", "-out_sharp", "-out_active",
      "-out_map_to_self", "-out_map_to", "-out_neighbors", "-out_isovert",
      "-write_isov_info", "-s", "-time", "-unknown"};

  PARAMETER get_parameter_token(const char * s)
  // convert string s into parameter token
  {
    string str(s);

    for (int i = 0; i < int(UNKNOWN_PARAM); i++)
      if (str == parameter_string[i])
        { return(PARAMETER(i)); }
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

  int get_option_int
  (const char * option, const char * value_string)
  {
    int x;
    if (!IJK::string2val(value_string, x)) {
      cerr << "Error in argument for option: " << option << endl;
      cerr << "Non-integer character in string: " << value_string << endl;
      exit(50);
    }

    return(x);
  }

  float get_option_float
  (const char * option, const char * value_string)
  {
    float x;
    if (!IJK::string2val(value_string, x)) {
      cerr << "Error in argument for option: " << option << endl;
      cerr << "Non-numertic character in string: " << value_string << endl;
      exit(50);
    }

    return(x);
  }

  /// Get string and convert to list of arguments.
  template <typename ETYPE>
  void get_option_multiple_arguments
  (const char * option, const char * value_string, std::vector<ETYPE> & v)
  {
    if (!IJK::string2vector(value_string, v)) {
      cerr << "Error in argument for option: " << option << endl;
      cerr << "Non-numeric character in string: " << value_string << endl;
      exit(50);
    }
  }

  void get_isovert_type
  (const char * option, const char * value0, const char * value1,
   ISOVERT_TYPE & isovert_type,
   MSHARP_ISOVERT_TYPE & isovert_type2)
  {
    string s0 = value0;
    if (s0 == "corner") {
      isovert_type = CORNER_ISOVERT;
    }
    else if (s0 == "edge") {
      isovert_type = EDGE_ISOVERT;
    }
    else if (s0 == "sharp") {
      isovert_type = SHARP_ISOVERT;
    }
    else if (s0 == "smooth") {
      isovert_type = SMOOTH_ISOVERT;
    }
    else if (s0 == "all") {
      isovert_type = ALL_ISOVERT;
    }
    else {
      cerr << "Usage error.  Illegal first argument " << value0
           << " to option " << option << ".";
      exit(220);
    }

    string s1 = value1;
    if (s1 == "selected") {
      isovert_type2 = SELECTED_ISOVERT;
    }
    else if (s1 == "uncovered") {
      isovert_type2 = UNCOVERED_ISOVERT;
    }
    else if (s1 == "all") {
      isovert_type2 = ACTIVE_ISOVERT;
    }
    else {
      cerr << "Usage error.  Illegal second argument " << value1
           << " to option " << option << ".";
      exit(221);
    }
  }

  void lindstrom_deprecated()
  {
    cerr << "*** Warning: Options -lindstrom -lindstrom2 -lindstrom_fast -no_lindstrom are all deprecated." << endl;
    cerr << "   Program mergesharp now always uses creates 3x3 matrix from the gradients"
         << "   and applies lindstrom's algorithm to this matrix."
         << endl;
  }

  // Set output neighbors flags.
  void set_output_neighbors_flags
  (const char * option, const char * value_string, INPUT_INFO & input_info)
  {
    std::vector<std::string> distance_string;

    input_info.flag_output_neighbors = true;

    IJK::string2vector(value_string, distance_string);

    for (int i = 0; i < distance_string.size(); i++) {

      if (distance_string[i] == "222") {
        input_info.flag_output_neighbors_222 = true;
        input_info.flag_output_neighbors_322 = true;
        input_info.flag_output_neighbors_332 = true;
      }
      else if (distance_string[i] == "22x") {
        input_info.flag_output_neighbors_222 = true;
        input_info.flag_output_neighbors_221 = true;
        input_info.flag_output_neighbors_220 = true;
        input_info.flag_output_neighbors_322 = true;
        input_info.flag_output_neighbors_321 = true;
        input_info.flag_output_neighbors_320 = true;
        input_info.flag_output_neighbors_332 = true;
        input_info.flag_output_neighbors_331 = true;
        input_info.flag_output_neighbors_330 = true;
      }
      else if (distance_string[i] == "321") {
        input_info.flag_output_neighbors_322 = true;
        input_info.flag_output_neighbors_321 = true;
        input_info.flag_output_neighbors_332 = true;
        input_info.flag_output_neighbors_331 = true;
      }
      else if (distance_string[i] == "32x") {
        input_info.flag_output_neighbors_322 = true;
        input_info.flag_output_neighbors_321 = true;
        input_info.flag_output_neighbors_320 = true;
        input_info.flag_output_neighbors_332 = true;
        input_info.flag_output_neighbors_331 = true;
        input_info.flag_output_neighbors_330 = true;
      }
      else if (distance_string[i] == "33x") {
        input_info.flag_output_neighbors_332 = true;
        input_info.flag_output_neighbors_331 = true;
        input_info.flag_output_neighbors_330 = true;
      }
      else {
        cerr << "Usage error.  Illegal argument " << value_string
             << " to option " << option << "." << endl;
        exit(227);
      }
    }
  }

  // Set flag in input_info based on param.
  // Return false if param is not valid or not a flag.
  bool set_input_info_flag
  (const PARAMETER param, INPUT_INFO & input_info)
  {
    switch(param) {

    case TRIMESH_PARAM:
      input_info.flag_convert_quad_to_tri = true;
      input_info.quad_tri_method = SPLIT_MAX_ANGLE;
      break;

    case UNIFORM_TRIMESH_PARAM:
      input_info.flag_convert_quad_to_tri = true;
      input_info.quad_tri_method = UNIFORM_TRI;
      break;

    case GRAD2HERMITE_PARAM:
      input_info.flag_grad2hermite = true;
      input_info.vertex_position_method = EDGEI_GRADIENT;
      input_info.is_vertex_position_method_set = true;
      break;

    case GRAD2HERMITE_INTERPOLATE_PARAM:
      input_info.vertex_position_method = EDGEI_INTERPOLATE;
      input_info.flag_grad2hermiteI = true;
      break;

    case SHARP_EDGEI_PARAM:
      input_info.use_sharp_edgeI = true;
      input_info.is_use_sharp_edgeI_set = true;
      break;

    case INTERPOLATE_EDGEI_PARAM:
      input_info.use_sharp_edgeI = false;
      input_info.is_use_sharp_edgeI_set = true;
      break;

    case LINF_PARAM:
      input_info.use_Linf_dist = true;
      break;

    case NO_LINF_PARAM:
      input_info.use_Linf_dist = false;
      break;

    case USE_LINDSTROM_PARAM:
      input_info.use_lindstrom =true;
      lindstrom_deprecated();
      break;

    case USE_LINDSTROM2_PARAM:
      input_info.use_lindstrom = true;
      input_info.use_lindstrom2 = true;
      lindstrom_deprecated();
      break;

    case USE_LINDSTROM_FAST:
      input_info.use_lindstrom = true;
      input_info.use_lindstrom_fast = true;
      lindstrom_deprecated();
      break;

    case NO_LINDSTROM_PARAM:
      input_info.use_lindstrom = false;
      input_info.use_lindstrom2 = false;
      input_info.use_lindstrom_fast = false;
      lindstrom_deprecated();
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

    case SPLIT_NON_MANIFOLD_PARAM:
      input_info.allow_multiple_iso_vertices = true;
      input_info.flag_split_non_manifold = true;
      break;

    case SELECT_SPLIT_PARAM:
      input_info.allow_multiple_iso_vertices = true;
      input_info.flag_select_split = true;
      break;

    case ALLOW_CONFLICT_PARAM:
      input_info.flag_allow_conflict = true;
      input_info.is_conflict_set = true;
      break;

    case CLAMP_CONFLICT_PARAM:
      input_info.flag_clamp_conflict = true;
      input_info.is_conflict_set = true;
      break;

    case CENTROID_CONFLICT_PARAM:
      input_info.flag_clamp_conflict = false;
      input_info.is_conflict_set = true;
      break;

    case MERGE_SHARP_PARAM:
      input_info.flag_merge_sharp = true;
      break;

    case NO_MERGE_SHARP_PARAM:
      input_info.flag_merge_sharp = false;
      break;

    case CLAMP_FAR_PARAM:
      input_info.flag_clamp_far = true;
      break;

    case CENTROID_FAR_PARAM:
      input_info.flag_clamp_far = false;
      break;

    case RECOMPUTE_ISOVERT:
      input_info.flag_recompute_isovert = true;
      break;

    case NO_RECOMPUTE_ISOVERT:
      input_info.flag_recompute_isovert = false;
      break;

    case RECOMPUTE_USING_ADJACENT:
      input_info.flag_recompute_using_adjacent = true;
      input_info.flag_recompute_changing_gradS_offset = false;
      break;

    case NO_RECOMPUTE_USING_ADJACENT:
      input_info.flag_recompute_using_adjacent = false;
      break;

    case CHECK_TRIANGLE_ANGLE:
      input_info.flag_check_triangle_angle = true;
      break;

    case NO_CHECK_TRIANGLE_ANGLE:
      input_info.flag_check_triangle_angle = false;
      break;

    case MAP_EXTENDED_PARAM:
      input_info.flag_map_extended = true;
      break;

    case NO_MAP_EXTENDED_PARAM:
      input_info.flag_map_extended = false;
      break;

    case COLLAPSE_TRIANGLES_PARAM:
      input_info.flag_collapse_triangles_with_small_angles = true;
      break;

    case NO_COLLAPSE_TRIANGLES_PARAM:
      input_info.flag_collapse_triangles_with_small_angles = false;
      break;

    case SELECT_MOD6_PARAM:
      input_info.flag_select_mod6 = true;
      break;

    case DIST2CENTER_PARAM:
      input_info.flag_dist2centroid = false;
      break;

    case DIST2CENTROID_PARAM:
      input_info.flag_dist2centroid = true;
      break;

    case NO_ROUND_PARAM:
      input_info.flag_round = false;
      break;

    case CHECK_DISK_PARAM:
      input_info.flag_check_disk = true;
      break;

    case NO_CHECK_DISK_PARAM:
      input_info.flag_check_disk = false;
      break;

    case MANIFOLD_PARAM:
      input_info.allow_multiple_iso_vertices = true;
      input_info.flag_split_non_manifold = true;
      input_info.flag_check_disk = true;
      break;

    case OFF_PARAM:
      input_info.output_format = OFF;
      break;

    case IV_PARAM:
      input_info.output_format = IV;
      break;

    case KEEPV_PARAM:
      input_info.flag_delete_isolated_vertices = false;
      break;

    case OUTPUT_PARAM_PARAM:
      input_info.flag_output_param = true;
      break;

    case OUTPUT_SELECTED_PARAM:
      input_info.flag_output_selected = true;
      break;

    case OUTPUT_SHARP_PARAM:
      input_info.flag_output_sharp = true;
      break;

    case OUTPUT_ACTIVE_PARAM:
      input_info.flag_output_active = true;
      break;

    case OUTPUT_MAP_TO_SELF_PARAM:
      input_info.flag_output_map_to_self = true;
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

    case WRITE_ISOV_INFO_PARAM:
      input_info.flag_store_isovert_info = true;
      break;

    case SILENT_PARAM:
      input_info.flag_silent = true;
      break;

    case TIME_PARAM:
      input_info.report_time_flag = true;
      break;

    default:
      return(false);
    }

    return(true);
  }

  // Set value in input_info based on param.
  // Return false if param is not valid or does not take a value.
  bool set_input_info_value
  (const PARAMETER param, const char * option_string, 
   const char * value_string, INPUT_INFO & input_info)
  {
    ANGLE_TYPE degrees;

    switch(param) {

    case SUBSAMPLE_PARAM:
      input_info.subsample_resolution =  
        get_option_int(option_string, value_string);
      input_info.flag_subsample = true;
      break;

    case GRADIENT_PARAM:
      input_info.gradient_filename = value_string;
      break;

    case NORMAL_PARAM:
      input_info.normal_filename = value_string;
      input_info.vertex_position_method = EDGEI_INPUT_DATA;
      input_info.is_vertex_position_method_set = true;
      break;

    case POSITION_PARAM:
    case POS_PARAM:
      set_vertex_position_method(value_string, input_info);

      input_info.is_vertex_position_method_set = true;
      break;

    case MAX_EIGEN_PARAM:
      input_info.max_small_eigenvalue = 
        get_option_float(option_string, value_string);
      break;

    case MAX_DIST_PARAM:
      input_info.max_dist = 
        get_option_float(option_string, value_string);
      break;

    case MAX_MAG_PARAM:
      input_info.max_small_magnitude = 
        get_option_float(option_string, value_string);
      break;

    case MIN_TRIANGLE_ANGLE_PARAM:
      degrees =
        get_option_float(option_string, value_string);
      input_info.SetMinTriangleAngle(degrees);
      break;

    case MIN_NORMAL_ANGLE_PARAM:
      degrees = 
        get_option_float(option_string, value_string);
      input_info.SetMinNormalAngle(degrees);
      break;

    case SNAP_DIST_PARAM:
      input_info.snap_dist = 
        get_option_float(option_string, value_string);
      break;

    case MAX_GRAD_DIST_PARAM:
      input_info.max_grad_dist =
        get_option_int(option_string, value_string);
      if (input_info.max_grad_dist > 1) 
        { input_info.use_large_neighborhood = true; }
      break;

    case GRAD_S_OFFSET_PARAM:
      input_info.grad_selection_cube_offset = 
        get_option_float(option_string, value_string);
      break;

    case MIN_GRAD_S_OFFSET_PARAM:
      input_info.min_grad_selection_cube_offset = 
        get_option_float(option_string, value_string);
      break;

    case MERGE_SHARP_LINF_THRES_PARAM:
      input_info.linf_dist_thresh_merge_sharp = 
        get_option_float(option_string, value_string);
      input_info.flag_merge_sharp = true;
      break;

    case ROUND_PARAM:
      input_info.flag_round = true;
      input_info.round_denominator =
        get_option_int(option_string, value_string);
      break;

    case MINC_PARAM:
      get_option_multiple_arguments
        (option_string, value_string, input_info.minc);
      break;

    case MAXC_PARAM:
      get_option_multiple_arguments
        (option_string, value_string, input_info.maxc);
	  break;

    case OUTPUT_MAP_TO_PARAM:
      input_info.flag_output_map_to = true;
      input_info.to_cube = 
        get_option_float(option_string, value_string);
      break;

    case OUTPUT_NEIGHBORS_PARAM:
      set_output_neighbors_flags(option_string, value_string, input_info);
      break;

    case OUTPUT_FILENAME_PARAM:
      input_info.output_filename = value_string;
      break;

    default:
      return(false);
    }

    return(true);
  }

  // Set value in input_info based on param.
  // Return false if param is not valid or does not take a value.
  bool set_input_info_value
  (const PARAMETER param, const char * option_string, 
   const char * value_string0, const char * value_string1,
   INPUT_INFO & input_info)
  {
    switch(param) {

    case OUTPUT_ISOVERT_PARAM:
      
      input_info.flag_output_isovert = true;
      get_isovert_type(option_string, value_string0, value_string1,
                       input_info.output_isovert_type,
                       input_info.output_isovert_type2);

      // Default output filename
      input_info.output_isovert_filename = "isovert.off";

      break;

    default:
      return(false);
    }

    return(true);
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

  string remove_off_suffix(const string & filename)
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
    if (suffix == "off") { fname = prefix; }
    else { fname = filename; }

    return(fname);
  }

}

// Parse the next option in the command line.
// Return false if no next option or parse fails.
bool MERGESHARP::parse_command_option
(const int argc, char **argv, const int iarg, int & next_arg,
 INPUT_INFO & input_info)
{
  next_arg = iarg;

  if (iarg >= argc) { return(false); }

  if (argv[iarg][0] != '-') { return(false); }

  PARAMETER param = get_parameter_token(argv[iarg]);

  if (param == UNKNOWN_PARAM) { return(false); }

  if (set_input_info_flag(param, input_info)) {
    next_arg = iarg+1;
    return(true);
  }

  if (iarg+1 >= argc) { return(false); }

  if (set_input_info_value(param, argv[iarg], argv[iarg+1], input_info)) {
    next_arg = iarg+2;
    return(true);
  }

  if (iarg+2 >= argc) { return(false); }

  if (set_input_info_value
      (param, argv[iarg], argv[iarg+1], argv[iarg+2], input_info)) {
    next_arg = iarg+3;
    return(true);
  }

  return(false);
}

// Parse the isovalue(s) and filename.
void MERGESHARP::parse_isovalue_and_filename
(const int argc, char **argv, const int iarg, INPUT_INFO & input_info)
{
  // remaining parameters should be list of isovalues followed
  // by input file name

  // check for more parameter tokens
  for (int j = iarg; j < argc; j++) {
    if (get_parameter_token(argv[j]) != UNKNOWN_PARAM) {
      // argv[iarg] is not an isovalue
      cerr << "Error. Illegal parameter: " << argv[iarg] << endl;
      usage_error(argv[0]);
    }
  }

  if (iarg+2 > argc) {
    cerr << "Error.  Missing input isovalue or input file name." << endl;
    usage_error(argv[0]);
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
      usage_error(argv[0]);
    };

    input_info.isovalue.push_back(value);
  }

  input_info.scalar_filename = argv[argc-1];
}

// Set input_info defaults.
void MERGESHARP::set_input_info_defaults(INPUT_INFO & input_info)
{
  if (!input_info.is_vertex_position_method_set) {
    input_info.vertex_position_method = GRADIENT_POSITIONING;
    input_info.SetGradSelectionMethod(GRAD_NS);
  }

  if (!input_info.is_conflict_set && input_info.flag_merge_sharp) {
    // Set merge_sharp defaults.
    input_info.flag_allow_conflict = true;
    input_info.flag_clamp_conflict = false;
  }

  if (input_info.vertex_position_method == EDGEI_GRADIENT) {
    input_info.use_sharp_edgeI = true;
  }

  if (input_info.vertex_position_method == EDGEI_INTERPOLATE) {
    input_info.use_sharp_edgeI = false;
  }
}

// Check input_info.
void MERGESHARP::check_input_info(const INPUT_INFO & input_info)
{
  if (input_info.is_use_sharp_edgeI_set) {
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

// Parse the command line.
void MERGESHARP::parse_command_line
(int argc, char **argv, INPUT_INFO & input_info)
{
  if (argc == 1) { usage_error(argv[0]); };

  int iarg = 1;
  while (iarg < argc) {

    int next_arg;
    if (!parse_command_option(argc, argv, iarg, next_arg, input_info) ||
        iarg == next_arg) { 

      if (argv[iarg][0] == '-') {
        PARAMETER param = get_parameter_token(argv[iarg]);
        if (param == HELP_PARAM) { help(argv[0]); }
      }
      break; 
    }

    iarg = next_arg;
  }

  parse_isovalue_and_filename(argc, argv, iarg, input_info);
  set_input_info_defaults(input_info);
  check_input_info(input_info);
}

// Check input information/flags.
bool MERGESHARP::check_input
(const INPUT_INFO & input_info,
 const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
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

void MERGESHARP::read_nrrd_file
(const char * input_filename, SHARPISO_SCALAR_GRID & scalar_grid,
 NRRD_INFO & nrrd_info)
{
  IJK::PROCEDURE_ERROR error("read_nrrd_file");

  IJK::GRID_NRRD_IN<int, int> nrrd_in;
  IJK::NRRD_DATA<int,int> nrrd_header;

  nrrd_in.ReadScalarGrid(input_filename, scalar_grid, nrrd_header, error);
  if (nrrd_in.ReadFailed()) { throw error; }

  if (scalar_grid.Dimension() < 1) {
    cerr << "Illegal scalar grid dimension.  Dimension must be at least 1." 
         << endl;
    exit(20);
  };

  std::vector<COORD_TYPE> grid_spacing;
  nrrd_header.GetSpacing(grid_spacing);

  nrrd_info.dimension = scalar_grid.Dimension();
  for (int d = 0; d < scalar_grid.Dimension(); d++) {
    nrrd_info.grid_spacing.push_back(grid_spacing[d]); 
    scalar_grid.SetSpacing(d, grid_spacing[d]);
  };
  
}

void MERGESHARP::read_nrrd_file
(const char * input_filename, GRADIENT_GRID & gradient_grid,
 NRRD_INFO & nrrd_info)
{
  IJK::PROCEDURE_ERROR error("read_nrrd_file");

  GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in_gradient;
  NRRD_DATA<int,AXIS_SIZE_TYPE> nrrd_header;

  nrrd_in_gradient.ReadVectorGrid
    (input_filename, gradient_grid, nrrd_header, error);
  if (nrrd_in_gradient.ReadFailed()) { throw error; }

  if (gradient_grid.Dimension() < 1) {
    cerr << "Illegal gradient grid dimension.  Dimension must be at least 1." 
         << endl;
    exit(20);
  };

  nrrd_info.dimension = gradient_grid.Dimension();

  std::vector<COORD_TYPE> grid_spacing;
  nrrd_header.GetSpacing(grid_spacing);

  nrrd_info.dimension = gradient_grid.Dimension();
  for (int d = 0; d < gradient_grid.Dimension(); d++) {
    nrrd_info.grid_spacing.push_back(grid_spacing[d+1]); 
    gradient_grid.SetSpacing(d, grid_spacing[d+1]);
  };
}

void MERGESHARP::read_nrrd_file
(const char * input_filename, SHARPISO_SCALAR_GRID & scalar_grid,
 NRRD_INFO & nrrd_info, IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  read_nrrd_file(input_filename, scalar_grid, nrrd_info);
  io_time.read_nrrd_time = wall_time.getElapsed();
}

// **************************************************
// READ OFF FILE
// **************************************************

void MERGESHARP::read_off_file
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
void MERGESHARP::read_off_file
(const char * input_filename,
 std::vector<COORD_TYPE> & coord, std::vector<GRADIENT_COORD_TYPE> & normal)
{
  std::vector<int> simplex_vert;

  read_off_file(input_filename, coord, normal, simplex_vert);
}


// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

void report_far_neighbors
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert);


void MERGESHARP::output_dual_isosurface
(const OUTPUT_INFO & output_info, const MERGESHARP_DATA & mergesharp_data,
 const DUAL_ISOSURFACE & dual_isosurface, const ISOVERT & isovert, 
 const MERGESHARP_INFO & mergesharp_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info3D
      (output_info, mergesharp_data, dual_isosurface, mergesharp_info);

    report_isovert_info(output_info, mergesharp_data.ScalarGrid(), isovert);

    if (output_info.flag_output_neighbors) {
      report_far_neighbors
        (output_info, mergesharp_data.ScalarGrid(), isovert);
    }
  }

  if (!output_info.nowrite_flag) {
    bool flag_reorder_quad_vertices = true;
    write_dual_mesh3D
      (output_info, dual_isosurface, flag_reorder_quad_vertices, io_time);

    if (output_info.flag_store_isovert_info) {
      write_isovert_info(output_info, mergesharp_info.sharpiso.vertex_info);
	}
  }

  if (output_info.flag_output_isovert) 
    { write_isovert(output_info, isovert); }
}


// **************************************************
// WRITE_DUAL_MESH
// **************************************************

// Write dual mesh.
void MERGESHARP::write_dual_mesh3D
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
void MERGESHARP::write_dual_mesh3D
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
void MERGESHARP::write_dual_mesh3D
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
void MERGESHARP::rescale_vertex_coord
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
    error.AddMessage
      ("Size of grid spacing array does not equal volume dimension.");
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
void MERGESHARP::rescale_vertex_coord
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
// WRITE ISOSURFACE VERTICES
// **************************************************

void add_coord_to_vector
(const COORD_TYPE coord[DIM3], std::vector<COORD_TYPE> & vcoord)
{
  NUM_TYPE k = vcoord.size();
  vcoord.resize(k+DIM3);
  copy_coord_3D(coord, &(vcoord.front())+k);
}

void MERGESHARP::write_isovert
(const OUTPUT_INFO & output_info, const ISOVERT & isovert)
{
  vector<COORD_TYPE> vertex_coord;
  vector<VERTEX_INDEX> simplex_vert;
  ofstream output_file;
  IJK::PROCEDURE_ERROR error("write_isovert");

  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

    if  (output_info.output_isovert_type2 == UNCOVERED_ISOVERT) {
      if ((isovert.NumEigenvalues(i) <= 1) &&
          (isovert.gcube_list[i].maps_to_cube == isovert.CubeIndex(i))) {
        add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord);
      }
    }
    else {
      if ((output_info.output_isovert_type2 == SELECTED_ISOVERT) &&
          (isovert.gcube_list[i].flag != SELECTED_GCUBE))
        { continue; }

      switch (output_info.output_isovert_type) {

      case CORNER_ISOVERT:
        if (isovert.NumEigenvalues(i) == 3)
          { add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord); }
        break;

      case EDGE_ISOVERT:
        if (isovert.NumEigenvalues(i) == 2)
          { add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord); }
        break;

      case SHARP_ISOVERT:
        if (isovert.NumEigenvalues(i) > 1)
          { add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord); }
        break;

      case SMOOTH_ISOVERT:
        if (isovert.NumEigenvalues(i) == 1)
          { add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord); }
        break;

      case ALL_ISOVERT:
        add_coord_to_vector(isovert.IsoVertCoord(i), vertex_coord);
      }
    }
  }

  if (output_info.output_isovert_filename == "") {
	  
    error.AddMessage("Programming error. Missing output filename.");
    throw error;
  }

  output_file.open(output_info.output_isovert_filename.c_str(), ios::out);
  if (!output_file.good()) {
    cerr << "Unable to open output file " 
         << output_info.output_isovert_filename << "." << endl;
    exit(65);
  };

  ijkoutOFF(output_file, DIM3, vertex_coord, simplex_vert);

  output_file.close();

  if (!output_info.flag_silent)
    cout << "Wrote isosurface vertices to file: " 
         << output_info.output_isovert_filename << endl;
}


// **************************************************
// ANALYZE ADJACENT
// **************************************************


bool find_corner(const ISOVERT & isovert, const vector<VERTEX_INDEX> & list,
                 VERTEX_INDEX & icubeB)
{
  icubeB = 0;
  for (int i = 0; i < list.size(); i++) {
    INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(list[i]);
    if (gcube_index == ISOVERT::NO_INDEX) { continue; }

    if (isovert.NumEigenvalues(gcube_index) == 3) {
      icubeB = list[i];
      return(true);
    }
  }

  return(false);
}

// Compute distance along each axis and then sort by decreasing magnitude.
void compute_dist
(const ISOVERT & isovert, const VERTEX_INDEX & icubeA, 
 const VERTEX_INDEX & icubeB, GRID_COORD_TYPE dist[DIM3])
{
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB);

  set_coord_3D(0, dist);
  if (gcubeA_index == ISOVERT::NO_INDEX || gcubeB_index == ISOVERT::NO_INDEX)
    { return; }

  subtract_coord_3D(isovert.gcube_list[gcubeA_index].cube_coord,
                    isovert.gcube_list[gcubeB_index].cube_coord,
                    dist);
  abs_coord_3D(dist, dist);

  if (dist[1] < dist[2]) { std:: swap(dist[1], dist[2]); }
  if (dist[0] < dist[1]) { std:: swap(dist[0], dist[1]); }
  if (dist[1] < dist[2]) { std:: swap(dist[1], dist[2]); }
}

bool is_dist_lt(const GRID_COORD_TYPE distA[DIM3],
                const GRID_COORD_TYPE distB[DIM3])
{
  if (distA[0] < distB[0]) { return(true); }
  if (distA[0] > distB[0]) { return(false); }

  // distA[0] == distB[0]
  if (distA[1] < distB[1]) { return(true); }
  if (distA[1] > distB[1]) { return(false); }

  // distA[0] == distB[0] and distA[1] == distB[1]
  if (distA[2] < distB[2]) { return(true); }

  return(false);
}

void find_closest(const ISOVERT & isovert, const VERTEX_INDEX icubeA,
                  const vector<VERTEX_INDEX> & list, VERTEX_INDEX & icubeB)
{
  GRID_COORD_TYPE distB[DIM3], distC[DIM3];

  if (list.size() == 0) { 
    icubeB = 0;
    return; 
  }

  icubeB = list[0];
  compute_dist(isovert, icubeA, icubeB, distB);

  for (int i = 1; i < list.size(); i++) {
    VERTEX_INDEX icubeC = list[i];
    compute_dist(isovert, icubeA, icubeC, distC);
    if (is_dist_lt(distC, distB)) {
      icubeB = icubeC;
      copy_coord_3D(distC, distB);
    }
  }
}

// Find closest cube in list, which is in opposite direction from icubeB
void find_closest_opposite
(const ISOVERT & isovert, const VERTEX_INDEX icubeA, 
 const VERTEX_INDEX icubeB,
 const vector<VERTEX_INDEX> & list, VERTEX_INDEX & icubeC,
 bool & flag_found)
{
  const INDEX_DIFF_TYPE gcubeA_index = isovert.GCubeIndex(icubeA);
  const INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(icubeB);
  GRID_COORD_TYPE distC[DIM3], distD[DIM3];
  GRID_BOX region(DIM3);

  flag_found = false;

  if (gcubeA_index == ISOVERT::NO_INDEX || gcubeB_index == ISOVERT::NO_INDEX)
    { return; }

  region.SetMinCoord(isovert.gcube_list[gcubeA_index].cube_coord);
  region.SetMaxCoord(isovert.gcube_list[gcubeA_index].cube_coord);
  region.Extend(isovert.gcube_list[gcubeB_index].cube_coord);

  icubeC = 0;
  if (list.size() == 0) { return; }

  for (int i = 0; i < list.size(); i++) {

    VERTEX_INDEX icubeD = list[i];
    INDEX_DIFF_TYPE gcubeD_index = isovert.GCubeIndex(icubeD);
    if (gcubeD_index == ISOVERT::NO_INDEX) { continue; }

    if (region.Contains(isovert.gcube_list[gcubeD_index].cube_coord))
      { continue; }

    compute_dist(isovert, icubeA, icubeD, distD);
    if (!flag_found || is_dist_lt(distD, distC)) {
      icubeC = icubeD;
      copy_coord_3D(distD, distC);
      flag_found = true;
    }
  }
}


// Return true if distance matches one of the output neighbor flags.
bool matches_output_neighbor_flag
(const GRID_COORD_TYPE dist[DIM3], const OUTPUT_INFO & output_info)
{
  if (dist[0] < 2) { return(false); }
  if (dist[1] < 2) { return(false); }

  if (dist[0] == 2 && dist[1] == 2) {

    if (dist[2] == 0) { 
      if (output_info.flag_output_neighbors_220) { return(true); }
    }
    else if (dist[2] == 1) {
      if (output_info.flag_output_neighbors_221) { return(true); }
    }
    else if (dist[2] == 2) {
      if (output_info.flag_output_neighbors_222) { return(true); }
    }
  }
  else if (dist[0] >= 3) {

    if (dist[1] == 2) {

      if (dist[2] == 0) {
        if (output_info.flag_output_neighbors_320) { return(true); }
      }
      else if (dist[2] == 1) {
        if (output_info.flag_output_neighbors_321) { return(true); }
      }
      else if (dist[2] == 2) {
        if (output_info.flag_output_neighbors_322) { return(true); }
      }
    }
    else if (dist[1] >= 3) {
      if (dist[2] == 0) {
        if (output_info.flag_output_neighbors_330) { return(true); }
      }
      else if (dist[2] == 1) {
        if (output_info.flag_output_neighbors_331) { return(true); }
      }
      else if (dist[2] >= 2) {
        if (output_info.flag_output_neighbors_332) { return(true); }
      }
    }
  }

  return(false);
}


void report_far_neighbors
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert)
{
  const int bin_width = output_info.bin_width;
  GRID_COORD_TYPE dist[DIM3];
  vector<VERTEX_INDEX> selected_list;
  vector<VERTEX_INDEX> list;

  BIN_GRID<VERTEX_INDEX> bin_grid;
  init_bin_grid(grid, bin_width, bin_grid);

  for (int i = 0; i < isovert.gcube_list.size(); i++) {
    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
      VERTEX_INDEX icubeA = isovert.CubeIndex(i);
      bin_grid_insert(grid, bin_width, icubeA, bin_grid);
    }
  }

  NUM_TYPE count = 0;

  cout << endl;
  cout << "Neighboring cube pairs which are far apart:" << endl;
  for (int i = 0; i < isovert.gcube_list.size(); i++) {

    if (isovert.gcube_list[i].flag != SELECTED_GCUBE) { continue; }

    VERTEX_INDEX icubeA = isovert.CubeIndex(i);

    get_selected(grid, icubeA, bin_grid, bin_width, selected_list);

    // Copy selected_list to list, removing icubeA.
    list.clear();
    for (int i = 0; i < selected_list.size(); i++) {
      if (selected_list[i] != icubeA) 
        { list.push_back(selected_list[i]); }
    }

    if (list.size() == 0) { continue; }

    VERTEX_INDEX icubeB, icubeC;
    bool flag = find_corner(isovert, list, icubeB);
    if (flag) {
      compute_dist(isovert, icubeA, icubeB, dist);
      if (dist[0] > 2) 
        { find_closest(isovert, icubeA, list, icubeB); }
    }
    else {
      find_closest(isovert, icubeA, list, icubeB); 
    }

    compute_dist(isovert, icubeA, icubeB, dist);

    if (matches_output_neighbor_flag(dist, output_info)) {
      grid.PrintIndexAndCoord
        (cout, "Cubes ", icubeA, " and (closest) ", icubeB, ".  Distance: ");
      print_coord3D(cout, dist);
      cout << endl;

      count++;
    }

    find_closest_opposite(isovert, icubeA, icubeB, list, icubeC, flag);
    if (flag) {

      compute_dist(isovert, icubeA, icubeC, dist);


    if (matches_output_neighbor_flag(dist, output_info)) {
        grid.PrintIndexAndCoord
          (cout, "Cubes ", icubeA, " and (opposite) ", icubeC, 
           ".  Distance: ");
        print_coord3D(cout, dist);
        cout << endl;

        count++;
      }
    }
  }

  cout << "Reported " << count << " pairs of neighboring cubes." << endl;
  cout << endl;
}


// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

void MERGESHARP::report_num_cubes
(const SHARPISO_GRID & full_scalar_grid, const INPUT_INFO & input_info,
 const MERGESHARP_DATA & mergesharp_data)
{
  const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
  const int num_cubes_in_mergesharp_data =
    mergesharp_data.ScalarGrid().ComputeNumCubes();

  if (!input_info.use_stdout && !input_info.flag_silent) {

    if (input_info.flag_subsample) {
      // subsampled grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_mergesharp_data << " subsampled grid cubes." << endl;
    }
    else if (input_info.flag_supersample) {
      // supersample grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_mergesharp_data << " supersampled grid cubes." << endl;
    }
    else {
      // use full_scalar_grid
      cout << num_grid_cubes << " grid cubes." << endl;
    }
  }

}

void MERGESHARP::report_mergesharp_param
(const MERGESHARP_PARAM & mergesharp_param)
{
  const VERTEX_POSITION_METHOD vpos_method = 
    mergesharp_param.vertex_position_method;

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
        mergesharp_param.GradSelectionMethod();
      get_grad_selection_string(grad_selection_method, s);

      cout << s << endl;
    }
    break;

  case EDGEI_INPUT_DATA:
    cout << "Hermite data" << endl;
    break;

  default:
    cout << "Unkown" << endl;
    break;
  }

  if (vpos_method != CUBECENTER && vpos_method != CENTROID_EDGE_ISO) {

    if (mergesharp_param.flag_merge_sharp) {
      cout << "Merge isosurface vertices near sharp edges and corners." 
           << endl;
      cout << "Merge distance (Linf) threshold: "
           << mergesharp_param.linf_dist_thresh_merge_sharp << endl;
    }

    cout << "Maximum small eigenvalue: "
         << mergesharp_param.max_small_eigenvalue << endl;
    cout << "Max (Linf) distance from cube to isosurface vertex: "
         << mergesharp_param.max_dist << endl;

    if (vpos_method != EDGEI_INPUT_DATA) {
      cout << "Max small gradient magnitude: "
           << mergesharp_param.max_small_magnitude << endl;
      cout << "Gradient selection cube offset: "
           << mergesharp_param.grad_selection_cube_offset << endl;
    }

    if (mergesharp_param.use_lindstrom) {
      cout << "Using Lindstrom formula." << endl;
    }

    if (mergesharp_param.flag_dist2centroid) {
      cout << "Using distance to centroid." << endl;
    }
    else {
      cout << "Using distance to center." << endl;
    }

  }

  if (mergesharp_param.flag_merge_sharp) {

    if (mergesharp_param.flag_allow_conflict) {
      cout << "Allow isosurface vertex conflicts." << endl;
    }
    else {
      if (mergesharp_param.flag_clamp_conflict) 
        {   cout <<"Do NOT allow conflicts, instead "<<endl;
			cout << "Resolve conflict by clamping coordinates." << endl; }
      else
        { 
			cout <<"Do NOT allow conflicts, instead "<<endl;
			cout << "Resolve conflict by reverting to centroid." << endl; }
    }
  }

  if (vpos_method != CUBECENTER && vpos_method != EDGEI_INPUT_DATA) {

    if (mergesharp_param.use_sharp_edgeI) {
      cout << "Using sharp formula to calculate intersections of isosurface and grid edges." << endl;
    }
    else {
      cout << "Using interpolation to calculate intersections of isosurface and grid edges." << endl;
    }
  }

  if (mergesharp_param.flag_merge_sharp) {
    if (mergesharp_param.flag_recompute_isovert) {
      cout << "Recompute isosurface vertex positions for unselected/uncovered cubes." << endl;
    }
  }

  if (mergesharp_param.allow_multiple_iso_vertices) {
    cout << "Allow multiple isosurface vertices per grid cube." << endl;

    if (mergesharp_param.flag_resolve_ambiguous_facets &&
        !mergesharp_param.flag_merge_sharp) {
      cout << "Resolve ambiguous facets." << endl;
    }
    else if (mergesharp_param.flag_separate_neg) {
      cout << "Separate negative vertices." << endl;
    }
    else {
      cout << "Separate positive vertices." << endl;
    }
  }
  else {
    cout << "One isosurface vertex per grid cube." << endl;
  }

  if (mergesharp_param.flag_delete_isolated_vertices) {
    cout << "Remove isosurface vertices which are not in any isosurface polygons." << endl;
  }
  else {
    cout << "Keep all isosurface vertices, including vertices which are not in any polygons." << endl;
  }

  cout << endl;
}


void MERGESHARP::report_iso_info3D
(const OUTPUT_INFO & output_info, const MERGESHARP_DATA & mergesharp_data,
 const DUAL_ISOSURFACE & dual_isosurface,
 const MERGESHARP_INFO & mergesharp_info)
{
  const char * indent4 = "    ";
  string grid_element_name = "cubes";

  VERTEX_INDEX numv = (dual_isosurface.vertex_coord.size())/DIM3;
  VERTEX_INDEX num_tri = 
    (dual_isosurface.tri_vert.size())/NUM_VERT_PER_TRI;
  VERTEX_INDEX num_quad = 
    (dual_isosurface.quad_vert.size())/NUM_VERT_PER_QUAD;
  VERTEX_INDEX num_grid_cubes = mergesharp_info.grid.num_cubes;
  VERTEX_INDEX num_non_empty_cubes = mergesharp_info.scalar.num_non_empty_cubes;

  float percent = 0.0;
  if (num_grid_cubes > 0)
    { percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
  int ipercent = int(100*percent);
  cout << "  Isovalue " << output_info.isovalue << ".  "
       << numv << " isosurface vertices.  " << endl;
  if (num_tri+num_quad == 0) {
    cout << "No isosurface polygons." << endl;
  }
  else {
    cout << "    ";
    if (num_tri > 0) {
      cout << num_tri << " isosurface triangles.  ";
    }
    if (num_quad > 0) {
      cout << num_quad << " isosurface quadrilaterals.";
    }
    cout << endl;
  }
  



  if (output_info.flag_output_alg_info) {
    VERTEX_POSITION_METHOD vpos_method = output_info.VertexPositionMethod();

    if (vpos_method == GRADIENT_POSITIONING ||
        vpos_method == EDGEI_INTERPOLATE ||
        vpos_method == EDGEI_GRADIENT ||
        vpos_method == EDGEI_INPUT_DATA) {

      cout << endl;
      if (!output_info.flag_merge_sharp) {
        cout << "  # of conflicts: "
             << mergesharp_info.sharpiso.num_conflicts << endl;
      }
      cout << "  # of sharp corners: "
           << mergesharp_info.sharpiso.num_sharp_corners << endl;
      cout << "  # of isosurface vertices on sharp edges: "
           << mergesharp_info.sharpiso.num_sharp_edges << endl;
      cout << "  # of smooth isosurface vertices: "
           << mergesharp_info.sharpiso.num_smooth_vertices << endl;
      if (!output_info.flag_merge_sharp && output_info.use_Linf_dist) {
        cout << "  # of vertices at min Linf distance to cube center: "
             << mergesharp_info.sharpiso.num_Linf_iso_vertex_locations << endl;
      }

      if (output_info.flag_merge_sharp) {
        cout << "  # of merged isosurface vertices: "
             << mergesharp_info.sharpiso.num_merged_iso_vertices << endl;
      }

      if (output_info.allow_multiple_iso_vertices) {

        if (!output_info.flag_merge_sharp) {
          if (output_info.flag_resolve_ambiguous_facets) {
            cout << "  # of non-ambiguous cubes: "
                 << mergesharp_info.sharpiso.num_cube_not_ambiguous << endl;
            cout << "  # of separate positive cubes: "
                 << mergesharp_info.sharpiso.num_cube_separate_pos << endl;
            cout << "  # of separate negative cubes: "
                 << mergesharp_info.sharpiso.num_cube_separate_neg << endl;
            cout << "  # of unresolved ambiguous cubes: "
                 << mergesharp_info.sharpiso.num_cube_unresolved_ambiguity 
                 << endl;
          }
        }

        cout << "  # of cubes with single isov: "
             << mergesharp_info.sharpiso.num_cube_single_isov << endl;
        cout << "  # of cubes with multi isov: "
             << mergesharp_info.sharpiso.num_cube_multi_isov << endl;

        if (output_info.flag_merge_sharp) { 

          if (output_info.flag_split_non_manifold) {
            cout << "  # of cubes changed to 2 isov to avoid non-manifold edges: "
                 << mergesharp_info.sharpiso.num_non_manifold_split
                 << endl;
          }

          if (output_info.flag_select_split) {
            cout << "  # of cubes changed in selecting isosurface patch splits: "
                 << mergesharp_info.sharpiso.num_1_2_change << endl;
          }

        }
      }

      if (output_info.flag_merge_sharp) {
        if (output_info.flag_check_disk) {
          cout << "  # of merges blocked by non-disk isosurface patches: "
               << mergesharp_info.sharpiso.num_non_disk_isopatches << endl;
        }
      }


      cout << endl;
    }

    if (vpos_method == CENTROID_EDGE_ISO) {
      if (output_info.allow_multiple_iso_vertices) {
        cout << endl;
        cout << "  # of cubes with single isov: "
             << mergesharp_info.sharpiso.num_cube_single_isov << endl;
        cout << "  # of cubes with multi isov: "
             << mergesharp_info.sharpiso.num_cube_multi_isov << endl;
        cout << endl;
      }
    }
  }
}

void output_yes_no(const char * s, const bool flag)
{
  if (flag) 
    { cout << s << "yes"; }
  else 
    { cout << s << "no"; }
}

/// Report information about isosurface vertices
void report_isovert_cube_info
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid, 
 const ISOVERT & isovert, const VERTEX_INDEX cube_index)
{
  const INDEX_DIFF_TYPE gcube_index = isovert.GCubeIndex(cube_index);
  const NUM_TYPE num_eigenvalues = 
    isovert.gcube_list[gcube_index].num_eigenvalues;
  string s;

  if (gcube_index == ISOVERT::NO_INDEX) { return; }

  cout << "Cube " << cube_index << ": ";
  ijkgrid_output_vertex_coord(cout, isovert.index_grid, cube_index);

  if (!grid.IsUnitSpacing()) {
    COORD_TYPE scaled_cube_coord[DIM3];
    grid.ComputeScaledCoord(cube_index, scaled_cube_coord);
    cout << "  ";
    IJK::print_coord3D(cout, scaled_cube_coord);
  }

  cout << "    Num eigenvalues: " << num_eigenvalues;
  cout << endl;
  cout << "    Isovert coord A: ";
  IJK::print_coord3D
    (cout, isovert.gcube_list[gcube_index].isovert_coord);
  cout << " B: ";
  IJK::print_coord3D
    (cout, isovert.gcube_list[gcube_index].isovert_coordB);
  cout << endl;

  if (num_eigenvalues == 1) {
    cout << "    Orthogonal direction: ";
    IJK::print_coord3D
      (cout, isovert.gcube_list[gcube_index].direction);
    cout << endl;
  }
  if (num_eigenvalues == 2) {
    cout << "    Sharp edge direction: ";
    IJK::print_coord3D
      (cout, isovert.gcube_list[gcube_index].direction);
    cout << endl;
  }

  convert2string(isovert.gcube_list[gcube_index].flag, s);
  cout << "    Dist (Linf): " 
       << isovert.gcube_list[gcube_index].linf_dist;
  cout << "  Type: " << s;
  cout << "  Maps to: " 
       << isovert.gcube_list[gcube_index].maps_to_cube
       << endl;

  output_yes_no("    Conflict? ", 
                isovert.gcube_list[gcube_index].flag_conflict);
  output_yes_no("  SVD coord far? ", 
                isovert.gcube_list[gcube_index].flag_far);
  output_yes_no("  Centroid coord? ", 
                isovert.gcube_list[gcube_index].flag_centroid_location);
  cout << endl;
  if (output_info.flag_recompute_using_adjacent) {
    output_yes_no("    Recomputed using adjacent? ",
                  isovert.gcube_list[gcube_index].flag_recomputed_using_adjacent);
    cout << endl;
  }
  if (isovert.gcube_list[gcube_index].cover_type != NO_ADJACENT_CUBE) {
    convert2string(isovert.gcube_list[gcube_index].cover_type, s);
    cout << "    Cover type: " << s << endl;
  }

  cout << "    ";
  if (isovert.gcube_list[gcube_index].flag_using_substitute_coord)
    { cout << "Substitute coord.  "; }
  if (isovert.gcube_list[gcube_index].flag_coord_from_other_cube) 
    { cout << "Coord from other cube.  "; }
  else if (isovert.gcube_list[gcube_index].flag_coord_from_vertex)
    { cout << "Coord from vertex.  "; }
  else if (isovert.gcube_list[gcube_index].flag_coord_from_edge)
    { cout << "Coord from edge.  "; }

  if (isovert.gcube_list[gcube_index].flag_ignore_mismatch) 
    { cout << "Mismatch ignored.  ";  }

  cout << "Isotable index: " 
       << int(isovert.gcube_list[gcube_index].table_index);
  cout << endl;
}

/// Report information about cubes containing corner or edge iso vertices.
void MERGESHARP::report_sharp_cubes
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert)
{
  std::vector<NUM_TYPE> sharp_gcube_list;

  // Get corner or edge cubes (in sorted order).
  get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

  cout << endl;
  cout << "Cubes containing corners or edges: " << endl;
  for (NUM_TYPE i = 0; i < sharp_gcube_list.size(); i++) {
    NUM_TYPE gcube_index = sharp_gcube_list[i];
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    report_isovert_cube_info(output_info, grid, isovert, cube_index);
  }
  cout << endl;
}

/// Report information about active cubes.
/// A cube is active if it has some bipolar edge.
void MERGESHARP::report_active_cubes
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert)
{
  cout << endl;
  cout << "Active cubes (containing bipolar edges): " << endl;
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    VERTEX_INDEX cube_index = isovert.CubeIndex(i);
    report_isovert_cube_info(output_info, grid, isovert, cube_index);
  }
  cout << endl;
}

/// Report information about cubes which map to themselves.
void MERGESHARP::report_cubes_which_map_to_self
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert)
{
  cout << endl;
  cout << "Cubes which map to themselves:" << endl;
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    VERTEX_INDEX cube_index = isovert.CubeIndex(i);
    if (isovert.gcube_list[i].maps_to_cube == cube_index) {
      report_isovert_cube_info(output_info, grid, isovert, cube_index);
    }
  }
  cout << endl;
}

/// Report information about cubes which map to a specific cube.
void MERGESHARP::report_cubes_which_map_to
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert, const VERTEX_INDEX & to_cube)
{
  cout << endl;
  grid.PrintIndexAndCoord
    (cout, "Cubes which map to cube ", to_cube, ":\n");
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    if (isovert.gcube_list[i].maps_to_cube == to_cube) {
      VERTEX_INDEX cube_index = isovert.CubeIndex(i);
      report_isovert_cube_info(output_info, grid, isovert, cube_index);
    }
  }
  cout << endl;
}


/// Report information about isosurface vertices
void MERGESHARP::report_isovert_info
(const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
 const ISOVERT & isovert)
{
  const VERTEX_POSITION_METHOD vpos_method = output_info.VertexPositionMethod();
  std::vector<NUM_TYPE> sharp_gcube_list;


  if (vpos_method != GRADIENT_POSITIONING &&
      vpos_method != EDGEI_INTERPOLATE &&
      vpos_method != EDGEI_GRADIENT &&
      vpos_method != EDGEI_INPUT_DATA) { return; }

  if (output_info.flag_output_selected) {

    // Get corner or edge cubes (in sorted order).
    get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

    cout << endl;
    cout << "Selected cubes containing corners or edges: " << endl;
    for (NUM_TYPE i = 0; i < sharp_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sharp_gcube_list[i];

      if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {
        VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
        report_isovert_cube_info(output_info, grid, isovert, cube_index);
      }
    }
    cout << endl;
  }

  if (output_info.flag_output_sharp) 
    { report_sharp_cubes(output_info, grid, isovert); }

  if (output_info.flag_output_active) 
    { report_active_cubes(output_info, grid, isovert); }

  if (output_info.flag_output_map_to_self) 
    { report_cubes_which_map_to_self(output_info, grid, isovert); }

  if (output_info.flag_output_map_to) {
    report_cubes_which_map_to
      (output_info, grid, isovert, output_info.to_cube); 
  }

}


// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

void MERGESHARP::report_mergesharp_time
(const INPUT_INFO & input_info, const MERGESHARP_TIME & mergesharp_time,
 const char * mesh_type_string)
{
  cout << "CPU time to run Marching Cubes: "
       << mergesharp_time.total << " seconds." << endl;

  if ((input_info.VertexPositionMethod() != CUBECENTER &&
       input_info.VertexPositionMethod() != CENTROID_EDGE_ISO) &&
      input_info.flag_merge_sharp) {

    cout << "    Time to position "
         << mesh_type_string << " vertices: "
         << mergesharp_time.position << " seconds." << endl;

    cout << "    Time to extract " << mesh_type_string << " triangles: "
         << mergesharp_time.extract << " seconds." << endl;

    cout << "    Time to merge sharp "
         << mesh_type_string << " vertices: "
         << mergesharp_time.merge_sharp << " seconds." << endl;
  }
  else {

    cout << "    Time to extract " << mesh_type_string << " triangles: "
         << mergesharp_time.extract << " seconds." << endl;

    cout << "    Time to merge identical "
         << mesh_type_string << " vertices: "
         << mergesharp_time.merge_identical << " seconds." << endl;

    cout << "    Time to position "
         << mesh_type_string << " vertices: "
         << mergesharp_time.position << " seconds." << endl;
  }

}


void MERGESHARP::report_time
(const INPUT_INFO & input_info, const IO_TIME & io_time,
 const MERGESHARP_TIME & mergesharp_time, const double total_elapsed_time)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * INTERVAL_VOLUME_STRING = "interval volume";
  const char * mesh_type_string = NULL;

  mesh_type_string = ISOSURFACE_STRING;

  cout << "Time to read file " << input_info.scalar_filename << ": "
       << io_time.read_nrrd_time << " seconds." << endl;

  cout << "Time to read " << mesh_type_string << " lookup tables: "
       << io_time.read_table_time << " seconds." << endl;

  report_mergesharp_time(input_info, mergesharp_time, mesh_type_string);
  if (!input_info.nowrite_flag) {
    cout << "Time to write "
         << mesh_type_string << ": "
         << io_time.write_time << " seconds." << endl;
  };
  cout << "Total elapsed time: " << total_elapsed_time
       << " seconds." << endl;
}

// **************************************************
// WRITE ISOSURFACE VERTEX INFORMATION TO FILE
// **************************************************

namespace {

  // Construct isovert info filename from "from_filename".
  void construct_isovert_info_filename
  (const std::string & from_filename, std::string & info_filename)
  {
    info_filename = remove_off_suffix(from_filename);
    info_filename += ".isov_info";
  }

}

void MERGESHARP::write_isovert_info
(const OUTPUT_INFO & output_info,
 const std::vector<DUAL_ISOVERT_INFO> & isovert_info)
{
  ofstream info_file;
  IJK::PROCEDURE_ERROR error("write_isovert_info");

  string info_filename;
  construct_isovert_info_filename
    (output_info.output_filename, info_filename);

  info_file.open(info_filename.c_str(), ios::out);
  if (!info_file.good()) {
    cerr << "Unable to open isovert info file " << info_filename << "." 
         << endl;
    exit(95);
  };

  for (NUM_TYPE i = 0; i < isovert_info.size(); i++) {
    info_file << i;
    info_file << " " << int(isovert_info[i].cube_index);
    info_file << " " << int(isovert_info[i].table_index);
    info_file << " " << int(isovert_info[i].patch_index);
    info_file << " " << int(isovert_info[i].num_eigenvalues);
    info_file << " " << int(isovert_info[i].flag_centroid_location);
    info_file << endl;
  }

  info_file.close();

  if (!output_info.flag_silent) {
    cout << "Wrote isosurface vertex info to file: " << info_filename << endl;

  }
}

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

// local namespace
namespace {


  void usage_msg(std::ostream & out, const char * command_path)
  {
    string command_name;
    string prefix, suffix;

    command_name = command_path;

    // remove path from file name
    split_string(command_path, PATH_DELIMITER, prefix, suffix);
    if (suffix != "") { command_name = suffix; }

    out << "Usage: " << command_name 
        << " [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
  }


  void options_msg()
  {
    cerr << "OPTIONS:" << endl;
    cerr << "  [-subsample S]" << endl;
    cerr << "  [-position {centroid|cube_center|gradC|gradN|gradCS|gradNS|gradXS"
         << endl
         << "              gradIE|gradIES|gradIEDir|gradCD|gradNIE|gradNIES|gradBIES"
         << endl
         << "              gradES|gradEC}]" << endl;
    cerr << "  [-gradient {gradient_nrrd_filename}]"
         << " [-normal {normal_off_filename}]" << endl;
    cerr << "  [-merge_sharp | -no_merge_sharp] [-merge_linf_th <D>]" << endl;
    cerr << "  [-grad2hermite | -grad2hermiteI]" << endl;
    cerr << "  [-manifold] [-select_split]" << endl;
    cerr << "  [-trimesh]" << endl;
    cerr << "  [-max_eigen {max}]" << endl;
    cerr << "  [-max_dist {D}] [-max_mag {M}] [-snap_dist {D}]" << endl;    
    cerr << "  [-gradS_offset {offset}] [-min_gradS_offset {offset}] [-max_grad_dist {D}]" << endl;
    cerr << "  [-min_triangle_angle {A}] [-min_normal_angle {A}]" << endl;
    cerr << "  [-sharp_edgeI | -interpolate_edgeI]" << endl;
    cerr << "  [-lindstrom]" << endl;
    cerr << "  [-single_isov | -multi_isov | -split_non_manifold]" << endl;
    cerr << "  [-sep_pos | -sep_neg | -resolve_ambig]" << endl;
    cerr << "  [-check_disk | -no_check_disk]" << endl;
    cerr << "  [-allow_conflict |-clamp_conflict | -centroid_conflict]"\
         << endl;
    cerr << "  [-clamp_far] [-centroid_far]" << endl;
    cerr << "  [-recompute_isovert | -no_recompute_isovert]"<<endl;
    cerr << "  [-check_triangle_angle | -no_check_triangle_angle]"<<endl;
    cerr << "  [-Linf | -no_Linf]" << endl;
    cerr << "  [-dist2center | -dist2centroid]" << endl;
    cerr << "  [-no_round | -round <n>]" << endl;
    cerr << "  [-map_extended]" <<endl;
    cerr << "  [-keepv]" << endl;
    cerr << "  [-off|-iv] [-o {output_filename}] [-stdout]"
         << endl;
    cerr << "  [-s] [-out_param] [-info]" << endl;
    cerr << "  [-out_selected] [-out_sharp] [-out_active] [-out_map_to_self]" 
         << endl;
    cerr << "  [-out_map_to {cube index}] [-out_neighbors \" 222|22x|32x|33x \"]" 
         << endl;
    cerr << "  [-out_isovert [corner|edge|sharp|smooth|all] [selected|all|uncovered]" << endl;
    cerr << "  [-write_isov_info] [-nowrite] [-time]" << endl;
    cerr << "  [-help]" << endl;
  }

}

void MERGESHARP::usage_error(const char * command_path)
{
  usage_msg(cerr, command_path);
  options_msg();
  exit(10);
}

void MERGESHARP::help(const char * command_path)
{
  usage_msg(cout, command_path);
  cout << endl;
  cout << "mergesharp - Merge sharp isosurface generation algorithm." << endl;
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
  cout << "  -position gradXS: Position isosurface vertices using svd"
       << endl;
  cout << "       on selected vertex gradients of cube and cube neighbors."
       << endl;
  cout << "       Neighbors include diagonal neighbors." << endl;
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
  cout << "  -position gradBIES: Position isosurface vertices using svd"
       << endl;
  cout << "       on selected gradients on boundary of zero gradient region."
       << endl;
  cout << "       Use only gradients from vertices incident on bipolar edges."
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
  cout << "       Use endpoint gradients to compute isosurface-edge intersections." << endl;
  cout << "  -gradient {gradient_nrrd_filename}: Read gradients from gradient nrrd file." << endl;
  cout << "  -normal {normal_off_filename}: Read edge-isosurface intersections"
       << endl
       << "      and normals from OFF file normal_off_filename." << endl;
  cout << "  -merge_sharp:   Merge vertices near sharp edges/corners." << endl;
  cout << "  -grad2hermite:  Convert gradient to hermite data." << endl;
  cout << "  -grad2hermiteI: Convert gradient to hermite data using linear interpolation." << endl;
  cout << "  -manifold:      Output is the embedding of a manifold." << endl
       << "     Equivalent to \"-check_disk -split_non_manifold\"." << endl
       << "     Note: The manifold may be self intersecting." << endl;
  cout << "  -select_split:  Select which cube has split isosurface vertices"
       << endl
       << "     where two adjacent cubes share an ambiguous facet." << endl;
  cout << "  -single_isov: Each intersected cube generates a single isosurface vertex." << endl;
  cout << "  -multi_isov:  An intersected cube may generate multiple isosurface vertices."  << endl;
  cout << "  -split_non_manifold:  Split vertices to avoid non-manifold edges."
       << endl;
  cout << "  -sep_pos:     Use dual isosurface table separating positive "
       << endl
       << "                grid vertices." << endl;
  cout << "  -sep_neg:     Use dual isosurface table separating negative "
       << endl
       << "                grid vertices." << endl;
  cout << "  -resolve_ambig:  Selectively resolve ambiguities." << endl
       << "       Note: Not all ambiguities may be resolved." << endl
       << "             Unresolved ambiguities may create non-manifold regions."
       << endl;
  cout << "  -merge_linf_th {D} : Do not select sharp vertices further"
       << endl
       << "            than Linf dist D from cube center." << endl;
  cerr << "  -max_eigen {max}: Set maximum small eigenvalue to max."
       << endl;
  cerr << "  -max_dist {D}:    Set max Linf distance from cube to isosurface vertex."
       << endl;
  cerr << "  -max_mag {max}:  Set maximum small gradient magnitude to max."
       << endl;
  cerr << "           Gradients with magnitude below max are ignored." << endl;
  cerr << "  -snap_dist {D}:  Snap points within distance D to cube."
       << endl;
  cerr << "  -max_grad_dist {D}: Max distance (integer) of gradients to cube. (Default 1.)" << endl;
  cerr << "  -gradS_offset {offset}: Set cube offset for gradient selection to offset."
       << endl;
  cerr << "           Offset measures distance from cube boundary." << endl;
  cerr << "  -min_triangle_angle {A}:  Avoid cube merging which creates triangles" << endl
       << "             with angle less than {A}." << endl;
  cerr << "  -min_normal_angle {A}:  Allow cube merging if angle between" << endl
       << "             normals of new and old triangle is less than {A}." 
       << endl;
  cout << "  -lindstrom:   Use Lindstrom's equation to compute sharp point."
       << endl;
  cout << "  -allow_conflict:  Allow more than one isosurface vertex in a cube."
       << endl;
  cout << "  -clamp_conflict:  Settle conflicts by clamping to cube."
       << " (Default.)"  << endl;
  cout << "  -centroid_conflict:  Settle conflicts by using centroid."
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
  cout << "  -recompute_isovert:    Recompute isosurface vertex locations"
       << endl
       << "             for unavailable cubes." << endl;
  cout << "  -no_recompute_isovert: Don't recompute isosurface vertex locations"
       << endl
       << "             for unavailable cubes "<<endl;
  cout << "  -dist2center:  Use distance to center in lindstrom." << endl;
  cout << "  -dist2centroid:  Use distance to centroid of isourface-edge"
       << endl
       << "                   intersections in lindstrom." << endl;
  cout << "  -Linf:     Use Linf metric to resolve conflicts." << endl;
  cout << "  -no_Linf:  Don't use Linf metric to resolve conflicts." << endl;
  cout << "  -no_round:  Don't round coordinates." << endl;
  cout << "  -round <n>: Round coordinates to nearest 1/n." << endl;
  cout << "              Suggest using n=16,32,64,... or 2^k for some k."
       << endl;
  cerr << "  -keepv:    Keep isosurface vertices.  Do not remove isosurface vertices"
       << endl
       << "             which do not lie in any isosurface polygon." << endl;
  cout << "  -check_disk: Check that merged vertices form a disk." << endl;
  cout << "  -no_check_disk: Skip disk check for merged vertices." << endl;
  cout << "  -trimesh:   Output triangle mesh." << endl;
  cout << "  -map_extended: Use the extended version of mapping to sharp vertices." << endl;
  cout << "  -off: Output in geomview OFF format. (Default.)" << endl;
  cout << "  -iv: Output in OpenInventor .iv format." << endl;
  cout << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
  cout << "  -stdout: Write isosurface to standard output." << endl;
  cout << "  -nowrite: Don't write isosurface." << endl;
  cout << "  -time: Output running time." << endl;
  cout << "  -s: Silent mode." << endl;
  cout << "  -out_param: Print mergesharp parameters." << endl;
  cout << "  -info: Print algorithm information." << endl;
  cout << "  -out_selected: Output list of selected cubes." << endl;
  cout << "  -out_sharp: Output list of cubes containing corner or edge vertices." << endl;
  cout << "  -out_active: Output list of active cubes." << endl;
  cout << "  -out_map_to_self: Output list of cubes which map to themselves." 
       << endl;
  cout << "  -out_map_to {cube index}: Output list of cubes which map to given cube index." << endl;
  cerr << "  -out_neighbors {\" 222 | 22x | 32x | 33x \"}:" << endl;
  cerr << "       Output pairs of selected neighboring cubes at given distances." << endl;
  cerr << "       Distances: 222 = (2,2,2). 22x = (2,2,x). 32x = (3,2,x). 33x = (3,3,x)." << endl;
  cerr << "       22x |: Output neighbors at distance (2,2,2) (or greater.)"
       << endl;
  cout << "  -out_isovert [arg1] [arg2]:  Output isosurface vertices to off file." << endl;
  cout << "  -write_isov_info:  Write isosurface vertex information to file."
       << endl;
  cout << "     File format: isosurface vertex index, cube index," << endl;
  cout << "       isosurface lookup table index, isosurface patch index," 
       << endl;
  cout << "       number of eigenvalues, centroid location flag." << endl;
  cout << "     If centroid location flag is 1, location is centroid" << endl;
  cout << "       of (grid edge)-isosurface intersections." << endl;
  cout << "  -help: Print this help message." << endl;
  exit(20);
}


// **************************************************
// CLASS IO_INFO
// **************************************************

/// IO information
void MERGESHARP::IO_INFO::Init()
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
  flag_output_selected = false;
  flag_output_sharp = false;
  flag_output_active = false;
  flag_output_map_to_self = false;
  flag_output_map_to = false;
  to_cube = 0;
  flag_output_isovert = false;
  flag_recompute_isovert = true; // recompute the isovert for unavailable cubes
  flag_check_triangle_angle = true;
  flag_output_isovert = false;
  grid_spacing.resize(3,1);

  flag_output_neighbors = false;
  flag_output_neighbors_222 = false;
  flag_output_neighbors_221 = false;
  flag_output_neighbors_220 = false;
  flag_output_neighbors_322 = false;
  flag_output_neighbors_321 = false;
  flag_output_neighbors_320 = false;
  flag_output_neighbors_332 = false;
  flag_output_neighbors_331 = false;
  flag_output_neighbors_330 = false;
}

void MERGESHARP::IO_INFO::Set(const IO_INFO & io_info)
{
  *this = io_info;
}

// **************************************************
// CLASS INPUT_INFO
// **************************************************

/// Clear input information
void MERGESHARP::INPUT_INFO::Init()
{
  Clear();
}

/// Clear input information
void MERGESHARP::INPUT_INFO::Clear()
{
  is_vertex_position_method_set = false;
  is_use_sharp_edgeI_set = false;
  is_conflict_set = false;
  isovalue.clear();
  isovalue_string.clear();
  isotable_directory = "";
}

// **************************************************
// class OUTPUT_INFO
// **************************************************

void MERGESHARP::OUTPUT_INFO::Init()
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
    string ofilename = 
      remove_nrrd_suffix(input_info.scalar_filename);

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

void MERGESHARP::set_mergesharp_data
(const INPUT_INFO & input_info, MERGESHARP_DATA & mergesharp_data, 
 MERGESHARP_TIME & mergesharp_time)
{
  PROCEDURE_ERROR error("set_mergesharp_data");

  if (!mergesharp_data.IsScalarGridSet()) {
    error.AddMessage
      ("Programming error. Scalar field must be set before set_mergesharp_data is called.");
    throw error;
  }

  // Set data structures in mergesharp_data
  mergesharp_data.Set(input_info);
}

void MERGESHARP::set_input_info
(const NRRD_INFO & nrrd_info, INPUT_INFO & input_info)
{
  input_info.grid_spacing.clear();
  for (int d = 0; d < nrrd_info.dimension; d++) {
    input_info.grid_spacing.push_back(nrrd_info.grid_spacing[d]);
  }
}

void MERGESHARP::set_output_info
(const INPUT_INFO & input_info,
 const int i, OUTPUT_INFO & output_info)
{

  // Set data structures in mergesharp_data
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

void MERGESHARP::set_color_alternating
(const SHARPISO_GRID & grid, const vector<VERTEX_INDEX> & cube_list,
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
void MERGESHARP::construct_gradient_filename
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

void MERGESHARP::OUTPUT_INFO::SetOutputTriMesh(const bool flag)
{
  flag_output_tri_mesh = flag;

  if (flag_output_tri_mesh)
    { num_vertices_per_isopoly = 3; }
  else
    { num_vertices_per_isopoly = 4; }
}

bool MERGESHARP::OUTPUT_INFO::OutputTriMesh() const
{
  return(flag_output_tri_mesh);
}

NUM_TYPE MERGESHARP::OUTPUT_INFO::NumVerticesPerIsopoly() const
{
  return(num_vertices_per_isopoly);
}

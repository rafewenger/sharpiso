/// \file mergesharpIO.h
/// IO classes and routines for mergesharp.

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

#ifndef _MERGESHARPIO_
#define _MERGESHARPIO_

#include <ctime>
#include <string>

#include "ijk.txx"
#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"
#include "sharpiso_eigen.h"
#include "ijkNrrd.h"

namespace MERGESHARP {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  typedef float COLOR_TYPE;           /// Color type.

  // **************************************************
  // NRRD INFORMATION
  // **************************************************

  /// NRRD information.
  class NRRD_INFO {

  public:
    int dimension;                        ///< Volume dimension.
    COORD_ARRAY grid_spacing;             ///< Grid spacing.

  public:
    NRRD_INFO();
    ~NRRD_INFO();

    void Clear();
  };
           
  // **************************************************
  // OUTPUT_FORMAT
  // **************************************************

  ///  Output format.
  typedef enum { OFF, IV } OUTPUT_FORMAT;

  typedef enum { SELECTED_ISOVERT, UNCOVERED_ISOVERT, ACTIVE_ISOVERT }
    MSHARP_ISOVERT_TYPE;


  // **************************************************
  // IO INFORMATION
  // **************************************************

  /// IO information
  class IO_INFO:public MERGESHARP_PARAM {

  protected:
    void Init();

  public:
    int dimension;
    COORD_ARRAY grid_spacing;
    const char * scalar_filename;       ///< Input scalar file name.
    const char * gradient_filename;     ///< Input gradient file name.

    /// Input edge-isosurface intersection normal file name.
    const char * normal_filename;       

    const char * output_filename;
    OUTPUT_FORMAT output_format;
    bool report_time_flag;
    bool use_stdout;
    bool nowrite_flag;
    bool flag_output_alg_info;    ///< Print algorithm information.
    bool flag_silent;
    bool flag_subsample;
    int subsample_resolution;
    bool flag_supersample;
    int supersample_resolution;
    bool flag_color_alternating; ///< Color simplices in alternating cubes
    int region_length;
    bool flag_output_param;      ///< Output algorithm parameters.
    bool flag_output_selected;   ///< Output information about selected cubes.
    bool flag_output_sharp;      ///< Output information about sharp cubes.
    bool flag_output_active;     ///< Output information about active cubes.
    bool flag_output_isovert;    ///< Output isosurface vertices to separate file.
    ISOVERT_TYPE output_isovert_type;
    MSHARP_ISOVERT_TYPE output_isovert_type2;
    std::string output_isovert_filename;
    std::vector<COORD_TYPE> minc;
    std::vector<COORD_TYPE> maxc;

  public:
    IO_INFO() { Init(); };
    ~IO_INFO() { Init(); };

    void Set(const IO_INFO & io_info);
  };


  // **************************************************
  // INPUT INFORMATION
  // **************************************************

  /// Input information
  class INPUT_INFO:public IO_INFO {

  protected:
    void Init();
    void Clear();

  public:
    bool is_vertex_position_method_set;
    bool is_use_sharp_edgeI_set;
    bool is_conflict_set;
    SCALAR_ARRAY isovalue;        ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    std::string isotable_directory;

    /// List of high resolution arguments,
    ///   e.g., "-highres {coord list}".
    std::vector<std::string> high_resolution_option;

  public:
    INPUT_INFO() { Init(); };
    ~INPUT_INFO() { Clear(); };
  };

  // **************************************************
  // OUTPUT INFORMATION
  // **************************************************

  /// Output information.
  class OUTPUT_INFO:public IO_INFO {

  protected:
    bool flag_output_tri_mesh;
    int num_vertices_per_isopoly;

    void Init();

  public:
    std::string output_filename;
    SCALAR_TYPE isovalue;
    int grow_factor;
    int shrink_factor;

    OUTPUT_INFO() { Init(); };
    ~OUTPUT_INFO() { Init(); };

    // Set functions
    void SetOutputTriMesh(const bool flag);

    // Get functions
    bool OutputTriMesh() const;
    NUM_TYPE NumVerticesPerIsopoly() const;
  };

  // **************************************************
  // TIMING FUNCTIONS/CLASSES
  // **************************************************

  /// Elapsed CPU time.
  class ELAPSED_CPU_TIME {

  protected:
    clock_t t;

  public:
    ELAPSED_CPU_TIME() { t = clock(); };

    clock_t getElapsed() {
      clock_t old_t = t;
      t = clock();
      return(t - old_t);
    };
  };

  /// Elapsed wall time.
  class ELAPSED_TIME {

  protected:
    time_t t;

  public:
    ELAPSED_TIME() { time(&t);  };

    double getElapsed() {
      time_t old_t = t;
      time(&t);
      return(difftime(t,old_t));
    };
  };

  /// IO time.
  struct IO_TIME {
    double read_table_time; ///< Wall time to read isosurface lookup table.
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };

  // **************************************************
  // PARSE COMMAND LINE
  // **************************************************

  /// Parse the next option in the command line.
  /// Return false if no next option or parse fails.
  bool parse_command_option
  (const int argc, char **argv, const int iarg, int & next_arg,
   INPUT_INFO & input_info);

  /// Parse the isovalue(s) and filename.
  void parse_isovalue_and_filename
  (const int argc, char **argv, const int iarg, INPUT_INFO & input_info);

  /// Parse the command line.
  /// Command line is control parameters, followed by one or more isovalues,
  ///   followed by input file name.
  void parse_command_line(int argc, char **argv, INPUT_INFO & input_info);

  /// Check input information in input_info
  bool check_input
  (const INPUT_INFO & input_info, 
   const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   IJK::ERROR & error);

  /// Set input_info defaults.
  void set_input_info_defaults(INPUT_INFO & input_info);

  /// Check input_info.
  void check_input_info(const INPUT_INFO & input_info);


  // **************************************************
  // READ NEARLY RAW RASTER DATA (nrrd) FILE
  // **************************************************

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
  (const char * input_filename, SHARPISO_SCALAR_GRID & scalar_grid, 
   NRRD_INFO & nrrd_info, IO_TIME & io_time);

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
  (const char * input_filename, SHARPISO_SCALAR_GRID & scalar_grid, 
   NRRD_INFO & nrrd_info);

  /// Read a nearly raw raster gradient data (nrrd) file.
  void read_nrrd_file
  (const char * input_filename, GRADIENT_GRID & gradient_grid, 
   NRRD_INFO & nrrd_info);

  // **************************************************
  // READ OFF FILE
  // **************************************************

  /// Read OFF file.
  void read_off_file
  (const char * input_filename,
   std::vector<COORD_TYPE> & coord, 
   std::vector<GRADIENT_COORD_TYPE> & normal, 
   std::vector<int> & simplex_vert);


  /// Read OFF file.  Ignore simplex vert.
  void read_off_file
  (const char * input_filename,
   std::vector<COORD_TYPE> & coord, std::vector<GRADIENT_COORD_TYPE> & normal);

  // **************************************************
  // OUTPUT ISOSURFACE
  // **************************************************

  /// Output dual isosurface.
  void output_dual_isosurface
  (const OUTPUT_INFO & output_info, const MERGESHARP_DATA & mergesharp_data,
   const DUAL_ISOSURFACE & dual_isosurface, const ISOVERT & isovert, 
   const MERGESHARP_INFO & mergesharp_info, IO_TIME & io_time);

  // **************************************************
  // RESCALE ROUTINES
  // **************************************************

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
  void rescale_vertex_coord
  (const OUTPUT_INFO & output_info, std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale vertex coordinates by grow and shrink factor and by grid_spacing.
  /// Precondition: grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
  (const int grow_factor, const int shrink_factor,
   const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);

  // **************************************************
  // WRITE_DUAL_MESH
  // **************************************************

  /// Write dual mesh.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] Array of triangle vertices.
  /// @param quad_vert[] Array of triangle vertices.
  /// @param flag_reorder_vertices If true, need to reorder quad vertices
  ///           for clockwise/counter-clockwise order.
  void write_dual_mesh3D
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & tri_vert,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const bool flag_reorder_quad_vertices);


  // Write dual mesh.
  // Return time to write data.
  void write_dual_mesh3D
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & tri_vert,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const bool flag_reorder_quad_vertices,
   IO_TIME & io_time);

  /// Write dual mesh.
  /// Version with DUAL_ISOSURFACE parameter.
  void write_dual_mesh3D
  (const OUTPUT_INFO & output_info,
   const DUAL_ISOSURFACE & dual_isosurface,
   const bool flag_reorder_quad_vertices,
   IO_TIME & io_time);

  // **************************************************
  // SET ROUTINES
  // **************************************************

  /// Set mergesharp_data based on input_info.
  /// Precondition: Scalar field in mergesharp_data must be set before
  ///   this routines is called.
  void set_mergesharp_data
  (const INPUT_INFO & input_info, MERGESHARP_DATA & mergesharp_data, 
   MERGESHARP_TIME & mergesharp_time);

  /// Copy nrrd_info into input_info.
  void set_input_info
  (const NRRD_INFO & nrrd_info, INPUT_INFO & input_info);

  /// Set output_info based on isotable, input_info and isovalue index i.
  void set_output_info
  (const INPUT_INFO & input_info, 
   const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
  (const SHARPISO_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
   COLOR_TYPE * color);

  // **************************************************
  // CONSTRUCT FILENAME
  // **************************************************

  /// Construct gradient filename from scalar filename.
  void construct_gradient_filename
  (const char * scalar_filename, std::string & gradient_filename);

  // **************************************************
  // REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

  /// Report number of grid cubes (and number of subsampled grid cubes.)
  void report_num_cubes
  (const SHARPISO_GRID & full_grid, const INPUT_INFO & input_info, 
   const MERGESHARP_DATA & mergesharp_data);

  /// Report mergesharp parameters.
  void report_mergesharp_param(const MERGESHARP_PARAM & mergesharp_param);

  void report_iso_info3D
  (const OUTPUT_INFO & output_info, const MERGESHARP_DATA & mergesharp_data,
   const DUAL_ISOSURFACE & dual_isosurface,
   const MERGESHARP_INFO & mergesharp_info);

  /// Report information about isosurface vertices
  void report_isovert_info
  (const OUTPUT_INFO & output_info, const SHARPISO_GRID & grid,
   const ISOVERT & isovert);

  /// Report information about cubes containing corner or edge iso vertices.
  void report_sharp_cubes
  (const SHARPISO_GRID & grid, const ISOVERT & isovert);

  /// Report information about active cubes.
  /// A cube is active if it has some bipolar edge.
  void report_active_cubes
  (const SHARPISO_GRID & grid, const ISOVERT & isovert);


  // **************************************************
  // REPORT TIMING INFORMATION
  // **************************************************

  void report_mergesharp_time
  (const INPUT_INFO & input_info, const MERGESHARP_TIME & mergesharp_time, 
   const char * mesh_type_string);

  void report_time
  (const INPUT_INFO & input_info, const IO_TIME & io_time, 
   const MERGESHARP_TIME & mergesharp_time, const double total_elapsed_time);

  // **************************************************
  // WRITE ISOSURFACE VERTICES OR VERTEX INFORMATION TO FILE
  // **************************************************

  void write_isovert_info
  (const OUTPUT_INFO & output_info,
   const std::vector<DUAL_ISOVERT_INFO> & isovert_info);

  void write_isovert
  (const OUTPUT_INFO & output_info, const ISOVERT & isovert);


  // **************************************************
  // USAGE/HELP MESSAGES
  // **************************************************

  void usage_error(const char * command_name);
  void help(const char * command_name);
}

#endif

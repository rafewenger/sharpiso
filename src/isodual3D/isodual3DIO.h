/// \file isodual3DIO.h
/// IO classes and routines for isodual3D.

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

#ifndef _ISODUAL3DIO_
#define _ISODUAL3DIO_

#include <ctime>
#include <string>

#include "ijk.txx"
#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"
#include "sharpiso_eigen.h"
#include "ijkNrrd.h"

namespace ISODUAL3D {

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

  typedef enum { OFF, IV } OUTPUT_FORMAT;   ///< Output format.

// **************************************************
// IO INFORMATION
// **************************************************
  /// IO information
  class IO_INFO:public ISODUAL_PARAM {

  protected:
    void Init();

  public:
    int dimension;
    COORD_ARRAY grid_spacing;
    char * scalar_filename;       ///< Input scalar file name.
    char * gradient_filename;     ///< Input gradient file name.
    char * output_filename;
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
    bool flag_color_alternating;  ///< Color simplices in alternating cubes
    int region_length;
    bool flag_output_param;

  public:
    IO_INFO() { Init(); };
    ~IO_INFO() { Init(); };

    void Set(const IO_INFO & io_info);
  };

  /// IO information
  class INPUT_INFO:public IO_INFO {

  protected:
    void Init();
    void Clear();

  public:
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

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, INPUT_INFO & input_info);

  /// Check input information in input_info
  bool check_input
    (const INPUT_INFO & input_info, 
     const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
     IJK::ERROR & error);

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, ISODUAL_SCALAR_GRID & scalar_grid, 
     NRRD_INFO & nrrd_info, IO_TIME & io_time);

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, ISODUAL_SCALAR_GRID & scalar_grid, 
     NRRD_INFO & nrrd_info);

  /// Read a nearly raw raster gradient data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, GRADIENT_GRID & gradient_grid, 
     NRRD_INFO & nrrd_info);

// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

  /// Output dual isosurface.
  void output_dual_isosurface
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

  /// Output isosurface of quadrilaterals.
  void output_quad_isosurface
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & quad_vert,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

  /// Output isosurface of triangles.
  void output_tri_isosurface
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & tri_vert,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

  void output_dual_isosurface_color
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

  void output_dual_isosurface_color
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

  void output_dual_isosurface_color_alternating
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const DUAL_ISOSURFACE & dual_isosurface,
     const ISODUAL_INFO & isodual_info, IO_TIME & io_time);

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

  /// Write dual isosurface quadrilateral mesh.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param quad_vert[] List of triangle vertices.
  ///        quad_vert[4*i+k] is k'th vertex of triangle i.
  ///        Quadrilateral vertices are listed in order:
  ///            Lower-Left, Lower-Right, Upper-Left, Upper-Right
  void write_dual_quad_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & quad_vert);

  void write_dual_quad_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time);

  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     IO_TIME & io_time);

  /// Write dual isosurface triangular mesh.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & quad_vert);

  /// Write dual isosurface triangular mesh.
  /// Record write time.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & tri_vert,
   IO_TIME & io_time);

// **************************************************
// SET ROUTINES
// **************************************************

  /// Set isodual_data based on input_info.
  /// Precondition: Scalar field in isodual_data must be set before
  ///   this routines is called.
  void set_isodual_data
    (const INPUT_INFO & input_info, ISODUAL_DATA & isodual_data, ISODUAL_TIME & isodual_time);

  /// Copy nrrd_info into input_info.
  void set_input_info
    (const NRRD_INFO & nrrd_info, INPUT_INFO & input_info);

  /// Set output_info based on isotable, input_info and isovalue index i.
  void set_output_info
    (const INPUT_INFO & input_info, 
     const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
    (const ISODUAL_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
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
    (const ISODUAL_GRID & full_grid, const INPUT_INFO & input_info, 
     const ISODUAL_DATA & isodual_data);

  /// Report isodual parameters.
  void report_isodual_param(const ISODUAL_PARAM & isodual_param);

  /// Report information about isosurface.
  void report_iso_info
    (const OUTPUT_INFO & output_info, const ISODUAL_DATA & isodual_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const ISODUAL_INFO & isodual_info);

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

  void report_isodual_time
    (const INPUT_INFO & input_info, const ISODUAL_TIME & isodual_time, 
     const char * mesh_type_string);

  void report_time
    (const INPUT_INFO & input_info, const IO_TIME & io_time, 
     const ISODUAL_TIME & isodual_time, const double total_elapsed_time);

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

  void usage_error();
  void help();
}

#endif

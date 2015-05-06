/// \file ijkscalarinfo.cxx
/// get info from scalar grid data
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008,2009,2011,2012 Rephael Wenger

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
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkcoord.txx"
#include "ijkcontourtree.txx"
#include "ijkprint.txx"
#include "ijkscalar_grid.txx"
#include "ijkinterpolate.txx"
#include "ijkscalarinfo.h"
#include "ijkNrrd.h"

using namespace std;
using namespace IJK;
using namespace IJKSCALARINFO;
using IJKDATATABLE::compute_num_buckets;
using IJKDATATABLE::compute_bucket;

// global variables
char * input_filename = NULL;
char * output_filename = NULL;
char * file_label = NULL;
bool flag_scalar_freq_table = false;
bool flag_isoarea_table = false;
bool flag_isoarea_edge_table = false;
bool flag_ivol_edge_table = false;
bool flag_ivol_table = false;
bool flag_total_gradient_table = false;
bool flag_mean_gradient_table = false;
bool flag_zero_gradient_table = false;
bool flag_num_components_table = false;
bool flag_box_dim = false;
bool flag_edge_dim = false;
bool flag_edge = false;
bool flag_cube = false;
bool flag_vertex_info = false;
bool flag_vertex_index = false;
bool flag_cube_info = false;
bool flag_plot = false;
bool flag_normalize = false;
bool flag_standard_deviation = false;  // if true, compute gradient standard deviation
bool flag_edge_span_freq = false;      // compute edge span frequencies
bool flag_cube_span_freq = false;      // compute cube span frequencies
INTERVAL_TYPE bucket_interval = 1.0;
bool flag_interval = false;
int num_buckets = 100;
bool flag_num_buckets = false;
bool flag_laplace = false;
bool flag_sum_from0 = false;
bool flag_subsample = false;
bool flag_silent_write = false; // if true, suppress message "Writing table..."
int subsample_resolution = 2;
int vertex_index = 0;
std::vector<GRID_COORD_TYPE> cube_coord;
std::vector<GRID_COORD_TYPE> vertex_coord;
int DEFAULT_TABLE_COLUMN_WIDTH = 8;
int DEFAULT_TABLE_PRECISION = 4;
int max_output_local_max_size = 10;
SCALAR_TYPE min_scalar = 0, max_scalar = 0;
bool is_min_scalar_set = false;
bool is_max_scalar_set = false;

// compute frequency routines
void compute_frequency_tables
(const SGRID & scalar_grid, const int scalar_type);
void compute_gradient_distribution_tables(const SGRID & scalar_grid);
template <class ITYPE>
NUM_TYPE compute_num_rows(const SGRID & scalar_grid, const ITYPE interval);
SCALAR_TYPE round_down_min_scalar
(const SCALAR_TYPE min_scalar, const SCALAR_TYPE max_scalar);
template <class TABLE_TYPE, class ITYPE>
void set_scalar_column(const ITYPE interval, TABLE_TYPE & table);
template <class TABLE_TYPE>
void compute_scalar_frequency
(const BUCKET_GRID & bucket_grid, TABLE_TYPE & table);
void measure_isoarea_using_edges
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void measure_isoarea_using_cubes
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void measure_ivol_using_edges
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 const INTERVAL_TYPE interval, DATA_TABLE & table);
void measure_ivol_using_cubes
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 const INTERVAL_TYPE interval, DATA_TABLE & table);
void compute_total_gradient_edge
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table);
void compute_total_gradient_cube
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void compute_total_gradient_Laplace
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table);
void compute_total_gradient
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table);
void compute_mean_gradient
(const DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> & total_gradient,
 const DATA_COLUMN<NUM_TYPE, NUM_TYPE> & iso_area,
 DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> & mean_gradient,
 const int num_rows);
void compute_gradient_standard_deviation_edge
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void compute_gradient_standard_deviation_cube
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void compute_zero_gradient
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void compute_num_components
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table);
void compute_box_dim
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table);
void compute_edge_dim
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table);
void compute_num_zero_cubes
(const SGRID & scalar_grid, NUM_TYPE & num_zero_cubes);

// compute other information
void compute_min_scalar_interval
(const SGRID & scalar_grid, SCALAR_TYPE & min_interval, bool & flag_uniform);

// output routines
void report_info(const SGRID & scalar_grid, const int scalar_type);
void output_grid_info(ostream & output, const SGRID & scalar_grid);
void output_vertex_info(ostream & output, const SGRID & scalar_grid,
                        const int vertex_index);
void output_cube_info(ostream & output, const SGRID & scalar_grid,
                      const std::vector<GRID_COORD_TYPE> & cube_coord);
template <class NTYPE_X, class XTYPE, class NTYPE_Y, class YTYPE>
void output_local_max
(ostream & out, const DATA_COLUMN<NTYPE_X,XTYPE> & x,
 const DATA_COLUMN<NTYPE_Y,YTYPE> & y);
template <class NTYPE_X, class XTYPE, class NTYPE_Y, class YTYPE,
          class NTYPE>
void output_data
(ostream & out, const DATA_COLUMN<NTYPE_X,XTYPE> & x,
 const DATA_COLUMN<NTYPE_Y,YTYPE> & y,
 const NTYPE istart, const NTYPE interval, const NTYPE precision);
template <class NTYPE, class DTYPE>
void output_average
(ostream & out, const AVERAGE_STAT<NTYPE,DTYPE> average, const char * label);
template <class NTYPE> NTYPE choose_interval(const NTYPE num_rows);
template <class TABLE_TYPE>
void write_table_gplt(const string & filename, const TABLE_TYPE & table);
template <class TABLE_TYPE>
void write_table_gplt(const string & filename_prefix, 
                      const string & filename_suffix,
                      const TABLE_TYPE & table);
template <class TABLE_TYPE>
void write_table_gplt(ofstream & ofile, const TABLE_TYPE & scalar_freq);
bool some_table_flag_is_true();

typedef enum
  {FREQ_PARAM, EDGE_PARAM, CUBE_PARAM, GRADIENT_PARAM, ZERO_GRADIENT_PARAM,
   STANDARD_DEVIATION_PARAM, COMPONENTS_PARAM,  ALL_PARAM, NONE_PARAM,
   IVOL_PARAM, ISO_PARAM, LAPLACE_PARAM, 
   EDGE_SPAN_PARAM, CUBE_SPAN_PARAM,
   BOX_DIM_PARAM, EDGE_DIM_PARAM,
   MIN_SCALAR_PARAM, MAX_SCALAR_PARAM,
   SUM_FROM0, SUM_TO0,
   PLOT_PARAM, NORMALIZE_PARAM, INTERVAL_PARAM, NUM_BUCKETS_PARAM, NUMB_PARAM,
   VERTEX_PARAM, VC_PARAM, CC_PARAM,
   SUBSAMPLE_PARAM, OUTPUT_FILENAME_PARAM, FILE_LABEL_PARAM,
   SILENT_WRITE_PARAM,
   HELP_PARAM, UNKNOWN_PARAM} PARAMETER;
const char * parameter_string[] = 
  {"-freq", "-edge", "-cube", "-gradient", "-zero_gradient", 
   "-deviation", "-components", "-all", "-none",
   "-ivol", "-iso", "-laplace", "-edge_span", "-cube_span",
   "-box_dim", "-edge_dim",
   "-min_scalar", "-max_scalar", "-sum_from0", "-sum_to0",
   "-plot", "-normalize", "-interval", "-num_buckets", "-numb",
   "-vertex", "-vc", "-cc", "-subsample", "-o", "-file_label", 
   "-silent_write", "-help", "-unknown"};

// check routines
template <class ITYPE>
bool check_interval(const ITYPE interval, ERROR & error);
template <class NUM_TYPE>
bool check_spread(const NUM_TYPE spread, ERROR & error);

// misc routines
void split_string(const string & s, const char c,
                  string & prefix, string & suffix);
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  IJK::ERROR error("ijkscalarinfo");

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    assert(input_filename != NULL);

    int dimension = 0;
    Nrrd *nin;
    int scalar_type = 0;

    /* get scalar field data from nrrd file */
    nin = nrrdNew();
    if (nrrdLoad(nin, input_filename, NULL)) {
      char *err = biffGetDone(NRRD);
      cerr << "Error reading: " << input_filename << endl;
      cerr << "  Error: " << err << endl;
      exit(35);
    };
    dimension = nin->dim;
    scalar_type = nin->type;

    if (dimension < 1) {
      cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
      exit(20);
    };

    IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);
    size_t size[NRRD_DIM_MAX];

    nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size); 

    for (int d = 0; d < dimension; d++) { axis_size[d] = size[d]; }

    SGRID scalar_grid(dimension, axis_size.PtrConst());
    nrrd2scalar(nin, scalar_grid.ScalarPtr());

    if (!flag_subsample) {
      report_info(scalar_grid, scalar_type);
    }
    else {
      // subsample grid
      SGRID subsampled_grid;
      subsampled_grid.Subsample(scalar_grid, subsample_resolution);

      report_info(subsampled_grid, scalar_type);
    }

    nrrdNuke(nin);
  }
  catch(IJK::ERROR error) {
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
// Report Information
// **************************************************

void report_info(const SGRID & scalar_grid, const int scalar_type)
{
  output_grid_info(cout, scalar_grid);

  if (flag_vertex_info && !flag_subsample) {
    if (!flag_vertex_index) 
      { vertex_index = scalar_grid.ComputeVertexIndex(vertex_coord); }

    output_vertex_info(cout, scalar_grid, vertex_index); 
  };

  if (flag_cube_info && !flag_subsample)
    { output_cube_info(cout, scalar_grid, cube_coord); };

  if (some_table_flag_is_true()) {

    SCALAR_TYPE min_scalar_interval;
    bool flag_uniform;
    compute_min_scalar_interval(scalar_grid, min_scalar_interval, flag_uniform);

    if (!flag_uniform) {
      compute_frequency_tables(scalar_grid, scalar_type);
    }
  }

  if (flag_edge_span_freq || flag_cube_span_freq) {
    compute_gradient_distribution_tables(scalar_grid);
  }

}


// **************************************************
// Compute Buckets
// **************************************************

/// Compute bucket containing each vertex.
/// Scalar values are divided into buckets.
/// bucket_grid.Scalar(iv0) = Index of bucket containing scalar value of vertex iv0.
template <class ITYPE>
void compute_vertex_buckets
(const SGRID & scalar_grid, const ITYPE interval, const SCALAR_TYPE min_scalar,
 const NUM_TYPE num_buckets, BUCKET_GRID & bucket_grid)
{
  const SCALAR_TYPE * scalar_ptr = scalar_grid.ScalarPtrConst();
  BUCKET_TYPE * bucket_ptr = bucket_grid.ScalarPtr();
  while (scalar_ptr != scalar_grid.End()) {
    const SCALAR_TYPE s = *scalar_ptr;
    NUM_TYPE ibucket = compute_bucket(s, min_scalar, interval, num_buckets);
    *bucket_ptr = ibucket;

    scalar_ptr++;
    bucket_ptr++;
  }

}

// **************************************************
// Compute Derivatives
// **************************************************

/// Compute second derivative.
/// Precondition: vertex is not on boundary of grid.
template <class RTYPE, class ATYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE SecondDerivative
(const ATYPE axis_index, const ITYPE * increment, const VTYPE iv,
 const STYPE * scalar)
{
  VTYPE iv_prev = iv - increment[axis_index];
  VTYPE iv_next = iv + increment[axis_index];
  RTYPE s_prev = scalar[iv_prev];
  RTYPE s = scalar[iv];
  RTYPE s_next = scalar[iv_next];
  RTYPE second_derivative = s_prev - 2*s + s_next;

  return(second_derivative);
}

/// Compute second derivative.
/// Precondition: vertex is not on boundary of grid.
template <class ATYPE, class ITYPE, class VTYPE>
inline float SecondDerivative
(const ATYPE axis_index, const ITYPE * increment, const VTYPE iv,
 const float * scalar)
// Precondition: vertex is not on boundary of grid
{
  return(SecondDerivative<float>(axis_index, increment, iv, scalar));
}

/// Compute Laplacian.
/// Precondition: vertex is not on boundary of grid.
template <class RTYPE, class DTYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE Laplace
(const DTYPE dimension, const ITYPE * increment, const VTYPE iv,
 const STYPE * scalar)
// Precondition: vertex is not on boundary of grid
{
  RTYPE result = 0;
  for (DTYPE d = 0; d < dimension; d++) {
    RTYPE second_derivative = 
      SecondDerivative<RTYPE>(d, increment, iv, scalar);
    result += second_derivative;
  }
  return(result);
}

/// Compute Laplacian.
/// Precondition: vertex is not on boundary of grid.
template <class DTYPE, class ITYPE, class VTYPE>
inline float Laplace
(const DTYPE dimension, const ITYPE * increment, const VTYPE iv,
 const float * scalar)
// Precondition: vertex is not on boundary of grid
{
  return(Laplace<float>(dimension, increment, iv, scalar));
}

/// Compute second derivative (boundary vertices.)
/// Precondition: Axis size is at least 3.
template <class RTYPE, class ATYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE SecondDerivativeBnd
(const ATYPE axis_index, const ITYPE * extended_increment,  const VTYPE iv, 
 const STYPE * scalar)
{
  ITYPE i0 = extended_increment[axis_index];
  ITYPE i1 = extended_increment[axis_index+1];

  VTYPE k = iv % i1;
  VTYPE jv = iv;
  if (k < i0) {
    jv = iv + i0;
  }
  else if (k >= i1 - i0) {
    jv = iv - i0;
  };

  return(SecondDerivative<RTYPE>(axis_index, extended_increment, jv, scalar));
}

/// Compute second derivative (boundary vertices.)
/// Precondition: Axis size is at least 3.
template <class ATYPE, class ITYPE, class VTYPE, class STYPE>
inline float SecondDerivativeBnd
(const ATYPE axis_index, const ITYPE * extended_increment,  const VTYPE iv, 
 const STYPE * scalar)
{
  return(SecondDerivativeBnd<float>
         (axis_index, extended_increment, iv, scalar));
}

/// Compute Laplacian (boundary vertices.)
/// Precondition: All axis sizes are at least 3.
template <class RTYPE, class DTYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE LaplaceBnd
(const DTYPE dimension, const ITYPE * extended_increment, const VTYPE iv,
 const STYPE * scalar)
{
  RTYPE result = 0;
  for (DTYPE d = 0; d < dimension; d++) {
    RTYPE second_derivative = 
      SecondDerivativeBnd<RTYPE>(d, extended_increment, iv, scalar);
    result += second_derivative;
  }
  return(result);
}

/// Compute Laplacian (boundary vertices.)
/// Precondition: All axis sizes are at least 3.
template <class DTYPE, class ITYPE, class VTYPE>
inline float LaplaceBnd
(const DTYPE dimension, const ITYPE * extended_increment, const VTYPE iv,
 const float * scalar)
{
  return(LaplaceBnd<float>(dimension, extended_increment, iv, scalar));
}

/// Compute second derivative using Richardson Extrapolation.
/// Precondition: vertex is not on boundary or adjacent to boundary.
template <class RTYPE, class ATYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE SecondDerivativeRichardsonExtrapolation
(const ATYPE axis_index, const ITYPE * increment, const VTYPE iv2,
 const STYPE * scalar)
{
  VTYPE inc = increment[axis_index];
  VTYPE inc2 = 2*inc;

  // s_i is scalar value of vertex v_i where order of vertices 
  //   along the axis is: v0, v1, v2, v3, v4
  RTYPE s0 = scalar[iv2-inc2];
  RTYPE s1 = scalar[iv2-inc];
  RTYPE s2 = scalar[iv2];
  RTYPE s3 = scalar[iv2+inc];
  RTYPE s4 = scalar[iv2+inc2];

  RTYPE sderiv0 = (s0 - 2*s2 + s4)/4.0;
  RTYPE sderiv1 = s1 - 2*s2 + s3;
  RTYPE second_derivative = (sderiv1*4.0 - sderiv0)/3.0;
  
  return(second_derivative);
}

/// Compute second derivative using Richardson Extrapolation.
/// Precondition: vertex is not on boundary or adjacent to boundary.
template <class ATYPE, class ITYPE, class VTYPE>
inline float SecondDerivativeRichardsonExtrapolation
(const ATYPE axis_index, const ITYPE * increment, const VTYPE iv,
 const float * scalar)
{
  return(SecondDerivativeRichardsonExtrapolation<float>
         (axis_index, increment, iv, scalar));
}

/// Compute Laplacian.
/// Precondition: vertex is not on boundary or adjacent to boundary.
template <class RTYPE, class DTYPE, class ITYPE, class VTYPE, class STYPE>
inline RTYPE LaplaceInterior2
(const DTYPE dimension, const ITYPE * increment, const VTYPE iv,
 const STYPE * scalar)
{
  RTYPE result = 0;
  for (DTYPE d = 0; d < dimension; d++) {
    RTYPE second_derivative = 
      SecondDerivativeRichardsonExtrapolation<RTYPE>(d, increment, iv, scalar);
    result += second_derivative;
  }
  return(result);
}

/// Compute Laplacian.
/// Precondition: vertex is not on boundary or adjacent to boundary.
template <class DTYPE, class ITYPE, class VTYPE>
inline float LaplaceInterior2
(const DTYPE dimension, const ITYPE * increment, const VTYPE iv,
 const float * scalar)
{
  return(LaplaceInterior2<float>(dimension, increment, iv, scalar));
}

// **************************************************
// Compute Frequencies
// **************************************************

void compute_frequency_tables
(const SGRID & scalar_grid, const int scalar_type)
{
  const NUM_TYPE dimension = scalar_grid.Dimension();
  IJK::PROCEDURE_ERROR error("compute_frequency_tables");

  SCALAR_TYPE min_scalar = scalar_grid.FindMinScalar();
  SCALAR_TYPE max_scalar = scalar_grid.FindMaxScalar();

  min_scalar = round_down_min_scalar(min_scalar, max_scalar);

  NUM_TYPE num_rows = 1;
  if (flag_interval) {

    if (!check_interval(bucket_interval, error)) { throw error; };
    num_rows = 
      compute_num_buckets<NUM_TYPE>(min_scalar, max_scalar, bucket_interval);
  }
  else if (!flag_num_buckets &&
           (scalar_type != nrrdTypeFloat &&
            scalar_type != nrrdTypeDouble)) {
    bucket_interval = 1;
    num_rows = compute_num_rows(scalar_grid, bucket_interval);
  }
  else {
    num_rows = num_buckets;

    if (num_buckets <= 1) { 
      bucket_interval = max_scalar - min_scalar+1;
    }
    else {
      bucket_interval = (max_scalar-min_scalar)/(num_buckets-1);
    }

    double min_floor = floor(double(min_scalar));
    if (min_scalar - bucket_interval < min_floor) {
      // snap min_scalar to min_floor;
      min_scalar = SCALAR_TYPE(min_floor);
    };

    double max_ceil = ceil(double(max_scalar));
    if (max_scalar + bucket_interval > max_ceil) {
      // snap max_scalar to max_ceil;
      max_scalar = SCALAR_TYPE(max_ceil);
    }

    if (num_buckets <= 1) { 
      bucket_interval = max_scalar - min_scalar+1;
    }
    else {
      bucket_interval = (max_scalar-min_scalar)/(num_buckets-1);
    }

  }

  DATA_TABLE table(num_rows, min_scalar);

  set_scalar_column(bucket_interval, table);

  BUCKET_GRID bucket_grid(dimension, scalar_grid.AxisSize());
  compute_vertex_buckets(scalar_grid, bucket_interval, min_scalar, num_rows, 
                         bucket_grid);

  if (flag_scalar_freq_table) 
    { compute_scalar_frequency(bucket_grid, table); };

  if (flag_isoarea_table || flag_isoarea_edge_table ||
      flag_mean_gradient_table || flag_edge_dim || flag_box_dim) {

    if (flag_edge || flag_edge_dim || 
        flag_isoarea_edge_table || flag_laplace)
      { measure_isoarea_using_edges(scalar_grid, bucket_grid, table); };

    if (flag_cube || flag_box_dim)
      { measure_isoarea_using_cubes(scalar_grid, bucket_grid, table); }
  }

  if (flag_ivol_table) {
    if (flag_edge || flag_ivol_edge_table)
      { measure_ivol_using_edges
          (scalar_grid, bucket_grid, bucket_interval, table); };
    if (flag_cube)
      { measure_ivol_using_cubes
          (scalar_grid, bucket_grid, bucket_interval, table); };
  }

  if (flag_total_gradient_table || flag_mean_gradient_table) {
    compute_total_gradient(scalar_grid, bucket_grid, table); 
  };


  if (flag_mean_gradient_table) {
    if (flag_edge) { 
      compute_mean_gradient
        (table.total_gradient_edge, table.isoarea_edge, 
         table.mean_gradient_edge, num_rows);

      if (flag_standard_deviation) {
        compute_gradient_standard_deviation_edge
          (scalar_grid, bucket_grid, table);
      }
    }

    if (flag_cube) { 
      compute_mean_gradient
        (table.total_gradient_cube, table.isoarea_cube, 
         table.mean_gradient_cube, num_rows);

      if (flag_standard_deviation) {
        compute_gradient_standard_deviation_cube
          (scalar_grid, bucket_grid, table);
      }
    }

    if (flag_laplace) { 
      compute_mean_gradient
        (table.total_gradient_laplace, table.isoarea_edge, 
         table.mean_gradient_laplace, num_rows);

      if (flag_standard_deviation) {
        /* NEED TO COMPLETE
         */
      }

    }

  }

  if (flag_edge_dim) 
    { compute_edge_dim(scalar_grid, bucket_grid, table); }

  if (flag_box_dim) 
    { compute_box_dim(scalar_grid, bucket_grid, table); }

  if (flag_zero_gradient_table)
    { compute_zero_gradient(scalar_grid, bucket_grid, table); }

  if (flag_num_components_table)
    { compute_num_components(scalar_grid, bucket_grid, table); }

  // Compute table sums
  table.ComputeSum();

  // Output local max
  if (flag_scalar_freq_table) 
    { output_local_max(cout, table.scalar, table.scalar_freq); }

  if (flag_isoarea_table) {
    if (flag_edge || flag_isoarea_edge_table) 
      { output_local_max(cout, table.scalar, table.isoarea_edge); }
    if (flag_cube)
      { output_local_max(cout, table.scalar, table.isoarea_cube); }
  }

  if (flag_ivol_table) {
    if (flag_edge || flag_ivol_edge_table) 
      { output_local_max(cout, table.scalar, table.ivol_edge); }
    if (flag_cube) 
      { output_local_max(cout, table.scalar, table.ivol_cube); }
  }

  if (flag_total_gradient_table) {
    if (flag_edge) 
      { output_local_max(cout, table.scalar, table.total_gradient_edge); };
    if (flag_cube) 
      { output_local_max(cout, table.scalar, table.total_gradient_cube); };
    if (flag_laplace) 
      { output_local_max(cout, table.scalar, table.total_gradient_laplace); };
  }

  if (flag_mean_gradient_table) {

    if (flag_edge) 
      { output_local_max(cout, table.scalar, table.mean_gradient_edge); }

    if (flag_cube) 
      { output_local_max(cout, table.scalar, table.mean_gradient_cube); }

    if (flag_laplace) 
      { output_local_max(cout, table.scalar, table.mean_gradient_laplace); }
  }

  if (flag_edge_dim) {
    int interval = choose_interval(table.NumRows());
    output_data(cout, table.scalar, table.edge_dim, interval, interval, 2);
  }

  if (flag_box_dim) {
    int interval = choose_interval(table.NumRows());
    output_data(cout, table.scalar, table.box_dim, interval, interval, 2);
  }

  if (flag_num_components_table) {
    output_local_max(cout, table.scalar, table.num_components);
  }

  cout << endl;

  if (flag_plot) {

    if (output_filename != NULL) {

      if (!flag_edge) { table.isoarea_edge.Hide(); };
      if (!flag_total_gradient_table) { 
        table.total_gradient_edge.Hide(); 
        table.total_gradient_cube.Hide(); 
        table.total_gradient_laplace.Hide(); 
      };

      write_table_gplt(output_filename, table);
    }
    else {

      string output_prefix;
      string output_suffix = ".gplt";

      // create output filename
      string fname = input_filename;

      // remove path from file name
      string prefix, suffix;
      split_string(fname, '/', prefix, suffix);
      if (suffix != "") { fname = suffix; }
      split_string(fname, '.', prefix, suffix);
      if (suffix == "nrrd" || suffix == "nhdr") {
        output_prefix = prefix;
      }
      else {
        output_prefix = input_filename;
      }

      if (file_label != NULL)
        { output_suffix = string(".") + file_label + output_suffix; }

      if (flag_scalar_freq_table) {
        table.HideAllExceptScalar();
        table.scalar_freq.Show();
        write_table_gplt(output_prefix, "scalar_freq"+output_suffix, table);
      };

      if (flag_isoarea_table) {

        if (flag_edge || flag_isoarea_edge_table) {
          table.HideAllExceptScalar();
          table.isoarea_edge.Show();
          write_table_gplt(output_prefix, "isoarea.edge"+output_suffix, table);
        }

        if (flag_cube) {
          table.HideAllExceptScalar();
          table.isoarea_cube.Show();
          write_table_gplt(output_prefix, "isoarea.cube"+output_suffix, table);
        }
      }

      if (flag_ivol_table) {

        if (flag_edge || flag_ivol_edge_table) {
          table.HideAllExceptScalar();
          table.ivol_edge.Show();
          write_table_gplt(output_prefix, "ivol.edge"+output_suffix, table);
        }

        if (flag_cube) {
          table.HideAllExceptScalar();
          table.ivol_cube.Show();
          write_table_gplt(output_prefix, "ivol.cube"+output_suffix, table);
        }
      }

      if (flag_total_gradient_table) {

        if (flag_edge) {
          table.HideAllExceptScalar();
          table.total_gradient_edge.Show();
          write_table_gplt(output_prefix, "total_gradient.edge"+output_suffix, table);
        }

        if (flag_cube) {
          table.HideAllExceptScalar();
          table.total_gradient_cube.Show();
          write_table_gplt(output_prefix, "total_gradient.cube"+output_suffix, table);
        }

        if (flag_laplace) {
          table.HideAllExceptScalar();
          table.total_gradient_laplace.Show();
          write_table_gplt(output_prefix, "total_gradient.laplace"+output_suffix, table);
        }
      };

      if (flag_mean_gradient_table) {

        if (flag_edge) {
          table.HideAllExceptScalar();
          table.mean_gradient_edge.Show();
          write_table_gplt(output_prefix, "mean_gradient.edge"+output_suffix, table);

          if (flag_standard_deviation) {
            table.HideAllExceptScalar();
            table.gradient_standard_deviation_edge.Show();
            write_table_gplt(output_prefix, "gradient_deviation.edge"+output_suffix, table);
          }
        }

        if (flag_cube) {
          table.HideAllExceptScalar();
          table.mean_gradient_cube.Show();
          write_table_gplt(output_prefix, "mean_gradient.cube"+output_suffix, table);

          if (flag_standard_deviation) {
            table.HideAllExceptScalar();
            table.gradient_standard_deviation_cube.Show();
            write_table_gplt(output_prefix, "gradient_deviation.cube"+output_suffix, table);
          }

        }

        if (flag_laplace) {
          table.HideAllExceptScalar();
          table.mean_gradient_laplace.Show();
          write_table_gplt(output_prefix, "mean_gradient.laplace"+output_suffix, table);
        }

      };

      if (flag_zero_gradient_table) {
        table.HideAllExceptScalar();
        table.zero_gradient_edge.Show();
        write_table_gplt(output_prefix, "zero_gradient.edge"+output_suffix, table);
      }

      if (flag_edge_dim) {
        table.HideAllExceptScalar();
        table.edge_dim.Show();
        write_table_gplt(output_prefix, "edge_dim"+output_suffix, table);
      }

      if (flag_box_dim) {
        table.HideAllExceptScalar();
        table.box_dim.Show();
        write_table_gplt(output_prefix, "box_dim"+output_suffix, table);
      }

      if (flag_num_components_table) {
        table.HideAllExceptScalar();
        table.num_components.Show();
        write_table_gplt(output_prefix, "num_components"+output_suffix, table);
      }
    }
  }

}


SCALAR_TYPE round_down_min_scalar
(const SCALAR_TYPE min_scalar, const SCALAR_TYPE max_scalar)
{
  SCALAR_TYPE min_scalar2 = min_scalar;
  SCALAR_TYPE span = max_scalar-min_scalar;

  if (span == 0) { return(min_scalar); };
  if (min_scalar == 0) { return(min_scalar); }

  if (long(min_scalar) == min_scalar) {
    // min_scalar is an integer
    return(min_scalar);
  }

  if (span >= 10) {
    min_scalar2 = SCALAR_TYPE(floor(float(min_scalar)));
  }
  else if (span >= 5) {
    min_scalar2 = SCALAR_TYPE(floor(float(2*min_scalar))/2);
  }
  else if (span >= 1) {
    min_scalar2 = SCALAR_TYPE(floor(float(10*min_scalar))/10);
  }
  return(min_scalar2);
}

template <class ITYPE>
NUM_TYPE compute_num_rows(const SGRID & scalar_grid, const ITYPE interval)
{
  IJK::PROCEDURE_ERROR error("compute_num_rows");

  if (!check_interval(interval, error)) { throw error; };

  if (scalar_grid.NumVertices() == 0) { return(0);  }
  else {

    const SCALAR_TYPE min_scalar = scalar_grid.FindMinScalar();
    const SCALAR_TYPE max_scalar = scalar_grid.FindMaxScalar();

    NUM_TYPE num_rows = 
      compute_num_buckets<NUM_TYPE>(min_scalar, max_scalar, interval);

    return(num_rows);
  };
}

template <class TABLE_TYPE, class ITYPE>
void set_scalar_column(const ITYPE interval, TABLE_TYPE & table)
{
  IJK::PROCEDURE_ERROR error("set_scalar_column");

  if (!check_interval(interval, error)) { throw error; };

  const SCALAR_TYPE min_scalar = table.MinScalar();

  table.scalar.Include();
  table.scalar.SetAtIntervals(min_scalar, interval);
}

template <class TABLE_TYPE>
void compute_scalar_frequency
(const BUCKET_GRID & bucket_grid, TABLE_TYPE & table)
{
  table.scalar_freq.Include();

  if (table.NumRows() == 0) { return; };

  table.scalar_freq.SetAll(0);

  const BUCKET_TYPE * bucket_ptr = bucket_grid.ScalarPtrConst();
  while (bucket_ptr != bucket_grid.End()) {
    NUM_TYPE irow = *bucket_ptr;
    table.scalar_freq.Increment(irow);
    bucket_ptr++;
  }

}

void measure_isoarea_using_edges
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("measure_isoarea_using_edges");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  IJK::ARRAY<VERTEX_INDEX> axis_increment( dimension );

  if (!table.scalar.IsIncluded()) {
    error.AddMessage("Programming error. Table scalar column is not set.");
    throw error;
  }

  table.isoarea_edge.Include();

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  ARRAY<NUM_TYPE> num_begin(num_rows, 0);
  ARRAY<NUM_TYPE> num_end(num_rows, 0);
  compute_increment(bucket_grid, &(axis_increment[0]));

  for (NUM_TYPE d = 0; d < dimension; d++) {

    NUM_TYPE numv_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d, numv_in_grid_facet);

    ARRAY<VERTEX_INDEX> vlist(numv_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d, false, vlist.Ptr());

    if (axis_size[d] <= 0) { return; }

    for (VERTEX_INDEX k = 0; k < numv_in_grid_facet; k++) {
      for (VERTEX_INDEX i = 0; i+1 < axis_size[d]; i++) {

        VERTEX_INDEX iv0 = vlist[k] + axis_increment[d]*i;
        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        if (s0 > s1) {
          swap(s0, s1);
          swap(irow0, irow1);
        }

        if (s0 > table.scalar[irow0]) { irow0++; };
        if (s1 <= table.scalar[irow1]) { irow1--; };

        if (irow0 <= irow1) {
          num_begin[irow0]++;
          num_end[irow1]++;
        }
      }
    }
    vlist.Free();
  }

  NUM_TYPE ncurrent = num_begin[0];
  table.isoarea_edge.Set(0, ncurrent);
  for (NUM_TYPE irow = 1; irow < num_rows; irow++) {
    ncurrent = ncurrent + num_begin[irow] - num_end[irow-1];
    table.isoarea_edge.Set(irow, ncurrent);
  }
}

void measure_ivol_using_edges
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 const INTERVAL_TYPE bucket_interval, DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("measure_ivol_using_edges");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();

  IJK::ARRAY<VERTEX_INDEX> axis_increment( dimension );
  
  table.ivol_edge.Include();
  table.ivol_edge.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(bucket_grid, (&axis_increment[0]));

  for (NUM_TYPE d = 0; d < dimension; d++) {

    NUM_TYPE numv_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d, numv_in_grid_facet);

    ARRAY<VERTEX_INDEX> vlist(numv_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d, false, vlist.Ptr());

    if (axis_size[d] <= 0) { return; }

    for (VERTEX_INDEX k = 0; k < numv_in_grid_facet; k++) {
      for (VERTEX_INDEX i = 0; i+1 < axis_size[d]; i++) {

        VERTEX_INDEX iv0 = vlist[k] + axis_increment[d]*i;
        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s1 < s0) {
          swap(irow0, irow1);
          swap(s1, s0);
        }

        // Divide by the dimension to compensate for the number
        //   of grid edges being about (dimension) x (number of grid cubes).
        WEIGHT_TYPE w = 1.0/dimension;
        if (s1-s0 > bucket_interval) 
          { w = double(bucket_interval)/(dimension*(s1-s0)); }
        for (int irow = irow0; irow <= irow1; irow++)
          { table.ivol_edge.Add(irow, w); }
      }
    }
  }
}

void compute_total_gradient_edge
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_total_gradient_edge");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  IJK::ARRAY<VERTEX_INDEX> axis_increment( dimension );

  table.total_gradient_edge.Include();
  table.total_gradient_edge.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(bucket_grid, &(axis_increment[0]));

  NUM_TYPE num_cubes_in_grid_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_grid_facet0);

  ARRAY<VERTEX_INDEX> vlist(num_cubes_in_grid_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, vlist.Ptr());

  if (axis_size[0] <= 0) { return; }

  for (VERTEX_INDEX k = 0; k < num_cubes_in_grid_facet0; k++) {
    for (VERTEX_INDEX iv0 = vlist[k]; iv0 < vlist[k]+axis_size[0]-1; iv0++) {

      for (NUM_TYPE d = 0; d < dimension; d++) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s0 > s1) {
          swap(s0, s1);
          swap(irow0, irow1);
        }
        GRADIENT_TYPE g = s1 - s0;

        if (s0 > table.scalar[irow0]) { irow0++; };
        if (s1 <= table.scalar[irow1]) { irow1--; };

        table.total_gradient_edge.Add(irow0, irow1, g);
      }
    }
  }
  vlist.Free();

  NUM_TYPE num_outer_vertices;
  compute_num_outer_vertices(dimension, axis_size, num_outer_vertices);
  ARRAY<VERTEX_INDEX> outer_vlist(num_outer_vertices);
  get_outer_grid_vertices(dimension, axis_size, outer_vlist.Ptr());
  //DEBUG
  //GRID_COORD_TYPE coord0[dimension];
  IJK::ARRAY<GRID_COORD_TYPE> coord0( dimension );
  for (VERTEX_INDEX k = 0; k < num_outer_vertices; k++) {
    VERTEX_INDEX iv0 = outer_vlist[k];

    bucket_grid.ComputeCoord(iv0, &(coord0[0]));

    NUM_TYPE irow0 = bucket_grid.Scalar(iv0);

    for (NUM_TYPE d = 0; d < dimension; d++) {
      if (coord0[d]+1 < axis_size[d]) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        GRADIENT_TYPE g = 1.0;
        if (s0 < s1) { g = s1 - s0; }
        else { g = s0 - s1; }

        table.total_gradient_edge.Add(irow0, irow1, g);
      }
    }
  }

}


void measure_isoarea_using_cubes
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("measure_isoarea_using_cubes");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  const NUM_TYPE numv = bucket_grid.NumVertices();
  IJK::ARRAY <VERTEX_INDEX> axis_increment (dimension);
  

  if (!table.scalar.IsIncluded()) {
    error.AddMessage("Programming error. Table scalar column is not set.");
    throw error;
  }

  table.isoarea_cube.Include();

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  ARRAY<NUM_TYPE> num_begin(num_rows, 0);
  ARRAY<NUM_TYPE> num_end(num_rows, 0);

  compute_increment(bucket_grid, &(axis_increment[0]));

  // *** Expensive and convoluted way to find vertex with min scalar ***
  // *** Should be replaced ***
  MINMAX_REGIONS<INFO_GRID,BUCKET_TYPE> minmax_bucket;
  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;

  minmax_bucket.ComputeMinMax(bucket_grid, 1);
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  for (NUM_TYPE icube = 0; icube < minmax_bucket.NumRegions(); icube++) {

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    NUM_TYPE irow0 = minmax_bucket.Min(icube);
    NUM_TYPE irow1 = minmax_bucket.Max(icube);

    if (s0 > table.scalar[irow0]) { irow0++; };
    if (s1 <= table.scalar[irow1]) { irow1--; };

    if (irow0 <= irow1) {
      num_begin[irow0]++;
      num_end[irow1]++;
    }
  }

  NUM_TYPE ncurrent = num_begin[0];
  table.isoarea_cube.Set(0, ncurrent);
  for (NUM_TYPE irow = 1; irow < num_rows; irow++) {
    ncurrent = ncurrent + num_begin[irow] - num_end[irow-1];
    table.isoarea_cube.Set(irow, ncurrent);
  }
}

void measure_ivol_using_cubes
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 const INTERVAL_TYPE bucket_interval, DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("measure_ivol_using_cubes");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  const NUM_TYPE numv = bucket_grid.NumVertices();
  IJK::ARRAY<VERTEX_INDEX> axis_increment ( dimension );
  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];

  table.ivol_cube.Include();

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(bucket_grid, &(axis_increment[0]));

  MINMAX_REGIONS<INFO_GRID,BUCKET_TYPE> minmax_bucket;
  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;

  minmax_bucket.ComputeMinMax(bucket_grid, 1);
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  table.ivol_cube.SetAll(0);

  for (NUM_TYPE icube = 0; icube < minmax_bucket.NumRegions(); icube++) {

    NUM_TYPE irow0 = minmax_bucket.Min(icube);
    NUM_TYPE irow1 = minmax_bucket.Max(icube);

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    WEIGHT_TYPE w = 1.0;
    if (s1-s0 > bucket_interval) 
      { w = double(bucket_interval)/(s1-s0); }    
    for (int irow = irow0; irow <= irow1; irow++)
      {	table.ivol_cube.Add(irow, w); }
  }
}

void compute_total_gradient_cube
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_total_gradient_cube");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  const NUM_TYPE numv = bucket_grid.NumVertices();

  table.total_gradient_cube.Include();
  table.total_gradient_cube.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  MINMAX_REGIONS<INFO_GRID,BUCKET_TYPE> minmax_bucket;
  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;

  minmax_bucket.ComputeMinMax(bucket_grid, 1);
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  for (NUM_TYPE icube = 0; icube < minmax_bucket.NumRegions(); icube++) {

    NUM_TYPE irow0 = minmax_bucket.Min(icube);
    NUM_TYPE irow1 = minmax_bucket.Max(icube);

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    GRADIENT_TYPE g = s1 - s0;

    if (s0 > table.scalar[irow0]) { irow0++; };
    if (s1 <= table.scalar[irow1]) { irow1--; };

    table.total_gradient_cube.Add(irow0, irow1, g);
  }
}

void compute_total_gradient
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table)
{
  if (flag_edge) {
    compute_total_gradient_edge(scalar_grid, bucket_grid, table);
  }

  if (flag_cube) {
    compute_total_gradient_cube(scalar_grid, bucket_grid, table);
  }

  if (flag_laplace) {
    compute_total_gradient_Laplace(scalar_grid, bucket_grid, table);
  }

}

void compute_total_gradient_Laplace
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_total_gradient_Laplace");
  const NUM_TYPE dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();

  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];
  IJK::ARRAY <VERTEX_INDEX> axis_increment ( dimension );
  

  table.total_gradient_laplace.Include();
  table.total_gradient_laplace.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(dimension, axis_size, &(axis_increment[0]));

  AXIS_SIZE_TYPE boundary_width = 2;
  if (dimension > 1) {

    NUM_TYPE num_vertices_in_grid_facet_interior;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, NUM_TYPE(0), boundary_width,
       num_vertices_in_grid_facet_interior);

    if (num_vertices_in_grid_facet_interior > 0) {

      IJK::ARRAY<VERTEX_INDEX> vlist(num_vertices_in_grid_facet_interior);
      get_vertices_in_grid_facet0
        (dimension, axis_size, boundary_width, vlist.Ptr());

      for (NUM_TYPE j = 0; j < num_vertices_in_grid_facet_interior; j++) {
        vlist[j] += boundary_width;

        for (AXIS_SIZE_TYPE i = boundary_width; 
             i+boundary_width < axis_size[0]; i++) {
          VERTEX_INDEX iv = vlist[j];

          GRADIENT_TYPE L = LaplaceInterior2
            (dimension, &(axis_increment[0]), iv, scalar);

          NUM_TYPE irow = bucket_grid.Scalar(iv);
          table.total_gradient_laplace.Add(irow, -L);
          vlist[j]++;
        }
      }
      vlist.Free();
    }
  }
  else {

    for (VERTEX_INDEX iv = boundary_width; 
         iv+boundary_width < axis_size[0]; iv++) {

      GRADIENT_TYPE L = LaplaceInterior2
        (dimension, &(axis_increment[0]), iv, scalar);

      NUM_TYPE irow = bucket_grid.Scalar(iv);
      table.total_gradient_laplace.Add(irow, -L);
    }
  }

  NUM_TYPE num_boundary_vertices;
  compute_num_boundary_grid_vertices
    (dimension, axis_size, boundary_width, num_boundary_vertices);

  IJK::ARRAY<VERTEX_INDEX> boundary_vlist(num_boundary_vertices);
  get_boundary_grid_vertices(dimension, axis_size, boundary_width, 
                             boundary_vlist.Ptr());
  //DEBUG
  //VERTEX_INDEX extended_increment[dimension+1];
  IJK::ARRAY<VERTEX_INDEX>  extended_increment ( dimension + 1 );
  compute_increment(dimension, axis_size, &(extended_increment[0]));
  extended_increment[dimension] = scalar_grid.NumVertices();
  
  for (NUM_TYPE j = 0; j < num_boundary_vertices; j++) {
    VERTEX_INDEX iv = boundary_vlist[j];
    GRADIENT_TYPE L = LaplaceBnd(dimension, &(extended_increment[0]), iv, scalar);
    NUM_TYPE irow = bucket_grid.Scalar(iv);

    table.total_gradient_laplace.Add(irow, -L);
  }

  if (flag_sum_from0) {
    // sum to 0
    for (NUM_TYPE k = 1; k < num_rows; k++) {
      GRADIENT_TYPE v = table.total_gradient_laplace[k-1];
      table.total_gradient_laplace.Add(k, v);
    }
  }
  else {
    // sum to 0
    for (NUM_TYPE k = 1; k < num_rows; k++) {
      GRADIENT_TYPE v = table.total_gradient_laplace[num_rows-k];
      table.total_gradient_laplace.Add(num_rows-k-1, v);
    }
  }
}


void compute_mean_gradient
(const DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> & total_gradient,
 const DATA_COLUMN<NUM_TYPE, NUM_TYPE> & iso_area,
 DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> & mean_gradient,
 const int num_rows)
{
  IJK::PROCEDURE_ERROR error("compute_mean_gradient");

  if (!iso_area.IsIncluded()) {
    error.AddMessage("Programming error. Isosurface area not computed.");
    throw error;
  }

  if (!total_gradient.IsIncluded()) {
    error.AddMessage("Programming error. Total gradient not computed.");
    throw error;
  }

  mean_gradient.Include();

  for (NUM_TYPE i = 0; i < num_rows; i++) {

    GRADIENT_TYPE v = total_gradient[i];
    NUM_TYPE n = iso_area[i];
    if (n > 0) {
      GRADIENT_TYPE x = v/n;
      mean_gradient.Set(i, x); 
    }
    else { mean_gradient.Set(i, 0); }
  };
  
}

void compute_gradient_standard_deviation_edge
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_gradient_standard_deviation_edge");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];
  IJK::ARRAY < VERTEX_INDEX > axis_increment ( dimension );
  table.gradient_standard_deviation_edge.Include();
  table.gradient_standard_deviation_edge.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(bucket_grid, &(axis_increment[0]));

  NUM_TYPE num_cubes_in_grid_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_grid_facet0);

  ARRAY<VERTEX_INDEX> vlist(num_cubes_in_grid_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, vlist.Ptr());

  if (axis_size[0] <= 0) { return; }

  for (VERTEX_INDEX k = 0; k < num_cubes_in_grid_facet0; k++) {
    for (VERTEX_INDEX iv0 = vlist[k]; iv0 < vlist[k]+axis_size[0]-1; iv0++) {

      for (NUM_TYPE d = 0; d < dimension; d++) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s0 > s1) {
          swap(s0, s1);
          swap(irow0, irow1);
        }
        GRADIENT_TYPE g = s1 - s0;

        if (s0 > table.scalar[irow0]) { irow0++; };
        if (s1 <= table.scalar[irow1]) { irow1--; };

        for (int i = irow0; i <= irow1; i++) {
          table.gradient_standard_deviation_edge.AddSquare
            (i, g-table.mean_gradient_edge[i]);
        }
      }
    }
  }
  vlist.Free();

  NUM_TYPE num_outer_vertices;
  compute_num_outer_vertices(dimension, axis_size, num_outer_vertices);
  ARRAY<VERTEX_INDEX> outer_vlist(num_outer_vertices);
  get_outer_grid_vertices(dimension, axis_size, outer_vlist.Ptr());
  //DEBUG
  //GRID_COORD_TYPE coord0[dimension];
  IJK::ARRAY<GRID_COORD_TYPE> coord0 (dimension);
  
  for (VERTEX_INDEX k = 0; k < num_outer_vertices; k++) {
    VERTEX_INDEX iv0 = outer_vlist[k];

    bucket_grid.ComputeCoord(iv0, &(coord0[0]));

    NUM_TYPE irow0 = bucket_grid.Scalar(iv0);

    for (NUM_TYPE d = 0; d < dimension; d++) {
      if (coord0[d]+1 < axis_size[d]) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        GRADIENT_TYPE g = 1.0;
        if (s0 < s1) { g = s1 - s0; }
        else { g = s0 - s1; }

        if (s0 > table.scalar[irow0]) { irow0++; };
        if (s1 <= table.scalar[irow1]) { irow1--; };

        for (int i = irow0; i <= irow1; i++) {
          table.gradient_standard_deviation_edge.AddSquare
            (i, g-table.mean_gradient_edge[i]); 
        }
      }
    }
  }

  for (int i = 0; i < num_rows; i++) {
    if (table.isoarea_edge[i] > 0) {
      GRADIENT_TYPE variance = 
        table.gradient_standard_deviation_edge[i]/table.isoarea_edge[i];
      GRADIENT_TYPE standard_dev = sqrt(variance);
      table.gradient_standard_deviation_edge.Set(i, standard_dev);
    }
    else {
      table.gradient_standard_deviation_edge.Set(i, 0);
    }
  }
}

void compute_gradient_standard_deviation_cube
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_gradient_standard_deviation_cube");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  const NUM_TYPE numv = bucket_grid.NumVertices();

  table.gradient_standard_deviation_cube.Include();
  table.gradient_standard_deviation_cube.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  MINMAX_REGIONS<INFO_GRID,BUCKET_TYPE> minmax_bucket;
  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;

  minmax_bucket.ComputeMinMax(bucket_grid, 1);
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  for (NUM_TYPE icube = 0; icube < minmax_bucket.NumRegions(); icube++) {

    NUM_TYPE irow0 = minmax_bucket.Min(icube);
    NUM_TYPE irow1 = minmax_bucket.Max(icube);

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    GRADIENT_TYPE g = s1 - s0;

    if (s0 > table.scalar[irow0]) { irow0++; };
    if (s1 <= table.scalar[irow1]) { irow1--; };

    for (int i = irow0; i < irow1; i++) {
      table.gradient_standard_deviation_cube.AddSquare
        (i, g-table.mean_gradient_cube[i]); 
    }
  }

  for (int i = 0; i < num_rows; i++) {
    if (table.isoarea_cube[i] > 0) {
      GRADIENT_TYPE variance = 
        table.gradient_standard_deviation_cube[i]/table.isoarea_cube[i];
      GRADIENT_TYPE standard_dev = sqrt(variance);
      table.gradient_standard_deviation_cube.Set(i, standard_dev);
    }
    else {
      table.gradient_standard_deviation_cube.Set(i, 0);
    }
  }

}

void compute_zero_gradient
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("compute_zero_gradient");
  const NUM_TYPE dimension = bucket_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = bucket_grid.AxisSize();
  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];
  IJK::ARRAY <VERTEX_INDEX> axis_increment ( dimension );

  table.zero_gradient_edge.Include();
  table.zero_gradient_edge.SetAll(0);

  const SCALAR_TYPE min_scalar = table.MinScalar();
  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  compute_increment(bucket_grid, &(axis_increment[0]));

  for (NUM_TYPE d = 0; d < dimension; d++) {

    NUM_TYPE numv_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d, numv_in_grid_facet);

    ARRAY<VERTEX_INDEX> vlist(numv_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d, false, vlist.Ptr());

    if (axis_size[d] <= 0) { return; }

    for (VERTEX_INDEX k = 0; k < numv_in_grid_facet; k++) {
      for (VERTEX_INDEX i = 0; i+1 < axis_size[d]; i++) {

        VERTEX_INDEX iv0 = vlist[k] + axis_increment[d]*i;
        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s0 == s1) {
          NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
          table.zero_gradient_edge.Increment(irow0);
        }
      }
    }
    vlist.Free();
  }

}

void compute_num_components
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid, 
 DATA_TABLE & table)
{
  using namespace IJK::CONTOUR_TREE;
  CONNECTIVITY join_connectivity = EDGE_CONNECT;
  CONNECTIVITY split_connectivity = CUBE_CONNECT;

  vector< CONTOUR_TREE_NODE< TREE_NODE_ID<int,TREE_NODE> > > contour_tree;

  const AUGMENT_FLAG augment_flag = AUGMENT_NONE;
  const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE;

  table.num_components.Include();
  table.num_components.SetAll(0);

  NUM_TYPE num_rows = table.NumRows();
  if (num_rows == 0) { return; };

  construct_contour_tree
    (scalar_grid, join_connectivity, split_connectivity,
     contour_tree, augment_flag, merge_ident_flag);

  ARRAY<NUM_TYPE> num_begin(num_rows, 0);
  ARRAY<NUM_TYPE> num_end(num_rows, 0);

  for (NUM_TYPE i = 0; i < contour_tree.size(); i++) {
    if (!contour_tree[i].IsDeleted() && !contour_tree[i].IsRoot()) {
      NUM_TYPE iv0 = contour_tree[i].Ident();
      NUM_TYPE iv1 = contour_tree[i].Parent()->Ident();
      SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
      SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

      NUM_TYPE irow0 = bucket_grid.Scalar(iv0);
      NUM_TYPE irow1 = bucket_grid.Scalar(iv1);

      if (s0 > s1) {
        swap(s0, s1);
        swap(irow0, irow1);
      }

      if (s0 > table.scalar[irow0]) { irow0++; };
      if (s1 <= table.scalar[irow1]) { irow1--; };

      if (irow0 <= irow1) {
        num_begin[irow0]++;
        num_end[irow1]++;
      }
    }
  }

  NUM_TYPE ncurrent = num_begin[0];
  table.num_components.Set(0, ncurrent);
  for (NUM_TYPE irow = 1; irow < num_rows; irow++) {
    ncurrent = ncurrent + num_begin[irow] - num_end[irow-1];
    table.num_components.Set(irow, ncurrent);
  }

}

/// Compute number of cubes with scalar value zero at all cube vertices.
void compute_num_zero_cubes
(const SGRID & scalar_grid, NUM_TYPE & num_zero_cubes)
{
  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;

  num_zero_cubes = 0;
  minmax_scalar.ComputeMinMax(scalar_grid, 1);
  for (NUM_TYPE icube = 0; icube < minmax_scalar.NumRegions(); icube++) {
    if (minmax_scalar.Min(icube) == 0 &&
        minmax_scalar.Max(icube) == 0) 
      { num_zero_cubes++; }
  }
}

// **************************************************
// Compute Gradient Distribution
// **************************************************

GRADIENT_TYPE compute_max_cube_span(const SGRID & scalar_grid)
{
  IJK::PROCEDURE_ERROR error("compute_max_cube_span");

  GRADIENT_TYPE max_cube_span = 0;

  if (scalar_grid.AxisSize(0) <= 0) { return(0); }

  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  for (NUM_TYPE icube = 0; icube < minmax_scalar.NumRegions(); icube++) {

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    GRADIENT_TYPE g = s1 - s0;

    if (max_cube_span < g) { max_cube_span = g; };
  }

  return(max_cube_span);
}

GRADIENT_TYPE compute_max_gradient_edge(const SGRID & scalar_grid)
{
  IJK::PROCEDURE_ERROR error("compute_max_gradient_edge");
  const NUM_TYPE dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];
  IJK::ARRAY<VERTEX_INDEX> axis_increment (dimension);

  GRADIENT_TYPE max_gradient = 0;

  compute_increment(scalar_grid, &(axis_increment[0]));

  NUM_TYPE num_cubes_in_grid_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_grid_facet0);

  ARRAY<VERTEX_INDEX> vlist(num_cubes_in_grid_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, vlist.Ptr());

  if (axis_size[0] <= 0) { return(0); }

  for (VERTEX_INDEX k = 0; k < num_cubes_in_grid_facet0; k++) {
    for (VERTEX_INDEX iv0 = vlist[k]; iv0 < vlist[k]+axis_size[0]-1; iv0++) {

      for (NUM_TYPE d = 0; d < dimension; d++) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s0 > s1) { swap(s0, s1);}
        GRADIENT_TYPE g = s1 - s0;

        if (max_gradient < g) { max_gradient = g; };
      }
    }
  }
  vlist.Free();

  NUM_TYPE num_outer_vertices;
  compute_num_outer_vertices(dimension, axis_size, num_outer_vertices);
  ARRAY<VERTEX_INDEX> outer_vlist(num_outer_vertices);
  get_outer_grid_vertices(dimension, axis_size, outer_vlist.Ptr());
  //DEBUG
  //GRID_COORD_TYPE coord0[dimension];
  IJK::ARRAY<GRID_COORD_TYPE> coord0 ( dimension );

  for (VERTEX_INDEX k = 0; k < num_outer_vertices; k++) {
    VERTEX_INDEX iv0 = outer_vlist[k];

    scalar_grid.ComputeCoord(iv0, &(coord0[0]));

    for (NUM_TYPE d = 0; d < dimension; d++) {
      if (coord0[d]+1 < axis_size[d]) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        GRADIENT_TYPE g = 1.0;
        if (s0 < s1) { g = s1 - s0; }
        else { g = s0 - s1; }

        if (max_gradient < g) { max_gradient = g; };
      }
    }
  }

  return(max_gradient);
}

template <class ITYPE>
NUM_TYPE compute_num_gradient_rows
(const SGRID & scalar_grid, const ITYPE interval, const GRADIENT_TYPE max_gradient)
{
  IJK::PROCEDURE_ERROR error("compute_num_gradient_rows");

  if (!check_interval(interval, error)) { throw error; };

  if (scalar_grid.NumVertices() == 0) { return(0);  }
  else {

    NUM_TYPE num_rows = 
      compute_num_buckets<NUM_TYPE>(0.0, max_gradient, interval);

    return(num_rows);
  };
}

template <class TABLE_TYPE, class ITYPE>
void set_gradient_column(const ITYPE interval, TABLE_TYPE & table)
{
  IJK::PROCEDURE_ERROR error("set_gradient_column");

  if (!check_interval(interval, error)) { throw error; };

  table.gradient.Include();
  table.gradient.SetAtIntervals(0, interval);
}

template <class TABLE_TYPE, class STAT_TYPE, class ITYPE>
void compute_edge_span_frequency
(const SGRID & scalar_grid, const ITYPE interval, 
 TABLE_TYPE & table, STAT_TYPE & statistics)
{
  const NUM_TYPE dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const NUM_TYPE num_rows = table.NumRows();
  //DEBUG
  //VERTEX_INDEX axis_increment[dimension];
  IJK::ARRAY <VERTEX_INDEX> axis_increment ( dimension );

  IJK::PROCEDURE_ERROR error("edge_span_frequency");

  table.edge_span_freq.Include();
  table.edge_span_freq.SetAll(0);
  statistics.average_edge_span.SetToZero();

  if (table.NumRows() == 0) { return; };

  compute_increment(scalar_grid, &(axis_increment[0]));

  NUM_TYPE num_cubes_in_grid_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_grid_facet0);

  ARRAY<VERTEX_INDEX> vlist(num_cubes_in_grid_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, vlist.Ptr());

  if (axis_size[0] <= 0) { return; }

  for (VERTEX_INDEX k = 0; k < num_cubes_in_grid_facet0; k++) {
    for (VERTEX_INDEX iv0 = vlist[k]; iv0 < vlist[k]+axis_size[0]-1; iv0++) {

      for (NUM_TYPE d = 0; d < dimension; d++) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];
        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        if (s0 > s1) { swap(s0, s1);}
        GRADIENT_TYPE g = s1 - s0;

        if (!is_min_scalar_set || min_scalar <= s1) {
          if (!is_max_scalar_set || max_scalar >= s0) {
            NUM_TYPE irow = compute_bucket(g, 0.0, interval, num_rows);
            table.edge_span_freq.Increment(irow);

            statistics.average_edge_span.Add(g);
            if (g == 0) { statistics.num_edges_with_span_zero++; }
          }
        }
      }
    }
  }
  vlist.Free();

  NUM_TYPE num_outer_vertices;
  compute_num_outer_vertices(dimension, axis_size, num_outer_vertices);
  ARRAY<VERTEX_INDEX> outer_vlist(num_outer_vertices);
  get_outer_grid_vertices(dimension, axis_size, outer_vlist.Ptr());

  //DEBUG
  //GRID_COORD_TYPE coord0[dimension];
  IJK::ARRAY <GRID_COORD_TYPE>  coord0 ( dimension );

  for (VERTEX_INDEX k = 0; k < num_outer_vertices; k++) {
    VERTEX_INDEX iv0 = outer_vlist[k];

    scalar_grid.ComputeCoord(iv0, &(coord0[0]));

    for (NUM_TYPE d = 0; d < dimension; d++) {
      if (coord0[d]+1 < axis_size[d]) {

        VERTEX_INDEX iv1 = iv0 + axis_increment[d];

        SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

        GRADIENT_TYPE g = 1.0;
        if (s0 < s1) { g = s1 - s0; }
        else { g = s0 - s1; }

        if (!is_min_scalar_set || min_scalar <= s1) {
          if (!is_max_scalar_set || max_scalar >= s0) {
            NUM_TYPE irow = compute_bucket(g, 0.0, interval, num_rows);
            table.edge_span_freq.Increment(irow);

            statistics.average_edge_span.Add(g);
          }
        }
      }
    }
  }
}

template <class TABLE_TYPE, class STAT_TYPE, class ITYPE>
void compute_cube_span_frequency
(const SGRID & scalar_grid, const ITYPE interval, 
 TABLE_TYPE & table, STAT_TYPE & statistics)
{
  const NUM_TYPE num_rows = table.NumRows();
  IJK::PROCEDURE_ERROR error("cube_span_frequency");

  table.cube_span_freq.Include();
  table.cube_span_freq.SetAll(0);
  statistics.average_cube_span.SetToZero();

  if (table.NumRows() == 0) { return; };
  if (scalar_grid.AxisSize(0) <= 0) { return; }

  MINMAX_REGIONS<INFO_GRID,SCALAR_TYPE> minmax_scalar;
  minmax_scalar.ComputeMinMax(scalar_grid, 1);

  for (NUM_TYPE icube = 0; icube < minmax_scalar.NumRegions(); icube++) {

    SCALAR_TYPE s0 = minmax_scalar.Min(icube);
    SCALAR_TYPE s1 = minmax_scalar.Max(icube);

    GRADIENT_TYPE g = s1 - s0;

    if (!is_min_scalar_set || min_scalar <= s1) {
      if (!is_max_scalar_set || max_scalar >= s0) {
        NUM_TYPE irow = compute_bucket(g, 0.0, interval, num_rows);
        table.cube_span_freq.Increment(irow);

        statistics.average_cube_span.Add(g);
        if (g == 0) { statistics.num_cubes_with_span_zero++; }
      }
    }
  }
}

void compute_gradient_distribution_tables(const SGRID & scalar_grid)
{
  const NUM_TYPE dimension = scalar_grid.Dimension();

  GRADIENT_TYPE max_cube_span = compute_max_cube_span(scalar_grid);

  NUM_TYPE num_rows = 1;
  if (flag_interval) {
    num_rows = compute_num_gradient_rows(scalar_grid, bucket_interval, max_cube_span);
  }
  else {
    num_rows = num_buckets;

    if (num_buckets <= 1) { 
      bucket_interval = max_cube_span+1;
    }
    else {
      bucket_interval = max_cube_span/(num_buckets-1);
    }

    double max_ceil = ceil(double(max_cube_span));
    if (max_cube_span + bucket_interval > max_ceil) {
      // snap max_cube_span to max_ceil;
      max_cube_span = GRADIENT_TYPE(max_cube_span);
    }

    if (num_buckets <= 1) { 
      bucket_interval = max_cube_span + 1;
    }
    else {
      bucket_interval = max_cube_span/(num_buckets-1);
    }

  }

  GRADIENT_TABLE gradient_table(num_rows);
  GRADIENT_STATISTICS gradient_statistics;

  set_gradient_column(bucket_interval, gradient_table);

  if (flag_edge_span_freq)
    { compute_edge_span_frequency
        (scalar_grid, bucket_interval, gradient_table, gradient_statistics); };

  if (flag_cube_span_freq)
    { compute_cube_span_frequency
        (scalar_grid, bucket_interval, gradient_table, gradient_statistics); };

  // Compute table sums
  gradient_table.ComputeSum();


  if (flag_edge_span_freq) {
    NUM_TYPE num_edges = scalar_grid.ComputeNumEdges();
    NUM_TYPE num_edges_with_span_zero = 
      gradient_statistics.num_edges_with_span_zero;
    output_average(cout, gradient_statistics.average_edge_span,
                     "edge span"); 
    cout << "Num edges with zero span: " << num_edges_with_span_zero;
    if (num_edges > 0) {
      cout << " ";
      print_percent(cout, num_edges_with_span_zero, num_edges);
    }
    cout << endl;
  }

  if (flag_cube_span_freq) {
    NUM_TYPE num_cubes = scalar_grid.ComputeNumCubes();
    NUM_TYPE num_cubes_with_span_zero = 
      gradient_statistics.num_cubes_with_span_zero;
    output_average(cout, gradient_statistics.average_cube_span,
                     "cube span"); 
    cout << "Num cubes with zero span: " << num_cubes_with_span_zero;
    if (num_cubes > 0) {
      cout << " ";
      print_percent(cout, num_cubes_with_span_zero, num_cubes);
    }
    cout << endl;
  }

  if (flag_plot) {

    if (output_filename != NULL) {
      write_table_gplt(output_filename, gradient_table);
    }
    else {

      string output_prefix;
      string output_suffix = ".gplt";

      // create output filename
      string fname = input_filename;

      // remove path from file name
      string prefix, suffix;
      split_string(fname, '/', prefix, suffix);
      if (suffix != "") { fname = suffix; }
      split_string(fname, '.', prefix, suffix);
      if (suffix == "nrrd" || suffix == "nhdr") {
        output_prefix = prefix;
      }
      else {
        output_prefix = input_filename;
      }

      if (file_label != NULL)
        { output_suffix = string(".") + file_label + output_suffix; }

      if (flag_edge_span_freq) {
        gradient_table.HideAllExceptGradient();
        gradient_table.edge_span_freq.Show();
        write_table_gplt(output_prefix, "edge_span_freq"+output_suffix, gradient_table);
      };

      if (flag_cube_span_freq) {
        gradient_table.HideAllExceptGradient();
        gradient_table.cube_span_freq.Show();
        write_table_gplt(output_prefix, "cube_span_freq"+output_suffix, gradient_table);
      };

    }
  }

}

// **************************************************
// Compute Fractal Edge & Box Span Dimension
// **************************************************

#ifndef log2
double log2(const double x)
{
  double y = log(x)/log(2.0);
  return(y);
}
#endif

void compute_edge_dim
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table)
{
  const NUM_TYPE dimension = bucket_grid.Dimension();
  IJK::PROCEDURE_ERROR error("compute_box_dim");

  if (!table.isoarea_edge.IsIncluded()) {
    error.AddMessage("Programming error. Edge based isoarea not computed.");
    throw error;
  }
  
  table.edge_dim.Include();

  SGRID subsampled_grid;
  subsampled_grid.Subsample(scalar_grid, 2);

  BUCKET_GRID subsampled_bucket_grid
    (dimension, subsampled_grid.AxisSize());

  compute_vertex_buckets
    (subsampled_grid, bucket_interval, table.MinScalar(),
     table.NumRows(), subsampled_bucket_grid);

  DATA_TABLE table2(table.NumRows(), table.MinScalar());
  set_scalar_column(bucket_interval, table2);

  measure_isoarea_using_edges
    (subsampled_grid, subsampled_bucket_grid, table2);

  for (NUM_TYPE irow = 0; irow < table.NumRows(); irow++) {
    NUM_TYPE m1 = table.isoarea_edge[irow];
    NUM_TYPE m2 = table2.isoarea_edge[irow];
    FRACTAL_DIM_TYPE edge_dim;
    if (m1 < 1) 
      { edge_dim = 0; }
    else if (m2 < 1) 
      { edge_dim = dimension; }
    else {
      edge_dim = log2(double(m1)) - log2(double(m2));
    }
    if (edge_dim > dimension) { edge_dim = dimension; };
    if (edge_dim < 0) { edge_dim = 0; };
    table.edge_dim.Set(irow, edge_dim);
  }
  
}

void compute_box_dim
(const SGRID & scalar_grid, const BUCKET_GRID & bucket_grid,
 DATA_TABLE & table)
{
  const NUM_TYPE dimension = bucket_grid.Dimension();
  IJK::PROCEDURE_ERROR error("compute_box_dim");

  if (!table.isoarea_cube.IsIncluded()) {
    error.AddMessage("Programming error. Cube based isoarea not computed.");
    throw error;
  }
  
  table.box_dim.Include();

  SGRID subsampled_grid;
  subsampled_grid.Subsample(scalar_grid, 2);

  BUCKET_GRID subsampled_bucket_grid
    (dimension, subsampled_grid.AxisSize());

  compute_vertex_buckets
    (subsampled_grid, bucket_interval, table.MinScalar(),
     table.NumRows(), subsampled_bucket_grid);

  DATA_TABLE table2(table.NumRows(), table.MinScalar());
  set_scalar_column(bucket_interval, table2);

  measure_isoarea_using_cubes
    (subsampled_grid, subsampled_bucket_grid, table2);

  for (NUM_TYPE irow = 0; irow < table.NumRows(); irow++) {
    NUM_TYPE m1 = table.isoarea_cube[irow];
    NUM_TYPE m2 = table2.isoarea_cube[irow];
    FRACTAL_DIM_TYPE box_dim;
    if (m1 < 1) 
      { box_dim = 0; }
    else if (m2 < 1) 
      { box_dim = dimension; }
    else {
      box_dim = log2(double(m1)) - log2(double(m2));
    }
    if (box_dim > dimension) { box_dim = dimension; };
    if (box_dim < 0) { box_dim = 0; };
    table.box_dim.Set(irow, box_dim);
  }
  
}

// **************************************************
// Compute Other Scalar Information
// **************************************************

/// Compute min difference between consecutive scalar values
void compute_min_scalar_interval
(const SGRID & scalar_grid, SCALAR_TYPE & min_interval, bool & flag_uniform)
{
  const NUM_TYPE num_vertices = scalar_grid.NumVertices();
  const NUM_TYPE num_zero = scalar_grid.CountScalar(0);
  const NUM_TYPE num_nonzero = num_vertices - num_zero;

  flag_uniform = true;
  min_interval = 0.0;

  if (num_nonzero == 0) // all vertices have scalar value zero
    { return; };

  // allocate space for all non-zero scalar values and one zero
  ARRAY<SCALAR_TYPE> scalar_list(num_nonzero+1);

  VERTEX_INDEX length = 0;
  if (num_zero > 0) {
    scalar_list[0] = 0;
    length++;
  }

  for (const SCALAR_TYPE * scalar_ptr = scalar_grid.ScalarPtrConst();
       scalar_ptr != scalar_grid.End(); scalar_ptr++) {
    if (*scalar_ptr != 0) {
      scalar_list[length] = *scalar_ptr;
      ++length;
    }
  }

  sort(scalar_list.Ptr(), scalar_list.Ptr()+length);
  SCALAR_TYPE s = scalar_list[0];
  for (int i = 0; i < length; i++) {
    SCALAR_TYPE s2 = scalar_list[i];
    if (s != s2) {
      if (flag_uniform == true) {
        min_interval = s2-s;
        flag_uniform = false;
      }
      else {
        SCALAR_TYPE diff = s2-s;
        if (min_interval > s2-s) { min_interval = diff; };
      }
      s = s2;
    }
  }

}


// **************************************************
// Output Information
// **************************************************

void output_grid_info(ostream & output, const SGRID & scalar_grid)
{
  if (flag_subsample) {
    output << "Subsampled grid." << endl;
    output << "Subsample resolution: " << subsample_resolution << endl;
    output << endl;
  }
  output << "Volume dimension: " << scalar_grid.Dimension() << endl;
  if (!flag_subsample) {
    output << "# vertices per grid axis:";
  }
  else {
    output << "# vertices per subsampled grid axis:";
  }
  for (NUM_TYPE d = 0; d < scalar_grid.Dimension(); d++) {
    output << "  " << scalar_grid.AxisSize(d);
  }
  output << endl;
  NUM_TYPE num_vertices = scalar_grid.NumVertices();
  NUM_TYPE num_edges = scalar_grid.ComputeNumEdges();
  NUM_TYPE num_cubes = scalar_grid.ComputeNumCubes();
  if (!flag_subsample) {
    output << "# grid vertices:  " << num_vertices << endl;
    output << "# grid edges:  " << num_edges << endl;
    output << "# grid cubes:  " << num_cubes << endl;
  }
  else {
    output << "# subsampled grid vertices:  " << num_vertices << endl;
    output << "# subsampled grid cubes:  " << num_cubes << endl;
  }
  if (!flag_subsample) { output << "Min/max scalar values:  "; }
  else { output << "Min/max subsampled scalar values:  "; }
  output << scalar_grid.FindMinScalar() << "  "
         << scalar_grid.FindMaxScalar() << endl;
  NUM_TYPE num_zero = scalar_grid.CountScalar(0);
  NUM_TYPE num_zero_cubes;
  compute_num_zero_cubes(scalar_grid, num_zero_cubes);
  if (!flag_subsample) {
    output << "# vertices with scalar value 0:  " << num_zero;
    if (num_vertices > 0) {
      output << " ";
      print_percent(output, num_zero, num_vertices);
    }
    output << endl;
    output << "# cubes with all vertex scalar values 0:  " << num_zero_cubes;
    if (num_cubes > 0) {
      output << " ";
      print_percent(output, num_zero_cubes, num_cubes);
    }
    output << endl;
  }
  else {
    output << "# subsampled vertices with scalar value 0:  " 
           << num_zero << endl;
    output << "# subsampled cubes with all vertex scalar values 0:  " 
           << num_zero_cubes << endl;
  }

  SCALAR_TYPE min_scalar_interval;
  bool flag_uniform;
  compute_min_scalar_interval(scalar_grid, min_scalar_interval, flag_uniform);
  if (flag_uniform) {
    if (num_vertices > 0) {
      if (!flag_subsample) {
        output << "All vertices have scalar value:  " 
               << scalar_grid.Scalar(0) << endl;
      }
      else {
        output << "All subsampled vertices have scalar value:  " 
               << scalar_grid.Scalar(0) << endl;
      }
    }
  }
  else {
    if (!flag_subsample) {
      output << "Min interval between distinct scalar values:  "
             << min_scalar_interval << endl;
    }
    else {
      output << "Min interval between distinct subsampled scalar values:  "
             << min_scalar_interval << endl;
    }
  }

  output << endl;
}

/// Output vertex information for vertex vertex_index.
void output_vertex_info(ostream & output, const SGRID & scalar_grid,
                        const int vertex_index)
{
  const int dimension = scalar_grid.Dimension();
  //DEBUG
  //GRID_COORD_TYPE coord[dimension];
  IJK::ARRAY < GRID_COORD_TYPE > coord ( dimension );

  if (vertex_index < scalar_grid.NumVertices()) {
    output << "Vertex: " << vertex_index << endl;
    output << "Scalar value: " << scalar_grid.Scalar(vertex_index)
           << endl;
    scalar_grid.ComputeCoord(vertex_index, &(coord[0]));
    output << "Coordinates:";
    for (int d = 0; d < dimension; d++) 
      { output << " " << coord[d]; }
    output << endl;
  }
  else {
    output << vertex_index << " is not a valid vertex index." << endl;
  }
}

/// Output cube information.
void output_cube_info(ostream & output, const SGRID & scalar_grid,
                      const std::vector<GRID_COORD_TYPE> & cube_coord)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  IJK::PROCEDURE_ERROR error("output_cube_info");

  if (cube_coord.size() != dimension) {
    error.AddMessage("Number of cube coordinates does not equal dimension.");
    throw error;
  }

  output << "Cube coordinates: ";
  for (int d = 0; d < dimension; d++) 
    { output << " " << cube_coord[d]; }
  output << endl;

  for (int d = 0; d < dimension; d++) {
    if (cube_coord[d] < 0 || cube_coord[d] >= scalar_grid.AxisSize(d)) {
      error.AddMessage("Illegal cube coordinate.");
      error.AddMessage("  cube_coord[", d, "] = ", cube_coord[d], ".");

      if (cube_coord[d] < 0) {
        error.AddMessage("Cube coordinate must be positive.");
      }
      else {
        error.AddMessage("Cube coordinate must less than ",
                         scalar_grid.AxisSize(d), ".");
      }
      throw error;
    }
  }

  const int numv = compute_num_cube_vertices(dimension);
  IJK::ARRAY<SCALAR_TYPE> cube_scalar(numv);

  VERTEX_INDEX iv0 = scalar_grid.ComputeVertexIndex(&(cube_coord[0]));

  output << "Cube vertices:";
  for (int k = 0; k < numv; k++)
    { output << " " << scalar_grid.CubeVertex(iv0, k); }
  output << endl;

  output << "Corner scalar values:";
  for (int k = 0; k < numv; k++) {
    VERTEX_INDEX iv = scalar_grid.CubeVertex(iv0, k);
    cube_scalar[k] = scalar_grid.Scalar(iv);
    output << " " << cube_scalar[k];
  }
  output << endl;

  if (dimension == 3) {
    IJK::ARRAY<int> parity_sign(numv);
    for (VERTEX_INDEX i = 0; i < numv; i++) {
      int p = compute_parity(dimension, i);
      parity_sign[i] = convert_parity_to_sign(p);
    }

    COORD_TYPE a;
    COORD_TYPE trilinear_midpoint[3];
    COORD_TYPE saddle_offset[3];
    COORD_TYPE saddle[3];
    int num_saddles;

    compute_saddle_3D
      (cube_scalar.PtrConst(), parity_sign.PtrConst(), a, 
       trilinear_midpoint, saddle_offset, num_saddles, 0.0001);

    output << "Num saddle points = " << num_saddles << endl;
    if (a != 0) {
      add_coord(dimension, &(cube_coord[0]), trilinear_midpoint,
                trilinear_midpoint);
      if (num_saddles == 2) {
        add_coord(dimension, trilinear_midpoint, saddle_offset, saddle);
        output << "Saddles = ("
               << saddle[0] << " " << saddle[1] << " " << saddle[2]
               << ") (";
        subtract_coord(dimension, trilinear_midpoint, saddle_offset, saddle);
        output << saddle[0] << " " << saddle[1] << " " << saddle[2]
               << ")";
        output << endl;
      }

      output << "trilinear midpoint = "
             << trilinear_midpoint[0] << " " << trilinear_midpoint[1]
             << " " << trilinear_midpoint[2] << endl;
    }
    else {
      output << "Degenerate case. a = " << a << endl;
    }
  }
}

/// Return true if irow is a local maximum.
template <class RTYPE, class NTYPE, class DTYPE,  class STYPE>
bool is_local_max
(const RTYPE irow, const DATA_COLUMN<NTYPE,DTYPE> & measure, 
 const STYPE spread)
{
  const NTYPE num_rows = measure.NumRows();

  RTYPE istart = 0;
  if (irow > spread) { istart = irow-spread; };
  for (RTYPE jrow = istart; jrow < irow; jrow++) {
    if (measure[irow] <= measure[jrow]) { return(false); }
  }

  for (RTYPE jrow = irow+1; jrow < num_rows; jrow++) {
    if (measure[irow] < measure[jrow]) { return(false); };
    RTYPE k = jrow-irow;
    if (measure[irow] > measure[jrow] && k >= spread) 
      { break; }
  }

  return(true);
}

/// Get list of local maxima.
/// Get only local maxima over range of values.
template <class NTYPE0, class DTYPE0, class NTYPE1, class DTYPE1, class STYPE>
void get_local_max
(const DATA_COLUMN<NTYPE0,DTYPE0> & scalar,
 const DATA_COLUMN<NTYPE1,DTYPE1> & measure,
 const STYPE spread, std::vector<DTYPE0> & list)
{
  const NTYPE0 num_rows = measure.NumRows();

  list.clear();

  for (NTYPE0 irow = 0; irow < num_rows; irow++) {
    if (is_local_max(irow, measure, spread)) {
      list.push_back(scalar[irow]);
    }
  }
}

/// Output x_i if y_i is a local maximum
/// Output only if y_i is local maxima over range of values.
template <class NTYPE_X, class XTYPE, class NTYPE_Y, class YTYPE,
          class STYPE>
void output_local_max
(ostream & out,
 const DATA_COLUMN<NTYPE_X,XTYPE> & x,
 const DATA_COLUMN<NTYPE_Y,YTYPE> & y,
 const STYPE max_output_size)
{
  IJK::PROCEDURE_ERROR error("output_local_max");
  const NTYPE_X num_rows = x.NumRows();

  if (!check_data_columns(x, y, error)) { throw error; }

  if (num_rows < 2) { return; };
  if (max_output_size < 1) { return; };

  vector<XTYPE> list;
  get_local_max(x, y, 1, list);

  NTYPE_X spread = 1;
  while (list.size() > max_output_size && spread <= list.size()) {
    if (spread < 5) { spread++; }
    else { spread = spread*2; };

    get_local_max(x, y, spread, list);
  }
  out << y.Label() << " local max";
  if (spread > 1) { out << " (partial list)"; }
  out << ":";

  for (typename vector<XTYPE>::const_iterator iter = list.begin(); 
       iter != list.end(); iter++) {
    out << "  " << *iter;
  }

  out << endl;
}

template <class NTYPE_X, class XTYPE, class NTYPE_Y, class YTYPE>
void output_local_max
(ostream & out, 
 const DATA_COLUMN<NTYPE_X,XTYPE> & x,
 const DATA_COLUMN<NTYPE_Y,YTYPE> & y)
{
  output_local_max(out, x, y, max_output_local_max_size);
}

/// Output some of the (x_i,y_i) values
template <class NTYPE_X, class XTYPE, class NTYPE_Y, class YTYPE,
          class NTYPE>
void output_data
(ostream & out,
 const DATA_COLUMN<NTYPE_X,XTYPE> & x,
 const DATA_COLUMN<NTYPE_Y,YTYPE> & y,
 const NTYPE istart, const NTYPE interval, const NTYPE precision)
{
  PROCEDURE_ERROR error("output_data");
  const NTYPE num_rows = x.NumRows();

  if (interval < 1) {
    error.AddMessage("Interval must be at least 1.");
    throw error;
  }

  if (istart < 0) {
    error.AddMessage("Start row must be at non-negative.");
    throw error;
  }

  out << y.Label() << ":";
  for (NTYPE irow = istart; irow < num_rows; irow += interval) {
    out << "  " << fixed << setprecision(0) << x[irow];
    out << " (" << fixed << setprecision(precision) << y[irow] << ")"; 
  }
  out << endl;
}

/// Output average
template <class NTYPE, class DTYPE>
void output_average
(ostream & out, const AVERAGE_STAT<NTYPE,DTYPE> average, const char * label)
{
  out << "Average " << label << ": " << setw(12) << left << average.Average()
      << "  (stdev: " << average.StandardDeviation() << ")"
      << endl;
}

// Choose the interval for output_data
template <class NTYPE>
NTYPE choose_interval(const NTYPE num_rows)
{
  NTYPE x = num_rows / 5;
  if (x < 1) { x = 1; };

  if (x <= 5) { return(x); }
  else if (x < 8) { return(5); }
  else if (x < 15) { return(10); }
  else if (x < 25) { return(20); }
  else if (x < 50) { return(25); }
  else if (x < 100) { return(50); }
  else if (x < 200) { return(100); }
  else if (x < 250) { return(200); }
  else if (x < 500) { return(250); }
  else if (x < 1000) { return(500); }
  else if (x < 2000) { return(1000); }
  else if (x < 2500) { return(2000); }
  else if (x < 5000) { return(2500); }
  else if (x < 10000) { return(5000); };
  
  return(10000);
}

// **************************************************
// Write Tables
// **************************************************

template <class TABLE_TYPE>
void write_table_gplt(const string & filename, const TABLE_TYPE & table)
{
  ofstream ofile(filename.c_str(), ios::out);
  if (!flag_silent_write) {
    cout << "Writing table: " << filename << endl;
  }
  write_table_gplt(ofile, table);
  ofile.close();
}

template <class TABLE_TYPE>
void write_table_gplt(const string & filename_prefix, 
                      const string & filename_suffix,
                      const TABLE_TYPE & table)
{
  string filename = filename_prefix + "." + filename_suffix;
  write_table_gplt(filename, table);
}

template <class TABLE_TYPE>
void write_table_gplt(ofstream & ofile, const TABLE_TYPE & table)
{

  ofile.precision(DEFAULT_TABLE_PRECISION);
  ofile.setf(ios::left);

  if (flag_normalize) {
    ofile << "# normalized values" << endl;
  }
  ofile << "#";
  table.WriteColumnLabels(ofile, "  ");
  ofile << endl;

  int width = DEFAULT_TABLE_COLUMN_WIDTH;
  if (flag_normalize) {
    double normalization_factor = 1.0/bucket_interval;
    table.WriteNormalizedColumnData(ofile, "  ", width, normalization_factor);
  }
  else {
    table.WriteColumnData(ofile, "  ", width);
  }
}

// **************************************************
// Check Routines
// **************************************************

template <class ITYPE>
bool check_interval(const ITYPE interval, ERROR & error)
{
  if (interval <= 0) {
    error.AddMessage("Programming error.  Interval must be positive.");
    return(false);
  }
  else { return(true); };
}

template <class NUM_TYPE>
bool check_spread(const NUM_TYPE spread, ERROR & error)
{
  if (spread <= 0) {
    error.AddMessage("Programming error.  Spread must be positive.");
    return(false);
  }
  else { return(true); };
}

// **************************************************
// Misc Routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void flag_all()
// set all measurement flags to true
{
  flag_scalar_freq_table = true;
  flag_edge = true;
  flag_cube = true;
  flag_laplace = true;
  flag_isoarea_table = true;
  flag_ivol_table = true;
  flag_total_gradient_table = true;
  flag_mean_gradient_table = true;
  flag_zero_gradient_table = true;
  flag_num_components_table = true;
  flag_edge_dim = true;
  flag_box_dim = true;
}

void set_table_flags_defaults()
// set missing table flags values to default values
{
  if (!some_table_flag_is_true()) {
    flag_scalar_freq_table = true;
    flag_isoarea_table = true;
    flag_ivol_table = true;
    flag_total_gradient_table = true;
    flag_mean_gradient_table = true;
    flag_zero_gradient_table = false;
    flag_num_components_table = false;
    flag_edge_dim = false;
    flag_box_dim = false;
  }

  if (!flag_edge && !flag_cube) {
    if (!flag_laplace) { flag_edge = true; }
    else {
      flag_isoarea_edge_table = true;
      flag_ivol_edge_table = true;
    }
  }
}

bool some_table_flag_is_true()
{
  if (flag_scalar_freq_table || flag_isoarea_table || flag_ivol_table ||
      flag_total_gradient_table || flag_mean_gradient_table || 
      flag_zero_gradient_table || flag_num_components_table ||
      flag_box_dim || flag_edge_dim || 
      flag_edge_span_freq || flag_cube_span_freq)
    { return(true); }
  else
    { return(false); }
}

void get_coord(const char * s, std::vector<GRID_COORD_TYPE> & coord)
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
    GRID_COORD_TYPE c;
    coord_string >> c;
    coord.push_back(c);
  }

  if (coord_string.fail() && !coord_string.eof()) {
    cerr << "Error reading cube coordinates: "
         << "\"" << s << "\"" << endl;
    cerr << "  Non-numeric character in coordinate string." << endl;
    exit(600);
  }
}

PARAMETER get_parameter_token(char * s)
// convert string s into parameter token
{
  for (int i = 0; i < int(UNKNOWN_PARAM); i++)
    if (strcmp(parameter_string[i], s) == 0)
      return(PARAMETER(i));
  return(UNKNOWN_PARAM);
}

void parse_command_line(int argc, char **argv)
{
  char * function_name = NULL;

  int iarg = 1;
  bool flag_none = false;
  while (iarg < argc && argv[iarg][0] == '-') {
    PARAMETER param = get_parameter_token(argv[iarg]);
    if (param == UNKNOWN_PARAM) break;

    switch(param) {

    case FREQ_PARAM:
      flag_scalar_freq_table = true;
      break;

    case EDGE_PARAM:
      flag_edge = true;
      break;

    case CUBE_PARAM:
      flag_cube = true;
      break;

    case GRADIENT_PARAM:
      flag_total_gradient_table = true;
      flag_mean_gradient_table = true;
      break;

    case STANDARD_DEVIATION_PARAM:
      flag_standard_deviation = true;
      break;

    case EDGE_SPAN_PARAM:
      flag_edge_span_freq = true;
      break;

    case CUBE_SPAN_PARAM:
      flag_cube_span_freq = true;
      break;

    case ZERO_GRADIENT_PARAM:
      flag_zero_gradient_table = true;
      break;

    case IVOL_PARAM:
      flag_ivol_table = true;
      break;

    case ISO_PARAM:
      flag_isoarea_table = true;
      break;

    case LAPLACE_PARAM:
      flag_laplace = true;
      break;

    case BOX_DIM_PARAM:
      flag_box_dim = true;
      break;

    case EDGE_DIM_PARAM:
      flag_edge_dim = true;
      break;

    case COMPONENTS_PARAM:
      flag_num_components_table = true;
      break;

    case VERTEX_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &vertex_index);
      flag_vertex_info = true;
      flag_vertex_index = true;
      break;

    case CC_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      get_coord(argv[iarg], cube_coord);
      flag_cube_info = true;
      flag_vertex_index = false;
      break;

    case VC_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      get_coord(argv[iarg], vertex_coord);
      flag_vertex_info = true;
      break;

    case INTERVAL_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &bucket_interval);
      flag_interval = true;
      break;

    case NUM_BUCKETS_PARAM:
    case NUMB_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &num_buckets);
      flag_num_buckets = true;
      break;

    case ALL_PARAM:
      flag_all();
      break;

    case NONE_PARAM:
      flag_none = true;
      break;

    case PLOT_PARAM:
      flag_plot = true;
      break;

    case NORMALIZE_PARAM:
      flag_normalize = true;
      break;

    case MIN_SCALAR_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &min_scalar);
      is_min_scalar_set = true;
      break;

    case MAX_SCALAR_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &max_scalar);
      is_max_scalar_set = true;
      break;

    case SUM_FROM0:
      flag_sum_from0 = true;
      break;

    case SUM_TO0:
      flag_sum_from0 = false;
      break;

    case SUBSAMPLE_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &subsample_resolution);
      flag_subsample = true;
      break;

    case OUTPUT_FILENAME_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      output_filename = argv[iarg];
      break;

    case FILE_LABEL_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      file_label = argv[iarg];
      break;

    case SILENT_WRITE_PARAM:
      flag_silent_write = true;
      break;

    case HELP_PARAM:
      help();
      break;

    default:
      cerr << "Error. Unknown option: " << argv[iarg] << endl;
      cerr << endl;
      usage_error();
      break;
    }

    iarg++;
  };

  if (argc != iarg+1) { usage_error(); };

  input_filename = argv[iarg];

  if (output_filename != NULL) { flag_plot = true; };

  if (output_filename != NULL && file_label != NULL) {
    cerr << endl;
    cerr << "*** Warning. File label is ignored when output file name is specified."
         << endl << endl;;
  }

  if (flag_num_buckets && flag_interval) {
    cerr << "Error.  Cannot use \"-interval\" with \"-num_buckets\" option."
         << endl;
    usage_error();
  }

  if (flag_vertex_info && flag_subsample) {
    cerr << "Error. Cannot use \"-vertex\" with \"-subsample\" option."
         << endl;
    usage_error();
  }

  if (flag_cube_info && flag_subsample) {
    cerr << "Error. Cannot use \"-cc\" with \"-subsample\" option."
         << endl;
    usage_error();
  }

  if (some_table_flag_is_true()) {
    if (flag_none) {
      cerr << "Error.  Cannot use \"-none\" parameter with \"-freq\", \"-edge\", \"-cube\"," << endl;
      cerr << "  \"-gradient\" or \"-all\" parameters." << endl;
      cerr << endl;
      usage_error();
    }

    set_table_flags_defaults();
  }
  else if (!flag_vertex_info && !flag_cube_info && !flag_none) {
    set_table_flags_defaults();
  }

}

void split_string(const string & s, const char c,
                  string & prefix, string & suffix)
// split string at last occurrence of character c into prefix and suffix
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

void usage_error()
{
  cerr << "Usage: ijkscalarinfo [OPTIONS] <input filename>"
       << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "-plot | -freq | -iso | -ivol | -gradient | -zero_gradient | -components"
       << endl;
  cerr << "-edge | -cube" << endl;
  cerr << "-edge_dim | -box_dim" << endl;
  cerr << "-edge_span | -cube_span" << endl;
  cerr << "-min_scalar <min> | -max_scalar <max>" << endl;
  cerr << "-all | -none | -normalize" << endl;
  cerr << "-laplace | -sum_from0 | -sum_to0" << endl;
  cerr << "-interval <I> | [-num_buckets|-numb] <n>" << endl;
  cerr << "-vertex <vertex_index>" << endl;
  cerr << "-cc \"cube coordinates\"" << endl;
  cerr << "-subsample <subsample resolution>" << endl;
  cerr << "-file_label <label> | -o <output filename>" << endl;
  cerr << "-help" << endl;
  exit(10);
}

void help()
{
  cerr << "Output information about a regular scalar grid." << endl;
  cerr << "  Input file format is Nearly Raw Raster Data (nrrd)." << endl;

  cerr << "Usage: ijkscalarinfo [OPTIONS] <input filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "-plot:      Create gnuplot (.gplt) files of each output measurement." 
       << endl;
  cerr << "-freq:      Output frequency of grid vertex scalar values. (Default.)" << endl;
  cerr << "-iso:       Output isosurface area. (Default.)" << endl;
  cerr << "-ivol:      Output volume of interval volume. (Default.)" << endl;
  cerr << "-gradient:  Output total and mean gradients. (Default.)" << endl;
  cerr << "-zero_gradient:  Output number of zero gradient edges." << endl;
  cerr << "-components: Output number of connected components." << endl;
  cerr << "-edge:      Output area/volume based on grid edges. (Default.)" 
       << endl;
  cerr << "-cube:      Output area/volume based on grid cubes." << endl;
  cerr << "-edge_dim:  Output fractal edge counting dimension." << endl
       << "            The fractal edge counting dimension is log2(m1)-log2(m2)" << endl
       << "            where m1 is the number of intersected edges in the original grid" << endl
       << "            and m2 is the number of intersected edges in a subsampled grid" << endl
       << "            with subsample resolution 2." << endl;
  cerr << "-box_dim:   Output fractal box counting dimension." << endl
       << "            The fractal box counting dimension is log2(m1)-log2(m2)" << endl
       << "            where m1 is the number of intersected cubes in the original grid" << endl
       << "            and m2 is the number of intersected cubes in a subsampled grid" << endl
       << "            with subsample resolution 2." << endl;
  cerr << "-edge_span: Output frequency plot of edge spans." << endl;
  cerr << "-cube_span: Output frequency plot of cube spans." << endl
       << "            If min_scalar is set, only edges/cubes with scalar value" << endl
       << "              at least the minimum are included in the edge/cube span frequencies." << endl
       << "            If max_scalar is set, only edges/cubes with scalar value" << endl
       << "              at least the maximum are included in the edge/cube span frequencies." << endl;
  cerr << "-min_scalar <min>: Minimum scalar values." << endl;
  cerr << "-max_scalar <max>: Maximum scalar values." << endl;
  cerr << "-laplace:   Compute gradient on isosurface by summing Laplacian." 
       << endl;
  cerr << "-sum_from0: Sum Laplacian from minimum value." << endl;
  cerr << "-sum_to0:   Sum Laplacian from maximum value." << endl;
  cerr << "-all:       Output all measurements." << endl;
  cerr << "-none:      Output no table measurements." << endl;
  cerr << "-normalize: Normalize frequencies by dividing by total. Normalized values are between 0.0 and 1.0." << endl;
  cerr << "-interval <I>:  Bucket interval. Must be positive floating point number." << endl;
  cerr << "-num_buckets <n>:  Number of buckets. Must be a positive integer." << endl;
  cerr << "-numb <n>:  Same as -num_buckets <n>." << endl;
  cerr << "-vertex <vertex_index>:  Output information about vertex <vertex_index>."
       << endl;
  cerr << "-cc \"<cube coordinates>\":  Output information about cube with coordinates <cube coordinates>." << endl;
  cerr << "                             <cube coordinates> must be enclosed between double apostrophes." << endl;
  cerr << "-subsample <n>:  Compute information for subsampled grid with subsample resolution <n>." << endl;
  cerr << "                 Cannot be used with \"-vertex\" or \"-cc\" options."
       << endl;
  cerr << "-o <output filename>:  Output measurements to file <output filename>." << endl;
  cerr << "   Note: Without the -o option, each measurement is placed in a separate file." << endl;
  cerr << "-file_label <label>:  Add <label> to all file names." << endl;
  cerr << "-help:      Print this help message." << endl;
  exit(10);
}

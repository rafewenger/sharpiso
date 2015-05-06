/// \file ijkscalarinfo.h
/// get info from scalar grid data
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008,2012 Rephael Wenger

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

#ifndef _IJKSCALARINFO_
#define _IJKSCALARINFO_

#include <iostream>
#include <string>

#include "ijk.txx"
#include "ijkdatatable.txx"
#include "ijkgrid.txx"

/// Types and data structures for ijkscalarinfo
namespace IJKSCALARINFO {

  using std::string;
  using std::ostream;
  using std::endl;

  typedef float SCALAR_TYPE;
  typedef int AXIS_SIZE_TYPE;
  typedef int VERTEX_INDEX;
  typedef int GRID_COORD_TYPE;
  typedef float COORD_TYPE;
  typedef int BUCKET_TYPE;
  typedef float INTERVAL_TYPE;
  typedef double GRADIENT_TYPE;
  typedef double WEIGHT_TYPE;
  typedef double FRACTAL_DIM_TYPE;
  typedef IJK::GRID_SIZE_TYPE NUM_TYPE;
  typedef IJK::GRID_PLUS<NUM_TYPE, AXIS_SIZE_TYPE, VERTEX_INDEX, NUM_TYPE> 
    INFO_GRID;
  typedef IJK::SCALAR_GRID<INFO_GRID, SCALAR_TYPE> SGRID;
  typedef IJK::SCALAR_GRID<INFO_GRID, BUCKET_TYPE> BUCKET_GRID;

  typedef enum { EDGE_GRADIENT, CUBE_GRADIENT, LAPLACE_GRADIENT } 
    GRADIENT_METHOD;

  class DATA_TABLE;
  class GRADIENT_TABLE;

  template <class NTYPE, class DTYPE>
    class DATA_COLUMN: public IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE> {

    protected:
      bool is_included;
      bool is_hidden;
      DTYPE sum;

      void Init() { is_included = false; is_hidden = false; sum = 0; };

    public:
    DATA_COLUMN<NTYPE,DTYPE>() :
      IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>() { Init(); };
    DATA_COLUMN<NTYPE,DTYPE>(const string & label) :
      IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>(label) { Init(); };

    // set functions
    void Include();
    void Hide() { is_hidden = true; };
    void Show() { is_hidden = false; };
    void ComputeSum();

    // get functions
    bool IsIncluded() const { return(is_included); };
    bool IsHidden() const { return(is_hidden); };
    DTYPE Sum() const { return(sum); };

    // write functions
    void WriteLabel(std::ostream & out, const string & separator) const;
    void WriteData(std::ostream & out, const string & separator, 
		   const NTYPE width, const NTYPE irow) const;
    void WriteNormalizedData
      (std::ostream & out, const string & separator, const NTYPE width,
       const double normalization_factor, const NTYPE irow) const;

    friend class DATA_TABLE;
    friend class GRADIENT_TABLE;
  };

  class DATA_TABLE:public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  protected:
    SCALAR_TYPE min_scalar;
    void Init(const NUM_TYPE num_rows, const SCALAR_TYPE min_scalar);

  public:
    DATA_COLUMN<NUM_TYPE, SCALAR_TYPE> scalar;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> scalar_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> isoarea_edge;
    DATA_COLUMN<NUM_TYPE, WEIGHT_TYPE> ivol_edge;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> isoarea_cube;
    DATA_COLUMN<NUM_TYPE, WEIGHT_TYPE> ivol_cube;
    DATA_COLUMN<NUM_TYPE, FRACTAL_DIM_TYPE> edge_dim;
    DATA_COLUMN<NUM_TYPE, FRACTAL_DIM_TYPE> box_dim;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> total_gradient_edge;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> total_gradient_cube;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> total_gradient_laplace;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> mean_gradient_edge;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> mean_gradient_cube;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> mean_gradient_laplace;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> gradient_standard_deviation_edge;
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> gradient_standard_deviation_cube;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> zero_gradient_edge;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> num_components;
    
    DATA_TABLE(const NUM_TYPE num_rows, const SCALAR_TYPE min_scalar):
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(num_rows),
      scalar("scalar"), scalar_freq("scalar-freq"),
      isoarea_edge("isoarea-edge"), ivol_edge("ivol-edge"),
      isoarea_cube("isoarea-cube"), ivol_cube("ivol-cube"),
      edge_dim("edge-dim"), box_dim("box-dim"),
      total_gradient_edge("total_gradient-edge"),
      total_gradient_cube("total_gradient-cube"),
      total_gradient_laplace("total_gradient-laplace"),
      mean_gradient_edge("mean_gradient-edge"), 
      mean_gradient_cube("mean_gradient-cube"), 
      mean_gradient_laplace("mean_gradient-laplace"), 
      gradient_standard_deviation_edge("gradient_standard_deviation-edge"), 
      gradient_standard_deviation_cube("gradient_standard_deviation-cube"), 
      zero_gradient_edge("zero_gradient-edge"),
      num_components("num_components")
      { Init(num_rows, min_scalar); };

    // set routines
    void HideAllExceptScalar();
    void ComputeSum();

    // get routines
    SCALAR_TYPE MinScalar() const { return(min_scalar); };

    // write routines
    void WriteColumnLabels(ostream & out, const string & separator) const;
    void WriteColumnData(ostream & out, const string & separator, const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (ostream & out, const string & separator, const NUM_TYPE width,
       const double normalization_factor) const;
  };

  /// Table of gradient distribution
  class GRADIENT_TABLE:public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  protected:
    void Init(const NUM_TYPE num_rows);
    
  public:
    DATA_COLUMN<NUM_TYPE, GRADIENT_TYPE> gradient;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> edge_span_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> cube_span_freq;

    GRADIENT_TABLE(const NUM_TYPE num_rows):
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(num_rows),
      gradient("gradient"), 
      edge_span_freq("edge_span_freq"),
      cube_span_freq("cube_span_freq")
      { Init(num_rows); };

    // set routines
    void HideAllExceptGradient();
    void ComputeSum();

    // write routines
    void WriteColumnLabels(ostream & out, const string & separator) const;
    void WriteColumnData(ostream & out, const string & separator, const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (ostream & out, const string & separator, const NUM_TYPE width,
       const double normalization_factor) const;
  };

  /// Average statistics
  template <class NTYPE, class DTYPE>
    class AVERAGE_STAT {

    protected:

    NTYPE num;                         ///< Number of terms added to sum.
    DTYPE sum;                         ///< Sum.
    DTYPE sum_of_squares;              ///< Sum of squares.

    public:
    
    AVERAGE_STAT();                    ///< Constructor.

    // Get functions
    NTYPE Num() const 
      { return(num); };
    DTYPE Sum() const 
      { return(sum); };
    DTYPE SumOfSquares() const 
      { return(sum_of_squares); };

    // Compute functions
    DTYPE Average() const;
    DTYPE StandardDeviation() const;

    // Set functions
    void SetToZero();                  ///< Set all values to zero.
    void Add(const DTYPE x);           ///< Add x to sum. Increment num.
  };

  /// Gradient statistics
  class GRADIENT_STATISTICS {

  public:
    AVERAGE_STAT<NUM_TYPE, GRADIENT_TYPE> average_edge_span;
    AVERAGE_STAT<NUM_TYPE, GRADIENT_TYPE> average_cube_span;
    NUM_TYPE num_edges_with_span_zero;
    NUM_TYPE num_cubes_with_span_zero;

    GRADIENT_STATISTICS() {
      num_edges_with_span_zero = 0;
      num_cubes_with_span_zero = 0;
    };
  };


// **************************************************
// DATA_COLUMN member functions
// **************************************************

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::Include()
    {
      is_included = true;

      if (this->NumRows() > 0 && this->data == NULL) 
	{ this->data = new DTYPE[this->NumRows()]; }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::ComputeSum()
    {
      if (IsIncluded()) 
        { sum = IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>::Sum(); }
      else
        { sum = 0; }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::WriteLabel
    (ostream & out, const string & separator) const
    {
      if (IsIncluded() && !IsHidden()) 
        { out << separator << this->Label(); }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::WriteData
    (ostream & out, const string & separator, const NTYPE width, const NTYPE irow) const
    {
      if (IsIncluded() && !IsHidden()) {
        out << separator;
        out.width(width);
        out << this->data[irow]; 
      }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::WriteNormalizedData
    (ostream & out, const string & separator, const NTYPE width,
     const double normalization_factor, const NTYPE irow) const
    {
      if (IsIncluded() && !IsHidden()) {
        double v = normalization_factor*double(this->data[irow])/Sum();
        out << separator;
        out.width(width);
        out << v; 
      }
    }

#define M_APPLY_TO_ALL_EXCEPT_SCALAR(F_X) \
  scalar_freq.F_X;			  \
  isoarea_edge.F_X;\
  ivol_edge.F_X;\
  isoarea_cube.F_X;\
  ivol_cube.F_X;\
  edge_dim.F_X;\
  box_dim.F_X;\
  total_gradient_edge.F_X;\
  total_gradient_cube.F_X;\
  total_gradient_laplace.F_X;\
  mean_gradient_edge.F_X;\
  mean_gradient_cube.F_X;\
  mean_gradient_laplace.F_X;\
  gradient_standard_deviation_edge.F_X;\
  gradient_standard_deviation_cube.F_X;\
  zero_gradient_edge.F_X;\
  num_components.F_X;

#define M_APPLY_TO_ALL(F_X) \
  scalar.F_X;\
  M_APPLY_TO_ALL_EXCEPT_SCALAR(F_X);

  void DATA_TABLE::Init
    (const NUM_TYPE num_rows, const SCALAR_TYPE min_scalar)
  {
    this->min_scalar = min_scalar;

    // set num_rows but do not allocate memory
    M_APPLY_TO_ALL(num_rows = num_rows);
  }

  void DATA_TABLE::HideAllExceptScalar()
  { 
    M_APPLY_TO_ALL_EXCEPT_SCALAR(Hide()); 
    scalar.Show();
  }

  void DATA_TABLE::ComputeSum()
  { M_APPLY_TO_ALL(ComputeSum()); }

  void DATA_TABLE::WriteColumnLabels(ostream & out, const string & separator) const
  {
    M_APPLY_TO_ALL(WriteLabel(out, separator));
  }

  void DATA_TABLE::WriteColumnData
    (ostream & out, const string & separator, const NUM_TYPE width) const
  {
    for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
      scalar.WriteData(out, "", width, irow);
      M_APPLY_TO_ALL_EXCEPT_SCALAR(WriteData(out, "  ", width, irow));
      out << endl;
    }
  }

  void DATA_TABLE::WriteNormalizedColumnData
    (ostream & out, const string & separator, const NUM_TYPE width,
     const double normalization_factor) const
  {
    for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
      scalar.WriteData(out, "", width, irow);
      M_APPLY_TO_ALL_EXCEPT_SCALAR
        (WriteNormalizedData(out, "  ", width, normalization_factor, irow));
      out << endl;
    }
  }

#undef M_APPLY_TO_ALL_EXCEPT_SCALAR
#undef M_APPLY_TO_ALL

// **************************************************
// GRADIENT_TABLE member functions
// **************************************************

  void GRADIENT_TABLE::Init(const NUM_TYPE num_rows)
  {
    // set num_rows but do not allocate memory
    gradient.num_rows = num_rows;
    edge_span_freq.num_rows = num_rows;
    cube_span_freq.num_rows = num_rows;
  }

  void GRADIENT_TABLE::HideAllExceptGradient()
  { 
    gradient.Show();
    edge_span_freq.Hide();
    cube_span_freq.Hide();
  }

  void GRADIENT_TABLE::ComputeSum()
  {
    edge_span_freq.ComputeSum();
    cube_span_freq.ComputeSum();
  }

  void GRADIENT_TABLE::WriteColumnLabels(ostream & out, const string & separator) const
  {
    edge_span_freq.WriteLabel(out, separator);
    cube_span_freq.WriteLabel(out, separator);
  }

  void GRADIENT_TABLE::WriteColumnData
    (ostream & out, const string & separator, const NUM_TYPE width) const
  {
    for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
      gradient.WriteData(out, "", width, irow);
      edge_span_freq.WriteData(out, "  ", width, irow);
      cube_span_freq.WriteData(out, "  ", width, irow);
      out << endl;
    }
  }

  void GRADIENT_TABLE::WriteNormalizedColumnData
    (ostream & out, const string & separator, const NUM_TYPE width,
     const double normalization_factor) const
  {
    for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
      gradient.WriteData(out, "", width, irow);
      edge_span_freq.WriteNormalizedData(out, "  ", width, normalization_factor, irow);
      cube_span_freq.WriteNormalizedData(out, "  ", width, normalization_factor, irow);
      out << endl;
    }
  }

// **************************************************
// STATISTICS member functions
// **************************************************

  template <class NTYPE, class DTYPE>
    AVERAGE_STAT<NTYPE,DTYPE>::AVERAGE_STAT()
    {
      SetToZero();
    };

  template <class NTYPE, class DTYPE>
    void AVERAGE_STAT<NTYPE,DTYPE>::SetToZero()
    {
      num = 0;
      sum = 0;
      sum_of_squares = 0;
    };

  template <class NTYPE, class DTYPE>
    void AVERAGE_STAT<NTYPE,DTYPE>::Add(const DTYPE x)
    {
      sum += x;
      sum_of_squares += (x*x);
      num++;
    };

  template <class NTYPE, class DTYPE>
    DTYPE AVERAGE_STAT<NTYPE,DTYPE>::Average() const
    {
      if (num > 0) { return(sum/num); }
      else { return(0); };
    };

  template <class NTYPE, class DTYPE>
    DTYPE AVERAGE_STAT<NTYPE,DTYPE>::StandardDeviation() const
    {
      if (num > 0) { 
	DTYPE average = Average();
	DTYPE average_sum_of_squares = SumOfSquares()/Num();
	DTYPE diff = average_sum_of_squares - average*average;

	if (diff <= 0) { 
	  // if numerical error causes diff to be negative, return 0.
	  return(0); 
	}

	DTYPE stdev = DTYPE(sqrt(double(diff)));
	return(stdev);
      }
      else { return(0); };
    };

};

#endif

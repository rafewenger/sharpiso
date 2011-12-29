/// \file ijksvector_grid.txx
/// ijk templates defining vector grid classes and functions.
/// Note: vector grid stores mathematical vectors, not C++ vectors.
/// - Version 0.1.1

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

#ifndef _IJKVECTOR_GRID_
#define _IJKVECTOR_GRID_

#include <algorithm>
#include <cmath>

#include "ijk.txx"
#include "ijkgrid.txx"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID_BASE
  // **************************************************

  /// Class storing vector grid information.
  /// Does not allocate or destroy any vector information.
  /// Note: VECTOR_GRID_BASE does not provide any way to allocate memory
  ///    for vector or set vector to point to existing memory.
  /// @tparam GRID_CLASS Base grid class.
  /// @tparam STYPE Scalar type. Vector coordinates are of type STYPE.
  /// @tparam LTYPE Length type.
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  class VECTOR_GRID_BASE:public GRID_CLASS {

  protected:
    LTYPE vector_length;   ///< Length of each vector.
    STYPE * vec;
    // Coordinate c of vector(x0,x1,x2,...) equals
    //     vec[c + vector_length*(x0 + x1*axis_size[0] + 
    //                            x2*axis_size[0]*axis_size[1] + ...)]

    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    typedef STYPE SCALAR_TYPE;
    typedef LTYPE LENGTH_TYPE;
    
  public:
    VECTOR_GRID_BASE() { vec = NULL; }; ///< Constructor
    /// Constructor
    /// Note: Constructor does not allocate memory for array vec[].
    VECTOR_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size,
                     const LTYPE vector_length);
    /// Destructor
    /// Note: Destructor does not free memory of array vec[].
    ~VECTOR_GRID_BASE() {};

    // copy constructor and assignment: NOT IMPLEMENTED
    VECTOR_GRID_BASE(const VECTOR_GRID_BASE & vector_grid); 
    const VECTOR_GRID_BASE & operator = (const VECTOR_GRID_BASE & right);

    // set functions
    template <typename VTYPE2, typename LTYPE2>
    void Set               ///< Set value of vertex \a iv, coord \a ic, to \a s.
    (const VTYPE2 iv, const LTYPE2 ic, const STYPE s) 
    { vec[iv*vector_length+ic] = s; };
    template <typename VTYPE2>
    void Set               ///< Set vector value of vertex \a iv to \a sv
    (const VTYPE2 iv, const STYPE * v);
    void SetAll(const STYPE s);    ///< Set all vector coord values to \a s.
    void SetAll(const STYPE * v);  ///< Set all vector values to \a v.
    void CopyVector        ///< Copy vector values from \a vector_grid.
    (const VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> & vector_grid);
    template <typename DTYPE2, typename ATYPE2, typename LTYPE2>
    void SetSize     ///< Set dimensions, axis sizes and vector length
    (const DTYPE2 dimension, const ATYPE2 * axis_size, 
     const LTYPE2 vector_length);

    /// Set vector values of vertices in subsampled grid to s.
    void SetSubsample(const ATYPE scale, const STYPE s);

    /// Set vector values of vertices in region to \a s.
    template <typename GTYPE>
    void SetRegion(const IJK::BOX<GTYPE> & box, 
                   const ATYPE scale, const STYPE s);

    /// Set vector values of vertices in list of regions to \a s.
    template <typename GTYPE>
    void SetRegion(const std::vector< IJK::BOX<GTYPE> > & box_list, 
                   const ATYPE scale, const STYPE s);

    // get functions
    LTYPE VectorLength() const { return(vector_length); };
    STYPE * VectorPtr() { return(vec); };
    STYPE * VectorPtr(const VTYPE iv) 
    { return(vec+iv*vector_length); };
    const STYPE * VectorPtrConst() const { return(vec); };
    const STYPE * VectorPtrConst(const VTYPE iv) const 
    { return(vec+iv*vector_length); };
    const STYPE * End() const 
    { return(vec+vector_length*this->NumVertices()); };
    STYPE Vector(const VTYPE iv, const LTYPE ic) const 
    { return(vec[ic+iv*vector_length]); };
    STYPE ComputeMagnitudeSquared(const VTYPE iv) const;
    STYPE ComputeMagnitude(const VTYPE iv) const;
  };

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID_ALLOC
  // **************************************************

  /// Class storing vector grid information
  /// Allocates and destroys array vec[]
  template <class VECTOR_GRID_BASE_CLASS>
  class VECTOR_GRID_ALLOC:public VECTOR_GRID_BASE_CLASS {

  protected:
    typedef typename VECTOR_GRID_BASE_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename VECTOR_GRID_BASE_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename VECTOR_GRID_BASE_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename VECTOR_GRID_BASE_CLASS::NUMBER_TYPE NTYPE;
    typedef typename VECTOR_GRID_BASE_CLASS::SCALAR_TYPE STYPE;
    typedef typename VECTOR_GRID_BASE_CLASS::LENGTH_TYPE LTYPE;

    void FreeAll();

    /// Copy (x0,x1,...,xd) in vector_grid to p*(x0,x1,...,xd) 
    ///   in current grid where p is supersample_period. 
    template <class GTYPE, class PTYPE>
    void SupersampleCopy
    (const GTYPE & vector_grid, const PTYPE supersample_period);

    /// Linear interpolate between supersampled vertices.
    template <typename PTYPE>
    void LinearInterpolate(const PTYPE supersample_period);

  public:
    VECTOR_GRID_ALLOC();
    VECTOR_GRID_ALLOC(const DTYPE dimension, const ATYPE * axis_size,
                      const LTYPE vector_length);
    ~VECTOR_GRID_ALLOC() { FreeAll(); };

    /// Copy vector grid
    template <typename GTYPE>
    void Copy(const GTYPE & vector_grid); 
    VECTOR_GRID_ALLOC                                   ///< Copy VECTOR_GRID
    (const VECTOR_GRID_BASE_CLASS & vector_grid2)
    { Copy(vector_grid2); }
    const VECTOR_GRID_ALLOC & operator =                ///< Copy VECTOR_GRID
    (const VECTOR_GRID_BASE_CLASS & right);

    // set functions
    template <typename DTYPE2, typename ATYPE2, typename LTYPE2>
    void SetSize     ///< Set dimensions, axis sizes and vector length
    (const DTYPE2 dimension, const ATYPE2 * axis_size, 
     const LTYPE2 vector_length);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2, typename LTYPE2>
    void SetSize       ///< Set dimensions, axis size and vector length
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2, 
     const LTYPE2 vector_length);
    template <typename GTYPE, typename LTYPE2, typename STYPE2>
    void SetSize       ///< Set dimensions, axis size and vector length
    (const VECTOR_GRID_BASE<GTYPE,LTYPE2,STYPE2> & grid2);


    /// Subsample \a vector_grid2. Resizes current grid.
    template <typename GTYPE, typename PTYPE>
    void Subsample   
    (const GTYPE & vector_grid2, const PTYPE subsample_period);

    /// Supersample \a vector_grid2.  Resizes current grid.
    template <typename GTYPE, typename PTYPE>
    void Supersample
    (const GTYPE & vector_grid2, const PTYPE supersample_period);

  };

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID
  // **************************************************

  /// Class storing vector grid information
  /// Default class using VECTOR_GRID_BASE
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  class VECTOR_GRID:
    public VECTOR_GRID_ALLOC< VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> >
  {
  protected:
    typedef VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> BCLASS;
    typedef typename BCLASS::DIMENSION_TYPE DTYPE;
    typedef typename BCLASS::AXIS_SIZE_TYPE ATYPE;

  public:

    VECTOR_GRID(){};
    VECTOR_GRID(const DTYPE dimension, const ATYPE * axis_size,
                const LTYPE vector_length):
      VECTOR_GRID_ALLOC<BCLASS>(dimension, axis_size, vector_length) {};
    ~VECTOR_GRID(){};
  };

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID_WRAPPER
  // **************************************************

  /// vector grid wrapper for vector array
  /// Does not allocate or destroy any vector information
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  class VECTOR_GRID_WRAPPER:
    public VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> {

  protected:
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    VECTOR_GRID_WRAPPER(const DTYPE dimension, const ATYPE * axis_size,
                        const LTYPE vector_length, STYPE * v);
    ~VECTOR_GRID_WRAPPER() { this->vec = NULL; };
    // Note: constructor and destructor do not allocate or free vector memory
  };

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID_BASE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  VECTOR_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size,
                   const LTYPE vector_length):
    GRID_CLASS(dimension,axis_size)
  { 
    vec = NULL; 
    this->vector_length = vector_length; 
  };

  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  template <typename DTYPE2, typename ATYPE2, typename LTYPE2>
  void VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size,
          const LTYPE2 vector_length)
  {
    GRID_CLASS::SetSize(dimension, axis_size);
    this->vector_length = vector_length;
  }

  /// Set vector value to v
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  template <typename VTYPE2>
  void VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  Set(const VTYPE2 iv, const STYPE * v)
  {
    STYPE * vec_iv = vec+iv*VectorLength();
    for (LTYPE ic = 0; ic < VectorLength(); ic++) 
      { vec_iv[ic] = v[ic]; }
  }

  /// Set all scalar values to s.
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  void VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  SetAll(const STYPE s)
  {
    for (STYPE * ptr = vec; ptr != End(); ptr++)
      { *ptr = s; }
  }

  /// Set all vectors to v
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  void VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  SetAll(const STYPE * v)
  {
    const LTYPE vector_length = VectorLength();
    LTYPE j = 0;
    
    for (VTYPE i = 0; i < vector_length*this->NumVertices(); i++) {
      vec[i] = v[j];
      j = (j+1)%vector_length;
    }
  }

  /// Copy vector values of scalar_grid to current grid.
  /// Precondition: vector_grid has same axis_size as current grid.
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  void VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::CopyVector
  (const VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> & vector_grid)
  {
    std::copy(vector_grid.VectorPtrConst(), vector_grid.End(),this->vec);
  }

  /// Compute vector magnitude squared
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  STYPE VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  ComputeMagnitudeSquared(const VTYPE iv) const
  {
    STYPE magnitude_squared = 0;
    const STYPE * v_ptr = VectorPtrConst(iv);
    for (DTYPE d = 0; d < this->VectorLength(); d++) {
      STYPE s = *v_ptr;
      magnitude_squared += s*s;
      v_ptr++;
    }

    return(magnitude_squared);
  }


  /// Compute vector magnitude
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  STYPE VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>::
  ComputeMagnitude(const VTYPE iv) const
  {
    STYPE magnitude = ComputeMagnitudeSquared(iv);
    if (magnitude > 0.0) 
      { magnitude = std::sqrt(magnitude); }

    return(magnitude);
  }

  // **************************************************
  // TEMPLATE CLASS VECTOR_GRID_ALLOC MEMBER FUNCTIONS
  // **************************************************

  template <typename VECTOR_GRID_BASE_CLASS>
  VECTOR_GRID_ALLOC<VECTOR_GRID_BASE_CLASS>::VECTOR_GRID_ALLOC()
    // constructor
  {
    this->vec = NULL;
  }

  template <typename VECTOR_GRID_BASE_CLASS>
  VECTOR_GRID_ALLOC<VECTOR_GRID_BASE_CLASS>::
  VECTOR_GRID_ALLOC(const DTYPE dimension, const ATYPE * axis_size, 
                    const LTYPE vector_length) :
    VECTOR_GRID_BASE_CLASS(dimension, axis_size, vector_length)
    // constructor
  {
    this->vec = new STYPE[this->NumVertices()*this->VectorLength()];
  }

  template <typename BASE_CLASS>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::FreeAll()
  {
    delete [] this->vec;
    this->vec = NULL;
    BASE_CLASS::FreeAll();
  }

  template <typename BASE_CLASS>
  template <typename DTYPE2, typename ATYPE2, typename LTYPE2>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size,
          const LTYPE2 vector_length)
  {
    NTYPE numv = this->NumVertices();
    LTYPE vlength = this->VectorLength();
    if (this->vec == NULL || !CompareSize(dimension, axis_size)) {
      BASE_CLASS::SetSize(dimension, axis_size, vector_length);
      if (this->vec == NULL || numv != this->NumVertices() ||
          vlength != this->VectorLength()) {
        if (this->vec != NULL) { delete [] this->vec; };
        this->vec = new STYPE[this->NumVertices()*this->VectorLength()];
      }
    }
  }

  template <typename BASE_CLASS>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2, typename LTYPE2>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2, const LTYPE2 vector_length)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize(), vector_length);
  }

  template <typename BASE_CLASS>
  template <typename GTYPE, typename LTYPE2, typename STYPE2>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::SetSize
  (const VECTOR_GRID_BASE<GTYPE,LTYPE2,STYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize(), grid2.VectorLength());
  }

  /// Copy vector grid
  template <typename BASE_CLASS>
  template <typename GTYPE>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::Copy(const GTYPE & vector_grid2)
  {
    this->SetSize(vector_grid2);
    CopyVector(vector_grid2);
  }

  template <typename BASE_CLASS>
  const VECTOR_GRID_ALLOC<BASE_CLASS> & 
  VECTOR_GRID_ALLOC<BASE_CLASS>::operator = 
  (const BASE_CLASS & right)
    // copy assignment	
  {
    if (&right != this) {         // avoid self-assignment
      Copy(right);
    }
  }

  template <typename BASE_CLASS>
  template <typename GCLASS, typename PTYPE>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::Subsample
  (const GCLASS & vector_grid, const PTYPE subsample_period)
  {
    const DTYPE dimension = vector_grid.Dimension();
    IJK::ARRAY<VTYPE> axis_increment(dimension);
    IJK::ARRAY<VTYPE> vprimary(dimension);
    IJK::ARRAY<ATYPE> subsampled_axis_size(dimension);
    IJK::PROCEDURE_ERROR error("VECTOR_GRID::Subsample");

    for (DTYPE d = 0; d < dimension; d++) {
      subsampled_axis_size[d] = 
        compute_subsample_size(vector_grid.AxisSize(d), subsample_period);
    }

    SetSize(dimension, subsampled_axis_size.PtrConst(),
            vector_grid.VectorLength());

    if (this->NumVertices() < 1) { return; };

    compute_increment(dimension, vector_grid.AxisSize(), 
                      axis_increment.Ptr());

    for (DTYPE d = 0; d < dimension; d++) { vprimary[d] = 0; };

    VTYPE iv = 0;
    while (vprimary[0] < vector_grid.NumVertices()) {
      VTYPE jv = vprimary[0];
      this->Set(iv, vector_grid.VectorPtrConst(jv));
      iv++;

      DTYPE d = 0;
      while (d+1 < dimension && 
             (vprimary[d] + subsample_period*axis_increment[d] >=
              vprimary[d+1] + axis_increment[d+1])) {
        d++;
      }

      // go to next subspace of dimension d
      vprimary[d] = vprimary[d] + subsample_period*axis_increment[d];

      for (DTYPE d2 = 0; d2 < d; d2++) 
        { vprimary[d2] = vprimary[d]; };
    }

    if (iv != this->NumVertices()) {
      error.AddMessage("Programming error.  Subsampled ", iv, 
                       " grid vertices.");
      error.AddMessage("Subsampled grid should have ", this->NumVertices(),
                       " vertices.");
      throw error;
    }

  }


  template <typename BASE_CLASS>
  template <typename GTYPE, typename PTYPE>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::Supersample
  (const GTYPE & vector_grid2, const PTYPE supersample_period)
  {
    const DTYPE dimension = vector_grid2.Dimension();
    IJK::ARRAY<ATYPE> supersampled_axis_size(dimension);
    IJK::PROCEDURE_ERROR error("SUPERSAMPLE_GRID::Supersample");

    for (DTYPE d = 0; d < dimension; d++) {
      supersampled_axis_size[d] = 
        compute_supersample_size(vector_grid2.AxisSize(d), supersample_period);
    }

    SetSize(dimension, supersampled_axis_size.PtrConst());

    if (this->NumVertices() < 1) { return; };

    SupersampleCopy(vector_grid2, supersample_period);
    LinearInterpolate(supersample_period);
  }

  template <typename BASE_CLASS>
  template <typename GTYPE, typename PTYPE>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::SupersampleCopy
  (const GTYPE & vector_grid2, const PTYPE supersample_period)
  {
    IJK::CONSTANT<ATYPE,ATYPE> period(supersample_period);
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    NTYPE numv0;
    compute_num_vertices_in_grid_facet
      (vector_grid2, 0, numv0);


    IJK::ARRAY<VTYPE> vlist0(numv0);
    IJK::ARRAY<VTYPE> vlist1(numv0);
  
    get_vertices_in_grid_facet0(vector_grid2, vlist0.Ptr());

    std::copy(this->AxisSize(), this->AxisSize()+this->Dimension(), 
              subgrid_axis_size.Ptr());
    subgrid_axis_size[0] = 1;

    subsample_subgrid_vertices
      (*this, 0, subgrid_axis_size.PtrConst(), period, vlist1.Ptr());

    for (VTYPE x0 = 0; x0 < vector_grid2.AxisSize(0); x0++) {
      VTYPE x1 = x0*supersample_period;
      for (VTYPE i = 0; i < numv0; i++) {
        VTYPE v0 = vlist0[i]+x0;
        VTYPE v1 = vlist1[i]+x1;
        STYPE s = STYPE(vector_grid2.Scalar(v0));
        this->vec[v1] = s;
      }
    }

  }

  template <typename BASE_CLASS>
  template <typename PTYPE>
  void VECTOR_GRID_ALLOC<BASE_CLASS>::LinearInterpolate
  (const PTYPE supersample_period)
  {
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::ARRAY<ATYPE> subsample_period(dimension);
    IJK::ARRAY<VTYPE> axis_increment(dimension);

    compute_increment(*this, axis_increment.Ptr());

    for (DTYPE d = 0; d < this->Dimension(); d++) {
      for (DTYPE j = 0; j < d; j++) { subsample_period[j] = 1; };
      for (DTYPE j = d; j < this->Dimension(); j++) 
        { subsample_period[j] = supersample_period; };
      for (DTYPE j = 0; j < this->Dimension(); j++) 
        { subgrid_axis_size[j] = AxisSize(j); };
      subgrid_axis_size[d] = 1;

      NTYPE numv;
      compute_subsample_size
        (this->Dimension(), subgrid_axis_size.PtrConst(), 
         subsample_period.PtrConst(), numv);

      IJK::ARRAY<VTYPE> vlist(numv);
      subsample_subgrid_vertices
        (*this, 0, subgrid_axis_size.PtrConst(), 
         subsample_period.PtrConst(), vlist.Ptr());

      for (VTYPE x = 0; x+1 < this->AxisSize(d); x += supersample_period) {

        VTYPE inc0 = x*axis_increment[d];
        VTYPE inc1 = inc0 + supersample_period*axis_increment[d];
        for (VTYPE i = 0; i < numv; i++) {
          VTYPE v0 = vlist[i] + inc0;
          VTYPE v1 = vlist[i] + inc1;

          VTYPE v2 = v0;
          for (VTYPE j = 1; j < supersample_period; j++) {
            v2 += axis_increment[d];
            STYPE s0 = this->vec[v0];
            STYPE s1 = this->vec[v1];
            this->vec[v2] = 
              linear_interpolate(s0, 0, s1, supersample_period, j);
          }
        }
      }
    }
  }

  // ******************************************************
  // TEMPLATE CLASS VECTOR_GRID_WRAPPER MEMBER FUNCTIONS
  // ******************************************************

  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  VECTOR_GRID_WRAPPER<GRID_CLASS,LTYPE,STYPE>::
  VECTOR_GRID_WRAPPER
  (const DTYPE dimension, const ATYPE * axis_size, 
   const LTYPE vector_length, STYPE * v):
    VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE>
  (dimension, axis_size, vector_length)
  { this->vec = v; }

  // **************************************************
  // TEMPLATE OUTPUT FUNCTIONS
  // **************************************************

  /// Output vector coord.  Template specialization for bool and char classes.
  template <class STYPE>
  inline void output_vector_coord(std::ostream & out, const STYPE s)
  { out << s; }

  /// Output vector coord.  Template specialization for bool class.
  template <>
  inline void output_vector_coord(std::ostream & out, const bool s)
  { out << int(s); }

  /// Output vector coord.  Template specialization for char class.
  template <>
  inline void output_vector_coord(std::ostream & out, const char s)
  { out << int(s); }

  /// Output vector coord.  Template specialization for unsigned char class.
  template <>
  inline void output_vector_coord(std::ostream & out, const unsigned char s)
  { out << int(s); }

  /// Output vector.
  template <typename LTYPE, typename STYPE>
  inline void output_vector
  (std::ostream & out, const LTYPE vector_length, const STYPE * v)
  {
    out << "(";
    for (LTYPE i = 0; i < vector_length; i++) {
      output_vector_coord(out, v[i]);
      if (i+1 < vector_length) { out << ","; }
    }
    out << ")";
  }

  /// Output vector grid (for debugging purposes)
  template <typename GRID_CLASS, typename LTYPE, typename STYPE>
  void output_vector_grid
  (std::ostream & out, 
   const VECTOR_GRID_BASE<GRID_CLASS,LTYPE,STYPE> & vector_grid)
  {
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

    const DTYPE dimension = vector_grid.Dimension();
    const ATYPE * axis_size = vector_grid.AxisSize();
    const LTYPE vector_length = vector_grid.VectorLength();

    using namespace std;
    
    if (vector_grid.NumVertices() < 1) { return; };

    if (vector_grid.Dimension() <= 1) {
      output_vector(out, vector_length, vector_grid.VectorPtrConst(0));
      for (NTYPE i = 1; i < vector_grid.NumVertices(); ++i) {
        out << " ";
        output_vector(out, vector_length, vector_grid.VectorPtrConst(i));
      }
      out << endl;
    }
    else {
      NTYPE num_vertices_in_ridge;
      compute_num_vertices_in_grid_ridge
        (dimension, axis_size, 0, 1, num_vertices_in_ridge);
      IJK::ARRAY<VTYPE> vlist(num_vertices_in_ridge);
      get_vertices_in_grid_ridge(dimension, axis_size, 0, 1, vlist.Ptr());

      for (NTYPE i = 0; i < num_vertices_in_ridge; i++) {
        VTYPE iv = vlist[i];
        for (ATYPE j1 = 0; j1 < axis_size[1]; j1++) {
          output_vector(out, vector_length, vector_grid.VectorPtrConst(iv));
          ++iv;
          for (ATYPE j0 = 1; j0 < axis_size[0]; j0++) {
            out << " ";
            output_vector(out, vector_length, vector_grid.VectorPtrConst(iv));
            ++iv;
          }
          out << endl;
        }
        if (num_vertices_in_ridge > 1) { out << endl; }
      }
    }
  }

}

#endif

/// \file ijkobject_grid.txx
/// ijk templates defining object grid classes and functions.
/// Objects are typically C++ classes.
/// - Version 0.1.1

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

#ifndef _IJKOBJECT_GRID_
#define _IJKOBJECT_GRID_

#include <algorithm>
#include <cmath>

#include "ijk.txx"
#include "ijkgrid.txx"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS OBJECT_GRID_BASE
  // **************************************************

  /// Class storing object grid information.
  /// Does not allocate or destroy any vector information.
  /// Note: OBJECT_GRID_BASE does not provide any way to allocate memory
  ///    for object or set object to point to existing memory.
  /// @tparam GRID_CLASS Base grid class.
  /// @tparam OTYPE Object type.
  template <typename GRID_CLASS, typename OTYPE>
  class OBJECT_GRID_BASE:public GRID_CLASS {

  protected:
    OTYPE * object;        ///< Array of objects.

    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VITYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    typedef OTYPE OBJECT_TYPE;

  public:
    OBJECT_GRID_BASE() {};       ///< Constructor
    /// Constructor
    /// Note: Constructor does not allocate memory for array object[].
    OBJECT_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size);

    /// Destructor
    /// Note: Destructor does not free memory of array object[].
    ~OBJECT_GRID_BASE() {};

    // copy constructor and assignment: NOT IMPLEMENTED
    OBJECT_GRID_BASE(const OBJECT_GRID_BASE & object_grid);
    const OBJECT_GRID_BASE & operator = (const OBJECT_GRID_BASE & right);

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize     ///< Set dimensions and axis sizes
    (const DTYPE2 dimension, const ATYPE2 * axis_size);

    // get functions

    /// Get const object at grid vertex iv.
    template <typename VITYPE2>
    const OTYPE & operator ()(const VITYPE2 iv) const
    { return(object[iv]); }

    /// Get object at grid vertex iv.  Object can be modified.
    template <typename VITYPE2>
    OTYPE & operator [] (const VITYPE2 iv)
    { return(object[iv]); }
  };

  // **************************************************
  // TEMPLATE CLASS OBJECT_GRID_ALLOC
  // **************************************************

  /// Class storing object grid information.
  /// Allocates and destroys array object.
  template <class OBJECT_GRID_BASE_CLASS>
  class OBJECT_GRID_ALLOC:public OBJECT_GRID_BASE_CLASS {

  protected:
    typedef typename OBJECT_GRID_BASE_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename OBJECT_GRID_BASE_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename OBJECT_GRID_BASE_CLASS::VERTEX_INDEX_TYPE VITYPE;
    typedef typename OBJECT_GRID_BASE_CLASS::NUMBER_TYPE NTYPE;
    typedef typename OBJECT_GRID_BASE_CLASS::OBJECT_TYPE OTYPE;

    void FreeAll();

  public:
    OBJECT_GRID_ALLOC();
    OBJECT_GRID_ALLOC(const DTYPE dimension, const ATYPE * axis_size);
    ~OBJECT_GRID_ALLOC() { FreeAll(); };

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize     ///< Set dimensions and axis sizes.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VITYPE2,
              typename NTYPE2>
    void SetSize       ///< Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VITYPE2,NTYPE2> & grid2);
    template <typename GTYPE, typename VCTYPE2>
    void SetSize       ///< Set dimensions and axis size.
    (const OBJECT_GRID_BASE<GTYPE,VCTYPE2> & grid2);
  };

  // **************************************************
  // TEMPLATE CLASS OBJECT_GRID
  // **************************************************

  /// Class storing object grid information
  /// Default class using OBJECT_GRID_BASE
  template <typename GRID_CLASS, typename OTYPE>
  class OBJECT_GRID:
    public OBJECT_GRID_ALLOC< OBJECT_GRID_BASE<GRID_CLASS,OTYPE> >
  {
  protected:
    typedef OBJECT_GRID_BASE<GRID_CLASS,OTYPE> BCLASS;
    typedef typename BCLASS::DIMENSION_TYPE DTYPE;
    typedef typename BCLASS::AXIS_SIZE_TYPE ATYPE;

  public:

    OBJECT_GRID(){};
    OBJECT_GRID(const DTYPE dimension, const ATYPE * axis_size):
      OBJECT_GRID_ALLOC<BCLASS>(dimension, axis_size) {};
    ~OBJECT_GRID(){};
  };

  // **************************************************
  // TEMPLATE CLASS OBJECT_GRID_BASE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor
  template <typename GRID_CLASS, typename OTYPE>
  OBJECT_GRID_BASE<GRID_CLASS,OTYPE>::
  OBJECT_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size):
    GRID_CLASS(dimension,axis_size)
  {
    object = NULL;
  };

  template <typename GRID_CLASS, typename OTYPE>
  template <typename DTYPE2, typename ATYPE2>
  void OBJECT_GRID_BASE<GRID_CLASS,OTYPE>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    GRID_CLASS::SetSize(dimension, axis_size);
  }


  // **************************************************
  // TEMPLATE CLASS OBJECT_GRID_ALLOC MEMBER FUNCTIONS
  // **************************************************

  template <typename OBJECT_GRID_BASE_CLASS>
  OBJECT_GRID_ALLOC<OBJECT_GRID_BASE_CLASS>::OBJECT_GRID_ALLOC()
    // constructor
  {
    this->object = NULL;
  }

  template <typename OBJECT_GRID_BASE_CLASS>
  OBJECT_GRID_ALLOC<OBJECT_GRID_BASE_CLASS>::
  OBJECT_GRID_ALLOC(const DTYPE dimension, const ATYPE * axis_size) :
    OBJECT_GRID_BASE_CLASS(dimension, axis_size)
    // constructor
  {
    this->object = new OTYPE[this->NumVertices()];
  }

  template <typename BASE_CLASS>
  void OBJECT_GRID_ALLOC<BASE_CLASS>::FreeAll()
  {
    delete [] this->object;
    this->object = NULL;
    BASE_CLASS::FreeAll();
  }

  template <typename BASE_CLASS>
  template <typename DTYPE2, typename ATYPE2>
  void OBJECT_GRID_ALLOC<BASE_CLASS>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    NTYPE numv = this->NumVertices();
    if (this->object == NULL || !this->CompareSize(dimension, axis_size)) {
      BASE_CLASS::SetSize(dimension, axis_size);
      if (this->object == NULL || numv != this->NumVertices()) {
        if (this->object != NULL) { delete [] this->object; };
        this->object = new OTYPE[this->NumVertices()];
      }
    }
  }

  template <typename BASE_CLASS>
  template <typename DTYPE2, typename ATYPE2, typename VITYPE2,
            typename NTYPE2>
  void OBJECT_GRID_ALLOC<BASE_CLASS>::SetSize
  (const GRID<DTYPE2,ATYPE2,VITYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  template <typename BASE_CLASS>
  template <typename GTYPE, typename OTYPE2>
  void OBJECT_GRID_ALLOC<BASE_CLASS>::SetSize
  (const OBJECT_GRID_BASE<GTYPE,OTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  // **************************************************
  // TEMPLATE OUTPUT FUNCTIONS
  // **************************************************

  /// Output object grid (for debugging purposes)
  /// @pre Requires that operator "ostream << (ostream,OBJECT)" is defined.
  template <typename GRID_CLASS, typename OTYPE>
  void output_object_grid
  (std::ostream & out,
   const OBJECT_GRID_BASE<GRID_CLASS,OTYPE> & object_grid)
  {
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VITYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

    const DTYPE dimension = object_grid.Dimension();
    const ATYPE * axis_size = object_grid.AxisSize();

    using namespace std;

    if (object_grid.NumVertices() < 1) { return; };

    if (object_grid.Dimension() <= 1) {
      out << object_grid(0);
      for (NTYPE i = 1; i < object_grid.NumVertices(); ++i) {
        out << " " << object_grid(i);
      }
      out << endl;
    }
    else {
      NTYPE num_vertices_in_ridge;
      compute_num_vertices_in_grid_ridge
        (dimension, axis_size, 0, 1, num_vertices_in_ridge);
      IJK::ARRAY<VITYPE> vlist(num_vertices_in_ridge);
      get_vertices_in_grid_ridge(dimension, axis_size, 0, 1, vlist.Ptr());

      for (NTYPE i = 0; i < num_vertices_in_ridge; i++) {
        VITYPE iv = vlist[i];
        for (ATYPE j1 = 0; j1 < axis_size[1]; j1++) {
          out << object_grid(iv);
          ++iv;
          for (ATYPE j0 = 1; j0 < axis_size[0]; j0++) {
            out << " " << object_grid(iv);
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

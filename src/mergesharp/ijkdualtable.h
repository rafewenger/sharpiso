/// \file ijkdualtable.h
/// Class containing dual lookup table of isosurface vertices.
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

#ifndef _IJKDUALTABLE_
#define _IJKDUALTABLE_

#include "ijk.txx"

/// Classes and routines for storing and manipulating 
///   dual isosurface lookup table.
namespace IJKDUALTABLE {

  // **************************************************
  // TYPES
  // **************************************************

  typedef int TABLE_INDEX;                    ///< Index of table entry.
  typedef unsigned char ISODUAL_VERTEX_INDEX; ///< Index of isosurface vertex.

  // Forward definition.
  class FIND_COMPONENT;
  
  // **************************************************
  // COMPUTE FUNCTIONS
  // **************************************************

  /// Compute complement index.
  template <typename NTYPE>
  inline TABLE_INDEX compute_complement
  (const TABLE_INDEX ival, const NTYPE num_table_entries)
  { return(num_table_entries-1-ival); }
    
  // **************************************************
  // ISODUAL TABLE
  // **************************************************

  /// Dual isosurface lookup table.
  /// Stores isosurface vertices and incident faces for each configuration 
  ///   of +/- labels at cube vertices.
  class ISODUAL_TABLE {

  protected:

    /// Entry in the dual isosurface lookup table.
    class ISODUAL_TABLE_ENTRY {

    public:
      int num_vertices;           ///< Number of dualiso vertices in cube.
      ISODUAL_TABLE_ENTRY();      ///< constructor
      ~ISODUAL_TABLE_ENTRY();     ///< destructor

      /// incident_isovertex[kf] = Isosurface vertex incident on face kf.
      ///       Face kf is dual to polytope edge kf.
      ISODUAL_VERTEX_INDEX * incident_isovertex;

      /// is_bipolar[ke] = True if polytope edge ke is bipolar.
      ///       Cube edge ke is dual to isosurface face kf.
      bool * is_bipolar;

      /// Allocate incident_isovert[] and is_bipolar[].
      void Allocate(const int num_poly_edges);

      bool Check(IJK::ERROR & error_msg) const;
      void FreeAll();                            // free all memory
    };


  public:

  /// Index of entry in isosurface lookup table.
  /// Define within ISODUAL_TABLE for use in templates.
  typedef IJKDUALTABLE::TABLE_INDEX TABLE_INDEX;    

  protected:

    int dimension;                  ///< Dimension
    int num_poly_vertices;          ///< Number of polytope vertices.
    int num_poly_edges;             ///< Number of polytope edges;
    ISODUAL_TABLE_ENTRY * entry;    ///< Array of dual isosurface table entries.
    long num_table_entries;         ///< Number of entries in table.

    /// Maximum number of vertices allowed for cube.
    int max_num_vertices; 

    /// True, if array num_table_entries[] is allocated.
    bool is_table_allocated;  

    /// Initialization routine.
    void Init(const int dimension);

 public:
    ISODUAL_TABLE();
    ISODUAL_TABLE(const int d);
    ~ISODUAL_TABLE();                ///< Destructor

    // Get functions.
    int Dimension() const            ///< Return dimension.
    { return(dimension); };
    int NumPolyVertices() const      ///< Return number of polytope vertices.
    { return(num_poly_vertices); };
    int NumPolyEdges() const         ///< Return number of polytope edges.
    { return(num_poly_edges); };

    /// Return number of lookup table entries.
    int NumTableEntries() const { return(num_table_entries); };

    /// Return complement of table index it
    int Complement(const int it) const
    { return(compute_complement(it, num_table_entries)); }

    /// Return number of vertices in isosurface patch for table entry \a it.
    int NumIsoVertices(const TABLE_INDEX it) const
    { return(entry[it].num_vertices); }; 

    /// Return index of isosurface vertex incident on face kf.
    /// Undefined if polytope edge k is not bipolar.
    /// @param it Index of table entry.
    /// @param kf Isosurface face kf, dual to polytope edge kf.
    ISODUAL_VERTEX_INDEX IncidentIsoVertex
    (const TABLE_INDEX it, const int kf) const
    { return(entry[it].incident_isovertex[kf]); };

    /// Return true if edge ke is bipolar.
    /// @param it Index of table entry.
    /// @param ke Polytope edge ke, dual to isosurface face ke.
    bool IsBipolar(const TABLE_INDEX it, const int ke) const
    { return(entry[it].is_bipolar[ke]); };


    /// Return maximum number of polytope vertices permitted in any table.
    /// Note: Even tables for polytopes of this size are probably impossible 
    ///   to compute/store.
    int MaxNumVertices() const { return(max_num_vertices); };

    /// Return true if table memory is allocated.
    bool IsTableAllocated() const
    { return(is_table_allocated); };

    // Set functions.
    void SetDimension(const int d);
    void SetNumPolyVertices(const int num_vertices);
    void SetNumPolyEdges(const int num_edges);

    /// Allocate table
    void SetNumTableEntries(const int num_table_entries);

    // Check functions
    bool CheckDimension(const int d) const;
    bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
    bool CheckTable(IJK::ERROR & error_msg) const;
    bool Check(IJK::ERROR & error_msg) const;

    virtual void FreeAll();                     /// Free all memory.
  };

  typedef ISODUAL_TABLE * ISODUAL_TABLE_PTR;

  // **************************************************
  // ISODUAL CUBE TABLE
  // **************************************************

  /// ISODUAL_TABLE based on cube.
  class ISODUAL_CUBE_TABLE:public ISODUAL_TABLE {

  protected:
    bool flag_separate_neg;                 ///< If true, separate negative vertices.

    /// If true, always separate two diagonally opposite
    ///   positive or negative vertices.
    bool flag_always_separate_opposite;     

    /// Create table entries.
    /// @param flag_opposite_vertices If true, always separate two
    ///        diagonally opposite positive or negative vertices.
    void CreateTableEntries
      (const bool flag_separate_neg, const bool flag_separate_opposite);

  public:
    ISODUAL_CUBE_TABLE() {};
    ISODUAL_CUBE_TABLE(const int dimension);
    ISODUAL_CUBE_TABLE
      (const int dimension, const bool flag_opposite_vertices);
    ISODUAL_CUBE_TABLE
      (const int dimension, const bool separate_neg,
       const bool flag_separate_opposite);

    // Set functions.
    void SetDimension(const int d);

    void Create(const int dimension);
    void Create(const int dimension, const bool flag_separate_opposite);
    void Create(const int dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);

    /// Undefine function.
    void SetNumTableEntries(const int num_table_entries);
  };


  // **************************************************
  // CLASS FIND_COMPONENT
  // **************************************************

  /// Find connected component among cube vertices.
  class FIND_COMPONENT {

  protected:
    int dimension;
    int num_cube_vertices;
    bool * vertex_flag;
    int * component;

  public:
    FIND_COMPONENT(const int dimension);
    ~FIND_COMPONENT();

    // set functions
    void SetVertexFlags(const TABLE_INDEX ival);
    void NegateVertexFlags();
    void ClearAll();

    // get functions
    int Dimension() const
    { return(dimension); }
    bool VertexFlag(const int i) const
    { return(vertex_flag[i]); }
    int Component(const int i) const
    { return(component[i]); }
    int NumCubeVertices() const
    { return(num_cube_vertices); }

    /// Search starting at vertex i.
    /// @pre icomp is not zero.
    void Search(const int i, const int icomp);

    /// Search facet starting at vertex i.
    /// @pre Facet kf contains vertex i.
    /// @pre icomp is not zero.
    void SearchFacet(const int kf, const int i, const int icomp);

    /// Compute number of components.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    int ComputeNumComponents
    (const int ientry, const bool flag_positive);

    /// Compute number of components in facet.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    int ComputeNumComponentsInFacet
    (const int ientry, const int kf, const bool flag_positive);
  };
  
  // **************************************************
  // UTILITY FUNCTIONS
  // **************************************************

  /// Calculate number of entries required in ISOSURFACE_TABLE.
  unsigned long calculate_num_entries(const int num_vert, const int num_colors);

  /// Convert integer to boolean flags.
  void convert2bool
    (const TABLE_INDEX ival, bool * flag, const unsigned int num_flags);

  /// Return true if ival represents two diagonally opposite ones.
  bool is_two_opposite_ones
    (const TABLE_INDEX ival, const int num_bits);

  /// Return true if ival represents two diagonally opposite zeros.
  bool is_two_opposite_zeros
    (const TABLE_INDEX ival, const int num_bits);

}

#endif

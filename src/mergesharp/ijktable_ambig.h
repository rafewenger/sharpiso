/// \file ijktabl_ambige.h
/// Class containing ambiguity information about isosurface lookup tables.
/// Version 0.3.1

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

#ifndef _IJKTABLE_AMBIG_
#define _IJKTABLE_AMBIG_

#include "ijk.txx"

#include "ijktable_poly.h"

/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKTABLE {

  /// Index of entry in ambig table.
  typedef int AMBIG_TABLE_INDEX;    

  // **************************************************
  // AMBIGUITY INFORMATION FOR ISOSURFACE TABLE
  // **************************************************

  /// Isosurface table with ambiguity information.
  class ISOSURFACE_TABLE_AMBIG_INFO {

  protected:
    long num_table_entries;         ///< Number of table entries
    bool * is_ambiguous;            ///< True for ambiguous configurations.
    FACET_INDEX * 
    num_ambiguous_facets;         ///< Number of ambiguous facts.
    FACET_SET * ambiguous_facet;    ///< k'th bit is 1 if facet k is ambiguous

    void Init();                    ///< Initialization routine.
    void Alloc(const long num_table_entries);  ///< Allocate memory.
    void FreeAll();                 ///< Free all memory.

    /// Compute ambiguity information.
    void ComputeAmbiguityInformation
    (const ISOSURFACE_TABLE_POLYHEDRON & poly);

  public:
    
  // constructors
    ISOSURFACE_TABLE_AMBIG_INFO() { Init(); };
    ~ISOSURFACE_TABLE_AMBIG_INFO();                // destructor

    // get functions
    bool IsAmbiguous(const AMBIG_TABLE_INDEX it) const
    { return(is_ambiguous[it]); };
    FACET_INDEX NumAmbiguousFacets(const AMBIG_TABLE_INDEX it) const
    { return(num_ambiguous_facets[it]); };
    // get functions
    bool IsFacetAmbiguous
      (const AMBIG_TABLE_INDEX it, const FACET_INDEX jf) const
    { return((ambiguous_facet[it] & ((1L) << jf)) != 0); };
    long NumTableEntries() const { return(num_table_entries); };

    /// Set ambiguity table.
    void SetAmbiguityTable(const ISOSURFACE_TABLE_POLYHEDRON & poly);

    /// Set cube ambiguity table.
    /// @param use_lex_order If true, use lexicographic order 
    ///        on edges and facets.
    void SetCubeAmbiguityTable
    (const int dimension, const bool use_lex_order);

  };

  // **************************************************
  // AMBIGUITY ROUTINES
  // **************************************************

  /// Return true if isosurface topology is ambiguous
  /// @param vertex_sign[i] = Sign of isosurface vertex i. (0 or 1);
  bool is_poly_ambiguous
  (const ISOSURFACE_TABLE_POLYHEDRON & poly, const int * vertex_sign);

  /// Return true if facet jf is ambiguous.
  /// @param vertex_sign[i] = Sign of isosurface vertex i.
  bool is_facet_ambiguous
  (const ISOSURFACE_TABLE_POLYHEDRON & poly,
   const FACET_INDEX jf, const int * vertex_sign);

  /// Return number of vertices connected by edges to iv with same sign as iv.
  /// @param vertex_sign[i] = Sign of isosurface vertex i.
  int compute_num_connected
  (const ISOSURFACE_TABLE_POLYHEDRON & poly,
   const int iv, const int * vertex_sign);

  /// Return number of vertices connected by edges in facet jf to vertex iv with same sign as iv.
  /// @param jf = Facet index.
  /// @param iv = Vertex index.  Precondition: Vertex iv is in facet jf.
  /// @param vertex_sign[i] = Sign of isosurface vertex i.
  int compute_num_connected_in_facet
  (const ISOSURFACE_TABLE_POLYHEDRON & poly,
   const FACET_INDEX jf, const int iv, const int * vertex_sign);

  /// Compute ambiguous facets.
  /// @param vertex_sign[i] = Sign of isosurface vertex i.
  void compute_ambiguous_facets
  (const ISOSURFACE_TABLE_POLYHEDRON & poly, const int * vertex_sign,
   FACET_SET & facet_set, FACET_INDEX & num_ambiguous_facets);

}

#endif

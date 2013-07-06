/// \file ijkdualtable_ambig.h
/// Class containing ambiguity information about dual isosurface lookup table
///   for cubes.
/// Version 0.0.1

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

#ifndef _IJKDUALTABLE_AMBIG_
#define _IJKDUALTABLE_AMBIG_

#include "ijk.txx"

#include "ijkdualtable.h"


/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKDUALTABLE {

  // **************************************************
  // TYPES
  // **************************************************

  /// Index of entry in ambig table.
  typedef int AMBIG_TABLE_INDEX;    

  typedef unsigned char FACET_INDEX;   ///< Index of facet.
  typedef int FACET;          ///< Bits representing vertices in facet.
  typedef int FACET_SET;      ///< Bits representing set of facets.

  // **************************************************
  // AMBIGUITY INFORMATION FOR ISODUAL CUBE TABLE
  // **************************************************

  /// Ambiguity information for isodual cube table.
  class ISODUAL_CUBE_TABLE_AMBIG_INFO {

  protected:
    int dimension;                  ///< Dimension.
    long num_table_entries;         ///< Number of table entries
    bool * is_ambiguous;            ///< True for ambiguous configurations.
    FACET_INDEX * 
    num_ambiguous_facets;           ///< Number of ambiguous facts.
    FACET_SET * ambiguous_facet;    ///< k'th bit is 1 if facet k is ambiguous

    void Init();                    ///< Initialization routine.
    void Alloc(const long num_table_entries);  ///< Allocate memory.
    void FreeAll();                 ///< Free all memory.

    /// Compute ambiguity information.
    void ComputeAmbiguityInformation();


  public:
    
  // constructors
    ISODUAL_CUBE_TABLE_AMBIG_INFO();
    ISODUAL_CUBE_TABLE_AMBIG_INFO(const int dimension);
    ~ISODUAL_CUBE_TABLE_AMBIG_INFO();                // destructor

    // get functions
    int Dimension() const { return(dimension); };
    bool IsAmbiguous(const AMBIG_TABLE_INDEX it) const
    { return(is_ambiguous[it]); };
    FACET_INDEX NumAmbiguousFacets(const AMBIG_TABLE_INDEX it) const
    { return(num_ambiguous_facets[it]); };
    // get functions
    bool IsFacetAmbiguous
      (const AMBIG_TABLE_INDEX it, const FACET_INDEX jf) const
    { return((ambiguous_facet[it] & ((1L) << jf)) != 0); };
    long NumTableEntries() const { return(num_table_entries); };
    FACET_SET AmbiguousFacetBits(const AMBIG_TABLE_INDEX it) const
    { return(ambiguous_facet[it]); };
                                 

    /// Set ambiguity table.
    void Set(const int dimension);
  };

  // **************************************************
  // AMBIGUITY ROUTINES
  // **************************************************

  bool is_cube_ambiguous
  (const AMBIG_TABLE_INDEX ientry, FIND_COMPONENT & find_component);

  bool is_cube_facet_ambiguous
  (const AMBIG_TABLE_INDEX ientry, const FACET_INDEX & kf, 
   FIND_COMPONENT & find_component);

  void compute_ambiguous_cube_facets
  (const AMBIG_TABLE_INDEX ientry, const FACET_INDEX num_facets,
   FACET_SET & facet_set,  FACET_INDEX & num_ambiguous_facets,
   FIND_COMPONENT & find_component);
}

#endif

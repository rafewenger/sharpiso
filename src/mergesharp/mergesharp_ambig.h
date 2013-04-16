/// \file mergesharp_ambig.h
/// Routines for handling ambiguous configurations.

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

#ifndef _MERGESHARP_AMBIG_
#define _MERGESHARP_AMBIG_

#include <string>
#include <vector>

#include "ijk.txx"
#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"

#include "ijkdualtable.h"
#include "ijktable_ambig.h"

/// mergesharp_ambig classes and routines.
namespace MERGESHARP {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  typedef unsigned char NUM_COMPONENTS_TYPE;

  class AMBIG_TABLE;

  // **************************************************
  // DETERMINE AMBIGUOUS FACETS
  // **************************************************

  /// Return true if facet is ambiguous.
  bool is_grid_facet_ambiguous
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX facet_v0,
   const NUM_TYPE orth_dir,
   const SCALAR_TYPE isovalue);

  /// Decide and return ambiguity status of facet ifacet.
  /// @param facet_ambig_status[] Array containing facet ambiguity.
  ///        facet_ambig_status[iv*DIM3+d] is the ambiguity of the facet
  ///          containing primary vertex v and orthodonal direction d.
  ///        Note: When iv is on the upper-rightmost grid boundary,
  ///          the facet iv*DIM3+d may not actually be contained in the grid.
  AMBIGUITY_STATUS decide_ambiguous_facet
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const AMBIG_TABLE & ambig_table,
   const SCALAR_TYPE isovalue,
   const VERTEX_INDEX facet_v0,
   const int facet_orth_dir,
   const MERGESHARP_PARAM & mergesharp_param);

  /// Propagate SEPARATE_POS and SEPARATE_NEG across cube facets.
  /// @param facet_ambig_status[] Array containing facet ambiguity.
  ///        facet_ambig_status[iv*DIM3+d] is the ambiguity of the facet
  ///          containing primary vertex v and orthodonal direction d.
  ///        Note: When iv is on the upper-rightmost grid boundary,
  ///          the facet iv*DIM3+d may not actually be contained in the grid.
  void propagate_sep
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const AMBIG_TABLE & ambig_table,
   std::vector<AMBIGUITY_TYPE> & facet_ambig_status);

  /// Set ambiguity of cubes in cube_list[].
  void set_cube_ambiguity
  (const SHARPISO_GRID & grid,
   const std::vector<ISO_VERTEX_INDEX> & cube_list,
   const std::vector<AMBIGUITY_TYPE> & facet_ambig_status,
   std::vector<AMBIGUITY_TYPE> & cube_ambig);

  /// Set ambiguity of cubes in cube_list[].
  void set_cube_ambiguity    
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & cube_list,
   const MERGESHARP_PARAM & mergesharp_param,
   std::vector<AMBIGUITY_TYPE> & cube_ambig);

  /// Set ambiguity information.
  void set_ambiguity_info
  (const std::vector<AMBIGUITY_TYPE> & cube_ambig,
   SHARPISO_INFO & sharpiso_info);


  // **************************************************
  // AMBIG_TABLE
  // **************************************************

  /// Cube ambiguity table.
  class AMBIG_TABLE:public IJKTABLE::ISOSURFACE_TABLE_AMBIG_INFO {

  protected:

    /// Number of connected components of positive vertices.
    NUM_COMPONENTS_TYPE * num_pos_components;

    /// Number of connected components of negative vertices.
    NUM_COMPONENTS_TYPE * num_neg_components;

    void Init();                    ///< Initialization routine.
    void Alloc(const long num_table_entries);  ///< Allocate memory.
    void FreeAll();                 ///< Free all memory.

    void ComputeNumConnectedComponents
    (const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly);

  public:

  // constructors 
    AMBIG_TABLE() { Init(); };
    ~AMBIG_TABLE();                // destructor

    /// Set cube ambiguity table.
    void SetCubeAmbiguityTable();

    // Get functions.
    NUM_COMPONENTS_TYPE NumPosComponents 
    (const IJKTABLE::AMBIG_TABLE_INDEX i) const
    { return(num_pos_components[i]); }

    NUM_COMPONENTS_TYPE NumNegComponents
    (const IJKTABLE::AMBIG_TABLE_INDEX i) const
    { return(num_neg_components[i]); }
  };


  // **************************************************
  // AMBIGUITY ROUTINES
  // **************************************************

  //// Return number of connected components of vertices with sign isign.
  /// @param vertex_sign[i] = Sign of isosurface vertex i.
  int compute_num_connected_components
  (const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly,
   const int sign, const int * vertex_sign);

}

#endif

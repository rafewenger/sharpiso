/// \file ijkdualtable_ambig.cxx
/// Class containing ambiguity information about dual isosurface lookup table
///   for cubes.
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

#include <ctype.h>
#include <limits>
#include <limits.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <set>
#include <vector>

#include "ijkcube.txx"
#include "ijkbits.txx"

#include "ijkdualtable_ambig.h"


using namespace IJK;
using namespace IJKDUALTABLE;
using namespace std;

// **************************************************
// AMBIGUITY INFORMATION FOR ISODUAL CUBE TABLE
// **************************************************

// Constructor
ISODUAL_CUBE_TABLE_AMBIG_INFO::ISODUAL_CUBE_TABLE_AMBIG_INFO()
{
  Init();
}

// Constructor
ISODUAL_CUBE_TABLE_AMBIG_INFO::ISODUAL_CUBE_TABLE_AMBIG_INFO
(const int dimension)
{
  Init();
  Set(dimension);
}

void ISODUAL_CUBE_TABLE_AMBIG_INFO::Init()
{
  num_table_entries = 0;
  is_ambiguous = NULL;
  num_ambiguous_facets = NULL;
  ambiguous_facet = NULL;
  dimension = 0;
}


/// Destructor.
ISODUAL_CUBE_TABLE_AMBIG_INFO::~ISODUAL_CUBE_TABLE_AMBIG_INFO()
{
  FreeAll();
}

/// Free all memory.
void ISODUAL_CUBE_TABLE_AMBIG_INFO::FreeAll()
{
  delete [] is_ambiguous;
  is_ambiguous = NULL;
  delete [] num_ambiguous_facets;
  num_ambiguous_facets = NULL;
  delete [] ambiguous_facet;
  ambiguous_facet = NULL;
  num_table_entries = 0;
}

/// Allocate memory.
/// Free any previously allocated memory.
void ISODUAL_CUBE_TABLE_AMBIG_INFO::
Alloc(const long num_table_entries)
{
  if (this->num_table_entries > 0) { FreeAll(); };

  is_ambiguous = new bool[num_table_entries];
  num_ambiguous_facets = new FACET_INDEX[num_table_entries];
  ambiguous_facet = new FACET_SET[num_table_entries];

  this->num_table_entries = num_table_entries;
}


/// Compute ambiguity information.
void ISODUAL_CUBE_TABLE_AMBIG_INFO::ComputeAmbiguityInformation()
{
  const int num_cube_facets = IJK::compute_num_cube_facets(dimension);
  FIND_COMPONENT find_component(dimension);
  IJK::PROCEDURE_ERROR error
    ("ISODUAL_CUBE_TABLE_AMBIG_INFO::ComputeAmbiguityInformation");

  for (AMBIG_TABLE_INDEX ientry = 0; ientry < num_table_entries; ientry++) {
    is_ambiguous[ientry] = is_cube_ambiguous(ientry, find_component);

    compute_ambiguous_cube_facets
      (ientry, num_cube_facets, ambiguous_facet[ientry], 
       num_ambiguous_facets[ientry], find_component);
  }
}

// Set ambiguity table.
void ISODUAL_CUBE_TABLE_AMBIG_INFO::Set(const int dimension)
{
  const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
  const int BASE2 = 2;
  IJK::PROCEDURE_ERROR error
    ("ISODUAL_CUBE_TABLE_AMBIG_INFO::SetAmbiguityTable");

  this->dimension = dimension;

  u_long num_table_entries;
  IJK::int_power(BASE2, num_cube_vertices, num_table_entries, error);

  // Allocate memory.
  Alloc(num_table_entries);

  ComputeAmbiguityInformation();
}

//**************************************************
// AMBIGUITY ROUTINES
//**************************************************

bool IJKDUALTABLE::is_cube_ambiguous
(const AMBIG_TABLE_INDEX ientry, FIND_COMPONENT & find_component)
{
  int num_pos_components = 
    find_component.ComputeNumComponents(ientry, true);
  int num_neg_components = 
    find_component.ComputeNumComponents(ientry, false);

  if (num_pos_components > 1 || num_neg_components > 1)
    { return(true); }
  else
    { return(false); }
}

bool IJKDUALTABLE::is_cube_facet_ambiguous
(const AMBIG_TABLE_INDEX ientry, const FACET_INDEX & kf, 
 FIND_COMPONENT & find_component)
{
  int num_pos_components = 
    find_component.ComputeNumComponentsInFacet(ientry, kf, true);
  int num_neg_components = 
    find_component.ComputeNumComponentsInFacet(ientry, kf, false);

  if (num_pos_components > 1 || num_neg_components > 1)
    { return(true); }
  else
    { return(false); }
}


void IJKDUALTABLE::compute_ambiguous_cube_facets
(const AMBIG_TABLE_INDEX ientry, const FACET_INDEX num_facets,
 FACET_SET & facet_set, FACET_INDEX & num_ambiguous_facets,
 FIND_COMPONENT & find_component)
{
  facet_set = 0;
  num_ambiguous_facets = 0;
  FACET_SET mask = FACET_SET(1);

  for (FACET_INDEX kf = 0; kf < num_facets; kf++) {
    if (is_cube_facet_ambiguous(ientry, kf, find_component)) {
      num_ambiguous_facets++;
      facet_set = (facet_set | mask);
    }
    mask = (mask << FACET_INDEX(1));
  }
}



/// \file ijktable_ambig.cxx
/// Class containing a table of isosurface patches in a given polyhedron.
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

#include "ijktable_ambig.h"

using namespace IJK;
using namespace IJKTABLE;
using namespace std;

//**************************************************
// AMBIGUITY INFORMATION FOR ISOSURFACE TABLE
//**************************************************

void ISOSURFACE_TABLE_AMBIG_INFO::Init()
{
  num_table_entries = 0;
  is_ambiguous = NULL;
  num_ambiguous_facets = NULL;
  ambiguous_facet = NULL;
}


/// Destructor.
ISOSURFACE_TABLE_AMBIG_INFO::~ISOSURFACE_TABLE_AMBIG_INFO()
{
  FreeAll();
}

/// Free all memory.
void ISOSURFACE_TABLE_AMBIG_INFO::FreeAll()
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
void ISOSURFACE_TABLE_AMBIG_INFO::
Alloc(const long num_table_entries)
{
  if (this->num_table_entries > 0) { FreeAll(); };

  is_ambiguous = new bool[num_table_entries];
  num_ambiguous_facets = new FACET_INDEX[num_table_entries];
  ambiguous_facet = new FACET_SET[num_table_entries];

  this->num_table_entries = num_table_entries;
}


/// Compute ambiguity information.
void ISOSURFACE_TABLE_AMBIG_INFO::ComputeAmbiguityInformation
(const ISOSURFACE_TABLE_POLYHEDRON & poly)
{
  const int poly_numv = poly.NumVertices();
  const int BASE2 = 2;
  IJK::ARRAY<int> vertex_sign(poly_numv);
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE_AMBIG_INFO::ComputeAmbiguityInformation");

  for (AMBIG_TABLE_INDEX i = 0; i < num_table_entries; i++) {
    convert2base(i, BASE2, vertex_sign.Ptr(), poly_numv, error);
    is_ambiguous[i] = 
      is_poly_ambiguous(poly, vertex_sign.PtrConst());

    compute_ambiguous_facets
      (poly, vertex_sign.PtrConst(),
       ambiguous_facet[i], num_ambiguous_facets[i]);
  }
}

// Set ambiguity table.
void ISOSURFACE_TABLE_AMBIG_INFO::SetAmbiguityTable
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly)
{
  const int poly_numv = poly.NumVertices();
  const int BASE2 = 2;
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE_AMBIG_INFO::SetAmbiguityTable");

  u_long num_table_entries;
  IJK::int_power(BASE2, poly_numv, num_table_entries, error);

  // Allocate memory.
  Alloc(num_table_entries);

  ComputeAmbiguityInformation(poly);
}

// Set cube ambiguity table.
// @param use_lex_order If true, use lexicographic order on edges and facets.
void ISOSURFACE_TABLE_AMBIG_INFO::SetCubeAmbiguityTable
(const int dimension, const bool use_lex_order)
{
  IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON cube(dimension);

  if (use_lex_order)
    { cube.GenCube(dimension); }
  else
    { cube.GenCubeOrderA(dimension); }

  SetAmbiguityTable(cube);
}

//**************************************************
// AMBIGUITY ROUTINES
//**************************************************

/// Return true if isosurface topology is ambiguous
/// @param vertex_sign[i] = Sign of isosurface vertex i. (0 or 1);
bool IJKTABLE::is_poly_ambiguous
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly, 
 const int * vertex_sign)
{
  const int numv = poly.NumVertices();

  int num0 = compute_num_connected(poly, 0, vertex_sign);
  int num1 = 0;
  for (int i = 0; i < numv; i++) {
    if (vertex_sign[i] != vertex_sign[0]) {
      num1 = compute_num_connected(poly, i, vertex_sign);
      break;
    }
  }

  if (num0 + num1 != numv) { return(true); }
  else { return(false); }
}

/// Return true if facet jf is ambiguous.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
bool IJKTABLE::is_facet_ambiguous
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly,
 const IJKTABLE::FACET_INDEX ifacet, const int * vertex_sign)
{
  const int numv = poly.NumVertices();

  int num0 = 0;
  int num1 = 0;

  for (int iv = 0; iv < numv; iv++) {
    if (poly.IsVertexInFacet(ifacet, iv) && vertex_sign[iv] == 0) {
      num0 = compute_num_connected_in_facet(poly, ifacet, iv, vertex_sign);
      break;
    }
  }

  for (int iv = 0; iv < numv; iv++) {
    if (poly.IsVertexInFacet(ifacet, iv) && vertex_sign[iv] == 1) {
      num1 = compute_num_connected_in_facet(poly, ifacet, iv, vertex_sign);
      break;
    }
  }

  if (num0 + num1 != poly.NumFacetVertices(ifacet)) {
    return(true);
  }
  else { return(false); }
}

/// Return number of vertices connected by edges to iv with same sign as iv.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
int IJKTABLE::compute_num_connected
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly,
 const int iv, const int * vertex_sign)
{
  const int numv = poly.NumVertices();
  const int nume = poly.NumEdges();
  const int isign = vertex_sign[iv];
  IJK::ARRAY<bool> visited(numv);

  for (int j = 0; j < numv; j++) { visited[j] = false; };
  visited[iv] = true;

  bool found_next = false;
  do {
    found_next = false;
    for (int k = 0; k < nume; k++) {
      int iv0 = poly.EdgeEndpoint(k, 0);
      int iv1 = poly.EdgeEndpoint(k, 1);
      if (vertex_sign[iv0] == isign && vertex_sign[iv1] == isign) {
        if (visited[iv0] && !visited[iv1]) {
          visited[iv1] = true;
          found_next = true;
        }
        else if (!visited[iv0] && visited[iv1]) {
          visited[iv0] = true;
          found_next = true;
        }
      }
    }
  } while(found_next);
  
  int count = 0;
  for (int j = 0; j < numv; j++) {
    if (visited[j]) { count++; };
  }

  return(count);
}

/// Compute ambiguous facets.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
void IJKTABLE::compute_ambiguous_facets
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const int * vertex_sign,
 FACET_SET & facet_set, FACET_INDEX & num_ambiguous_facets)
{
  const int numf = poly.NumFacets();
 
  num_ambiguous_facets = 0;
  facet_set = 0;

  for (int ifacet = 0; ifacet < numf; ifacet++) {
    if (is_facet_ambiguous(poly, ifacet, vertex_sign)) {
      facet_set = facet_set | ((1L) << ifacet);
      num_ambiguous_facets++;
    }
  }
}

/// Return number of vertices connected by edges in facet jf to vertex iv with same sign as iv.
/// @param jf = Facet index.
/// @param iv = Vertex index.  Precondition: Vertex iv is in facet jf.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
int IJKTABLE::compute_num_connected_in_facet
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly,
 const FACET_INDEX jf, const int iv, const int * vertex_sign)
{
  const int numv = poly.NumVertices();
  const int nume = poly.NumEdges();
  const int isign = vertex_sign[iv];
  IJK::ARRAY<bool> visited(numv);

  for (int j = 0; j < numv; j++) { visited[j] = false; };
  visited[iv] = true;

  bool found_next = false;
  do {
    found_next = false;
    for (int k = 0; k < nume; k++) {
      int iv0 = poly.EdgeEndpoint(k, 0);
      int iv1 = poly.EdgeEndpoint(k, 1);
      if (vertex_sign[iv0] == isign && vertex_sign[iv1] == isign) {
        if (poly.IsVertexInFacet(jf, iv0) && poly.IsVertexInFacet(jf, iv1)) {
          if (visited[iv0] && !visited[iv1]) {
            visited[iv1] = true;
            found_next = true;
          }
          else if (!visited[iv0] && visited[iv1]) {
            visited[iv0] = true;
            found_next = true;
          }
        }
      }
    }
  } while(found_next);
  
  int count = 0;
  for (int j = 0; j < numv; j++) {
    if (visited[j]) { count++; };
  }

  return(count);
}

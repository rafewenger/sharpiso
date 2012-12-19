/// \file ijkdualtable.cxx
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

#include <cctype>
#include <limits>
#include <limits.h>
#include <sstream>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <set>
#include <vector>

#include "ijkcube.txx"
#include "ijkbits.txx"

#include "ijkdualtable.h"

using namespace IJK;
using namespace IJKDUALTABLE;
using namespace std;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif



// **************************************************
// ISODUAL_TABLE
// **************************************************

ISODUAL_TABLE::ISODUAL_TABLE_ENTRY::ISODUAL_TABLE_ENTRY()
  // constructor
{
  num_vertices = 0;
  incident_isovertex = NULL;
  is_bipolar = NULL;
}

ISODUAL_TABLE::ISODUAL_TABLE_ENTRY::~ISODUAL_TABLE_ENTRY()
  // destructor
{
  FreeAll();
}

void ISODUAL_TABLE::ISODUAL_TABLE_ENTRY::Allocate
(const int num_poly_edges)
{
  if (incident_isovertex != NULL) {
    delete [] incident_isovertex; 
    incident_isovertex = NULL;
  }

  if (is_bipolar != NULL) {
    delete [] is_bipolar;
    is_bipolar = NULL;
  }

  incident_isovertex = new ISODUAL_VERTEX_INDEX[num_poly_edges];
  is_bipolar = new bool[num_poly_edges];
}

bool ISODUAL_TABLE::ISODUAL_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_vertices < 0) {
    error_msg.AddMessage
      ("Error.  Dual isosurface table entry contains negative number of isosurface vertices.");
    return(false);
  }

  if (incident_isovertex == NULL) {
    error_msg.AddMessage("Memory for incident isosurface vertices not allocated.");
    return(false);
  }

  if (is_bipolar == NULL) {
    error_msg.AddMessage("Memory for bipolar edge flags not allocated.");
    return(false);
  }

  return(true);
}

// Free all memory.
void ISODUAL_TABLE::ISODUAL_TABLE_ENTRY::FreeAll()
{
  delete [] incident_isovertex;
  incident_isovertex = NULL;
  delete [] is_bipolar;
  is_bipolar = NULL;
  num_vertices = 0;
}

// default constructor. dimension = 3
ISODUAL_TABLE::ISODUAL_TABLE()
{
  Init(3);
}

// Constructor.
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
ISODUAL_TABLE::ISODUAL_TABLE(const int d)
{
  Init(d);
}

/// Initialize
void ISODUAL_TABLE::Init(const int dimension)
{
  const char * procname = "ISODUAL_TABLE::Init";

  max_num_vertices = 20;
  // Note: Even tables for polytopes of this size are probably impossible 
  //   to compute/store

  num_table_entries = 0;
  entry = NULL;
  is_table_allocated = false;

  SetDimension(dimension);
}

// destructor
ISODUAL_TABLE::~ISODUAL_TABLE()
{
  FreeAll();
}

// Set dimension.
void ISODUAL_TABLE::SetDimension(const int dimension)
{
  const char * procname = "ISODUAL_TABLE::SetDimension";

  if (IsTableAllocated()) {
    throw PROCEDURE_ERROR
      (procname, "Dual isosurface table already allocated.");
  }

  this->dimension = dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polytope dimension.");
}

// Set number of polytope vertices.
void ISODUAL_TABLE::SetNumPolyVertices(const int num_vertices)
{
  const char * procname = "ISODUAL_TABLE::SetNumVertices";

  if (IsTableAllocated()) {
    throw PROCEDURE_ERROR
      (procname, "Dual isosurface table already allocated.");
  }

  num_poly_vertices = num_vertices;
}

// Set number of polytope edges
void ISODUAL_TABLE::SetNumPolyEdges(const int num_edges)
{
  const char * procname = "ISODUAL_TABLE::SetNumEdges";

  if (IsTableAllocated()) {
    throw PROCEDURE_ERROR
      (procname, "Dual isosurface table already allocated.");
  }

  num_poly_edges = num_edges;
}

// Allocate table
void ISODUAL_TABLE::SetNumTableEntries(const int num_table_entries)
{
  const char * procname = "ISODUAL_TABLE::SetNumTableEntries";

  if (entry != NULL) delete [] entry;
  entry = NULL;
  this->num_table_entries = 0;
  is_table_allocated = false;

  entry = new ISODUAL_TABLE_ENTRY[num_table_entries];
  if (entry == NULL && num_table_entries > 0)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for dual isosurface table.");

  for (int i = 0; i < num_table_entries; i++) {
    entry[i].Allocate(NumPolyVertices());
  }

  this->num_table_entries = num_table_entries;
  is_table_allocated = true;
}


// check dimension
bool ISODUAL_TABLE::CheckDimension(const int d) const
{
  if (d < 1)
    return(false);
  else
    return(true);
}

// check table
bool ISODUAL_TABLE::CheckTable(ERROR & error_msg) const
{
  if (NumPolyVertices() >= LONG_BIT) {
    error_msg.AddMessage("Too many polytope vertices");
    return(false);
  }

  if (NumPolyVertices() > MaxNumVertices()) {
    error_msg.AddMessage("Too many polytope vertices");
    return(false);
  }

  if (NumPolyVertices() < 1) {
    error_msg.AddMessage("Polytope must have at least one vertex.");
    return(false);
  }

  if (entry == NULL) {
    error_msg.AddMessage("Memory for dual isosurface table not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
                           it, ".");
      return(false);
    }

  return(true);
}

bool ISODUAL_TABLE::Check(ERROR & error_msg) const
{
  if (!CheckTable(error_msg)) return(false);
  return(true);
}

void ISODUAL_TABLE::FreeAll()
  // free all memory
{
  if (entry != NULL) {
    for (int i = 0; i < num_table_entries; i++)
      { entry[i].FreeAll(); }
    delete [] entry;
    entry = NULL;
  };
  num_table_entries = 0;
  is_table_allocated = false;
}

// **************************************************
// ISODUAL CUBE TABLE
// **************************************************

// Constructor.
ISODUAL_CUBE_TABLE::ISODUAL_CUBE_TABLE(const int dimension)
{
  Create(dimension);
}

// Constructor.
ISODUAL_CUBE_TABLE::ISODUAL_CUBE_TABLE
(const int dimension, const bool flag_opposite_vertices)
{
  Create(dimension, flag_opposite_vertices);
}

// Constructor.
ISODUAL_CUBE_TABLE::ISODUAL_CUBE_TABLE
(const int dimension, const bool flag_separate_neg,
 const bool flag_opposite_vertices)
{
  Create(dimension, flag_separate_neg, flag_opposite_vertices);
}

// Set dimension.
void ISODUAL_CUBE_TABLE::SetDimension(const int dimension)
{
  ISODUAL_TABLE::SetDimension(dimension);

  int num_vertices = compute_num_cube_vertices(dimension);
  this->SetNumPolyVertices(num_vertices);

  int num_edges = compute_num_cube_edges(dimension);
  this->SetNumPolyEdges(num_edges);
}

// Create table.
void ISODUAL_CUBE_TABLE::Create
(const int dimension, const bool flag_separate_neg,
 const bool flag_separate_opposite)
{
  SetDimension(dimension);

  TABLE_INDEX n = calculate_num_entries(NumPolyVertices(), 2);
  ISODUAL_TABLE::SetNumTableEntries(n);
  CreateTableEntries(flag_separate_neg, flag_separate_opposite);
}

// Create table.
void ISODUAL_CUBE_TABLE::Create
(const int dimension, const bool flag_separate_opposite)
{
  Create(dimension, true, flag_separate_opposite);
}

// Create table.
void ISODUAL_CUBE_TABLE::Create(const int dimension)
{
  Create(dimension, true);
}

// Create table entries.
// @param flag_separate_neg  If true, separate negative vertices.
// @param flag_separate_opposite If true, always separate two diagonally 
//        opposite negative or positive vertices 
void ISODUAL_CUBE_TABLE::CreateTableEntries
(const bool flag_separate_neg, const bool flag_separate_opposite)
{

  this->flag_separate_neg = flag_separate_neg;
  this->flag_always_separate_opposite = flag_separate_opposite;

  FIND_COMPONENT find_component(Dimension());
  IJK::CUBE_FACE_INFO<int, int, int> cube(Dimension());

  bool flag_separate_pos = (!flag_separate_neg);
  for (TABLE_INDEX ientry = 0; ientry < NumTableEntries(); ientry++) {

    find_component.ClearAll();
    find_component.SetVertexFlags(ientry);

    if (flag_separate_opposite) {
      if (flag_separate_neg) {
        if (is_two_opposite_ones(ientry, NumPolyVertices()))
          { find_component.NegateVertexFlags(); };
      }
      else {
        if (is_two_opposite_zeros(ientry, NumPolyVertices()))
          { find_component.NegateVertexFlags(); };
      }
    }


    int num_components(0);
    for (int i = 0; i < NumPolyVertices(); i++) {
      if (find_component.Component(i) == 0) {

        if (find_component.VertexFlag(i) == flag_separate_pos) {
          num_components++;
          find_component.Search(i, num_components);
        }
      }
    }

    entry[ientry].num_vertices = num_components;

    for (int ie = 0; ie < cube.NumEdges(); ie++) {
      int iv0 = cube.EdgeEndpoint(ie, 0);
      int iv1 = cube.EdgeEndpoint(ie, 1);

      if (find_component.VertexFlag(iv0) == find_component.VertexFlag(iv1)) {
        entry[ientry].is_bipolar[ie] = false;
        entry[ientry].incident_isovertex[ie] = 0;
      }
      else {
        entry[ientry].is_bipolar[ie] = true;
        int icomp = find_component.Component(iv0);
        if (find_component.VertexFlag(iv1) == flag_separate_pos) {
          // Vertex iv1 is negative.
          icomp = find_component.Component(iv1);
        }
        entry[ientry].incident_isovertex[ie] = icomp-1;
      }
    }

  }

}


// **************************************************
// CLASS FIND_COMPONENT
// **************************************************

FIND_COMPONENT::FIND_COMPONENT(const int dimension)
{
  this->dimension = dimension;
  num_cube_vertices = 1L << Dimension();

  vertex_flag = new bool[num_cube_vertices];
  component = new int[num_cube_vertices];
}

FIND_COMPONENT::~FIND_COMPONENT()
{
  if (component != NULL) { delete [] component; }
  component = NULL;

  if (vertex_flag != NULL) { delete [] vertex_flag; }
  vertex_flag = NULL;
}

void FIND_COMPONENT::ClearAll()
{
  for (int i = 0; i < NumCubeVertices(); i++) {
    vertex_flag[i] = false;
    component[i] = 0;
  }

}

void FIND_COMPONENT::SetVertexFlags(const TABLE_INDEX ival)
{
  convert2bool(ival, vertex_flag, num_cube_vertices);
}

void FIND_COMPONENT::NegateVertexFlags()
{
  for (int i = 0; i < num_cube_vertices; i++) 
    { vertex_flag[i] = (!vertex_flag[i]); }
}


void FIND_COMPONENT::Search(const int i, const int icomp)
{
  std::vector<int> vlist;

  vlist.reserve(num_cube_vertices);

  if (icomp == 0) {
    IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
    error.AddMessage("Programming error. Component number cannot be zero.");
    throw error;
  }

  bool flag = vertex_flag[i];
  vlist.push_back(i);
  component[i] = icomp;
  while(vlist.size() > 0) {
    int j = vlist.back();
    vlist.pop_back();
    for (int d = 0; d < dimension; d++) {
      int j2;
      compute_cube_vertex_neighbor(j, d, j2);
      if (vertex_flag[j2] == flag && component[j2] == 0) {
        vlist.push_back(j2);
        component[j2] = icomp;
      }
    }
  }
}


void FIND_COMPONENT::SearchFacet(const int kf, const int i, const int icomp)
{
  std::vector<int> vlist;

  vlist.reserve(num_cube_vertices);

  if (icomp == 0) {
    IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
    error.AddMessage("Programming error. Component number cannot be zero.");
    throw error;
  }

  if (!IJK::cube_facet_contains(dimension, kf, i)) {
    IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
    error.AddMessage("Programming error. Facet ", kf,
                     " does not contain vertex ", i, ".");
    throw error;
  }

  bool flag = vertex_flag[i];
  vlist.push_back(i);
  component[i] = icomp;
  while(vlist.size() > 0) {
    int j = vlist.back();
    vlist.pop_back();
    for (int d = 0; d < dimension; d++) {
      int j2;
      compute_cube_vertex_neighbor(j, d, j2);
      if (vertex_flag[j2] == flag && component[j2] == 0 &&
          IJK::cube_facet_contains(dimension, kf, i)) {
        vlist.push_back(j2);
        component[j2] = icomp;
      }
    }
  }
}

// Compute number of components of vertices.
int FIND_COMPONENT::ComputeNumComponents
(const int ientry, const bool flag_positive)
{
  ClearAll();
  SetVertexFlags(ientry);
  if (!flag_positive) { NegateVertexFlags(); }

  int num_components(0);
  for (int i = 0; i < NumCubeVertices(); i++) {
    if (Component(i) == 0) {

      if (VertexFlag(i)) {
        num_components++;
        Search(i, num_components);
      }
    }
  }

  return(num_components);
}

// Compute number of components of vertices in facet.
int FIND_COMPONENT::ComputeNumComponentsInFacet
(const int ientry, const int kf, const bool flag_positive)
{
  ClearAll();
  SetVertexFlags(ientry);
  if (!flag_positive) { NegateVertexFlags(); }

  int num_components(0);
  for (int i = 0; i < NumCubeVertices(); i++) {
    if (Component(i) == 0 &&
        IJK::cube_facet_contains(dimension, kf, i)) {

      if (VertexFlag(i)) {
        num_components++;
        SearchFacet(kf, i, num_components);
      }
    }
  }

  return(num_components);
}

// **************************************************
// UTILITY FUNCTIONS
// **************************************************

// calculate number of entries required in ISODUAL_TABLE
unsigned long IJKDUALTABLE::calculate_num_entries
(const int num_vert, const int num_colors)
{
  const char * procname = "calculate_num_entries";

  unsigned long num_table_entries = 0;

  if (num_colors < 1)
    throw PROCEDURE_ERROR(procname, "Number of colors must be positive.");

  const unsigned long max2 = ULONG_MAX/num_colors;

  num_table_entries = 1;
  for (int iv = 0; iv < num_vert; iv++) {
    if (num_table_entries > max2)
      throw PROCEDURE_ERROR(procname, "Number of entries is too large.");

    num_table_entries = num_table_entries * num_colors;
  };

  return(num_table_entries);
}
// convert integer bits to boolean flags
void IJKDUALTABLE::convert2bool
(const TABLE_INDEX ival, bool * flag, const unsigned int num_flags)
{
  const int base = 2;

  unsigned long jval = ival;
  for (int i = 0; i < int(num_flags); i++) {
    flag[i] = bool(jval % base);
    jval = jval/base;
  };

  if (jval != 0) {
    PROCEDURE_ERROR error("convert2bool");
    error.AddMessage("Error converting ", ival, " to base 2.");
    error.AddMessage("Output has more than ", num_flags, " digits.");

    throw error;
  };
}

/// Return true if ival represents two opposite ones.
bool IJKDUALTABLE::is_two_opposite_ones
(const TABLE_INDEX ival, const int num_bits)
{
  int num_zeros, num_ones;
  IJK::count_bits(ival, num_bits, num_zeros, num_ones);
 
  if (equals_reverse_bits(ival, num_bits) && (num_ones == 2))
    { return(true); }
  else
    { return(false); }
}

/// Return true if ival represents two opposite zeros.
bool IJKDUALTABLE::is_two_opposite_zeros
(const TABLE_INDEX ival, const int num_bits)
{
  int num_zeros, num_ones;
  IJK::count_bits(ival, num_bits, num_zeros, num_ones);

   if (equals_reverse_bits(ival, num_bits) && (num_zeros == 2))
    { return(true); }
  else
    { return(false); }
}

/// \file isodual3D_ambig.cxx
/// Routines for handling ambiguous configurations.
/// Version 0.0.1

/*
IJK: Isosurface Jeneration Kode
Copyright (C) 2012-2013 Rephael Wenger

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


#include "isodual3D_ambig.h"
#include "sharpiso_feature.h"

#include "ijkbits.txx"
#include "ijkisopoly.txx"

#include <queue>

using namespace IJK;
using namespace ISODUAL3D;
using namespace std;


// **************************************************
// DETERMINE AMBIGUOUS FACETS
// **************************************************

/// Return true if facet is ambiguous.
bool ISODUAL3D::is_grid_facet_ambiguous
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX facet_v0,
 const NUM_TYPE orth_dir,
 const SCALAR_TYPE isovalue)
{
  const NUM_TYPE dim1 = (orth_dir+1)%DIM3;
  const NUM_TYPE dim2 = (orth_dir+2)%DIM3;

  const VERTEX_INDEX v1 = scalar_grid.NextVertex(facet_v0, dim1);

  const SCALAR_TYPE s0 = scalar_grid.Scalar(facet_v0);
  const SCALAR_TYPE s1 = scalar_grid.Scalar(v1);

  if (s0 < isovalue) {
    if (s1 < isovalue) { return(false); }

    const VERTEX_INDEX v2 = scalar_grid.NextVertex(facet_v0, dim2);
    const SCALAR_TYPE s2 = scalar_grid.Scalar(v2);
    if (s2 < isovalue) { return(false); }

    const VERTEX_INDEX v3 = scalar_grid.NextVertex(v1, dim2);
    const SCALAR_TYPE s3 = scalar_grid.Scalar(v3);
    if (s3 >= isovalue) { return(false); }

    return(true);
  }
  else {
    if (s1 >= isovalue) { return(false); }

    const VERTEX_INDEX v2 = scalar_grid.NextVertex(facet_v0, dim2);
    const SCALAR_TYPE s2 = scalar_grid.Scalar(v2);
    if (s2 >= isovalue) { return(false); }

    const VERTEX_INDEX v3 = scalar_grid.NextVertex(v1, dim2);
    const SCALAR_TYPE s3 = scalar_grid.Scalar(v3);
    if (s3 < isovalue) { return(false); }

    return(true);
  }
}

/// Decide and return ambiguity status of facet ifacet.
AMBIGUITY_STATUS ISODUAL3D::decide_ambiguous_facet
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const AMBIG_TABLE & ambig_table,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX facet_v0,
 const int facet_orth_dir,
 const ISODUAL_PARAM & isodual_param)
{
  NUM_TYPE sharp_vertex_location;
  static COORD_TYPE coord[DIM3];

  if (!is_grid_facet_ambiguous
      (scalar_grid, facet_v0, facet_orth_dir, isovalue)) 
    { return(NOT_AMBIGUOUS); }

  scalar_grid.ComputeCoord(facet_v0, coord);
  if (coord[facet_orth_dir] <= 0 || 
      coord[facet_orth_dir]+1 >= scalar_grid.AxisSize(facet_orth_dir))
    // boundary facet
    { return(UNDECIDED_AMBIGUITY); }

  svd_compute_sharp_vertex_near_facet
    (scalar_grid, gradient_grid, facet_v0, facet_orth_dir, isovalue,
     isodual_param, sharp_vertex_location);

  if (sharp_vertex_location == 0) 
    { return(UNDECIDED_AMBIGUITY); }

  VERTEX_INDEX icube;
  if (sharp_vertex_location > 0) {
    icube = facet_v0;
  }
  else {
    icube = scalar_grid.PrevVertex(facet_v0, facet_orth_dir);
  }

  IJKTABLE::AMBIG_TABLE_INDEX it;
  compute_isotable_index
    (scalar_grid.ScalarPtrConst(), isovalue, icube,
     scalar_grid.CubeVertexIncrement(), NUM_CUBE_VERTICES3D, it);

  if (ambig_table.NumPosComponents(it) == 1) 
    { return(SEPARATE_POS); }
  else if (ambig_table.NumNegComponents(it) == 1) 
    { return(SEPARATE_NEG); }
  else 
    { return(UNDECIDED_AMBIGUITY); }

}

namespace {

  // Return true if cube has facet with a given ambig_status.
  bool cube_has_facet_with_ambig_status
  (const SHARPISO_GRID & grid,
   const VERTEX_INDEX icube,
   const AMBIGUITY_STATUS ambig_status,
   const std::vector<AMBIGUITY_TYPE> & facet_ambig_status)
  {
    for (int orth_dir = 0; orth_dir < DIM3; orth_dir++) {
      if (facet_ambig_status[icube*DIM3+orth_dir] == ambig_status) 
        { return(true); }

      VERTEX_INDEX iv1 = grid.NextVertex(icube, orth_dir);

      if (facet_ambig_status[iv1*DIM3+orth_dir] == ambig_status) 
        { return(true); }
    }

    return(false);
  }

  // Return true if some facet of two cubes have a given ambig_status.
  // @pre facet is internal grid facet.
  bool two_cubes_have_facet_with_ambig_status
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX facet_v0,
   const int facet_orth_dir,
   const AMBIGUITY_STATUS ambig_status,
   const std::vector<AMBIGUITY_TYPE> & facet_ambig_status)
  {
    const VERTEX_INDEX icube1 = facet_v0;

    if (cube_has_facet_with_ambig_status
        (scalar_grid, icube1, ambig_status, facet_ambig_status)) 
      { return(true); }

    const VERTEX_INDEX icube0 = 
      scalar_grid.PrevVertex(icube1, facet_orth_dir);

    if (cube_has_facet_with_ambig_status
        (scalar_grid, icube0, ambig_status, facet_ambig_status)) 
      { return(true); }

    return(false);
  }

  void set_cube_facets_undecided_ambig
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX icube,
   const AMBIGUITY_TYPE ambig_status,
   std::vector<AMBIGUITY_TYPE> & facet_ambig_status,
   std::queue<VERTEX_INDEX> facet_list)
  {
    for (int orth_dir = 0; orth_dir < DIM3; orth_dir++) {
      VERTEX_INDEX ifacet = icube*DIM3+orth_dir;
      if (facet_ambig_status[ifacet] == UNDECIDED_AMBIGUITY) {
        facet_ambig_status[ifacet] = ambig_status; 
        facet_list.push(ifacet);
      }

      VERTEX_INDEX iv1 = scalar_grid.NextVertex(icube, orth_dir);
      ifacet = iv1*DIM3+orth_dir;

      if (facet_ambig_status[ifacet] == UNDECIDED_AMBIGUITY) {
        facet_ambig_status[ifacet] = ambig_status; 
        facet_list.push(ifacet);
      }
    }
  }

  void set_two_cube_facets_undecided_ambig
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX facet_v0,
   const int facet_orth_dir,
   const AMBIGUITY_STATUS ambig_status,
   std::vector<AMBIGUITY_TYPE> & facet_ambig_status,
   std::queue<VERTEX_INDEX> & facet_list)
  {
    const VERTEX_INDEX icube1 = facet_v0;

    set_cube_facets_undecided_ambig
      (scalar_grid, icube1, ambig_status, facet_ambig_status, facet_list);

    const VERTEX_INDEX icube0 = 
      scalar_grid.PrevVertex(icube1, facet_orth_dir);

    set_cube_facets_undecided_ambig
      (scalar_grid, icube0, ambig_status, facet_ambig_status, facet_list);
  }

};


/// Propagate SEPARATE_POS and SEPARATE_NEG across cube facets.
/// @param facet_ambig_status[] Array containing facet ambiguity.
///        facet_ambig_status[iv*DIM3+d] is the ambiguity of the facet
///          containing primary vertex v and orthodonal direction d.
///        Note: When iv is on the upper-rightmost grid boundary,
///          the facet iv*DIM3+d may not actually be contained in the grid.
void ISODUAL3D::propagate_sep
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const AMBIG_TABLE & ambig_table,
 std::vector<AMBIGUITY_TYPE> & facet_ambig_status)
{
  std::queue<VERTEX_INDEX> facet_list;
  static COORD_TYPE coord[DIM3];

  for (int i = 0; i < facet_ambig_status.size(); i++) {
    if (facet_ambig_status[i] == SEPARATE_POS ||
        facet_ambig_status[i] == SEPARATE_NEG) 
      { facet_list.push(i); }
  }

  while (facet_list.size() != 0) {
    VERTEX_INDEX ifacet = facet_list.front();
    facet_list.pop();

    VERTEX_INDEX facet_v0 = ifacet/DIM3;
    VERTEX_INDEX facet_orth_dir = ifacet%DIM3;

    scalar_grid.ComputeCoord(facet_v0, coord);
    if (coord[facet_orth_dir] != 0 &&
        coord[facet_orth_dir]+1 != scalar_grid.AxisSize(facet_orth_dir)) {
      // Internal facet.
      
      AMBIGUITY_STATUS fambig = AMBIGUITY_STATUS(facet_ambig_status[ifacet]);
      AMBIGUITY_STATUS fambig_complement = SEPARATE_POS;

      if (fambig == SEPARATE_POS) 
        { fambig_complement = SEPARATE_NEG; }

      if (!two_cubes_have_facet_with_ambig_status
          (scalar_grid, facet_v0, facet_orth_dir, fambig_complement, 
           facet_ambig_status)) {

        set_two_cube_facets_undecided_ambig
          (scalar_grid, facet_v0, facet_orth_dir, fambig, 
           facet_ambig_status, facet_list);
      }
    }
  }
}


// Set ambiguity of cubes in cube_list[].
void ISODUAL3D::set_cube_ambiguity
(const SHARPISO_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & cube_list,
 const std::vector<AMBIGUITY_TYPE> & facet_ambig_status,
 std::vector<AMBIGUITY_TYPE> & cube_ambig)
{
  cube_ambig.resize(cube_list.size());

  for (int i = 0; i < cube_list.size(); i++) {
    bool flag_pos = 
      cube_has_facet_with_ambig_status
      (grid, cube_list[i], SEPARATE_POS, facet_ambig_status);
    bool flag_neg = 
      cube_has_facet_with_ambig_status
      (grid, cube_list[i], SEPARATE_NEG, facet_ambig_status);

    if (flag_pos) {
      if (flag_neg) 
        { cube_ambig[i] = CONFLICTING_SEPARATION; }
      else 
        { cube_ambig[i] = SEPARATE_POS; }
    }
    else if (flag_neg) 
      { cube_ambig[i] = SEPARATE_NEG; }
    else 
      { cube_ambig[i] = NOT_AMBIGUOUS; }

  }
}

/// Set ambiguity of cubes in cube_list[].
void ISODUAL3D::set_cube_ambiguity    
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & cube_list,
 const ISODUAL_PARAM & isodual_param,
 std::vector<AMBIGUITY_TYPE> & cube_ambig)
{
  const VERTEX_INDEX num_gridv = scalar_grid.NumVertices();
  const VERTEX_INDEX num_gridf = DIM3 * num_gridv;
  std::vector<AMBIGUITY_TYPE> 
    facet_ambig_status(num_gridf, AMBIGUITY_TYPE(AMBIGUITY_NOT_SET));
  AMBIG_TABLE ambig_table;

  ambig_table.SetCubeAmbiguityTable();

  // Set facet_ambig_status[].
  for (VERTEX_INDEX i = 0; i < cube_list.size(); i++) {
    VERTEX_INDEX iv0 = cube_list[i];

    for (int facet_orth_dir = 0; facet_orth_dir < DIM3; facet_orth_dir++) {
      VERTEX_INDEX ifacet0 = iv0*DIM3 + facet_orth_dir;
      if (facet_ambig_status[ifacet0] == AMBIGUITY_TYPE(AMBIGUITY_NOT_SET)) {
        AMBIGUITY_STATUS ambig_status = 
          decide_ambiguous_facet(scalar_grid, gradient_grid, ambig_table,
                                 isovalue, iv0, facet_orth_dir, isodual_param);
        facet_ambig_status[ifacet0] = AMBIGUITY_TYPE(ambig_status);
      }

      VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, facet_orth_dir);
      VERTEX_INDEX ifacet1 = iv1*DIM3 + facet_orth_dir;
      if (facet_ambig_status[ifacet1] == AMBIGUITY_TYPE(AMBIGUITY_NOT_SET)) {
        AMBIGUITY_STATUS ambig_status = 
          decide_ambiguous_facet(scalar_grid, gradient_grid, ambig_table,
                                 isovalue, iv1, facet_orth_dir, isodual_param);
        facet_ambig_status[ifacet1] = AMBIGUITY_TYPE(ambig_status);
      }
    }
  }

  propagate_sep(scalar_grid, ambig_table, facet_ambig_status);

  set_cube_ambiguity(scalar_grid, cube_list, facet_ambig_status, cube_ambig);
}

// Set ambiguity information.
void ISODUAL3D::set_ambiguity_info
(const std::vector<AMBIGUITY_TYPE> & cube_ambig,
 SHARPISO_INFO & sharpiso_info)
{
  for (ISO_VERTEX_INDEX i = 0; i < cube_ambig.size(); i++) {

    if (cube_ambig[i] == NOT_AMBIGUOUS) {
      sharpiso_info.num_cube_not_ambiguous++;
    }
    else if (cube_ambig[i] == SEPARATE_POS) {
      sharpiso_info.num_cube_separate_pos++;
    }
    else if (cube_ambig[i] == SEPARATE_NEG) {
      sharpiso_info.num_cube_separate_neg++;
    }
    else if (cube_ambig[i] == CONFLICTING_SEPARATION) {
      sharpiso_info.num_cube_unresolved_ambiguity++;
    }

  }
}


// **************************************************
// AMBIG_TABLE
// **************************************************

void ISODUAL3D::AMBIG_TABLE::Init()
{
  num_pos_components = NULL;
  num_neg_components = NULL;
}

/// Destructor.
ISODUAL3D::AMBIG_TABLE::~AMBIG_TABLE()
{
  FreeAll();
}

/// Free all memory.
void ISODUAL3D::AMBIG_TABLE::FreeAll()
{
  delete [] num_pos_components;
  num_pos_components = NULL;
  delete [] num_neg_components;
  num_neg_components = NULL;

  ISOSURFACE_TABLE_AMBIG_INFO::FreeAll();
}

/// Allocate memory.
/// Free any previously allocated memory.
void ISODUAL3D::AMBIG_TABLE::
Alloc(const long num_table_entries)
{
  if (this->num_table_entries > 0) { FreeAll(); };

  ISOSURFACE_TABLE_AMBIG_INFO::Alloc(num_table_entries);

  num_pos_components = new NUM_COMPONENTS_TYPE[num_table_entries];
  num_neg_components = new NUM_COMPONENTS_TYPE[num_table_entries];
}

// Set cube ambiguity table.
void ISODUAL3D::AMBIG_TABLE::SetCubeAmbiguityTable()
{
  IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON cube(DIM3);
  const int BASE2 = 2;
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE_AMBIG_INFO::SetCubeAmbiguityTable");

  cube.GenCubeOrderA(DIM3);

  u_long num_table_entries;
  IJK::int_power(BASE2, NUM_CUBE_VERTICES3D, num_table_entries, error);

  // Allocate memory.
  Alloc(num_table_entries);

  ComputeAmbiguityInformation(cube);
  ComputeNumConnectedComponents(cube);
}

// Compute number of connected components.
void ISODUAL3D::AMBIG_TABLE::
ComputeNumConnectedComponents(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly)
{
  const int poly_numv = poly.NumVertices();
  const int BASE2 = 2;
  IJK::ARRAY<int> vertex_sign(poly_numv);
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE_AMBIG_INFO::ComputeNumConnectedComponents");

  for (IJKTABLE::AMBIG_TABLE_INDEX i = 0; i < NumTableEntries(); i++) {
    convert2base(i, BASE2, vertex_sign.Ptr(), poly_numv, error);
    num_pos_components[i] = 
      compute_num_connected_components(poly, 1, vertex_sign.PtrConst());
    num_neg_components[i] = 
      compute_num_connected_components(poly, 0, vertex_sign.PtrConst());
  }
}


/// Return number of connected components of vertices with sign isign.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
int ISODUAL3D::compute_num_connected_components
(const IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON & poly,
 const int isign, const int * vertex_sign)
{
  const int numv = poly.NumVertices();
  const int nume = poly.NumEdges();
  IJK::ARRAY<bool> visited(numv);

  for (int j = 0; j < numv; j++) { visited[j] = false; };

  int num_components = 0;
  for (int j = 0; j < numv; j++) {

    if (visited[j] == false && vertex_sign[j] == isign) {

      num_components++;
      visited[j] = true;

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
    }
  }
  
  return(num_components);
}




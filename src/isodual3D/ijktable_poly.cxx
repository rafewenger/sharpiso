/// \file ijktable_poly.cxx
/// Class containing isosurface table polyhedron.
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

#include "ijktable_poly.h"

#include "ijkcube.txx"

using namespace IJK;
using namespace IJKTABLE;
using namespace std;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif

//**************************************************
// ISOSURFACE_TABLE_POLYHEDRON
//**************************************************

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON(const int d)
  // constructor
{
  Init();
  dimension = d;
}

ISOSURFACE_TABLE_POLYHEDRON::~ISOSURFACE_TABLE_POLYHEDRON()
  // destructor
{
  FreeAll();
}

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON
(const ISOSURFACE_TABLE_POLYHEDRON & init)
  // copy
{
  Init();

  *this = init;
}

const ISOSURFACE_TABLE_POLYHEDRON & ISOSURFACE_TABLE_POLYHEDRON::operator = 
(const ISOSURFACE_TABLE_POLYHEDRON & right)  
  // assign 
{
  if (&right != this) {         // avoid self-assignment
    FreeAll();
    dimension = right.Dimension();
    SetSize(right.NumVertices(), right.NumEdges(), right.NumFacets());

    // copy vertices
    for (int iv = 0; iv < NumVertices(); iv++)
      for (int ic = 0; ic < Dimension(); ic++)
        SetVertexCoord(iv, ic, right.VertexCoord(iv, ic));

    // copy edges
    for (int ie = 0; ie < NumEdges(); ie++)
      SetEdge(ie, right.EdgeEndpoint(ie, 0), right.EdgeEndpoint(ie, 1));

    /// **** NOTE: NEED TO MODIFY ***

    // copy facets
    for (int jf = 0; jf < NumFacets(); jf++) {
      facet[jf] = right.Facet(jf);
      int num_fv = right.NumFacetVertices(jf);
      SetNumFacetVertices(jf, num_fv);

      for (int k = 0; k < num_fv; k++) {
        int iv = right.FacetVertex(jf, k);
        SetFacetVertex(jf, k, iv);
      }
    }

  };

  return *this;
}

/// Free facet arrays.
void ISOSURFACE_TABLE_POLYHEDRON::FreeFacets()
{
  if (facet_vertex_list != NULL) {
    for (int jf = 0; jf < num_facets; jf++) {
      delete[] facet_vertex_list[jf];
      facet_vertex_list[jf] = NULL;
    }
  }
  delete [] facet_vertex_list;
  facet_vertex_list = NULL;

  delete [] num_facet_vertices;
  num_facet_vertices = NULL;
  delete [] facet;
  facet = NULL;
  num_facets = 0;
}

/// Free all memory.
void ISOSURFACE_TABLE_POLYHEDRON::FreeAll()
{
  FreeFacets();

  num_vertices = 0;
  num_edges = 0;
  delete [] vertex_coord;
  vertex_coord = NULL;
  delete [] edge_endpoint;
  edge_endpoint = NULL;
}

void ISOSURFACE_TABLE_POLYHEDRON::Init()
  // initialize
{
  dimension = 0;
  num_vertices = 0;
  num_edges = 0;
  num_facets = 0;
  vertex_coord = NULL;
  edge_endpoint = NULL;
  facet = NULL;
  num_facet_vertices = NULL;
  facet_vertex_list = NULL;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetDimension(const int d)
  // set polyhedron dimension
{
  FreeAll();
  num_vertices = num_edges = 0;

  dimension = d;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumVertices(const int numv)
  // set number of vertices
  // Must be called before setting polyhedron vertices
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumVertices";

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polyhedron dimension.");

  FreeAll();
  num_vertices = num_edges = 0;

  if (numv == 0)
    throw PROCEDURE_ERROR(procname, "Number of vertices must be non-zero.");

  // Note that even if numv <= LONG_BIT, there may not be enough 
  //   memory to store the isosurface table.
  if (numv > LONG_BIT)
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron vertices is too large.");

  num_vertices = numv;
  vertex_coord = new int[num_vertices*Dimension()];
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumEdges(const int nume)
  // set number of edges
  // Must be called before setting polyhedron edges
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumEdges";

  delete [] edge_endpoint;
  edge_endpoint = NULL;
  num_edges = 0;

  if (!CheckDimension()) 
    throw PROCEDURE_ERROR(procname, "Illegal dimension.");

  if (NumVertices() == 0)
    throw PROCEDURE_ERROR
      (procname, "Number of vertices must be set before number of edges.");

  if (nume < 1)
    throw PROCEDURE_ERROR(procname, "Number of edges must be non-zero.");

  if (nume > std::numeric_limits<EDGE_INDEX>::max())
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron edges is too large.");

  num_edges = nume;
  edge_endpoint = new int[num_edges*2];
}


void ISOSURFACE_TABLE_POLYHEDRON::SetNumFacets(const int numf)
  // set number of facets
  // Must be called before setting polyhedron facets
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetNumFacets";

  FreeFacets();

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal dimension.");

  if (NumVertices() == 0)
    throw PROCEDURE_ERROR
      (procname, "Number of vertices must be set before number of facets.");

  if (numf < 1)
    throw PROCEDURE_ERROR(procname, "Number of facets must be non-zero.");

  if (numf > std::numeric_limits<FACET_INDEX>::max())
    throw PROCEDURE_ERROR
      (procname, "Number of polyhedron facets is too large.");

  num_facets = numf;
  facet = new FACET[numf];
  num_facet_vertices = new int[numf];
  facet_vertex_list = new int *[numf];

  if (facet == NULL || num_facet_vertices == NULL || 
      facet_vertex_list == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for list of facets.");

  // initialize each facet to 0
  for (int jf = 0; jf < numf; jf++) {
    facet[jf] = 0;
    num_facet_vertices[jf] = 0;
    facet_vertex_list[jf] = NULL;
  };

}

/// Set number of vertices in facet \a jf.
void ISOSURFACE_TABLE_POLYHEDRON::SetNumFacetVertices
(const FACET_INDEX jf, const int numv)
{
  const char * procname = 
    "ISOSURFACE_TABLE_POLYHEDRON::SetNumFacetVertices";

  int jf_int = int(jf);
  if (jf_int < 0 || jf_int >= NumFacets())
    throw PROCEDURE_ERROR(procname, "Illegal facet index.");


  if (facet_vertex_list[jf] != NULL) {
    delete [] facet_vertex_list[jf];
    facet_vertex_list[jf] = NULL;
  }

  num_facet_vertices[jf] = numv;
  facet_vertex_list[jf] = new int[numv];
}

void ISOSURFACE_TABLE_POLYHEDRON::SetVertexCoord
(const int iv, const int ic, const int coord)
  // set polyhedron vertex coordinate
  // iv = vertex index.  In range [0..NumVertices()-1].
  // ic = coordinate index. In range [0..Dimension()-1].
  // coord = coordinate.  Must be even.
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetVertexCoord";

  if (iv < 0 || iv >= NumVertices())
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");

  if (ic < 0 || ic >= Dimension())
    throw PROCEDURE_ERROR(procname, "Illegal vertex coordinate index.");

  if (coord%2 != 0)
    throw PROCEDURE_ERROR
      (procname, "Illegal vertex coordinate.  Vertex coordinate must be even.");

  if (vertex_coord == NULL)
    throw PROCEDURE_ERROR
      (procname, "Vertex coordinate memory not allocated.");

  vertex_coord[iv*Dimension() + ic] = coord;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetEdge
(const EDGE_INDEX ie, const int iv0, const int iv1)
  // set polyhedron edge coordinate
  // ie = edge index.  In range [0..NumEdges()-1].
  // iv0 = vertex 0 index.  In range [0..NumVertices()-1].
  // iv1 = vertex 1 index.  In range [0..NumVertices()-1].
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetEdge";

  int ie_int = int(ie);
  if (ie_int < 0 || ie_int >= NumEdges())
    throw PROCEDURE_ERROR(procname, "Illegal edge index.");

  if (iv0 < 0 || iv0 > NumVertices() ||
      iv1 < 0 || iv1 > NumVertices())
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");

  if (edge_endpoint == NULL)
    throw PROCEDURE_ERROR(procname, "Edge endpoint memory not allocated.");

  edge_endpoint[2*int(ie)] = iv0;
  edge_endpoint[2*int(ie)+1] = iv1;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetFacetVertex
(const FACET_INDEX jf, const int k, const int iv)
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::SetFacetVertex";

  int jf_int = int(jf);
  if (jf_int < 0 || jf_int >= NumFacets())
    throw PROCEDURE_ERROR(procname, "Illegal facet index.");

  if (k < 0 || k >= NumFacetVertices(jf)) 
    throw PROCEDURE_ERROR
      (procname, "Illegal index to list of facet vertices.");

  if (iv < 0 || iv >= NumVertices())
    throw PROCEDURE_ERROR(procname, "Illegal vertex index.");

  facet_vertex_list[jf_int][k] = iv;

  long mask = 1L << iv;
  facet[jf_int] = facet[jf_int] | mask;
}

int ISOSURFACE_TABLE_POLYHEDRON::MidpointCoord
(const EDGE_INDEX ie, const int ic) const
  // return ic'th coordinate of midpoint of edge ie
  // Note: vertice coordinates must all be even so midpoint coordinate
  //   is an integer
{
  int iv0 = EdgeEndpoint(ie, 0);
  int iv1 = EdgeEndpoint(ie, 1);
  int coord0 = VertexCoord(iv0, ic);
  int coord1 = VertexCoord(iv1, ic);
  
  return((coord0+coord1)/2);
}


void ISOSURFACE_TABLE_POLYHEDRON::GenCube(const int cube_dimension)
// generate a cube
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenCube";

  dimension = cube_dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal cube dimension.");

  const int numv = 1L << Dimension();
  const int nume = (numv * Dimension()) / 2;
  const int numf = 2*Dimension();

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices(); iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
        coord = 2;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };

  // generate edges in lexicographic order
  int ie = 0;
  long control = 0;
  while (ie < NumEdges()) {

    // find first 0 bit in control
    int ic = 0;
    long mask = 1L;
    while ((mask & control) != 0) {
      mask = mask << 1;
      ic++;
    };

    // find start vertex by stripping ic bits from control
    int start = control;
    start = start >> ic;
    start = start << ic;

    for (int i = 0; i < (1L << ic); i++) {
      int iv0 = start + i;
      int iv1 = iv0 + (1L << ic);
      SetEdge(ie, iv0, iv1);
      ie++;
    };

    control++;
  };

  if (control+1 != (1L << Dimension()))
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  // generate facets
  const int num_vertices_per_facet = NumVertices()/2;

  for (int jf = 0; jf < numf; jf++) 
    { SetNumFacetVertices(jf, num_vertices_per_facet); };

  long mask = 1L;
  for (int ic = 0; ic < Dimension(); ic++) {
    int jf0 = 2*ic;
    int jf1 = 2*ic+1;
    int k0 = 0;
    int k1 = 0;

    for (int iv = 0; iv < NumVertices(); iv++) {
      long bit = mask & iv;
      if (bit == 0) {
        SetFacetVertex(jf0, k0, iv);
        k0++;
      }
      else {
        SetFacetVertex(jf1, k1, iv);
        k1++;
      };
    };

    if (k0 != NumFacetVertices(jf0) || k1 != NumFacetVertices(jf1))
      throw PROCEDURE_ERROR
        (procname, "Programming error in facet generation.");

    mask = mask << 1;
  };

}

// Generate a cube.
// Use order of facets and edges given by template class CUBE_FACE_INFO.
void ISOSURFACE_TABLE_POLYHEDRON::GenCubeOrderA(const int cube_dimension)
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenCubeOrderA";

  dimension = cube_dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal cube dimension.");

  CUBE_FACE_INFO<int,int,int> cube(cube_dimension);

  const int numv = cube.NumVertices();
  const int nume = cube.NumEdges();
  const int numf = cube.NumFacets();
  const int num_facet_vertices = cube.NumFacetVertices();

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices(); iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
        coord = 2;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };


  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    int iv0 = cube.EdgeEndpoint(ie, 0);
    int iv1 = cube.EdgeEndpoint(ie, 1);
    SetEdge(ie, iv0, iv1);
  }

  for (int ifacet = 0; ifacet < cube.NumFacets(); ifacet++) {
    SetNumFacetVertices(ifacet, num_facet_vertices);

    for (int k = 0; k < cube.NumFacetVertices(); k++) {
      int iv = cube.FacetVertex(ifacet, k);
      SetFacetVertex(ifacet, k, iv);
    }
  }

}


void ISOSURFACE_TABLE_POLYHEDRON::GenSimplex(const int simplex_dimension)
  // generate a simplex
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenSimplex";

  dimension = simplex_dimension;

  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal simplex dimension.");

  const int numv = Dimension() + 1;
  const int nume = (numv * Dimension()) / 2;
  const int numf = Dimension() + 1;

  SetSize(numv, nume, numf);

  // initialize all vertex coordinates to 0
  for (int iv = 0; iv < NumVertices(); iv++)
    for (int ic = 0; ic < Dimension(); ic++)
      SetVertexCoord(iv, ic, 0);

  // set vertex coordinates
  int ic = 0;
  for (int jv = 1; jv < NumVertices(); jv++) {
    SetVertexCoord(jv, ic, 2);
    ic++;
  };

  // generate edges in lexicographic order
  int ie = 0;
  for (int iv0 = 0; iv0 < NumVertices()-1; iv0++)
    for (int iv1 = iv0+1; iv1 < NumVertices(); iv1++) {
      SetEdge(ie, iv0, iv1);
      ie++;
    };

  if (ie != nume)
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  // generate facets
  const int num_vertices_per_facet = Dimension();

  for (int jf = 0; jf < numf; jf++) 
    { SetNumFacetVertices(jf, num_vertices_per_facet); };

  for (int jf = 0; jf < numf; jf++) {

    int k= 0;

    for (int jv = 0; jv < num_vertices_per_facet; jv++) {
      int iv = (jf + jv) % numv;
      SetFacetVertex(jf, k, iv);
      k++;
    };

    if (k != NumFacetVertices(jf))
      throw PROCEDURE_ERROR
        (procname, "Programming error in facet generation.");
  };
}


void ISOSURFACE_TABLE_POLYHEDRON::GenPyramid(const int pyramid_dimension)
  // generate a pyramid
{
  const char * procname = "ISOSURFACE_TABLE_POLYHEDRON::GenPyramid";

  dimension = pyramid_dimension;

  if (!CheckDimension() || dimension < 2)
    throw PROCEDURE_ERROR(procname, "Illegal pyramid dimension.");

  int numv = (1L << (Dimension()-1)) + 1;
  int num_base_edges = (numv-1) * (Dimension()-1)/2;
  int nume = num_base_edges + (numv-1);
  int numf = 2*Dimension()-1;

  int apex = numv-1;

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices()-1; iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
        coord = 4;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };

  // Set coordinate of pyramid
  for (int ic = 0; ic < Dimension(); ic++) {
    SetVertexCoord(apex, ic, 2);
  }

  // generate edges in lexicographic order
  int ie = 0;
  long control = 0;
  while (ie < num_base_edges) {

    // find first 0 bit in control
    int ic = 0;
    long mask = 1L;
    while ((mask & control) != 0) {
      mask = mask << 1;
      ic++;
    };

    // find start vertex by stripping ic bits from control
    int start = control;
    start = start >> ic;
    start = start << ic;

    for (int i = 0; i < (1L << ic); i++) {
      int iv0 = start + i;
      int iv1 = iv0 + (1L << ic);
      SetEdge(ie, iv0, iv1);
      ie++;
    };

    control++;
  };

  if (control+1 != (1L << (Dimension()-1)))
    throw PROCEDURE_ERROR(procname, "Programming error in edge generation.");

  for (int iv = 0; iv < NumVertices()-1; iv++) {
    SetEdge(num_base_edges+iv, iv, apex);
  }

  // generate facets containing apex
  int num_vertices_per_facet = 1+(numv-1)/2;

  for (int jf = 0; jf+1 < numf; jf++) 
    { SetNumFacetVertices(jf, num_vertices_per_facet); };

  long mask = 1L;
  for (int ic = 0; ic < Dimension()-1; ic++) {
    int jf0 = 2*ic;
    int jf1 = 2*ic+1;

    SetFacetVertex(jf0, 0, apex);
    SetFacetVertex(jf1, 0, apex);

    int k0 = 1;
    int k1 = 1;
    for (int iv = 0; iv < NumVertices()-1; iv++) {
      long bit = mask & iv;
      if (bit == 0) {
        SetFacetVertex(jf0, k0, iv);
        k0++;
      }
      else {
        SetFacetVertex(jf1, k1, iv);
        k1++;
      };
    };

    if (k0 != NumFacetVertices(jf0) || k1 != NumFacetVertices(jf1))
      throw PROCEDURE_ERROR
        (procname, "Programming error in facet generation.");

    mask = (mask << 1L);
  };

  // generate base facet
  int base_facet = numf-1;
  int num_vertices_in_base_facet = numv-1;;
  SetNumFacetVertices(base_facet, num_vertices_in_base_facet);


  for (int iv = 0; iv < NumVertices()-1; iv++) {
    SetFacetVertex(base_facet, iv, iv);
  };

}

bool ISOSURFACE_TABLE_POLYHEDRON::CheckDimension() const
  // check dimension
{
  if (dimension < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE_POLYHEDRON::Check(ERROR & error_msg) const
  // check polyhedron
{
  if (!CheckDimension()) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal polyhedron dimension ",
                         Dimension(), ".");
    return(false);
  }

  if (NumVertices() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal number of vertices.");
    return(false);
  };

  if (NumEdges() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Illegal number of edges.");
    return(false);
  };

  if (vertex_coord == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for vertex coordinate list not allocated.");
    return(false);
  };

  if (edge_endpoint == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for edge endpoint list not allocated.");
    return(false);
  }

  for (int iv = 0; iv < NumVertices(); iv++) {
    for (int ic = 0; ic < Dimension(); ic++) {
      if ((VertexCoord(iv, ic) % 2) != 0) {
        error_msg.ClearAll();
        error_msg.AddMessage("Vertex coordinates must be even integers.");
        return(false);
      };
    };
  };

  for (int ie = 0; ie < NumEdges(); ie++) {
    for (int ip = 0; ip < 2; ip++) {
      int iv = EdgeEndpoint(ie, ip);
      if (iv < 0 || iv >= NumVertices()) {
        error_msg.ClearAll();
        error_msg.AddMessage("Illegal edge endpoint ", iv,
                             " for edge ", ie, ".");
        return(false);
      };
    };
  };

  if (NumFacets() > 0) {
    if (facet == NULL) {
      error_msg.ClearAll();
      error_msg.AddMessage("Memory for facet list not allocated.");
      return(false);
    };
  };

  return(true);
}

// **************************************************
// GENERATE POLYTOPE
// **************************************************

void IJKTABLE::generate_prism
(const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
 ISOSURFACE_TABLE_POLYHEDRON & prism)
  // generate a prism with base base_polyhedron
  // first numv vertices have last coordinate = 0
  // last numv vertices have last coordinate = 2
  // first nume edges connect first numv vertices
  // second nume edges connect second numv vertices
  // third numv edges connect first to second set of vertices
  // (numv = # vertices in base_polyhedron; nume = # edges in base_polyhedron)
{
  int dim = base_polyhedron.Dimension();
  int numc = dim;
  int numv = base_polyhedron.NumVertices();
  int nume = base_polyhedron.NumEdges();
  int numf = base_polyhedron.NumFacets();
  int prism_dim = dim + 1;
  int prism_numc = prism_dim;
  int prism_lastc = prism_numc - 1;
  int prism_numv = numv * 2;
  int prism_nume = nume * 2 + numv;
  int prism_numf = 2 + numf;
  prism.SetDimension(prism_dim);
  prism.SetSize(prism_numv, prism_nume, prism_numf);

  // set prism vertex coord
  for (int iv = 0; iv < numv; iv++) {
    for (int ic = 0; ic < prism_lastc; ic++) {
      int coord = base_polyhedron.VertexCoord(iv, ic);
      prism.SetVertexCoord(iv, ic, coord);
      prism.SetVertexCoord(iv+numv, ic, coord);
    };
    prism.SetVertexCoord(iv, prism_lastc, 0);
    prism.SetVertexCoord(iv+numv, prism_lastc, 2);
  };

  // set edges
  for (int ie = 0; ie < base_polyhedron.NumEdges(); ie++) {
    int iv0 = base_polyhedron.EdgeEndpoint(ie, 0);
    int iv1 = base_polyhedron.EdgeEndpoint(ie, 1);
    prism.SetEdge(ie, iv0, iv1);
    prism.SetEdge(ie+nume, iv0+numv, iv1+numv);
  };

  for (int iv = 0; iv < base_polyhedron.NumVertices(); iv++) {
    int ie = 2*nume + iv;
    prism.SetEdge(ie, iv, iv+numv);
  };

  // set facets
  prism.SetNumFacetVertices(0, numv);
  prism.SetNumFacetVertices(1, numv);

  for (int iv = 0; iv < numv; iv++) {
    // set two base facets
    prism.SetFacetVertex(0, iv, iv);
    prism.SetFacetVertex(1, iv, iv+numv);
  };

  for (int jf = 0; jf < numf; jf++) {
    // set prism facet corresponding to original facet jf
    int prism_jf = 2 + jf;

    int base_num_fv = base_polyhedron.NumFacetVertices(jf);
    prism.SetNumFacetVertices(prism_jf, 2*base_num_fv);

    for (int k = 0; k < base_num_fv; k++) {
      int iv = base_polyhedron.FacetVertex(jf, k);
      prism.SetFacetVertex(prism_jf, k, iv);
      prism.SetFacetVertex(prism_jf, k+base_num_fv, iv+numv);
    }
  };


}

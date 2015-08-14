/// \file shrec_extract.cxx
/// Subroutines for extracting dual isosurface mesh

/*
  Copyright (C) 2011-2015 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "ijkgrid_macros.h"
#include "ijkisopoly.txx"
#include "ijktime.txx"

#include "shrec_extract.h"

using namespace IJK;
using namespace SHREC;


// **************************************************
// EXTRACT ISOPOLY
// **************************************************

/// Extract dual isosurface polytopes.
/// Returns list of isosurface polytope vertices.
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
void SHREC::extract_dual_isopoly
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 SHREC_INFO & shrec_info)
{
  shrec_info.time.extract = 0;

  clock_t t0 = std::clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, shrec_info.time.extract);
}


/// Extract dual isosurface polytopes.
/// Returns list of isosurface polytope vertices.
/// Return locations of isosurface vertices on each facet.
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
void SHREC::extract_dual_isopoly
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 SHREC_INFO & shrec_info)
{
  shrec_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, shrec_info.time.extract);
}

/// Extract dual isosurface polytopes from list of edges.
/// Returns list of isosurface polytope vertices.
/// Return locations of isosurface vertices on each facet.
/// @param edge_list = List of edges. Polytopes are dual to edges.
/// @pre Each edge is an internal grid edge.
void SHREC::extract_dual_isopoly_from_list
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<EDGE_INDEX> & edge_list,
 std::vector<ISO_VERTEX_INDEX> & iso_poly,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex)
{
  const int dimension = scalar_grid.Dimension();

  // initialize output
  iso_poly.clear();
  facet_vertex.clear();

  for (NUM_TYPE i = 0; i < edge_list.size(); i++) {
    EDGE_INDEX edge_index = edge_list[i];
    VERTEX_INDEX iend0 = edge_index/dimension;
    int edge_dir = edge_index%dimension;

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex);
  }
}


// **************************************************
// MAP TO ISOPOLY VERTICES
// **************************************************

// Map cube indices in iso_poly_vert to isosurface vertices
//   in isovert.gcube_list.
void SHREC::map_isopoly_vert
(const ISOVERT & isovert, std::vector<ISO_VERTEX_INDEX> & iso_poly_vert)
{
  for (VERTEX_INDEX i = 0; i < iso_poly_vert.size(); i++) {
    VERTEX_INDEX icube = iso_poly_vert[i];
    VERTEX_INDEX gcube_index = isovert.GCubeIndex(icube);
    iso_poly_vert[i] = gcube_index;
  }
}


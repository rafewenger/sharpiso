/// \file isodual3D_extract.cxx
/// Subroutines for extracting dual isosurface mesh

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

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

#include "ijkgrid_macros.h"
#include "isodual3D_extract.h"
#include "ijktime.txx"

using namespace IJK;
using namespace ISODUAL3D;

namespace {
  void extract_dual_isopatch_around_edge
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iend0, const VERTEX_INDEX iend1, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & iso_poly);

  void extract_dual_isopatch_around_edge_reverse_orient
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iend0, const VERTEX_INDEX iend1, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & iso_poly);
}

/// Extract isosurface polytopes
/// Returns list representing isosurface polytopes
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
void ISODUAL3D::extract_dual_isopoly
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 ISODUAL_INFO & isodual_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_dual_isopoly");

  isodual_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);


    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopatch_around_edge
        (scalar_grid, isovalue, iend0, iend1, edge_dir, iso_poly);
    }
    else if(is_end0_positive && !is_end1_positive) {
      extract_dual_isopatch_around_edge_reverse_orient
        (scalar_grid, isovalue, iend0, iend1, edge_dir, iso_poly);
    }

  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, isodual_info.time.extract);
}

namespace {

  /// Extract dual isosurface patch around an edge.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  /// @pre scalar_grid.NumCubeVertices() > 0.
  void extract_dual_isopatch_around_edge
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iend0, const VERTEX_INDEX iend1, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & iso_poly)
  {
    const int num_facet_vertices = scalar_grid.NumFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VERTEX_INDEX iv0 = 
      iend0 - scalar_grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (VERTEX_INDEX k = 0; k < num_facet_vertices; k++) 
      { iso_poly.push_back(scalar_grid.FacetVertex(iv0, edge_dir, k)); }
  }

  /// Extract dual isosurface patch around an edge, reverse orientation.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  /// @pre scalar_grid.NumCubeVertices() > 0.
  void extract_dual_isopatch_around_edge_reverse_orient
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iend0, const VERTEX_INDEX iend1, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & iso_poly)
  {
    const int num_facet_vertices = scalar_grid.NumFacetVertices();
    const int half_num_facet_vertices = num_facet_vertices/2;

    if (num_facet_vertices == 0) { return; };

    VERTEX_INDEX iv0 = 
      iend0 - scalar_grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    // Add last half_num_facet_vertices vertices first.
    for (VERTEX_INDEX k = half_num_facet_vertices; k < num_facet_vertices; k++) 
      { iso_poly.push_back(scalar_grid.FacetVertex(iv0, edge_dir, k)); }

    for (VERTEX_INDEX k = 0; k < half_num_facet_vertices; k++) 
      { iso_poly.push_back(scalar_grid.FacetVertex(iv0, edge_dir, k)); }
  }

}


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
#include "ijkisopoly.txx"
#include "ijktime.txx"

#include "isodual3D_extract.h"

using namespace IJK;
using namespace ISODUAL3D;


/// Extract dual isosurface polytopes.
/// Returns list of isosurface polytope vertices.
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
  isodual_info.time.extract = 0;

  clock_t t0 = std::clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, isodual_info.time.extract);
}


/// Extract dual isosurface polytopes.
/// Returns list of isosurface polytope vertices.
/// Return locations of isosurface vertices on each facet.
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
/// @param iso_poly[] = vector of isosurface polygope vertices
///   iso_simplices[numv_per_poly*ip+k] = 
///     cube containing k'th vertex of polytope ip.
void ISODUAL3D::extract_dual_isopoly
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
 std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 ISODUAL_INFO & isodual_info)
{
  isodual_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  iso_poly.clear();

  if (scalar_grid.NumCubeVertices() < 1) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

    extract_dual_isopoly_around_bipolar_edge
      (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex);
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, isodual_info.time.extract);
}


/// Extract dual isosurface polytopes using isovert data structure.
/// @param scalar_grid = scalar grid data
/// @param isovalue = isosurface scalar value
void ISODUAL3D::extract_dual_isopoly
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, 
 const ISOVERT & isovert,
 DUAL_ISOSURFACE & dual_isosurface,
 ISODUAL_INFO & isodual_info)
{
  IJK::PROCEDURE_ERROR error("extract_dual_isopoly");


  isodual_info.time.extract = 0;

  clock_t t0 = clock();

  // initialize output
  dual_isosurface.tri_vert.clear();
  dual_isosurface.quad_vert.clear();

  extract_dual_isopoly(scalar_grid, isovalue, dual_isosurface.quad_vert, 
                       isodual_info);

  // Map cube indices in quad_vert to cube locations in isovert.gcube_list()
  for (VERTEX_INDEX i = 0; i < dual_isosurface.quad_vert.size(); i++) {
    VERTEX_INDEX icube = dual_isosurface.quad_vert[i];
    VERTEX_INDEX gcube_index = isovert.sharp_ind_grid.Scalar(icube);

    if (gcube_index < 0 || gcube_index >= isovert.gcube_list.size()) {
      error.AddMessage("Programming error.  Inconsistency between quad_vert and isovert.");
      error.AddMessage("  quad_vert[", i, "] = ", icube, " but");
      error.AddMessage("  isovert.sharp_ind_grid[", icube, "] = ",
                       gcube_index, ".");
      throw error;
    }

    // Replace quad_vert[i].
    dual_isosurface.quad_vert[i] = gcube_index;
  }

  clock_t t1 = clock();
  clock2seconds(t1-t0, isodual_info.time.extract);
}

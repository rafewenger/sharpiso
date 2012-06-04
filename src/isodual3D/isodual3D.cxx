
/// \file isodual3D.cxx
/// Marching cubes/hypercubes isosurface generation
/// Version 0.0.1

/*
IJK: Isosurface Jeneration Kode
Copyright (C) 2011,2012 Rephael Wenger

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


#include "ijktime.txx"
#include "ijkisopoly.txx"

#include "isodual3D.h"
#include "isodual3D_ambig.h"
#include "isodual3D_extract.h"
#include "isodual3D_position.h"

using namespace IJK;
using namespace ISODUAL3D;
using namespace std;


// **************************************************
// DUAL CONTOURING (HYPERCUBES)
// **************************************************

/// Dual Contouring Algorithm.
void ISODUAL3D::dual_contouring
(const ISODUAL_DATA & isodual_data, const SCALAR_TYPE isovalue,
 DUAL_ISOSURFACE & dual_isosurface, ISODUAL_INFO & isodual_info)
{
  const int dimension = isodual_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = isodual_data.ScalarGrid().AxisSize();
  float merge_time = 0.0;
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t_start = clock();

  if (!isodual_data.Check(error)) { throw error; };

  dual_isosurface.Clear();
  isodual_info.time.Clear();

  ISO_MERGE_DATA merge_data(dimension, axis_size);

  if (isodual_data.IsGradientGridSet() &&
    isodual_data.VertexPositionMethod() == GRADIENT_POSITIONING
    || isodual_data.VertexPositionMethod() == EDGE_SIMPLE
    || isodual_data.VertexPositionMethod() == EDGE_COMPLEX) {
      dual_contouring
        (isodual_data.ScalarGrid(), isodual_data.GradientGrid(),
        isovalue, isodual_data,
        dual_isosurface.isopoly_vert, dual_isosurface.vertex_coord,
        merge_data, isodual_info);
  }
  else {
    dual_contouring
      (isodual_data.ScalarGrid(), isovalue,
      isodual_data.VertexPositionMethod(),
      dual_isosurface.isopoly_vert, dual_isosurface.vertex_coord,
      merge_data, isodual_info);
  }

  // store times
  clock_t t_end = clock();
  clock2seconds(t_end-t_start, isodual_info.time.total);
}


// **************************************************
// DUAL CONTOURING
// **************************************************

void ISODUAL3D::dual_contouring
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_POSITION_METHOD vertex_position_method,
 std::vector<VERTEX_INDEX> & isopoly_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
 // extract isosurface using Dual Contouring algorithm
 // returns list of isosurface simplex vertices
 //   and list of isosurface vertex coordinates
 // scalar_grid = scalar grid data
 // isovalue = isosurface scalar value
 // vertex_position_method = vertex position method
 // isopoly_vert[] = list of isosurface polytope vertices
 // vertex_coord[] = list of isosurface vertex coordinates
 //   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
 // merge_data = internal data structure for merging identical edges
 // isodual_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0 = clock();

  isopoly_vert.clear();
  vertex_coord.clear();
  isodual_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isopoly;
  extract_dual_isopoly
    (scalar_grid, isovalue, isopoly, isodual_info);
  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical(isopoly, iso_vlist, isopoly_vert, merge_data);
  clock_t t2 = clock();

  if (vertex_position_method == CUBECENTER) {
    position_dual_isovertices_cube_center
      (scalar_grid, iso_vlist, vertex_coord);
  }

  else {
    // default
    position_dual_isovertices_centroid
      (scalar_grid, isovalue, iso_vlist, vertex_coord);
  }
  clock_t t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}

void ISODUAL3D::dual_contouring
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 std::vector<VERTEX_INDEX> & isopoly_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
 // extract isosurface using Dual Contouring algorithm
 // returns list of isosurface simplex vertices
 //   and list of isosurface vertex coordinates
 // scalar_grid = scalar grid data
 // isovalue = isosurface scalar value
 // isopoly_vert[] = list of isosurface polytope vertices
 // vertex_coord[] = list of isosurface vertex coordinates
 //   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
 // merge_data = internal data structure for merging identical edges
 // isodual_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_POSITION_METHOD vertex_position_method =
    isodual_param.vertex_position_method;
  const bool use_selected_gradients =
    isodual_param.use_selected_gradients;
  const bool use_only_cube_gradients =
    isodual_param.use_only_cube_gradients;
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isodual_param.grad_selection_cube_offset;
  const bool allow_multiple_iso_vertices =
    isodual_param.allow_multiple_iso_vertices;
  const bool flag_separate_neg = isodual_param.flag_separate_neg;
  const bool flag_resolve_ambiguous_facets =
    isodual_param.flag_resolve_ambiguous_facets;
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0, t1, t2, t3;
  t0 = clock();

  isopoly_vert.clear();
  vertex_coord.clear();
  isodual_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isopoly;

  if (allow_multiple_iso_vertices) {

    bool flag_separate_opposite(true);
    IJKDUALTABLE::ISODUAL_CUBE_TABLE 
      isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

    std::vector<FACET_VERTEX_INDEX> facet_vertex;
    extract_dual_isopoly
      (scalar_grid, isovalue, isopoly, facet_vertex, isodual_info);

    t1 = clock();

    std::vector<ISO_VERTEX_INDEX> cube_list;
    std::vector<ISO_VERTEX_INDEX> isopoly_cube;      
    merge_identical(isopoly, cube_list, isopoly_cube, merge_data);

    std::vector<ISO_VERTEX_INDEX> iso_vlist_cube;
    std::vector<FACET_VERTEX_INDEX> iso_vlist_patch;

    if (flag_resolve_ambiguous_facets) {

      std::vector<AMBIGUITY_TYPE> cube_ambig(cube_list.size());
      std::vector<AMBIGUITY_TYPE> iso_vlist_cube_ambig;

      set_cube_ambiguity(scalar_grid, gradient_grid, isovalue,
                         cube_list, isodual_param, cube_ambig);

      set_ambiguity_info(cube_ambig, isodual_info.sharpiso);

      split_dual_isovert
        (scalar_grid, isodual_table, isovalue, 
         cube_list, cube_ambig, isopoly_cube, facet_vertex, 
         iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig,
         isopoly_vert);

      t2 = clock();

      if (vertex_position_method == GRADIENT_POSITIONING) {
        position_dual_isovertices_using_gradients
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig,
           vertex_coord, isodual_info.sharpiso);
      }
      /* NOT YET IMPLEMENTED
      else if (vertex_position_method == EDGE_COMPLEX) {
        // Position using SVD on grid edge-isosurface intersections.
        // Select endpoint gradient which determines edge-isosurface intersection.
        position_dual_isovertices_edgeI_sharp_gradients
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, iso_vlist_ambig, vertex_coord);
      }
      */
      else {
        error.AddMessage("Programming error. Positioning method error.");
        error.AddMessage
          ("  Positioning does not allow multiple isosurface vertices in a cube.");
        throw error;
      }

    }
    else {
      IJK::split_dual_isovert
        (scalar_grid, isodual_table, isovalue, 
         cube_list, isopoly_cube, facet_vertex,
         iso_vlist_cube, iso_vlist_patch, isopoly_vert);

      t2 = clock();

      if (vertex_position_method == GRADIENT_POSITIONING) {
        position_dual_isovertices_using_gradients
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, vertex_coord, isodual_info.sharpiso);
      }
      else if (vertex_position_method == EDGE_COMPLEX) {
        // Position using SVD on grid edge-isosurface intersections.
        // Select endpoint gradient which determines edge-isosurface intersection.
        position_dual_isovertices_edgeI_sharp_gradients
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, vertex_coord);
      }
      else {
        error.AddMessage("Programming error. Positioning method error.");
        error.AddMessage
          ("  Positioning does not allow multiple isosurface vertices in a cube.");
        throw error;
      }
    }

  }
  else {


    extract_dual_isopoly
      (scalar_grid, isovalue, isopoly, isodual_info);
    t1 = clock();

    std::vector<ISO_VERTEX_INDEX> iso_vlist;
    merge_identical(isopoly, iso_vlist, isopoly_vert, merge_data);
    t2 = clock();

    if (vertex_position_method == GRADIENT_POSITIONING) {
      position_dual_isovertices_using_gradients
        (scalar_grid, gradient_grid, isovalue, isodual_param,
         iso_vlist, vertex_coord, isodual_info.sharpiso);
    }
    else if (vertex_position_method == EDGE_SIMPLE) {
      position_dual_isovertices_edgeI_interpolate_gradients
        (scalar_grid, gradient_grid, isovalue,
         isodual_param, iso_vlist, vertex_coord);
    }
    else if (vertex_position_method == EDGE_COMPLEX) {
      // Position using SVD on grid edge-isosurface intersections.
      // Select endpoint gradient which determines edge-isosurface intersection.
      position_dual_isovertices_edgeI_sharp_gradients
        (scalar_grid, gradient_grid, isovalue,
         isodual_param, iso_vlist, vertex_coord);
    }
    else {
      // default
      position_dual_isovertices_centroid
        (scalar_grid, isovalue, iso_vlist, vertex_coord);
    }

    if (isodual_param.flag_reposition) {
      reposition_dual_isovertices
        (scalar_grid, gradient_grid, isovalue, isodual_param,
         iso_vlist, &(vertex_coord.front()), isodual_info.sharpiso);
    }

  }

  t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}


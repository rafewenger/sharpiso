/// \file isodual3D.cxx
/// Marching cubes/hypercubes isosurface generation
/// Version 0.0.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2013 Rephael Wenger

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


#include "ijkisopoly.txx"
#include "ijklist.txx"
#include "ijkmesh.txx"
#include "ijktime.txx"

#include "isodual3D.h"
#include "isodual3D_ambig.h"
#include "isodual3D_decimate.h"
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
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t_start = clock();

  if (!isodual_data.Check(error)) { throw error; };

  dual_isosurface.Clear();
  isodual_info.time.Clear();

  ISO_MERGE_DATA merge_data(dimension, axis_size);

  if (isodual_data.IsGradientGridSet() &&
      isodual_data.VertexPositionMethod() == GRADIENT_POSITIONING
      || isodual_data.VertexPositionMethod() == EDGEI_INTERPOLATE
      || isodual_data.VertexPositionMethod() == EDGEI_GRADIENT) {

    if (isodual_data.flag_merge_sharp) {
      ISOVERT isovert;
      dual_contouring_merge_sharp_from_grad
        (isodual_data.ScalarGrid(), isodual_data.GradientGrid(),
         isovalue, isodual_data, dual_isosurface, isovert,
         isodual_info);
    }
    else {
      dual_contouring_sharp
        (isodual_data.ScalarGrid(), isodual_data.GradientGrid(),
         isovalue, isodual_data, dual_isosurface,
         merge_data, isodual_info);
    }
  }
  else if (isodual_data.AreEdgeISet() &&
           isodual_data.VertexPositionMethod() == EDGEI_INPUT_DATA) {

    ISOVERT isovert;
    dual_contouring_merge_sharp_from_hermite
      (isodual_data.ScalarGrid(), 
       isodual_data.EdgeICoord(), isodual_data.EdgeINormalCoord(),
       isovalue, isodual_data, dual_isosurface, isovert,
       isodual_info);
  }
  else {
    dual_contouring
      (isodual_data.ScalarGrid(), isovalue, isodual_data,
       dual_isosurface.quad_vert, dual_isosurface.vertex_coord,
       merge_data, isodual_info);
  }

  // store times
  clock_t t_end = clock();
  clock2seconds(t_end-t_start, isodual_info.time.total);
}


// **************************************************
// DUAL CONTOURING USING SCALAR DATA
// **************************************************

void ISODUAL3D::dual_contouring
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
// extract isosurface using Dual Contouring algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// vertex_position_method = vertex position method
// isoquad_vert[] = list of isosurface quadrilateral vertices
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// isodual_info = information about running time and grid cubes and edges
{
  if (isodual_param.VertexPositionMethod() == CUBECENTER) {
    dual_contouring_cube_center
      (scalar_grid, isovalue, isoquad_vert, vertex_coord, 
       merge_data, isodual_info);
  }
  else {
    // Default: Position iso vertices at centroid.
    if (isodual_param.allow_multiple_iso_vertices) {
      dual_contouring_centroid_multiv
        (scalar_grid, isovalue, isodual_param.flag_separate_neg,
         isoquad_vert, vertex_coord, merge_data, isodual_info);
    }
    else {
      dual_contouring_centroid
        (scalar_grid, isovalue, isoquad_vert, vertex_coord, 
         merge_data, isodual_info);
    }
  }
}

/// Dual contouring algorithm.
/// Position isosurface vertices at cube centers.
void ISODUAL3D::dual_contouring_cube_center
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
{
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0 = clock();

  isoquad_vert.clear();
  vertex_coord.clear();
  isodual_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
  extract_dual_isopoly
    (scalar_grid, isovalue, isoquad_vert2, isodual_info);
  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
  clock_t t2 = clock();

  position_dual_isovertices_cube_center
    (scalar_grid, iso_vlist, vertex_coord);

  clock_t t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge_identical);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
void ISODUAL3D::dual_contouring_centroid
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
{
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0 = clock();

  isoquad_vert.clear();
  vertex_coord.clear();
  isodual_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isoquad_vert2;

  extract_dual_isopoly
    (scalar_grid, isovalue, isoquad_vert2, isodual_info);
  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
  clock_t t2 = clock();

  position_dual_isovertices_centroid
    (scalar_grid, isovalue, iso_vlist, vertex_coord);

  clock_t t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge_identical);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
/// Allow multiple isosurface vertices in a grid cube.
void ISODUAL3D::dual_contouring_centroid_multiv
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const bool flag_separate_neg,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
{
  const int dimension = scalar_grid.Dimension();
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0 = clock();

  isoquad_vert.clear();
  vertex_coord.clear();
  isodual_info.time.Clear();

  // Create dual isosurface lookup table.
  bool flag_separate_opposite(true);
  IJKDUALTABLE::ISODUAL_CUBE_TABLE 
    isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

  std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
  std::vector<FACET_VERTEX_INDEX> facet_vertex;
  extract_dual_isopoly
    (scalar_grid, isovalue, isoquad_vert2, facet_vertex, isodual_info);

  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> cube_list;
  std::vector<ISO_VERTEX_INDEX> isoquad_cube;      
  merge_identical(isoquad_vert2, cube_list, isoquad_cube, merge_data);

  std::vector<DUAL_ISOVERT> iso_vlist;
  VERTEX_INDEX num_split;
  IJK::split_dual_isovert
    (scalar_grid, isodual_table, isovalue, 
     cube_list, isoquad_cube, facet_vertex, 
     iso_vlist, isoquad_vert, num_split);
  isodual_info.sharpiso.num_cube_multi_isov = num_split;
  isodual_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;


  clock_t t2 = clock();

  position_dual_isovertices_centroid_multi
    (scalar_grid, isodual_table, isovalue, iso_vlist, vertex_coord);

  clock_t t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge_identical);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}

// **************************************************
// DUAL CONTOURING USING SCALAR & GRADIENT DATA
// **************************************************

// Extract isosurface using Dual Contouring algorithm
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use gradients to place isosurface vertices on sharp features. 
void ISODUAL3D::dual_contouring_sharp
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 DUAL_ISOSURFACE & dual_isosurface,
 MERGE_DATA & merge_data,
 ISODUAL_INFO & isodual_info)
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

  dual_isosurface.Clear();
  isodual_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isoquad_vert2;

  if (allow_multiple_iso_vertices) {

    bool flag_separate_opposite(true);
    IJKDUALTABLE::ISODUAL_CUBE_TABLE 
      isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

    std::vector<FACET_VERTEX_INDEX> facet_vertex;
    extract_dual_isopoly
      (scalar_grid, isovalue, isoquad_vert2, facet_vertex, isodual_info);

    t1 = clock();

    std::vector<ISO_VERTEX_INDEX> cube_list;
    std::vector<ISO_VERTEX_INDEX> isoquad_cube;      
    merge_identical(isoquad_vert2, cube_list, isoquad_cube, merge_data);


    if (flag_resolve_ambiguous_facets) {

      std::vector<ISO_VERTEX_INDEX> iso_vlist_cube;
      std::vector<FACET_VERTEX_INDEX> iso_vlist_patch;
      std::vector<VERTEX_PAIR> edge_list;
      std::vector<AMBIGUITY_TYPE> cube_ambig(cube_list.size());
      std::vector<AMBIGUITY_TYPE> iso_vlist_cube_ambig;

      set_cube_ambiguity(scalar_grid, gradient_grid, isovalue,
                         cube_list, isodual_param, cube_ambig);

      set_ambiguity_info(cube_ambig, isodual_info.sharpiso);

      VERTEX_INDEX num_split;
      split_dual_isovert_ambig
        (scalar_grid, isodual_table, isovalue, 
         cube_list, cube_ambig, isoquad_cube, facet_vertex, 
         iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig,
         dual_isosurface.quad_vert, num_split);
      isodual_info.sharpiso.num_cube_multi_isov = num_split;
      isodual_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;

      t2 = clock();

      if (vertex_position_method == GRADIENT_POSITIONING) {

        position_dual_isovertices_using_gradients
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig,
           dual_isosurface.vertex_coord, isodual_info.sharpiso);
      }
      else if (vertex_position_method == EDGEI_INTERPOLATE ||
               vertex_position_method == EDGEI_GRADIENT) {

        position_dual_isovertices_edgeI
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig, 
           vertex_position_method, dual_isosurface.vertex_coord,
           isodual_info.sharpiso);
      }
      else {
        error.AddMessage("Programming error. Positioning method error.");
        error.AddMessage
          ("  Positioning does not allow resolving ambiguities in a cube.");
        throw error;
      }

    }
    else {
      VERTEX_INDEX num_split;
      std::vector<DUAL_ISOVERT> iso_vlist;
      IJK::split_dual_isovert
        (scalar_grid, isodual_table, isovalue, 
         cube_list, isoquad_cube, facet_vertex,
         iso_vlist, dual_isosurface.quad_vert, num_split);

      isodual_info.sharpiso.num_cube_multi_isov = num_split;
      isodual_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;

      t2 = clock();

      if (vertex_position_method == GRADIENT_POSITIONING) {

        position_dual_isovertices_using_gradients_multi
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist, dual_isosurface.vertex_coord, isodual_info.sharpiso);
      }
      else if (vertex_position_method == EDGEI_INTERPOLATE ||
               vertex_position_method == EDGEI_GRADIENT) {

        position_dual_isovertices_edgeI_multi
          (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
           iso_vlist, vertex_position_method,
           dual_isosurface.vertex_coord, isodual_info.sharpiso);
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
      (scalar_grid, isovalue, isoquad_vert2, isodual_info);
    t1 = clock();

    std::vector<ISO_VERTEX_INDEX> iso_vlist;
    merge_identical(isoquad_vert2, iso_vlist, 
                    dual_isosurface.quad_vert, merge_data);
    t2 = clock();

    if (vertex_position_method == GRADIENT_POSITIONING) {
      position_dual_isovertices_using_gradients
        (scalar_grid, gradient_grid, isovalue, isodual_param,
         iso_vlist, dual_isosurface.vertex_coord, isodual_info.sharpiso);
    }
    else if (vertex_position_method == EDGEI_INTERPOLATE ||
             vertex_position_method == EDGEI_GRADIENT) {

      // Position using SVD on grid edge-isosurface intersections.
      // Select endpoint gradient which determines edge-isosurface intersection.
      position_dual_isovertices_edgeI
        (scalar_grid, gradient_grid, isovalue, isodual_param, iso_vlist, 
         vertex_position_method, dual_isosurface.vertex_coord,
         isodual_info.sharpiso);
    }
    else {
      // default
      position_dual_isovertices_centroid
        (scalar_grid, isovalue, iso_vlist, dual_isosurface.vertex_coord);
    }

  }

  t3 = clock();

  // store times
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, isodual_info.time.merge_identical);
  clock2seconds(t3-t2, isodual_info.time.position);
  clock2seconds(t3-t0, isodual_info.time.total);
}


// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use gradients to place isosurface vertices on sharp features. 
void ISODUAL3D::dual_contouring_merge_sharp_from_grad
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 DUAL_ISOSURFACE & dual_isosurface,
 ISOVERT & isovert,
 ISODUAL_INFO & isodual_info)
{
  ISOVERT_INFO isovert_info;
  PROCEDURE_ERROR error("dual_contouring");

  if (!gradient_grid.Check
      (scalar_grid, "gradient grid", "scalar grid", error))
    { throw error; }

  clock_t t0, t1, t2, t3, t4;

  dual_isosurface.Clear();
  isodual_info.time.Clear();

  t0 = clock();
    
  compute_dual_isovert
    (scalar_grid, gradient_grid, isovalue, isodual_param, isovert);

  t1 = clock();

  select_sharp_isovert(scalar_grid, isovalue, isodual_param, isovert);

  t2 = clock();

	if (isodual_param.flag_recompute_isovert)
	{
		recompute_isovert_positions
      (scalar_grid, gradient_grid, isovalue, isodual_param, isovert);
	}

  count_vertices(isovert, isovert_info);

  t3 = clock();

  isodual_info.time.merge_sharp = 0;
  dual_contouring_merge_sharp
    (scalar_grid, isovalue, isodual_param, dual_isosurface, isovert,
     isodual_info, isovert_info);

  t4 = clock();

  // store times
  float seconds;
  clock2seconds(t1-t0+t3-t2, isodual_info.time.position);
  clock2seconds(t2-t1, seconds);
  isodual_info.time.merge_sharp += seconds;
  clock2seconds(t4-t0, isodual_info.time.total);
}


// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use isosurface-edge intersections and normals (hermite data)
//   to place isosurface vertices on sharp features. 
void ISODUAL3D::dual_contouring_merge_sharp_from_hermite
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 DUAL_ISOSURFACE & dual_isosurface,
 ISOVERT & isovert,
 ISODUAL_INFO & isodual_info)
{
  ISOVERT_INFO isovert_info;
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t0, t1, t2, t3, t4;

  dual_isosurface.Clear();
  isodual_info.time.Clear();

  t0 = clock();

  compute_dual_isovert
    (scalar_grid, edgeI_coord, edgeI_normal_coord, 
     isovalue, isodual_param, isovert);

  t1 = clock();

  select_sharp_isovert(scalar_grid, isovalue, isodual_param, isovert);

  t2 = clock();

	if (isodual_param.flag_recompute_isovert)
	{
		recompute_isovert_positions
      (scalar_grid, edgeI_coord, isovalue, isodual_param, isovert);
	}

  count_vertices(isovert, isovert_info);

  t3 = clock();

  isodual_info.time.merge_sharp = 0;
  dual_contouring_merge_sharp
    (scalar_grid, isovalue, isodual_param, dual_isosurface, isovert,
     isodual_info, isovert_info);

  t4 = clock();

  // store times
  float seconds;
  clock2seconds(t1-t0+t3-t2, isodual_info.time.position);
  clock2seconds(t2-t1, seconds);
  isodual_info.time.merge_sharp += seconds;
  clock2seconds(t4-t0, isodual_info.time.total);
}

// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void ISODUAL3D::dual_contouring_merge_sharp
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 DUAL_ISOSURFACE & dual_isosurface,
 ISOVERT & isovert,
 ISODUAL_INFO & isodual_info,
 ISOVERT_INFO & isovert_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = isodual_param.flag_separate_neg;
  const bool allow_multiple_iso_vertices =
    isodual_param.allow_multiple_iso_vertices;
  std::vector<VERTEX_INDEX> quad_vert;
  clock_t t0, t1, t2;

  t0 = clock();

  std::vector<DUAL_ISOVERT> iso_vlist;
  std::vector<NUM_TYPE> new_isovert_index;
  std::vector<bool> flag_keep;

  if (allow_multiple_iso_vertices) {

    bool flag_separate_opposite(true);
    IJKDUALTABLE::ISODUAL_CUBE_TABLE 
      isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

    std::vector<ISO_VERTEX_INDEX> isoquad_cube;
    std::vector<FACET_VERTEX_INDEX> facet_vertex;

    extract_dual_isopoly
      (scalar_grid, isovalue, isoquad_cube, facet_vertex, isodual_info);

    map_isopoly_vert(isovert, isoquad_cube);
    t1 = clock();

    full_split_dual_isovert
      (scalar_grid, isodual_table, isovalue,
       isovert, isoquad_cube, facet_vertex, isodual_param,
       iso_vlist, quad_vert, isodual_info.sharpiso);

    merge_sharp_iso_vertices_multi
      (scalar_grid, isodual_table, isovalue, iso_vlist, isovert, 
       isodual_param, quad_vert, isodual_info.sharpiso);

    IJK::get_non_degenerate_quad_btlr
      (quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

    position_dual_isovertices_multi
      (scalar_grid, isodual_table, isovalue, isovert,
       iso_vlist, dual_isosurface.vertex_coord);
  }
  else {

    extract_dual_isopoly(scalar_grid, isovalue, quad_vert, isodual_info);

    map_isopoly_vert(isovert, quad_vert);
    t1 = clock();

    merge_sharp_iso_vertices
      (scalar_grid, isovalue, isovert, isodual_param,
       quad_vert, isodual_info.sharpiso);

    IJK::get_non_degenerate_quad_btlr
      (quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

    copy_isovert_positions
      (isovert.gcube_list, dual_isosurface.vertex_coord);

  }

  if (isodual_param.flag_delete_isolated_vertices) {
    IJK::delete_unreferenced_vertices_two_lists
      (dimension, dual_isosurface.vertex_coord,
       dual_isosurface.tri_vert, dual_isosurface.quad_vert,
       new_isovert_index, flag_keep);
  }

  t2 = clock();

  // Set isodual_info
  isodual_info.sharpiso.num_sharp_corners = isovert_info.num_sharp_corners;
  isodual_info.sharpiso.num_sharp_edges = isovert_info.num_sharp_edges;
  isodual_info.sharpiso.num_smooth_vertices = 
    isovert_info.num_smooth_vertices;

  if (isodual_param.flag_store_isovert_info) {
    if (isodual_param.flag_delete_isolated_vertices) {
      std::vector<DUAL_ISOVERT> iso_vlist2;
      delete_vertices(iso_vlist, new_isovert_index, flag_keep, iso_vlist2);

      set_isovert_info(iso_vlist2, isovert.gcube_list, 
                       isodual_info.sharpiso.vertex_info);
    }
    else {
      set_isovert_info(iso_vlist, isovert.gcube_list, 
                       isodual_info.sharpiso.vertex_info);
    }
  };

  // store times
  float seconds;
  clock2seconds(t1-t0, isodual_info.time.extract);
  clock2seconds(t2-t1, seconds);
  isodual_info.time.merge_sharp += seconds;
}


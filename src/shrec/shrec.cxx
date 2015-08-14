/// \file shrec.cxx
/// SHarp REConstruction of isosurfaces using dual contouring.
/// Version 0.1.0


/*
Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or any later version.

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

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

#include "shrec.h"
#include "shrec_ambig.h"
#include "shrec_select.h"
#include "shrec_merge.h"
#include "shrec_extract.h"
#include "shrec_position.h"
#include "sharpiso_intersect.h"


using namespace IJK;
using namespace SHREC;
using namespace std;


// **************************************************
// DUAL CONTOURING
// **************************************************

/// Dual Contouring Algorithm.
void SHREC::dual_contouring
(const SHREC_DATA & shrec_data, const SCALAR_TYPE isovalue,
 DUAL_ISOSURFACE & dual_isosurface, ISOVERT & isovert, 
 SHREC_INFO & shrec_info)
{
	const int dimension = shrec_data.ScalarGrid().Dimension();
	const AXIS_SIZE_TYPE * axis_size = shrec_data.ScalarGrid().AxisSize();
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t_start = clock();

	if (!shrec_data.Check(error)) { throw error; };

	dual_isosurface.Clear();
	shrec_info.time.Clear();

	ISO_MERGE_DATA merge_data(dimension, axis_size);

	if (shrec_data.IsGradientGridSet() &&
		(shrec_data.flag_grad2hermite || shrec_data.flag_grad2hermiteI)) {
			const GRADIENT_COORD_TYPE max_small_magnitude 
				= shrec_data.max_small_magnitude;

			std::vector<COORD_TYPE> edgeI_coord;
			std::vector<GRADIENT_COORD_TYPE> edgeI_normal_coord;

			if (shrec_data.flag_grad2hermiteI) {
				compute_all_edgeI_linear_interpolate
					(shrec_data.ScalarGrid(), shrec_data.GradientGrid(),
					isovalue, max_small_magnitude, edgeI_coord, edgeI_normal_coord);
			}
			else {
				compute_all_edgeI
					(shrec_data.ScalarGrid(), shrec_data.GradientGrid(),
					isovalue, max_small_magnitude, edgeI_coord, edgeI_normal_coord);
			}

			dual_contouring_merge_sharp_from_hermite
				(shrec_data.ScalarGrid(), edgeI_coord, edgeI_normal_coord,
				isovalue, shrec_data, dual_isosurface, isovert,
				shrec_info);
	}
	else if (shrec_data.IsGradientGridSet() &&
		(shrec_data.VertexPositionMethod() == GRADIENT_POSITIONING
		|| shrec_data.VertexPositionMethod() == EDGEI_INTERPOLATE
		|| shrec_data.VertexPositionMethod() == EDGEI_GRADIENT)) {

			if (shrec_data.flag_merge) {
				dual_contouring_merge_sharp_from_grad
					(shrec_data.ScalarGrid(), shrec_data.GradientGrid(),
					isovalue, shrec_data, dual_isosurface, isovert,
					shrec_info);
			}
			else {
				dual_contouring_sharp_from_grad
					(shrec_data.ScalarGrid(), shrec_data.GradientGrid(),
					isovalue, shrec_data, dual_isosurface,
					isovert, shrec_info);
			}
	}
	else if (shrec_data.AreEdgeISet() &&
		shrec_data.VertexPositionMethod() == EDGEI_INPUT_DATA) {

			dual_contouring_merge_sharp_from_hermite
				(shrec_data.ScalarGrid(), 
				shrec_data.EdgeICoord(), shrec_data.EdgeINormalCoord(),
				isovalue, shrec_data, dual_isosurface, isovert,
				shrec_info);
	}
	else {
		dual_contouring
			(shrec_data.ScalarGrid(), isovalue, shrec_data,
			dual_isosurface.quad_vert, dual_isosurface.vertex_coord,
			merge_data, shrec_info);
	}

	// store times
	clock_t t_end = clock();
	clock2seconds(t_end-t_start, shrec_info.time.total);
}

/// Dual Contouring Algorithm.
void SHREC::dual_contouring
(const SHREC_DATA & shrec_data, const SCALAR_TYPE isovalue,
 DUAL_ISOSURFACE & dual_isosurface, SHREC_INFO & shrec_info)
{
  ISOVERT isovert;

  dual_contouring
    (shrec_data, isovalue, dual_isosurface, isovert, shrec_info);
}


// **************************************************
// DUAL CONTOURING USING SCALAR DATA
// **************************************************

void SHREC::dual_contouring
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	SHREC_INFO & shrec_info)
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
	// shrec_info = information about running time and grid cubes and edges
{
	if (shrec_param.VertexPositionMethod() == CUBECENTER) {
		dual_contouring_cube_center
			(scalar_grid, isovalue, isoquad_vert, vertex_coord, 
			merge_data, shrec_info);
	}
	else {
		// Default: Position iso vertices at centroid.
		if (shrec_param.allow_multiple_iso_vertices) {
			dual_contouring_centroid_multiv
				(scalar_grid, isovalue, shrec_param.flag_separate_neg,
				isoquad_vert, vertex_coord, merge_data, shrec_info);
		}
		else {

			dual_contouring_centroid
				(scalar_grid, isovalue, isoquad_vert, vertex_coord, 
				merge_data, shrec_info);
		}
	}
}

/// Dual contouring algorithm.
/// Position isosurface vertices at cube centers.
void SHREC::dual_contouring_cube_center
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	SHREC_INFO & shrec_info)
{
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	shrec_info.time.Clear();

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, shrec_info);
	clock_t t1 = clock();

	std::vector<ISO_VERTEX_INDEX> iso_vlist;
	merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
	clock_t t2 = clock();

	position_dual_isovertices_cube_center
		(scalar_grid, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, shrec_info.time.merge_identical);
	clock2seconds(t3-t2, shrec_info.time.position);
	clock2seconds(t3-t0, shrec_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
void SHREC::dual_contouring_centroid
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	SHREC_INFO & shrec_info)
{
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	shrec_info.time.Clear();

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;

	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, shrec_info);
	clock_t t1 = clock();

	std::vector<ISO_VERTEX_INDEX> iso_vlist;
	merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
	clock_t t2 = clock();

	position_dual_isovertices_centroid
		(scalar_grid, isovalue, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, shrec_info.time.merge_identical);
	clock2seconds(t3-t2, shrec_info.time.position);
	clock2seconds(t3-t0, shrec_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
/// Allow multiple isosurface vertices in a grid cube.
void SHREC::dual_contouring_centroid_multiv
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const bool flag_separate_neg,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	SHREC_INFO & shrec_info)
{
	const int dimension = scalar_grid.Dimension();
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	shrec_info.time.Clear();

	// Create dual isosurface lookup table.
	bool flag_separate_opposite(true);

	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;
	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, facet_vertex, shrec_info);

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
	shrec_info.sharpiso.num_cube_multi_isov = num_split;
	shrec_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;


	clock_t t2 = clock();
	position_dual_isovertices_centroid_multi
		(scalar_grid, isodual_table, isovalue, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, shrec_info.time.merge_identical);
	clock2seconds(t3-t2, shrec_info.time.position);
	clock2seconds(t3-t0, shrec_info.time.total);
}


// **************************************************
// DUAL CONTOURING USING SCALAR & GRADIENT DATA
// **************************************************

// Extract dual contouring isosurface.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use gradients to place isosurface vertices on sharp features. 
void SHREC::dual_contouring_sharp_from_grad
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring_sharp_from_grid");

	if (!gradient_grid.Check
		(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	shrec_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, gradient_grid, isovalue, shrec_param, 
     shrec_param.vertex_position_method, isovert);

	select_non_smooth(isovert);

	t1 = clock();

	if (shrec_param.allow_multiple_iso_vertices) {

		if (shrec_param.flag_resolve_ambiguous_facets) {

			std::vector<VERTEX_INDEX> cube_list;
			std::vector<AMBIGUITY_TYPE> cube_ambig;

			get_cube_list(isovert, cube_list);
			set_cube_ambiguity(scalar_grid, gradient_grid, isovalue,
				cube_list, shrec_param, cube_ambig);
			set_ambiguity_info(cube_ambig, shrec_info.sharpiso);

			dual_contouring_extract_isopoly_multi
				(scalar_grid, isovalue, shrec_param, dual_isosurface, isovert,
				cube_ambig, shrec_info, isovert_info);
		}
		else {
			dual_contouring_extract_isopoly_multi
				(scalar_grid, isovalue, shrec_param, dual_isosurface, isovert,
				shrec_info, isovert_info);
		}
	}
	else {
		dual_contouring_extract_isopoly
			(scalar_grid, isovalue, shrec_param, dual_isosurface, isovert,
			shrec_info, isovert_info);
	}

	t2 = clock();

	// Set shrec_info
	count_vertices(isovert, isovert_info);
	shrec_info.sharpiso.Set(isovert_info);

	// store times
	float seconds;
	clock2seconds(t1-t0, seconds);
	shrec_info.time.position += seconds;
	clock2seconds(t2-t0, shrec_info.time.total);
}

// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void SHREC::dual_contouring_extract_isopoly
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;

	extract_dual_isopoly(scalar_grid, isovalue, 
		dual_isosurface.quad_vert, shrec_info);

	map_isopoly_vert(isovert, dual_isosurface.quad_vert);

	t1 = clock();

	copy_isovert_positions
		(isovert.gcube_list, dual_isosurface.vertex_coord);

	t2 = clock();

	if (shrec_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			shrec_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.position += seconds;
}

// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void SHREC::dual_contouring_extract_isopoly_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = shrec_param.flag_separate_neg;
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;

	bool flag_separate_opposite(true);
	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

	std::vector<VERTEX_INDEX> isoquad_cube;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;

	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_cube, facet_vertex, shrec_info);

	map_isopoly_vert(isovert, isoquad_cube);

	full_split_dual_isovert
		(scalar_grid, isodual_table, isovalue,
		isovert, isoquad_cube, facet_vertex, shrec_param,
		iso_vlist, dual_isosurface.quad_vert, shrec_info.sharpiso);

	t1 = clock();

	position_dual_isovertices_multi
		(scalar_grid, isodual_table, isovalue, isovert,
		iso_vlist, dual_isosurface.vertex_coord);

	t2 = clock();

	if (shrec_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			shrec_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.position += seconds;
}


// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// Resolve ambiguous facets.
// @pre isovert contains isovert locations.
void SHREC::dual_contouring_extract_isopoly_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	const std::vector<AMBIGUITY_TYPE> & cube_ambig,
	SHREC_INFO & shrec_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = shrec_param.flag_separate_neg;
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;

	bool flag_separate_opposite(true);
	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

	std::vector<VERTEX_INDEX> isoquad_cube;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;
	std::vector<VERTEX_INDEX> cube_list;

	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_cube, facet_vertex, shrec_info);

	map_isopoly_vert(isovert, isoquad_cube);

	get_cube_list(isovert, cube_list);

	VERTEX_INDEX num_split;
	split_dual_isovert_ambig
		(scalar_grid, isodual_table, isovalue, 
		cube_list, cube_ambig, isoquad_cube, facet_vertex, 
		iso_vlist, dual_isosurface.quad_vert, num_split);
	shrec_info.sharpiso.num_cube_multi_isov = num_split;
	shrec_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;

	t1 = clock();

	position_dual_isovertices_multi
		(scalar_grid, isodual_table, isovalue, isovert,
		iso_vlist, dual_isosurface.vertex_coord);

	t2 = clock();

	if (shrec_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			shrec_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.position += seconds;
}


// **************************************************
// MERGE SHARP
// **************************************************
/**
// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use gradients to place isosurface vertices on sharp features. 
*/
void SHREC::dual_contouring_merge_sharp_from_grad
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring");

	if (!gradient_grid.Check
		(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	shrec_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, gradient_grid, isovalue, shrec_param, 
     shrec_param.vertex_position_method, isovert);

	t1 = clock();

  if (shrec_param.flag_select_mod6) {
    select_sharp_isovert_mod6
      (scalar_grid, gradient_grid, isovalue, shrec_param, isovert);
  }
  else if (shrec_param.flag_select_mod3) {
    select_sharp_isovert_mod3
      (scalar_grid, gradient_grid, isovalue, shrec_param, isovert);
  }
  else {
    select_sharp_isovert
      (scalar_grid, gradient_grid, isovalue, shrec_param, isovert);
  }

	t2 = clock();

	if (shrec_param.flag_recompute_isovert) {
		recompute_isovert_positions
			(scalar_grid, gradient_grid, isovalue, shrec_param, isovert);
	}

	count_vertices(isovert, isovert_info);

	t3 = clock();

	shrec_info.time.merge_sharp = 0;
	dual_contouring_merge_sharp
		(scalar_grid, isovalue, shrec_param, dual_isosurface, isovert,
		shrec_info, isovert_info);

	t4 = clock();

	// store times
	float seconds;
	clock2seconds(t1-t0+t3-t2, shrec_info.time.position);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.merge_sharp += seconds;
	clock2seconds(t4-t0, shrec_info.time.total);
}


// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use isosurface-edge intersections and normals (hermite data)
//   to place isosurface vertices on sharp features. 
void SHREC::dual_contouring_merge_sharp_from_hermite
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	shrec_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, edgeI_coord, edgeI_normal_coord, 
		isovalue, shrec_param, isovert);

	t1 = clock();

	select_sharp_isovert(scalar_grid, isovalue, shrec_param, isovert);

	t2 = clock();

	if (shrec_param.flag_recompute_isovert) {
		recompute_isovert_positions
			(scalar_grid, edgeI_coord, isovalue, shrec_param, isovert);
	}

	count_vertices(isovert, isovert_info);

	t3 = clock();

	shrec_info.time.merge_sharp = 0;
	dual_contouring_merge_sharp
		(scalar_grid, isovalue, shrec_param, dual_isosurface, isovert,
		shrec_info, isovert_info);

	t4 = clock();

	// store times
	float seconds;
	clock2seconds(t1-t0+t3-t2, shrec_info.time.position);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.merge_sharp += seconds;
	clock2seconds(t4-t0, shrec_info.time.total);
}

// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void SHREC::dual_contouring_merge_sharp
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHREC_PARAM & shrec_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	SHREC_INFO & shrec_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = shrec_param.flag_separate_neg;
	const bool flag_split_non_manifold = shrec_param.flag_split_non_manifold;
	const bool flag_select_split = shrec_param.flag_select_split;
	const bool allow_multiple_iso_vertices =
		shrec_param.allow_multiple_iso_vertices;
	std::vector<VERTEX_INDEX> quad_vert;
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;
	std::vector<NUM_TYPE> new_isovert_index;
	std::vector<bool> flag_keep;

	if (allow_multiple_iso_vertices) {

		const NUM_TYPE num_gcube = isovert.gcube_list.size();
		std::vector<IJKDUALTABLE::TABLE_INDEX> table_index(num_gcube);
		std::vector<ISO_VERTEX_INDEX> cube_list(num_gcube);

		for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) 
		{ cube_list[i] = isovert.gcube_list[i].cube_index; }

		bool flag_separate_opposite(true);
		IJKDUALTABLE::ISODUAL_CUBE_TABLE 
			isodual_table(dimension, flag_separate_neg, flag_separate_opposite);
    IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO ambig_info(dimension);

		std::vector<ISO_VERTEX_INDEX> isoquad_cube;
		std::vector<FACET_VERTEX_INDEX> facet_vertex;

		extract_dual_isopoly
			(scalar_grid, isovalue, isoquad_cube, facet_vertex, shrec_info);

		map_isopoly_vert(isovert, isoquad_cube);
		t1 = clock();

		compute_cube_isotable_index
			(scalar_grid, isodual_table, isovalue, cube_list, table_index);

		if (flag_split_non_manifold || flag_select_split) {

			if (flag_split_non_manifold) {
				int num_non_manifold_split;
				IJK::split_non_manifold_isov_pairs
					(scalar_grid, isodual_table, ambig_info, cube_list,
					table_index, num_non_manifold_split);
				shrec_info.sharpiso.num_non_manifold_split = 
					num_non_manifold_split;
			}

			if (flag_select_split) {
				int num_1_2_change;
				IJK::select_split_1_2_ambig
					(scalar_grid, isodual_table, ambig_info, isovalue, cube_list,
					table_index, num_1_2_change);
				shrec_info.sharpiso.num_1_2_change =  num_1_2_change;
			}
		}

		int num_split;
		split_dual_isovert
			(isodual_table, cube_list, table_index, 
			isoquad_cube, facet_vertex, iso_vlist, quad_vert, num_split);

		shrec_info.sharpiso.num_cube_multi_isov = num_split;
		shrec_info.sharpiso.num_cube_single_isov = num_gcube - num_split;

		store_table_index(table_index, isovert.gcube_list);

		merge_sharp_iso_vertices_multi
			(scalar_grid, isodual_table, ambig_info, isovalue, iso_vlist, isovert, 
			shrec_param, quad_vert, shrec_info.sharpiso);

		IJK::get_non_degenerate_quad_btlr
			(quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

    if (shrec_param.flag_recompute_using_adjacent) {
      recompute_using_adjacent(scalar_grid, isovert);
    }

		position_merged_dual_isovertices_multi
			(scalar_grid, isodual_table, isovalue, isovert,
			iso_vlist, dual_isosurface.vertex_coord);
	}
	else {

		extract_dual_isopoly(scalar_grid, isovalue, quad_vert, shrec_info);

		map_isopoly_vert(isovert, quad_vert);
		t1 = clock();

		merge_sharp_iso_vertices
			(scalar_grid, isovalue, isovert, shrec_param,
			quad_vert, shrec_info.sharpiso);

		IJK::get_non_degenerate_quad_btlr
			(quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

		copy_isovert_positions
			(isovert.gcube_list, dual_isosurface.vertex_coord);

	}

	if (shrec_param.flag_delete_isolated_vertices) {
		IJK::delete_unreferenced_vertices_two_lists
			(dimension, dual_isosurface.vertex_coord,
			dual_isosurface.tri_vert, dual_isosurface.quad_vert,
			new_isovert_index, flag_keep);
	}

	t2 = clock();

	// Set shrec_info
	shrec_info.sharpiso.num_sharp_corners = isovert_info.num_sharp_corners;
	shrec_info.sharpiso.num_sharp_edges = isovert_info.num_sharp_edges;
	shrec_info.sharpiso.num_smooth_vertices = 
		isovert_info.num_smooth_vertices;

	if (shrec_param.flag_store_isovert_info) {
		if (shrec_param.flag_delete_isolated_vertices) {
			std::vector<DUAL_ISOVERT> iso_vlist2;
			delete_vertices(iso_vlist, new_isovert_index, flag_keep, iso_vlist2);

			set_isovert_info(iso_vlist2, isovert.gcube_list, 
				shrec_info.sharpiso.vertex_info);
		}
		else {
			set_isovert_info(iso_vlist, isovert.gcube_list, 
				shrec_info.sharpiso.vertex_info);
		}
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, shrec_info.time.extract);
	clock2seconds(t2-t1, seconds);
	shrec_info.time.merge_sharp += seconds;
}

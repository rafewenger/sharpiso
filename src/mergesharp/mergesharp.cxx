/// \file mergesharp.cxx
/// Marching cubes/hypercubes isosurface generation
/// Version 0.0.1

/*
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

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

#include "mergesharp.h"
#include "mergesharp_ambig.h"
#include "mergesharp_merge.h"
#include "mergesharp_extract.h"
#include "mergesharp_position.h"
#include "sharpiso_intersect.h"


using namespace IJK;
using namespace MERGESHARP;
using namespace std;


// **************************************************
// DUAL CONTOURING
// **************************************************

/// Dual Contouring Algorithm.
void MERGESHARP::dual_contouring
	(const MERGESHARP_DATA & mergesharp_data, const SCALAR_TYPE isovalue,
	DUAL_ISOSURFACE & dual_isosurface, MERGESHARP_INFO & mergesharp_info)
{
	const int dimension = mergesharp_data.ScalarGrid().Dimension();
	const AXIS_SIZE_TYPE * axis_size = mergesharp_data.ScalarGrid().AxisSize();
	ISOVERT isovert;
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t_start = clock();

	if (!mergesharp_data.Check(error)) { throw error; };

	dual_isosurface.Clear();
	mergesharp_info.time.Clear();

	ISO_MERGE_DATA merge_data(dimension, axis_size);

	if (mergesharp_data.IsGradientGridSet() &&
		(mergesharp_data.flag_grad2hermite || mergesharp_data.flag_grad2hermiteI)) {
			const GRADIENT_COORD_TYPE max_small_magnitude 
				= mergesharp_data.max_small_magnitude;

			std::vector<COORD_TYPE> edgeI_coord;
			std::vector<GRADIENT_COORD_TYPE> edgeI_normal_coord;

			if (mergesharp_data.flag_grad2hermiteI) {
				compute_all_edgeI_linear_interpolate
					(mergesharp_data.ScalarGrid(), mergesharp_data.GradientGrid(),
					isovalue, max_small_magnitude, edgeI_coord, edgeI_normal_coord);
			}
			else {
				compute_all_edgeI
					(mergesharp_data.ScalarGrid(), mergesharp_data.GradientGrid(),
					isovalue, max_small_magnitude, edgeI_coord, edgeI_normal_coord);
			}

			dual_contouring_merge_sharp_from_hermite
				(mergesharp_data.ScalarGrid(), edgeI_coord, edgeI_normal_coord,
				isovalue, mergesharp_data, dual_isosurface, isovert,
				mergesharp_info);
	}
	else if (mergesharp_data.IsGradientGridSet() &&
		(mergesharp_data.VertexPositionMethod() == GRADIENT_POSITIONING
		|| mergesharp_data.VertexPositionMethod() == EDGEI_INTERPOLATE
		|| mergesharp_data.VertexPositionMethod() == EDGEI_GRADIENT)) {

			if (mergesharp_data.flag_merge_sharp) {
				dual_contouring_merge_sharp_from_grad
					(mergesharp_data.ScalarGrid(), mergesharp_data.GradientGrid(),
					isovalue, mergesharp_data, dual_isosurface, isovert,
					mergesharp_info);
			}
			else {
				dual_contouring_sharp_from_grad
					(mergesharp_data.ScalarGrid(), mergesharp_data.GradientGrid(),
					isovalue, mergesharp_data, dual_isosurface,
					isovert, mergesharp_info);
			}
	}
	else if (mergesharp_data.AreEdgeISet() &&
		mergesharp_data.VertexPositionMethod() == EDGEI_INPUT_DATA) {

			dual_contouring_merge_sharp_from_hermite
				(mergesharp_data.ScalarGrid(), 
				mergesharp_data.EdgeICoord(), mergesharp_data.EdgeINormalCoord(),
				isovalue, mergesharp_data, dual_isosurface, isovert,
				mergesharp_info);
	}
	else {
		dual_contouring
			(mergesharp_data.ScalarGrid(), isovalue, mergesharp_data,
			dual_isosurface.quad_vert, dual_isosurface.vertex_coord,
			merge_data, mergesharp_info);
	}

	// store times
	clock_t t_end = clock();
	clock2seconds(t_end-t_start, mergesharp_info.time.total);
}


// **************************************************
// DUAL CONTOURING USING SCALAR DATA
// **************************************************

void MERGESHARP::dual_contouring
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	MERGESHARP_INFO & mergesharp_info)
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
	// mergesharp_info = information about running time and grid cubes and edges
{
	if (mergesharp_param.VertexPositionMethod() == CUBECENTER) {
		dual_contouring_cube_center
			(scalar_grid, isovalue, isoquad_vert, vertex_coord, 
			merge_data, mergesharp_info);
	}
	else {
		// Default: Position iso vertices at centroid.
		if (mergesharp_param.allow_multiple_iso_vertices) {
			cout <<"dual contouring centroid "<<endl;
			dual_contouring_centroid_multiv
				(scalar_grid, isovalue, mergesharp_param.flag_separate_neg,
				isoquad_vert, vertex_coord, merge_data, mergesharp_info);
			cout <<"return form dual contoruing centroid "<<endl;
		}
		else {

			dual_contouring_centroid
				(scalar_grid, isovalue, isoquad_vert, vertex_coord, 
				merge_data, mergesharp_info);
		}
	}
}

/// Dual contouring algorithm.
/// Position isosurface vertices at cube centers.
void MERGESHARP::dual_contouring_cube_center
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	MERGESHARP_INFO & mergesharp_info)
{
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	mergesharp_info.time.Clear();

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, mergesharp_info);
	clock_t t1 = clock();

	std::vector<ISO_VERTEX_INDEX> iso_vlist;
	merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
	clock_t t2 = clock();

	position_dual_isovertices_cube_center
		(scalar_grid, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, mergesharp_info.time.merge_identical);
	clock2seconds(t3-t2, mergesharp_info.time.position);
	clock2seconds(t3-t0, mergesharp_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
void MERGESHARP::dual_contouring_centroid
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	MERGESHARP_INFO & mergesharp_info)
{
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	mergesharp_info.time.Clear();

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;

	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, mergesharp_info);
	clock_t t1 = clock();

	std::vector<ISO_VERTEX_INDEX> iso_vlist;
	merge_identical(isoquad_vert2, iso_vlist, isoquad_vert, merge_data);
	clock_t t2 = clock();

	position_dual_isovertices_centroid
		(scalar_grid, isovalue, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, mergesharp_info.time.merge_identical);
	clock2seconds(t3-t2, mergesharp_info.time.position);
	clock2seconds(t3-t0, mergesharp_info.time.total);
}

/// Dual contouring algorithm.
/// Position isosurface vertices at centroid of isosurface-edge intersections.
/// Allow multiple isosurface vertices in a grid cube.
void MERGESHARP::dual_contouring_centroid_multiv
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const bool flag_separate_neg,
	std::vector<VERTEX_INDEX> & isoquad_vert,
	std::vector<COORD_TYPE> & vertex_coord,
	MERGE_DATA & merge_data,
	MERGESHARP_INFO & mergesharp_info)
{
	const int dimension = scalar_grid.Dimension();
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0 = clock();

	isoquad_vert.clear();
	vertex_coord.clear();
	mergesharp_info.time.Clear();

	// Create dual isosurface lookup table.
	bool flag_separate_opposite(true);

	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

	std::vector<ISO_VERTEX_INDEX> isoquad_vert2;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;
	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_vert2, facet_vertex, mergesharp_info);

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
	mergesharp_info.sharpiso.num_cube_multi_isov = num_split;
	mergesharp_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;


	clock_t t2 = clock();
	position_dual_isovertices_centroid_multi
		(scalar_grid, isodual_table, isovalue, iso_vlist, vertex_coord);

	clock_t t3 = clock();

	// store times
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, mergesharp_info.time.merge_identical);
	clock2seconds(t3-t2, mergesharp_info.time.position);
	clock2seconds(t3-t0, mergesharp_info.time.total);
}


// **************************************************
// DUAL CONTOURING USING SCALAR & GRADIENT DATA
// **************************************************

// Extract dual contouring isosurface.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use gradients to place isosurface vertices on sharp features. 
void MERGESHARP::dual_contouring_sharp_from_grad
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring_sharp_from_grid");

	if (!gradient_grid.Check
		(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	mergesharp_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, gradient_grid, isovalue, mergesharp_param, 
		mergesharp_param.vertex_position_method, isovert, isovert_info);

	select_non_smooth(isovert);

	t1 = clock();

	if (mergesharp_param.allow_multiple_iso_vertices) {

		if (mergesharp_param.flag_resolve_ambiguous_facets) {

			std::vector<VERTEX_INDEX> cube_list;
			std::vector<AMBIGUITY_TYPE> cube_ambig;

			get_cube_list(isovert, cube_list);
			set_cube_ambiguity(scalar_grid, gradient_grid, isovalue,
				cube_list, mergesharp_param, cube_ambig);
			set_ambiguity_info(cube_ambig, mergesharp_info.sharpiso);

			dual_contouring_extract_isopoly_multi
				(scalar_grid, isovalue, mergesharp_param, dual_isosurface, isovert,
				cube_ambig, mergesharp_info, isovert_info);
		}
		else {
			dual_contouring_extract_isopoly_multi
				(scalar_grid, isovalue, mergesharp_param, dual_isosurface, isovert,
				mergesharp_info, isovert_info);
		}
	}
	else {
		dual_contouring_extract_isopoly
			(scalar_grid, isovalue, mergesharp_param, dual_isosurface, isovert,
			mergesharp_info, isovert_info);
	}

	t2 = clock();

	// Set mergesharp_info
	count_vertices(isovert, isovert_info);
	mergesharp_info.sharpiso.Set(isovert_info);

	// store times
	float seconds;
	clock2seconds(t1-t0, seconds);
	mergesharp_info.time.position += seconds;
	clock2seconds(t2-t0, mergesharp_info.time.total);
}

// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void MERGESHARP::dual_contouring_extract_isopoly
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;

	extract_dual_isopoly(scalar_grid, isovalue, 
		dual_isosurface.quad_vert, mergesharp_info);

	map_isopoly_vert(isovert, dual_isosurface.quad_vert);

	t1 = clock();

	copy_isovert_positions
		(isovert.gcube_list, dual_isosurface.vertex_coord);

	t2 = clock();

	if (mergesharp_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			mergesharp_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.position += seconds;
}

// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void MERGESHARP::dual_contouring_extract_isopoly_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = mergesharp_param.flag_separate_neg;
	clock_t t0, t1, t2;

	t0 = clock();

	std::vector<DUAL_ISOVERT> iso_vlist;

	bool flag_separate_opposite(true);
	IJKDUALTABLE::ISODUAL_CUBE_TABLE 
		isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

	std::vector<VERTEX_INDEX> isoquad_cube;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;

	extract_dual_isopoly
		(scalar_grid, isovalue, isoquad_cube, facet_vertex, mergesharp_info);

	map_isopoly_vert(isovert, isoquad_cube);

	full_split_dual_isovert
		(scalar_grid, isodual_table, isovalue,
		isovert, isoquad_cube, facet_vertex, mergesharp_param,
		iso_vlist, dual_isosurface.quad_vert, mergesharp_info.sharpiso);

	t1 = clock();

	position_dual_isovertices_multi
		(scalar_grid, isodual_table, isovalue, isovert,
		iso_vlist, dual_isosurface.vertex_coord);

	t2 = clock();

	if (mergesharp_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			mergesharp_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.position += seconds;
}


// Extract dual contouring isosurface.
// Returns list of isosurface quad vertices
//   and list of isosurface vertex coordinates.
// Resolve ambiguous facets.
// @pre isovert contains isovert locations.
void MERGESHARP::dual_contouring_extract_isopoly_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	const std::vector<AMBIGUITY_TYPE> & cube_ambig,
	MERGESHARP_INFO & mergesharp_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = mergesharp_param.flag_separate_neg;
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
		(scalar_grid, isovalue, isoquad_cube, facet_vertex, mergesharp_info);

	map_isopoly_vert(isovert, isoquad_cube);

	get_cube_list(isovert, cube_list);

	VERTEX_INDEX num_split;
	split_dual_isovert_ambig
		(scalar_grid, isodual_table, isovalue, 
		cube_list, cube_ambig, isoquad_cube, facet_vertex, 
		iso_vlist, dual_isosurface.quad_vert, num_split);
	mergesharp_info.sharpiso.num_cube_multi_isov = num_split;
	mergesharp_info.sharpiso.num_cube_single_isov = cube_list.size() - num_split;

	t1 = clock();

	position_dual_isovertices_multi
		(scalar_grid, isodual_table, isovalue, isovert,
		iso_vlist, dual_isosurface.vertex_coord);

	t2 = clock();

	if (mergesharp_param.flag_store_isovert_info) {
		set_isovert_info(iso_vlist, isovert.gcube_list, 
			mergesharp_info.sharpiso.vertex_info);
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.position += seconds;
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
void MERGESHARP::dual_contouring_merge_sharp_from_grad
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring");

	if (!gradient_grid.Check
		(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	mergesharp_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, gradient_grid, isovalue, mergesharp_param, 
		mergesharp_param.vertex_position_method, isovert, isovert_info);

	t1 = clock();

	select_sharp_isovert(scalar_grid, isovalue, mergesharp_param, isovert);

	t2 = clock();

	if (mergesharp_param.flag_recompute_isovert) {

		//recompute_isovert_positions
		//(scalar_grid, gradient_grid, isovalue, mergesharp_param, isovert, isovert_info);
		// No longer computing covered point location, shifting to centroid, shift to old call 
		recompute_isovert_positions
			(scalar_grid, gradient_grid, isovalue, mergesharp_param, isovert);

	}
	//DEBUG
	/*COORD_TYPE coord1[DIM3]={0.0,0.0,0.0};
	std::cout <<"size "<<isovert.gcube_list.size() << std::endl;
	for (int i = 0; i < isovert.gcube_list.size(); i++)
	{
	 using namespace std;
	 cout <<"cube index "<< isovert.gcube_list[i].cube_index  <<" ";
	 scalar_grid.ComputeCoord(isovert.gcube_list[i].cube_index , coord1);
	 cout <<"{"<<coord1[0]<<","<<coord1[1]<<","<<coord1[2]<<"}";
	 cout <<"all coord "<< isovert.gcube_list[i].isovert_coord[0]
	 <<" "<<isovert.gcube_list[i].isovert_coord[1]
	 <<" "<<isovert.gcube_list[i].isovert_coord[2]
	 <<" eigen val large "<< int(isovert.gcube_list[i].num_eigenvalues)
	  <<" flag "<< isovert.gcube_list[i].flag << endl;
	 scalar_grid.ComputeScaledCoord(isovert.gcube_list[i].cube_index , coord1);
	
	 cout <<"cc "<< coord1[0]-scalar_grid.Spacing(0)/2.0<<" "
	  <<coord1[1]-scalar_grid.Spacing(1)/2.0 <<" "
	  <<coord1[2]-scalar_grid.Spacing(2)/2.0 <<"\n";
	}*/

	count_vertices(isovert, isovert_info);

	t3 = clock();

	mergesharp_info.time.merge_sharp = 0;
	dual_contouring_merge_sharp
		(scalar_grid, isovalue, mergesharp_param, dual_isosurface, isovert,
		mergesharp_info, isovert_info);

	t4 = clock();

	// store times
	float seconds;
	clock2seconds(t1-t0+t3-t2, mergesharp_info.time.position);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.merge_sharp += seconds;
	clock2seconds(t4-t0, mergesharp_info.time.total);
}


// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// Use isosurface-edge intersections and normals (hermite data)
//   to place isosurface vertices on sharp features. 
void MERGESHARP::dual_contouring_merge_sharp_from_hermite
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info)
{
	ISOVERT_INFO isovert_info;
	PROCEDURE_ERROR error("dual_contouring");

	clock_t t0, t1, t2, t3, t4;

	dual_isosurface.Clear();
	mergesharp_info.time.Clear();

	t0 = clock();

	compute_dual_isovert
		(scalar_grid, edgeI_coord, edgeI_normal_coord, 
		isovalue, mergesharp_param, isovert);

	t1 = clock();

	select_sharp_isovert(scalar_grid, isovalue, mergesharp_param, isovert);

	t2 = clock();

	if (mergesharp_param.flag_recompute_isovert) {
		recompute_isovert_positions
			(scalar_grid, edgeI_coord, isovalue, mergesharp_param, isovert);
	}

	count_vertices(isovert, isovert_info);

	t3 = clock();

	mergesharp_info.time.merge_sharp = 0;
	dual_contouring_merge_sharp
		(scalar_grid, isovalue, mergesharp_param, dual_isosurface, isovert,
		mergesharp_info, isovert_info);

	t4 = clock();

	// store times
	float seconds;
	clock2seconds(t1-t0+t3-t2, mergesharp_info.time.position);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.merge_sharp += seconds;
	clock2seconds(t4-t0, mergesharp_info.time.total);
}

// Extract dual contouring isosurface by merging grid cubes
//   around sharp vertices.
// Returns list of isosurface triangle and quad vertices
//   and list of isosurface vertex coordinates.
// @pre isovert contains isovert locations.
void MERGESHARP::dual_contouring_merge_sharp
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const MERGESHARP_PARAM & mergesharp_param,
	DUAL_ISOSURFACE & dual_isosurface,
	ISOVERT & isovert,
	MERGESHARP_INFO & mergesharp_info,
	ISOVERT_INFO & isovert_info)
{
	const int dimension = scalar_grid.Dimension();
	const bool flag_separate_neg = mergesharp_param.flag_separate_neg;
	const bool flag_split_non_manifold = mergesharp_param.flag_split_non_manifold;
	const bool flag_select_split = mergesharp_param.flag_select_split;
	const bool allow_multiple_iso_vertices =
		mergesharp_param.allow_multiple_iso_vertices;
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

		std::vector<ISO_VERTEX_INDEX> isoquad_cube;
		std::vector<FACET_VERTEX_INDEX> facet_vertex;

		extract_dual_isopoly
			(scalar_grid, isovalue, isoquad_cube, facet_vertex, mergesharp_info);

		map_isopoly_vert(isovert, isoquad_cube);
		t1 = clock();

		compute_cube_isotable_index
			(scalar_grid, isodual_table, isovalue, cube_list, table_index);

		if (flag_split_non_manifold || flag_select_split) {
			IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO ambig_info(dimension);

			if (flag_split_non_manifold) {
				int num_non_manifold_split;
				IJK::split_non_manifold_isov_pairs
					(scalar_grid, isodual_table, ambig_info, cube_list,
					table_index, num_non_manifold_split);
				mergesharp_info.sharpiso.num_non_manifold_split = 
					num_non_manifold_split;
			}

			if (flag_select_split) {
				int num_1_2_change;
				IJK::select_split_1_2_ambig
					(scalar_grid, isodual_table, ambig_info, isovalue, cube_list,
					table_index, num_1_2_change);
				mergesharp_info.sharpiso.num_1_2_change =  num_1_2_change;
			}
		}

		int num_split;
		split_dual_isovert
			(isodual_table, cube_list, table_index, 
			isoquad_cube, facet_vertex, iso_vlist, quad_vert, num_split);

		mergesharp_info.sharpiso.num_cube_multi_isov = num_split;
		mergesharp_info.sharpiso.num_cube_single_isov = num_gcube - num_split;

		store_table_index(table_index, isovert.gcube_list);

		merge_sharp_iso_vertices_multi
			(scalar_grid, isodual_table, isovalue, iso_vlist, isovert, 
			mergesharp_param, quad_vert, mergesharp_info.sharpiso);

		IJK::get_non_degenerate_quad_btlr
			(quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

		position_merged_dual_isovertices_multi
			(scalar_grid, isodual_table, isovalue, isovert,
			iso_vlist, dual_isosurface.vertex_coord);
	}
	else {

		extract_dual_isopoly(scalar_grid, isovalue, quad_vert, mergesharp_info);

		map_isopoly_vert(isovert, quad_vert);
		t1 = clock();

		merge_sharp_iso_vertices
			(scalar_grid, isovalue, isovert, mergesharp_param,
			quad_vert, mergesharp_info.sharpiso);

		IJK::get_non_degenerate_quad_btlr
			(quad_vert, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

		copy_isovert_positions
			(isovert.gcube_list, dual_isosurface.vertex_coord);

	}

	if (mergesharp_param.flag_delete_isolated_vertices) {
		IJK::delete_unreferenced_vertices_two_lists
			(dimension, dual_isosurface.vertex_coord,
			dual_isosurface.tri_vert, dual_isosurface.quad_vert,
			new_isovert_index, flag_keep);
	}

	t2 = clock();

	// Set mergesharp_info
	mergesharp_info.sharpiso.num_sharp_corners = isovert_info.num_sharp_corners;
	mergesharp_info.sharpiso.num_sharp_edges = isovert_info.num_sharp_edges;
	mergesharp_info.sharpiso.num_smooth_vertices = 
		isovert_info.num_smooth_vertices;

	if (mergesharp_param.flag_store_isovert_info) {
		if (mergesharp_param.flag_delete_isolated_vertices) {
			std::vector<DUAL_ISOVERT> iso_vlist2;
			delete_vertices(iso_vlist, new_isovert_index, flag_keep, iso_vlist2);

			set_isovert_info(iso_vlist2, isovert.gcube_list, 
				mergesharp_info.sharpiso.vertex_info);
		}
		else {
			set_isovert_info(iso_vlist, isovert.gcube_list, 
				mergesharp_info.sharpiso.vertex_info);
		}
	};

	// store times
	float seconds;
	clock2seconds(t1-t0, mergesharp_info.time.extract);
	clock2seconds(t2-t1, seconds);
	mergesharp_info.time.merge_sharp += seconds;
}

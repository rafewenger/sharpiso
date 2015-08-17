/// \file sharpiso_get_gradients.cxx
/// Get gradients in cube and cube neighbors.
/// Version 0.1.0

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

#include <iomanip>

#include "sharpiso_get_gradients.h"
#include "sharpiso_intersect.h"

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_scalar.txx"


// **************************************************
// GET GRADIENTS
// **************************************************

// local namespace
namespace {

	using namespace SHARPISO;

	inline void add_gradient
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv,
		std::vector<COORD_TYPE> & point_coord,
		std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		std::vector<SCALAR_TYPE> & scalar,
		NUM_TYPE & num_gradients)
	{
		NUM_TYPE ic = point_coord.size();
		point_coord.resize(ic+DIM3);
		gradient_grid.ComputeScaledCoord(iv, &(point_coord[ic]));

		gradient_coord.resize(ic+DIM3);
		std::copy(gradient_grid.VectorPtrConst(iv),
			gradient_grid.VectorPtrConst(iv)+DIM3,
			&(gradient_coord[ic]));

		scalar.push_back(scalar_grid.Scalar(iv));

		num_gradients++;
	}

	inline void add_gradient
		(const std::vector<COORD_TYPE> & edgeI_coord,
		const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
		const NUM_TYPE j,
		const SCALAR_TYPE s,
		std::vector<COORD_TYPE> & point_coord,
		std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		std::vector<SCALAR_TYPE> & scalar,
		NUM_TYPE & num_gradients)
	{
		NUM_TYPE ic = point_coord.size();
		point_coord.resize(ic+DIM3);
		std::copy(edgeI_coord.begin()+j*DIM3, edgeI_coord.begin()+j*DIM3+DIM3,
			point_coord.begin()+ic);

		gradient_coord.resize(ic+DIM3);
		std::copy(edgeI_normal_coord.begin()+j*DIM3, 
			edgeI_normal_coord.begin()+j*DIM3+DIM3,
			gradient_coord.begin()+ic);

		scalar.push_back(s);

		num_gradients++;
	}

	inline void add_large_gradient
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv,
		const GRADIENT_COORD_TYPE max_small_mag_squared,
		std::vector<COORD_TYPE> & point_coord,
		std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		std::vector<SCALAR_TYPE> & scalar,
		NUM_TYPE & num_gradients)
	{
		GRADIENT_COORD_TYPE magnitude_squared =
			gradient_grid.ComputeMagnitudeSquared(iv);

		if (magnitude_squared > max_small_mag_squared) {
			add_gradient(scalar_grid, gradient_grid, iv,
				point_coord, gradient_coord, scalar, num_gradients);
		}
	}

	inline void add_large_gradient
		(const std::vector<COORD_TYPE> & edgeI_coord,
		const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
		const NUM_TYPE j,
		const SCALAR_TYPE s,
		const GRADIENT_COORD_TYPE max_small_mag_squared,
		std::vector<COORD_TYPE> & point_coord,
		std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		std::vector<SCALAR_TYPE> & scalar,
		NUM_TYPE & num_gradients)
	{
		GRADIENT_COORD_TYPE magnitude_squared;

		IJK::compute_inner_product
			(DIM3, &(edgeI_normal_coord[j*DIM3]), &(edgeI_normal_coord[j*DIM3]),
			magnitude_squared);

		if (magnitude_squared > max_small_mag_squared) {
			add_gradient
				(edgeI_coord, edgeI_normal_coord, j, s,
				point_coord, gradient_coord, scalar, num_gradients);
		}
	}

	inline void add_selected_gradient
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX iv,
		const GRID_COORD_TYPE * cube_coord,
		const GRADIENT_COORD_TYPE max_small_mag_squared,
		const SCALAR_TYPE isovalue,
		const OFFSET_VOXEL & voxel,
		std::vector<COORD_TYPE> & point_coord,
		std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		std::vector<SCALAR_TYPE> & scalar,
		NUM_TYPE & num_gradients)
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		// static so not reallocated at each call
		static GRID_COORD_TYPE vertex_coord[DIM3];
		static SIGNED_COORD_TYPE coord[DIM3];

		GRADIENT_COORD_TYPE magnitude_squared =
			gradient_grid.ComputeMagnitudeSquared(iv);

		if (magnitude_squared > max_small_mag_squared) {

			gradient_grid.ComputeScaledCoord(iv, vertex_coord);
			for (DTYPE d = 0; d < DIM3; d++) {
				coord[d] = vertex_coord[d] - cube_coord[d] + 
					voxel.OffsetFactor()*scalar_grid.Spacing(d); 
			}
			const GRADIENT_COORD_TYPE * vertex_gradient_coord =
				gradient_grid.VectorPtrConst(iv);
			SCALAR_TYPE s = scalar_grid.Scalar(iv);

			if (iso_intersects_cube
				(voxel, coord, vertex_gradient_coord, s, isovalue)) {
					add_gradient(scalar_grid, gradient_grid, iv,
						point_coord, gradient_coord, scalar, num_gradients);
			}
		}
	}

	/// Delete unflagged list elements of list[] where flag[i] is true.
	/// Alters list[] and list_length.
	template <typename T>
	void get_flagged_list_elements
		(const bool flag[], T list[], NUM_TYPE & list_length)
	{
		NUM_TYPE new_list_length = 0;
		for (NUM_TYPE i = 0; i < list_length; i++) {
			if (flag[i]) {
				list[new_list_length] = list[i];
				new_list_length++;
			}
		}

		list_length = new_list_length;
	}

	/// Get selected vertices
	void get_selected_vertices
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const GET_GRADIENTS_PARAM & sharpiso_param,
		const OFFSET_VOXEL & voxel,
		const int num_vertices,
		VERTEX_INDEX vertex_list[],
		NUM_TYPE & num_selected)
	{
		const GRADIENT_COORD_TYPE max_small_mag = 
			sharpiso_param.max_small_magnitude;
		const GRADIENT_COORD_TYPE max_small_mag_squared = 
			max_small_mag * max_small_mag;

		// static so not reallocated at each call
		static GRID_COORD_TYPE cube_coord[DIM3];

		IJK::ARRAY<bool> vertex_flag(num_vertices, true);

		deselect_vertices_with_small_gradients
			(gradient_grid, vertex_list, num_vertices, max_small_mag_squared,
			vertex_flag.Ptr());

		if (sharpiso_param.use_selected_gradients &&
			!sharpiso_param.select_based_on_grad_dir) {

				scalar_grid.ComputeCoord(cube_index, cube_coord);

				deselect_vertices_based_on_isoplanes
					(scalar_grid, gradient_grid, cube_coord, voxel,
					isovalue, vertex_list, num_vertices, vertex_flag.Ptr());
		}

		num_selected = num_vertices;
		get_flagged_list_elements
			(vertex_flag.PtrConst(), vertex_list, num_selected);
	}

	/// Get selected vertices.
	/// Resizes vector vertex_list.
	void get_selected_vertices
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const GET_GRADIENTS_PARAM & sharpiso_param,
		const OFFSET_VOXEL & voxel,
		std::vector<VERTEX_INDEX> & vertex_list)
	{
		int num_selected;

		get_selected_vertices
			(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param, voxel,
			vertex_list.size(), &(vertex_list[0]), num_selected);
		vertex_list.resize(num_selected);
	}

	/// Select vertices with large gradient magnitudes.
	/// Returns vertex_list[] with selected vertices in first num_selected positions.
	void select_vertices_with_large_gradient_magnitudes
		(const GRADIENT_GRID_BASE & gradient_grid,
		const GRADIENT_COORD_TYPE max_small_mag_squared,
		const NUM_TYPE num_vertices,
		VERTEX_INDEX vertex_list[],
		NUM_TYPE & num_selected)
	{
		num_selected = 0;
		for (NUM_TYPE i = 0; i < num_vertices; i++) {

			VERTEX_INDEX iv = vertex_list[i];
			GRADIENT_COORD_TYPE magnitude_squared =
				gradient_grid.ComputeMagnitudeSquared(iv);

			if (magnitude_squared > max_small_mag_squared) {
				vertex_list[num_selected] = vertex_list[i];
				num_selected++;
			}
		}

	}

	/// Get cube vertices indicating selected gradients.
	/// @param sharpiso_param Determines which gradients are selected.
	/// @pre No duplicates allowed.
	void get_cube_vertices_with_selected_gradients
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const GET_GRADIENTS_PARAM & sharpiso_param,
		const OFFSET_VOXEL & voxel,
		VERTEX_INDEX cube_vertex_list[NUM_CUBE_VERTICES3D],
		NUM_TYPE & num_selected)
	{
		// NOTE: cube_vertex_list is an array.
		// static so not reallocated at each call
		static GRID_COORD_TYPE cube_coord[DIM3];

		NUM_TYPE num_vertices(0);
		IJK::ARRAY<bool> vertex_flag(NUM_CUBE_VERTICES3D, true);

		if (sharpiso_param.use_intersected_edge_endpoint_gradients) {

			if (sharpiso_param.select_based_on_grad_dir) {
				get_intersected_cube_edge_endpoints_select
					(scalar_grid, gradient_grid, cube_index, isovalue, 
					cube_vertex_list, num_vertices);
				num_selected = num_vertices;
			}
			else {
				get_intersected_cube_edge_endpoints
					(scalar_grid, cube_index, isovalue, cube_vertex_list, num_vertices);
			}
		}
		else if (sharpiso_param.use_gradients_determining_edge_intersections) {
			get_cube_vertices_determining_edge_intersections
				(scalar_grid, gradient_grid, cube_index, isovalue, 
				cube_vertex_list, num_vertices);
		}
		else {
			get_cube_vertices(scalar_grid, cube_index, cube_vertex_list);
			num_vertices = NUM_CUBE_VERTICES3D;
		}

		get_selected_vertices
			(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
			voxel, num_vertices, cube_vertex_list, num_selected);
	}

}

// Get selected grid vertex gradients.
/// @pre Size of vertex_flag[] is at least size of vertex_list[].
void SHARPISO::get_selected_vertex_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
	const bool vertex_flag[],
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	num_gradients = 0;

	for (NUM_TYPE i = 0; i < num_vertices; i++) {

		if (vertex_flag[i]) {

			VERTEX_INDEX iv = vertex_list[i];
			NUM_TYPE ic = point_coord.size();
			point_coord.resize(ic+DIM3);
			gradient_grid.ComputeScaledCoord(iv, &(point_coord[ic]));

			gradient_coord.resize(ic+DIM3);
			std::copy(gradient_grid.VectorPtrConst(iv),
				gradient_grid.VectorPtrConst(iv)+DIM3,
				&(gradient_coord[ic]));

			scalar.push_back(scalar_grid.Scalar(iv));

			num_gradients++;
		}
	}
}

/// Get selected grid vertex gradients.
/// std::vector variation.
void SHARPISO::get_selected_vertex_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const std::vector<VERTEX_INDEX> & vertex_list,
	const bool vertex_flag[],
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	const NUM_TYPE num_vertices = vertex_list.size();

	if (num_vertices > 0) {

		get_selected_vertex_gradients
			(scalar_grid, gradient_grid, &(vertex_list[0]), num_vertices,
			vertex_flag, point_coord, gradient_coord, scalar, 
			num_gradients);
	}
	else {
		num_gradients = 0;
	}
}

/// Get grid vertex gradients.
/// @pre point_coord[] is preallocated to size at least num_vertices*DIM3.
/// @pre gradient_coord[] is preallocated to size at least num_vertices*DIM3.
/// @pre scalar[] is preallocated to size at least num_vertices.
void SHARPISO::get_vertex_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
	COORD_TYPE point_coord[],
	GRADIENT_COORD_TYPE gradient_coord[],
	SCALAR_TYPE scalar[])
{
	for (NUM_TYPE i = 0; i < num_vertices; i++) {

		VERTEX_INDEX iv = vertex_list[i];
		gradient_grid.ComputeScaledCoord(iv, point_coord+i*DIM3);

		std::copy(gradient_grid.VectorPtrConst(iv),
			gradient_grid.VectorPtrConst(iv)+DIM3,
			gradient_coord+i*DIM3);

		scalar[i] = scalar_grid.Scalar(iv);
	}
}

/// Get grid vertex gradients.
void SHARPISO::get_vertex_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE>  & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar)
{
	point_coord.resize(num_vertices*DIM3);
	gradient_coord.resize(num_vertices*DIM3);
	scalar.resize(num_vertices);

	get_vertex_gradients
		(scalar_grid, gradient_grid, vertex_list, num_vertices,
		&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]));
}

/// Get grid vertex gradients.
void SHARPISO::get_vertex_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const std::vector<VERTEX_INDEX> & vertex_list,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE>  & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar)
{
	const NUM_TYPE num_vertices = vertex_list.size();

	point_coord.resize(num_vertices*DIM3);
	gradient_coord.resize(num_vertices*DIM3);
	scalar.resize(num_vertices);

	get_vertex_gradients
		(scalar_grid, gradient_grid, &(vertex_list[0]), num_vertices,
		&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]));
}

// Get all 8 cube gradients
void SHARPISO::get_cube_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3],
	SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D])
{

	for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) {
		VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
		scalar[k] = scalar_grid.Scalar(iv);
		IJK::copy_coord(DIM3, gradient_grid.VectorPtrConst(iv),
			gradient_coord+k*DIM3);
	}
}


/// Get gradients.
/// @param sharpiso_param Determines which gradients are selected.
/// @param flag_sort_gradients If true, sort gradients.  
///        Overrides flag_sort_gradients in sharpiso_param.
void SHARPISO::get_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const SCALAR_TYPE isovalue,
	const GET_GRADIENTS_PARAM & sharpiso_param,
	const OFFSET_VOXEL & voxel,
	const bool flag_sort_gradients,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	//DEBUG
	using namespace std;

	if (sharpiso_param.use_only_cube_gradients && 
		!sharpiso_param.allow_duplicates) 
	{

		// NOTE: cube_vertex_list is an array.
		VERTEX_INDEX cube_vertex_list[NUM_CUBE_VERTICES3D];

		get_cube_vertices_with_selected_gradients
			(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
			voxel, cube_vertex_list, num_gradients);

		if (flag_sort_gradients) {
			sort_vertices_by_isoplane_dist2cc
				(scalar_grid, gradient_grid, isovalue, cube_index,
				cube_vertex_list, num_gradients);
		}

		get_vertex_gradients
			(scalar_grid, gradient_grid, cube_vertex_list, num_gradients, 
			point_coord, gradient_coord, scalar);
	}
	else 
	{

		// NOTE: vertex_list is a C++ vector.
		std::vector<VERTEX_INDEX> vertex_list;
		//cube gradients
		if (sharpiso_param.use_only_cube_gradients) 
		{

			get_cube_vertices_determining_edgeI_allow_duplicates
				(scalar_grid, gradient_grid, cube_index, isovalue, 
				vertex_list);
		}
		//using large neighborhood
		else if (sharpiso_param.use_large_neighborhood) 
		{

			if (sharpiso_param.use_zero_grad_boundary) {

        if (sharpiso_param.use_new_version) {
          get_ie_endpoints_in_large_neighborhood
            (scalar_grid, gradient_grid, isovalue, cube_index,
             sharpiso_param, vertex_list);
        }
        else {
          get_ie_endpoints_in_large_neighborhood_old_version
            (scalar_grid, gradient_grid, isovalue, cube_index,
             sharpiso_param, vertex_list);
        }
			}
			else {

				get_vertices_with_large_gradient_magnitudes
					(scalar_grid, gradient_grid, cube_index, sharpiso_param, vertex_list);
			}
		}//large neighborhood end.
		//use intersected neighbors.
		else if (sharpiso_param.use_intersected_edge_endpoint_gradients) 
		{

			get_intersected_cube_neighbor_edge_endpoints  
				(scalar_grid, cube_index, isovalue, vertex_list);
		}
		else 
		{

			get_cube_neighbor_vertices(scalar_grid, cube_index, vertex_list);
		}

		get_selected_vertices
			(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
			voxel, vertex_list);

		if (flag_sort_gradients) {
			sort_vertices_by_isoplane_dist2cc
				(scalar_grid, gradient_grid, isovalue, cube_index, vertex_list);
		}

		get_vertex_gradients
			(scalar_grid, gradient_grid, vertex_list,
			point_coord, gradient_coord, scalar);
		num_gradients = vertex_list.size();
	}

	//DEBUG
	/* 
	using namespace std;
	cout <<"\ngradients: "<<endl;
	cout <<"num large "<< gradient_coord.size()/3<<" "<< num_gradients<<endl;

	if (num_gradients >0)
	{
	for (int i = 0; i < gradient_coord.size()/3; i++ )
	{

	float gradient_magnitude = 0.0;
	IJK::compute_magnitude(DIM3,&(gradient_coord[DIM3 * i]), gradient_magnitude );
	cout <<"["<<setw(7)<<point_coord[3*i+0] 
	<<" "<<setw(7)<<point_coord[3*i+1] 
	<<" "<<setw(4)<<point_coord[3*i+2]<<"] ["; 
	cout <<setw(3)<<point_coord[3*i+0]/scalar_grid.Spacing(0) 
	<<" "<<setw(3)<<point_coord[3*i+1]/scalar_grid.Spacing(1) 
	<<" "<<setw(3)<<point_coord[3*i+2]/scalar_grid.Spacing(2)<<"]";

	cout <<" s "<<scalar[i];
	cout <<setw(12)<<" [ "<<gradient_coord[3*i+0]/gradient_magnitude 
	<<" "<<setw(12)<<gradient_coord[3*i+1]/gradient_magnitude 
	<<" "<<setw(12)<<gradient_coord[3*i+2]/gradient_magnitude <<"]\n";
	}
	cout <<"\n";
	}*/

}

/// Get gradients from two cubes sharing a facet.
/// Used in getting gradients around a facet.
/// @param sharpiso_param Determines which gradients are selected.
void SHARPISO::get_two_cube_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX facet_v0,
	const NUM_TYPE orth_dir,
	const SCALAR_TYPE isovalue,
	const GET_GRADIENTS_PARAM & sharpiso_param,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	const GRADIENT_COORD_TYPE max_small_mag = 
		sharpiso_param.max_small_magnitude;
	const GRADIENT_COORD_TYPE max_small_mag_squared = 
		max_small_mag * max_small_mag;

	VERTEX_INDEX vertex_list[NUM_TWO_CUBE_VERTICES3D];
	NUM_TYPE num_vertices;

	get_intersected_two_cube_edge_endpoints
		(scalar_grid, facet_v0, orth_dir, isovalue,
		vertex_list, num_vertices);

	select_vertices_with_large_gradient_magnitudes
		(gradient_grid, max_small_mag_squared, 
		num_vertices, vertex_list, num_gradients);

	get_vertex_gradients
		(scalar_grid, gradient_grid, vertex_list, num_gradients,
		point_coord, gradient_coord, scalar);
}

void SHARPISO::get_large_cube_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GRADIENT_COORD_TYPE max_small_mag,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	IJK::PROCEDURE_ERROR error("get_large_cube_gradients");

	// Initialize num_gradients
	num_gradients = 0;

	for (NUM_TYPE k = 0; k < scalar_grid.NumCubeVertices(); k++) {
		VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
		add_large_gradient
			(scalar_grid, gradient_grid, iv, max_small_mag_squared,
			point_coord, gradient_coord, scalar, num_gradients);
	}

}

void SHARPISO::get_large_cube_neighbor_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GRADIENT_COORD_TYPE max_small_mag,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
	IJK::PROCEDURE_ERROR error("get_large_cube_neighbor_gradients");

	// Initialize num_gradients
	num_gradients = 0;

	scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

	get_large_cube_gradients
		(scalar_grid, gradient_grid, cube_index, max_small_mag,
		point_coord, gradient_coord, scalar, num_gradients);

	for (DTYPE d = 0; d < DIM3; d++) {

		if (cube_coord[d] > 0) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
				add_large_gradient
					(scalar_grid, gradient_grid, iv0, max_small_mag_squared,
					point_coord, gradient_coord, scalar, num_gradients);
			}

		}

		if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv2 = iv1 + 2*scalar_grid.AxisIncrement(d);
				add_large_gradient
					(scalar_grid, gradient_grid, iv2, max_small_mag_squared,
					point_coord, gradient_coord, scalar, num_gradients);
			}

		}

	}

}

void SHARPISO::get_selected_cube_neighbor_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
	const SCALAR_TYPE isovalue,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients,
	const OFFSET_VOXEL & voxel)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	IJK::ARRAY<GRID_COORD_TYPE> offset_111(DIM3, 1);
	IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
	IJK::ARRAY<GRID_COORD_TYPE> vertex_coord(DIM3);
	IJK::ARRAY<COORD_TYPE> coord(DIM3);
	IJK::ARRAY<COORD_TYPE> cube_diagonal_coord(DIM3*NUM_CUBE_VERTICES3D);
	IJK::PROCEDURE_ERROR error("get_large_cube_neighbor_gradients");

	// Initialize num_gradients
	num_gradients = 0;

	scalar_grid.ComputeScaledCoord(cube_index, cube_coord.Ptr());

	select_cube_gradients_based_on_isoplanes
		(scalar_grid, gradient_grid, cube_index, max_small_mag, isovalue,
		point_coord, gradient_coord, scalar, num_gradients, voxel);

	for (DTYPE d = 0; d < DIM3; d++) {

		if (cube_coord[d] > 0) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);

				add_selected_gradient
					(scalar_grid, gradient_grid, iv0, cube_coord.PtrConst(),
					max_small_mag_squared, isovalue, voxel,
					point_coord, gradient_coord, scalar, num_gradients);
			}

		}

		if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv2 = iv1 + 2*scalar_grid.AxisIncrement(d);

				add_selected_gradient
					(scalar_grid, gradient_grid, iv2, cube_coord.PtrConst(),
					max_small_mag_squared, isovalue, voxel,
					point_coord, gradient_coord, scalar, num_gradients);
			}

		}
	}

}

namespace {
	int axis_size_222[DIM3] = { 2, 2, 2 };
	SHARPISO_GRID grid_222(DIM3, axis_size_222);

	/// Set corner_flag[icorner] to true for any endpoint 
	///   of an intersected cube edge.
	/// Does NOT initialize corner_flag[] to false.
	void flag_intersected_cube_edge_endpoints
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		bool corner_flag[NUM_CUBE_VERTICES3D])
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		const DTYPE dimension = scalar_grid.Dimension();

		for (DTYPE d = 0; d < dimension; d++) {
			for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
				VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);
				SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
				SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

				if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

					VERTEX_INDEX icorner0 = grid_222.FacetVertex(0, d, k);
					VERTEX_INDEX icorner1 = grid_222.NextVertex(icorner0, d);

					corner_flag[icorner0] = true;
					corner_flag[icorner1] = true;
				}
			}
		}

	}

	/// Set corner_flag[icorner] to true for any endpoint 
	///   of an intersected edge of two cubes sharing a facet.
	/// Does NOT initialize corner_flag[] to false.
	void flag_intersected_two_cube_edge_endpoints
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX facet_v0,
		const NUM_TYPE orth_dir,
		const SCALAR_TYPE isovalue,
		bool vertex_flag[NUM_TWO_CUBE_VERTICES3D])
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		const DTYPE dimension = scalar_grid.Dimension();

		VERTEX_INDEX iv0 = scalar_grid.PrevVertex(facet_v0, orth_dir);

		// *** Unnecessarily recalculates edge intersections on middle facet. ***
		flag_intersected_cube_edge_endpoints
			(scalar_grid, iv0, isovalue, vertex_flag);
		flag_intersected_cube_edge_endpoints
			(scalar_grid, facet_v0, isovalue, vertex_flag+NUM_CUBE_FACET_VERTICES3D);
	}

	/// Return true if isoplane intersects ray from point 0.
	/// @param ray_dir Ray direction. 1 or -1.
	bool does_isoplane_intersect_ray
		(const SCALAR_TYPE s0,
		const GRADIENT_COORD_TYPE g0,
		const SCALAR_TYPE s1,
		const GRADIENT_COORD_TYPE ray_dir)
	{
		if (s0 < s1) {
			if (g0*ray_dir > 0) { return(true); }
		}
		else {
			if (g0*ray_dir < 0) { return(true); }
		}
		return(false);
	}

	void flag_intersected_cube_edge_endpoints_select
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		bool corner_flag[NUM_CUBE_VERTICES3D])
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		const DTYPE dimension = scalar_grid.Dimension();

		for (VERTEX_INDEX i = 0; i < NUM_CUBE_VERTICES3D; i++)
		{ corner_flag[i] = false; }

		for (DTYPE d = 0; d < dimension; d++) {
			for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
				VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);
				SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
				SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

				if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

					VERTEX_INDEX icorner0 = grid_222.FacetVertex(0, d, k);
					VERTEX_INDEX icorner1 = grid_222.NextVertex(icorner0, d);

					GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, d);
					GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, d);

					if (does_isoplane_intersect_ray(s0, g0, s1, 1))
					{ corner_flag[icorner0] = true; }
					if (does_isoplane_intersect_ray(s1, g1, s0, -1))
					{ corner_flag[icorner1] = true; }
				}
			}
		}

	}

}

void SHARPISO::get_intersected_edge_endpoint_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
	const SCALAR_TYPE isovalue,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	bool corner_flag[NUM_CUBE_VERTICES3D];

	IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);

	flag_intersected_cube_edge_endpoints
		(scalar_grid, cube_index, isovalue, corner_flag);

	for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
		if (corner_flag[j]) {
			// grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
			VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
			VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

			add_large_gradient
				(scalar_grid, gradient_grid, iv, max_small_mag_squared, 
				point_coord, gradient_coord, scalar, num_gradients);
		}
	}

}


namespace {

	// Return true if t0 should be selected.
	// @pre 0 <= t0 <= 1 and 0 <= t1 <= 1.
	bool select_t0
		(const SCALAR_TYPE s0, const SCALAR_TYPE s1,
		const GRADIENT_COORD_TYPE g0, const GRADIENT_COORD_TYPE g1,
		const COORD_TYPE t0, const COORD_TYPE t1)
	{
		if (s0 < s1) {
			if (g0 > g1) {
				if (t0 < t1) { return(false); }
			}
			else {
				// g0 <= g1
				if (t0 > t1) { return(false); }
			}
		}
		else if (s0 > s1) {
			if (g0 > g1) {
				if (t0 > t1) { return(false); }
			}
			else {
				// g0 <= g1
				if (t0 < t1) { return(false); }
			}
		}

		return(true);
	}

	void flag_cube_gradients_determining_edge_intersections
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		bool corner_flag[NUM_CUBE_VERTICES3D])
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		const DTYPE dimension = scalar_grid.Dimension();

		for (VERTEX_INDEX i = 0; i < NUM_CUBE_VERTICES3D; i++)
		{ corner_flag[i] = false; }

		for (DTYPE d = 0; d < dimension; d++) {
			for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
				VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);

				if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

					VERTEX_INDEX icorner0 = grid_222.FacetVertex(0, d, k);
					VERTEX_INDEX icorner1 = grid_222.NextVertex(icorner0, d);

					VERTEX_INDEX iv2;
					get_vertex_determining_edge_intersection
						(scalar_grid, gradient_grid, isovalue, iv0, iv1, d, iv2);

					if (iv2 == iv0) { corner_flag[icorner0] = true; }
					else { corner_flag[icorner1] = true; }
				}
			}
		}

	}

}

/// Get gradients of vertices which determine edge isosurface intersections.
void SHARPISO::get_gradients_determining_edge_intersections
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
	const SCALAR_TYPE isovalue,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	bool corner_flag[NUM_CUBE_VERTICES3D];

	flag_cube_gradients_determining_edge_intersections
		(scalar_grid, gradient_grid, cube_index, isovalue, corner_flag);

	for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
		if (corner_flag[j]) {
			// grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
			VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
			VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

			add_large_gradient
				(scalar_grid, gradient_grid, iv, max_small_mag_squared, 
				point_coord, gradient_coord, scalar, num_gradients);
		}
	}

}

namespace {

	int axis_size_444[DIM3] = { 4, 4, 4 };
	SHARPISO_GRID grid_444(DIM3, axis_size_444);

	void add_cube_flags_to_grid_444_flags
		(const bool cube_corner_flag[NUM_CUBE_VERTICES3D],
		const VERTEX_INDEX iw0,
		bool * vertex_flag)
	{
		for (VERTEX_INDEX k = 0; k < NUM_CUBE_VERTICES3D; k++)
			if (cube_corner_flag[k]) {
				VERTEX_INDEX iw = grid_444.CubeVertex(iw0, k);
				vertex_flag[iw] = true;
			}
	}

	void flag_intersected_neighbor_edge_endpoints
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		bool * vertex_flag)
	{
		typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

		const DTYPE dimension = scalar_grid.Dimension();
		IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
		bool corner_flag[NUM_CUBE_VERTICES3D];

		// Index of vertex in grid_444 with coordinates (1,1,1).
		const VERTEX_INDEX iw_111 = 21;

		for (VERTEX_INDEX i = 0; i < grid_444.NumVertices(); i++)
		{ vertex_flag[i] = false; }

		IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);
		flag_intersected_cube_edge_endpoints
			(scalar_grid, cube_index, isovalue, corner_flag);
		add_cube_flags_to_grid_444_flags(corner_flag, iw_111, vertex_flag);

		scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

		for (DTYPE d = 0; d < DIM3; d++) {

			if (cube_coord[d] > 0) {

				VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, d);
				VERTEX_INDEX iw0 = grid_444.PrevVertex(iw_111, d);

				// *** Unnecessarily recalculates some edge intersections.
				IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);
				flag_intersected_cube_edge_endpoints
					(scalar_grid, iv0, isovalue, corner_flag);
				add_cube_flags_to_grid_444_flags(corner_flag, iw0, vertex_flag);
			}

			if (cube_coord[d]+2 < scalar_grid.AxisSize(d)) {
				VERTEX_INDEX iv0 = scalar_grid.NextVertex(cube_index, d);
				VERTEX_INDEX iw0 = grid_444.NextVertex(iw_111, d);

				// *** Unnecessarily recalculates some edge intersections.
				IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);
				flag_intersected_cube_edge_endpoints
					(scalar_grid, iv0, isovalue, corner_flag);
				add_cube_flags_to_grid_444_flags(corner_flag, iw0, vertex_flag);

			}
		}


	}

}

void SHARPISO::get_intersected_neighbor_edge_endpoint_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const GRADIENT_COORD_TYPE max_small_mag,
	const SCALAR_TYPE isovalue,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;

	IJK::ARRAY<bool> vertex_flag(grid_444.NumVertices());

	GRID_COORD_TYPE cube_coord[DIM3];
	GRID_COORD_TYPE coord_inc[DIM3];
	GRID_COORD_TYPE vcoord[DIM3];
	GRID_COORD_TYPE coord_111[DIM3] = { 1, 1, 1 };

	flag_intersected_neighbor_edge_endpoints
		(scalar_grid, cube_index, isovalue, vertex_flag.Ptr());

	scalar_grid.ComputeCoord(cube_index, cube_coord);

	for (VERTEX_INDEX j = 0; j < grid_444.NumVertices(); j++) {
		if (vertex_flag[j]) {

			// *** SLOW COMPUTATION ***
			grid_444.ComputeCoord(j, coord_inc);
			IJK::add_coord_3D(cube_coord, coord_inc, vcoord);
			IJK::subtract_coord_3D(vcoord, coord_111, vcoord);
			VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(vcoord);

			add_large_gradient
				(scalar_grid, gradient_grid, iv, max_small_mag_squared, 
				point_coord, gradient_coord, scalar, num_gradients);
		}
	}

}

/// Get gradients from list of edge-isosurface intersections.
/// @param sharpiso_param Determines which gradients are selected.
/// @param flag_sort_gradients If true, sort gradients.  
///        Overrides flag_sort_gradients in sharpiso_param.
void SHARPISO::get_gradients_from_list
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
	const SHARPISO_EDGE_INDEX_GRID & edge_index,
	const VERTEX_INDEX cube_index,
	const SCALAR_TYPE isovalue,
	const GET_GRADIENTS_PARAM & sharpiso_param,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients)
{
	const GRADIENT_COORD_TYPE max_small_mag =
		sharpiso_param.max_small_magnitude;
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag*max_small_mag;
	IJK::PROCEDURE_ERROR error("get_gradients_from_list");

	for (NUM_TYPE edge_dir = 0; edge_dir < DIM3; edge_dir++) {

		for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
			VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, edge_dir, k);
			VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, edge_dir);

			if (IJK::is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {
				INDEX_DIFF_TYPE j = edge_index.Vector(iv0, edge_dir);

				if (j < 0) {
					error.AddMessage
						("Error.  Missing edge-isosurface intersection for edge (",
						iv0, ",", iv1, ").");
					throw error;
				}

				add_large_gradient
					(edgeI_coord, edgeI_normal_coord, j, isovalue, max_small_mag_squared,
					point_coord, gradient_coord, scalar, num_gradients);
			}
		}
	}
}


// **************************************************
// GET VERTICES
// **************************************************

// Get all 8 cube vertices.
void SHARPISO::get_cube_vertices
	(const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
	VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D])
{
	for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) 
	{ vertex_list[k] = grid.CubeVertex(cube_index, k); }
}


// Get vertices at endpoints of cube edges which intersect the isosurface.
void SHARPISO::get_intersected_cube_edge_endpoints
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
	const SCALAR_TYPE isovalue, 
	VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	bool corner_flag[NUM_CUBE_VERTICES3D];

	// Initialize.
	num_vertices = 0;

	IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);

	flag_intersected_cube_edge_endpoints
		(scalar_grid, cube_index, isovalue, corner_flag);

	for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
		if (corner_flag[j]) {
			// grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
			VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
			VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

			vertex_list[num_vertices] = iv;
			num_vertices++;
		}
	}

}

// Get vertices at endpoints of cube edges which intersect the isosurface.
// Select endpoints whose isoplane intersects the ray 
//   pointing to the opposing endpoint.
void SHARPISO::get_intersected_cube_edge_endpoints_select
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const SCALAR_TYPE isovalue, 
	VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	bool corner_flag[NUM_CUBE_VERTICES3D];

	// Initialize.
	num_vertices = 0;

	IJK::set_c_array(NUM_CUBE_VERTICES3D, 0, corner_flag);

	flag_intersected_cube_edge_endpoints_select
		(scalar_grid, gradient_grid, cube_index, isovalue, corner_flag);

	for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
		if (corner_flag[j]) {
			// grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
			VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
			VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

			vertex_list[num_vertices] = iv;
			num_vertices++;
		}
	}

}


// Get cube vertices which determining intersection of isosurface and edges.
void SHARPISO::get_cube_vertices_determining_edge_intersections
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
	const GRADIENT_GRID_BASE & gradient_grid, const VERTEX_INDEX cube_index,
	const SCALAR_TYPE isovalue, 
	VERTEX_INDEX vertex_list[NUM_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	bool corner_flag[NUM_CUBE_VERTICES3D];

	// Initialize.
	num_vertices = 0;

	flag_cube_gradients_determining_edge_intersections
		(scalar_grid, gradient_grid, cube_index, isovalue, corner_flag);

	for (VERTEX_INDEX j = 0; j < NUM_CUBE_VERTICES3D; j++) {
		if (corner_flag[j]) {
			// grid_222.CubeVertex(0,j) will probably be j, but no guarantees.
			VERTEX_INDEX icorner = grid_222.CubeVertex(0, j);
			VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, icorner);

			vertex_list[num_vertices] = iv;
			num_vertices++;
		}
	}

}


/// Get cube vertices determining the intersection of isosurface and edges.
/// Allow duplicate gradients.
void SHARPISO::get_cube_vertices_determining_edgeI_allow_duplicates
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index, const SCALAR_TYPE isovalue,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();

	for (DTYPE d = 0; d < dimension; d++) {
		for (VERTEX_INDEX k = 0; k < scalar_grid.NumFacetVertices(); k++) {
			VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, k);
			VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);

			if (is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {

				VERTEX_INDEX iv2;
				get_vertex_determining_edge_intersection
					(scalar_grid, gradient_grid, isovalue, iv0, iv1, d, iv2);

				vertex_list.push_back(iv2);
			}
		}
	}

}

/// Get vertices of cube and cube neighbors.
void SHARPISO::get_cube_neighbor_vertices
	(const SHARPISO_GRID & grid, const VERTEX_INDEX cube_index,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	// static so not reallocated at each call
	static GRID_COORD_TYPE cube_coord[DIM3];

	vertex_list.resize(NUM_CUBE_VERTICES3D);
	get_cube_vertices(grid, cube_index, &vertex_list[0]);

	grid.ComputeCoord(cube_index, cube_coord);

	for (DTYPE d = 0; d < DIM3; d++) {

		if (cube_coord[d] > 0) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv0 = grid.PrevVertex(iv1, d);
				vertex_list.push_back(iv0);
			}

		}

		if (cube_coord[d]+2 < grid.AxisSize(d)) {

			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv1 = grid.FacetVertex(cube_index, d, k);
				VERTEX_INDEX iv2 = iv1 + 2*grid.AxisIncrement(d);
				vertex_list.push_back(iv2);
			}

		}

	}

}

// Get vertices at endpoints of cube and neighbor edges 
//   which intersect the isosurface.
void SHARPISO::get_intersected_cube_neighbor_edge_endpoints
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
	const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & vertex_list)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	IJK::ARRAY<bool> vertex_flag (grid_444.NumVertices());

	GRID_COORD_TYPE cube_coord[DIM3];
	GRID_COORD_TYPE coord_inc[DIM3];
	GRID_COORD_TYPE vcoord[DIM3];
	GRID_COORD_TYPE coord_111[DIM3] = { 1, 1, 1 };

	flag_intersected_neighbor_edge_endpoints
		(scalar_grid, cube_index, isovalue, vertex_flag.Ptr());

	scalar_grid.ComputeCoord(cube_index, cube_coord);

	for (VERTEX_INDEX j = 0; j < grid_444.NumVertices(); j++) {
		if (vertex_flag[j]) {

			// *** SLOW COMPUTATION ***
			grid_444.ComputeCoord(j, coord_inc);
			IJK::add_coord_3D(cube_coord, coord_inc, vcoord);
			IJK::subtract_coord_3D(vcoord, coord_111, vcoord);
			VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(vcoord);
			vertex_list.push_back(iv);
		}
	}

}

// Get vertex whose gradient determines intersection of isosurface 
//    and edge (iv0,iv1)
void SHARPISO::get_vertex_determining_edge_intersection
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
	VERTEX_INDEX & iv2)
{
	const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
	const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);
	const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
	const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
	COORD_TYPE t0, t1;

	SHARPISO::compute_edge_intersection(s0, g0, isovalue, t0);
	SHARPISO::compute_edge_intersection(s1, -g1, isovalue, t1);
	t1 = 1-t1;

	iv2 = iv0;  // default

	if (s0 == isovalue) { return; }

	if (s1 == isovalue) { 
		iv2 = iv1; 
		return;
	}

	if (0 <= t1 && t1 <= 1) {
		if (t0 < 0 || t0 > 1) { iv2 = iv1; }
		else {
			if (!select_t0(s0, s1, g0, g1, t0, t1))
			{ iv2 = iv1; }
		}
	}

}

// Compute point on edge where gradient changes.
// @param zero_tolerance No division by numbers less than or equal 
//        to zero_tolerance.
// @pre zero_tolerance must be non-negative.
void SHARPISO::compute_gradient_change_on_edge
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
	const GRADIENT_COORD_TYPE zero_tolerance,
	bool & flag_no_split,
	SCALAR_TYPE & t_split, 
	SCALAR_TYPE & s_split)
{
	const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
	const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);

	const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
	const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
	const GRADIENT_COORD_TYPE gdiff = g0 - g1;

	flag_no_split = false;
	if (fabs(gdiff) <= zero_tolerance) { 
		flag_no_split = true;
		return; 
	}

	t_split = (s1-s0-g1)/gdiff;
	s_split = g0*t_split + s0;
}


// Get vertices at endpoints of facet edges which intersect the isosurface.
// Two cubes share a facet.
// @param facet_v0 Index of primary (lowest-left) facet vertex.
// @param orth_dir Direction orthogonal to facet.
// @pre Facet is in interior of grid.
void SHARPISO::get_intersected_two_cube_edge_endpoints
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX facet_v0,
	const int orth_dir, const SCALAR_TYPE isovalue, 
	VERTEX_INDEX vertex_list[NUM_TWO_CUBE_VERTICES3D], NUM_TYPE & num_vertices)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	bool vertex_flag[NUM_TWO_CUBE_VERTICES3D];

	// Initialize.
	num_vertices = 0;

	IJK::set_c_array(NUM_TWO_CUBE_VERTICES3D, 0, vertex_flag);

	flag_intersected_two_cube_edge_endpoints
		(scalar_grid, facet_v0, orth_dir, isovalue, vertex_flag);

	for (VERTEX_INDEX j = 0; j < NUM_TWO_CUBE_VERTICES3D; j++) {
		if (vertex_flag[j]) {
			VERTEX_INDEX jf = j/NUM_CUBE_FACET_VERTICES3D;
			VERTEX_INDEX k = j%NUM_CUBE_FACET_VERTICES3D;
			VERTEX_INDEX v1 = facet_v0;
			if (jf == 0) {
				v1 = scalar_grid.PrevVertex(facet_v0, orth_dir);
			}
			else if (jf == 2) {
				v1 = scalar_grid.NextVertex(facet_v0, orth_dir);
			}
			VERTEX_INDEX iv = scalar_grid.FacetVertex(v1, orth_dir, k);

			vertex_list[num_vertices] = iv;
			num_vertices++;
		}
	}

}

// Get cube vertices with large gradients.
void SHARPISO::get_cube_vertices_with_large_gradients
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GRADIENT_COORD_TYPE max_small_magnitude,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	for (NUM_TYPE i = 0; i < NUM_CUBE_VERTICES3D; i++) {
		VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, i);

		if (gradient_grid.IsMagnitudeGT(iv, max_small_magnitude)) 
		{ vertex_list.push_back(iv); }
	}
}

// Get vertices with large gradients magnitudes.
void SHARPISO::get_vertices_with_large_gradient_magnitudes
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GRADIENT_COORD_TYPE max_small_magnitude,
	const NUM_TYPE max_dist,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	COORD_TYPE cube_coord[DIM3];

	vertex_list.clear();

	get_cube_vertices_with_large_gradients
		(scalar_grid, gradient_grid, cube_index, max_small_magnitude,
		vertex_list);

	scalar_grid.ComputeCoord(cube_index, cube_coord);
	for (NUM_TYPE d = 0; d < DIM3; d++) {
		for (int i = 0; i < NUM_CUBE_FACET_VERTICES3D; i++) {
			VERTEX_INDEX iv0 = scalar_grid.FacetVertex(cube_index, d, i);
			VERTEX_INDEX iv = iv0;
			NUM_TYPE k = 0;
			while (k+1 < cube_coord[d] && k < max_dist) {
				iv = scalar_grid.PrevVertex(iv, d);
				k++;

				if (gradient_grid.IsMagnitudeGT(iv, max_small_magnitude)) {
					vertex_list.push_back(iv);
					// Don't get any other vertices in this direction.
					break;
				}
			}

			iv = scalar_grid.NextVertex(iv0, d);
			k = 0;
			while (k+cube_coord[d]+1 < scalar_grid.AxisSize(d) && k < max_dist) {
				iv = scalar_grid.NextVertex(iv, d);
				k++;

				if (gradient_grid.IsMagnitudeGT(iv, max_small_magnitude)) {
					vertex_list.push_back(iv);
					// Don't get any other vertices in this direction.
					break;
				}
			}
		}
	}
}

// Get vertices with large gradients magnitudes.
void SHARPISO::get_vertices_with_large_gradient_magnitudes
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GET_GRADIENTS_PARAM & gradient_param,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	const GRADIENT_COORD_TYPE max_small_magnitude =
		gradient_param.max_small_magnitude;
	const NUM_TYPE max_grad_dist = gradient_param.max_grad_dist;
	GRID_COORD_TYPE cube_coord[DIM3];

	get_vertices_with_large_gradient_magnitudes
		(scalar_grid, gradient_grid, cube_index, max_small_magnitude,
		max_grad_dist, vertex_list);

	if (gradient_param.use_diagonal_neighbors) {

		scalar_grid.ComputeCoord(cube_index, cube_coord);

		for (int zcoef = 0; zcoef < 2; zcoef++)
			for (int ycoef = 0; ycoef < 2; ycoef++)
				for (int xcoef = 0; xcoef < 2; xcoef++) {

					VERTEX_INDEX vertex_increment =
						(2*xcoef-1)*scalar_grid.AxisIncrement(0) +
						(2*ycoef-1)*scalar_grid.AxisIncrement(1) +
						(2*zcoef-1)*scalar_grid.AxisIncrement(2);
					VERTEX_INDEX iv0 =  cube_index + 
						xcoef*scalar_grid.AxisIncrement(0) +
						ycoef*scalar_grid.AxisIncrement(1) + 
						zcoef*scalar_grid.AxisIncrement(2);
					NUM_TYPE kmax = max_grad_dist;
					kmax = std::min(kmax, cube_coord[0]+xcoef*max_grad_dist);
					kmax = std::min(kmax, cube_coord[1]+ycoef*max_grad_dist);
					kmax = std::min(kmax, cube_coord[2]+zcoef*max_grad_dist);
					GRID_COORD_TYPE xdiff = scalar_grid.AxisSize(0)-cube_coord[0]-1;
					kmax = std::min(kmax, xdiff+(1-xcoef)*max_grad_dist);
					GRID_COORD_TYPE ydiff = scalar_grid.AxisSize(1)-cube_coord[1]-1;
					kmax = std::min(kmax, ydiff+(1-ycoef)*max_grad_dist);
					GRID_COORD_TYPE zdiff = scalar_grid.AxisSize(2)-cube_coord[2]-1;
					kmax = std::min(kmax, zdiff+(1-zcoef)*max_grad_dist);

					NUM_TYPE k = 0;
					while (k < kmax) {

						VERTEX_INDEX iv = iv0 + vertex_increment;
						k++;

						if (gradient_grid.IsMagnitudeGT(iv, max_small_magnitude)) {
							vertex_list.push_back(iv);
							// Don't get any other vertices in this direction.
							break;
						}
					}
				}
	}

}


// local namespace
namespace {

	/// Return true if iv is adjacent to a visited vertex.
	/// @param boundary_bits Boundary bits for iv.
	bool is_adjacent_to_visited
		(const SHARPISO_BOOL_GRID_BASE & visited,
		const VERTEX_INDEX kv, const long boundary_bits)
	{
		typedef SHARPISO_BOOL_GRID_BASE::DIMENSION_TYPE DTYPE;
		const DTYPE dimension = visited.Dimension();

		for (DTYPE d = 0; d < dimension; d++) {
			long mask = (1L << (2*d));
			long bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = visited.PrevVertex(kv, d);
				if (visited.Scalar(kv2)) 
				{ return(true); }
			}

			mask = (1L << (2*d+1));
			bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = visited.NextVertex(kv, d);
				if (visited.Scalar(kv2)) 
				{ return(true); }
			}
		}

		return(false);
	}

	/// Return true if iv is adjacent to a bipolar edge
	/// @param boundary_bits Boundary bits for iv.
	bool is_adjacent_to_bipolar_edge
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const VERTEX_INDEX kv, const long boundary_bits)
	{
		typedef SHARPISO_SCALAR_GRID_BASE::DIMENSION_TYPE DTYPE;
		const DTYPE dimension = scalar_grid.Dimension();

		for (DTYPE d = 0; d < dimension; d++) {
			long mask = (1L << (2*d));
			long bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = scalar_grid.PrevVertex(kv, d);
				if (IJK::is_gt_min_le_max(scalar_grid, kv, kv2, isovalue))
				{ return(true); }
			}

			mask = (1L << (2*d+1));
			bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = scalar_grid.NextVertex(kv, d);
				if (IJK::is_gt_min_le_max(scalar_grid, kv, kv2, isovalue))
				{ return(true); }
			}
		}

		return(false);
	}

};

// Get intersected edge endpoints in large neighborhood.
// Old, deprecated version.
void SHARPISO::get_ie_endpoints_in_large_neighborhood_old_version
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX cube_index,
	const GET_GRADIENTS_PARAM & gradient_param,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	typedef SHARPISO_SCALAR_GRID_BASE::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	const NUM_TYPE max_grad_dist = gradient_param.max_grad_dist;
	const GRADIENT_COORD_TYPE max_small_magnitude =
		gradient_param.max_small_magnitude;
	VERTEX_INDEX region_iv0, subgrid_cube_index;
	IJK::ARRAY<AXIS_SIZE_TYPE> region_axis_size(dimension);
	std::vector<VERTEX_INDEX> vlist2;
	long boundary_bits;

	vertex_list.clear();

	get_cube_vertices_with_large_gradients
		(scalar_grid, gradient_grid, cube_index, max_small_magnitude,
		vertex_list);

	IJK::compute_region_around_cube
		(cube_index, dimension, scalar_grid.AxisSize(), max_grad_dist,
		region_iv0, region_axis_size.Ptr());

	SHARPISO_INDEX_GRID subgrid(dimension, region_axis_size.PtrConst());
	subgrid.SetToVertexIndices(scalar_grid, region_iv0);

	// Locate subgrid_cube_index
	subgrid_cube_index = 0;
	for (VERTEX_INDEX iv = 0; iv < subgrid.NumVertices(); iv++) {
		if (subgrid.Scalar(iv) == cube_index)
		{ subgrid_cube_index = iv; }
	}

	SHARPISO_BOOL_GRID visited;
	visited.SetSize(subgrid);

	visited.SetAll(false);

	for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) {
		VERTEX_INDEX kv = subgrid.CubeVertex(subgrid_cube_index, k);
		visited.Set(kv, true);
		vlist2.push_back(kv);
	}

	while (vlist2.size() != 0) {
		VERTEX_INDEX kv = vlist2.back();
		vlist2.pop_back();

		subgrid.ComputeBoundaryBits(kv, boundary_bits);

		for (DTYPE d = 0; d < dimension; d++) {

			long mask = (1L << (2*d));
			long bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = subgrid.PrevVertex(kv, d);
				if (!visited.Scalar(kv2)) {
					VERTEX_INDEX iv2 = subgrid.Scalar(kv2);
					if (!gradient_grid.IsMagnitudeGT(iv2, max_small_magnitude)) {
						visited.Set(kv2, true);
						vlist2.push_back(kv2);
					}
				}
			}

			mask = (1L << (2*d+1));
			bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = subgrid.NextVertex(kv, d);
				if (!visited.Scalar(kv2)) {
					VERTEX_INDEX iv2 = subgrid.Scalar(kv2);
					if (!gradient_grid.IsMagnitudeGT(iv2, max_small_magnitude)) {
						visited.Set(kv2, true);
						vlist2.push_back(kv2);
					}
				}
			}
		}
	}

	SHARPISO_SCALAR_GRID scalar_subgrid;
	scalar_subgrid.SetSize(subgrid);

	scalar_subgrid.CopyRegion
		(scalar_grid, region_iv0, scalar_subgrid.AxisSize(), 0);

	for (VERTEX_INDEX kv = 0; kv < subgrid.NumVertices(); kv++) {
		if (!visited.Scalar(kv)) {
			subgrid.ComputeBoundaryBits(kv, boundary_bits);

			if (is_adjacent_to_visited(visited, kv, boundary_bits) &&
				is_adjacent_to_bipolar_edge
				(scalar_subgrid, isovalue, kv, boundary_bits)) {

					VERTEX_INDEX iv = subgrid.Scalar(kv);
					vertex_list.push_back(iv);
			}
		}
	}

}

// Get intersected edge endpoints in large neighborhood.
// Only extends boundary through vertices adjacent to the isosurface.
void SHARPISO::get_ie_endpoints_in_large_neighborhood
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX cube_index,
	const GET_GRADIENTS_PARAM & gradient_param,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	typedef SHARPISO_SCALAR_GRID_BASE::DIMENSION_TYPE DTYPE;

	const DTYPE dimension = scalar_grid.Dimension();
	const NUM_TYPE max_grad_dist = gradient_param.max_grad_dist;
	const GRADIENT_COORD_TYPE max_small_magnitude =
		gradient_param.max_small_magnitude;
	VERTEX_INDEX region_iv0, subgrid_cube_index;
	IJK::ARRAY<AXIS_SIZE_TYPE> region_axis_size(dimension);
	std::vector<VERTEX_INDEX> vlist2;
	long boundary_bits, boundary_bits2;

	vertex_list.clear();

	get_cube_vertices_with_large_gradients
		(scalar_grid, gradient_grid, cube_index, max_small_magnitude,
		vertex_list);

	IJK::compute_region_around_cube
		(cube_index, dimension, scalar_grid.AxisSize(), max_grad_dist,
		region_iv0, region_axis_size.Ptr());

	SHARPISO_INDEX_GRID subgrid(dimension, region_axis_size.PtrConst());
	subgrid.SetToVertexIndices(scalar_grid, region_iv0);

	// Locate subgrid_cube_index
	subgrid_cube_index = 0;
	for (VERTEX_INDEX iv = 0; iv < subgrid.NumVertices(); iv++) {
		if (subgrid.Scalar(iv) == cube_index)
		{ subgrid_cube_index = iv; }
	}

	SHARPISO_BOOL_GRID visited;
	visited.SetSize(subgrid);

	visited.SetAll(false);

	for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) {
		VERTEX_INDEX kv = subgrid.CubeVertex(subgrid_cube_index, k);
		visited.Set(kv, true);
		vlist2.push_back(kv);
	}

	SHARPISO_SCALAR_GRID scalar_subgrid;
	scalar_subgrid.SetSize(subgrid);

	scalar_subgrid.CopyRegion
		(scalar_grid, region_iv0, scalar_subgrid.AxisSize(), 0);

	while (vlist2.size() != 0) {
		VERTEX_INDEX kv = vlist2.back();
		vlist2.pop_back();

		subgrid.ComputeBoundaryBits(kv, boundary_bits);

		for (DTYPE d = 0; d < dimension; d++) {

			long mask = (1L << (2*d));
			long bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = subgrid.PrevVertex(kv, d);
				if (!visited.Scalar(kv2)) {
          subgrid.ComputeBoundaryBits(kv2, boundary_bits2);
          if (is_adjacent_to_bipolar_edge
              (scalar_subgrid, isovalue, kv2, boundary_bits2)) {

            VERTEX_INDEX iv2 = subgrid.Scalar(kv2);
            if (!gradient_grid.IsMagnitudeGT(iv2, max_small_magnitude)) {
              visited.Set(kv2, true);
              vlist2.push_back(kv2);
            }
          }
				}
			}

			mask = (1L << (2*d+1));
			bit = (boundary_bits & mask);
			if (bit == 0) {
				VERTEX_INDEX kv2 = subgrid.NextVertex(kv, d);
        subgrid.ComputeBoundaryBits(kv2, boundary_bits2);
				if (!visited.Scalar(kv2)) {
          subgrid.ComputeBoundaryBits(kv2, boundary_bits2);
          if (is_adjacent_to_bipolar_edge
              (scalar_subgrid, isovalue, kv2, boundary_bits2)) {

            VERTEX_INDEX iv2 = subgrid.Scalar(kv2);
            if (!gradient_grid.IsMagnitudeGT(iv2, max_small_magnitude)) {
              visited.Set(kv2, true);
              vlist2.push_back(kv2);
            }
					}
				}
			}
		}
	}

	for (VERTEX_INDEX kv = 0; kv < subgrid.NumVertices(); kv++) {
		if (!visited.Scalar(kv)) {
			subgrid.ComputeBoundaryBits(kv, boundary_bits);

			if (is_adjacent_to_visited(visited, kv, boundary_bits) &&
				is_adjacent_to_bipolar_edge
				(scalar_subgrid, isovalue, kv, boundary_bits)) {

					VERTEX_INDEX iv = subgrid.Scalar(kv);
					vertex_list.push_back(iv);
			}
		}
	}

}


// **************************************************
// SORT VERTICES
// **************************************************

/// Sort vertices based on the distance of the isoplane to point pcoord[].
void SHARPISO::sort_vertices_by_isoplane
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const COORD_TYPE pcoord[DIM3],
	VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices)
{
	IJK::ARRAY<COORD_TYPE> isoplane_dist(num_vertices);
	COORD_TYPE vcoord[DIM3];

	for (NUM_TYPE i = 0; i < num_vertices; i++) {
		VERTEX_INDEX iv = vertex_list[i];
		scalar_grid.ComputeCoord(iv, vcoord);

		compute_distance_to_gfield_plane
			(gradient_grid.VectorPtrConst(iv), vcoord, scalar_grid.Scalar(iv),
			pcoord, isovalue, isoplane_dist[i]);
	}

	// insertion sort
	for (NUM_TYPE i = 1; i < num_vertices; i++) {
		NUM_TYPE j = i;
		VERTEX_INDEX jv = vertex_list[j];
		COORD_TYPE jdist = isoplane_dist[j];
		while (j > 0 && jdist < isoplane_dist[j-1]) {
			vertex_list[j] = vertex_list[j-1];
			isoplane_dist[j] = isoplane_dist[j-1];
			j--;
		}
		vertex_list[j] = jv;
		isoplane_dist[j] = jdist;
	}

}

/// Sort vertices based on the distance of the isoplane to point pcoord[].
void SHARPISO::sort_vertices_by_isoplane
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const COORD_TYPE pcoord[DIM3],
	std::vector<VERTEX_INDEX> & vertex_list)
{
	sort_vertices_by_isoplane
		(scalar_grid, gradient_grid, isovalue, pcoord, 
		&(vertex_list[0]), vertex_list.size());
}

/// Sort vertices based on the distance of the isoplane distance to cube center.
void SHARPISO::sort_vertices_by_isoplane_dist2cc
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX cube_index,
	VERTEX_INDEX vertex_list[], const NUM_TYPE num_vertices)
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	// static so not reallocated at each call
	static COORD_TYPE cube_center[DIM3];

	scalar_grid.ComputeCoord(cube_index, cube_center);

	for (DTYPE d = 0; d < DIM3; d++)
	{ cube_center[d] += 0.5; }

	sort_vertices_by_isoplane
		(scalar_grid, gradient_grid, isovalue, cube_center, 
		vertex_list, num_vertices);
}

/// Sort vertices based on the distance of the isoplane distance to cube center.
void SHARPISO::sort_vertices_by_isoplane_dist2cc
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX cube_index,
	std::vector<VERTEX_INDEX> & vertex_list)
{
	sort_vertices_by_isoplane_dist2cc
		(scalar_grid, gradient_grid, isovalue, cube_index,
		&(vertex_list[0]), vertex_list.size());
}

// **************************************************
// SELECTION FUNCTIONS
// **************************************************

// Select gradients at cube vertices.
// Select large gradients which give a level set intersecting the cube.
void SHARPISO::select_cube_gradients_based_on_isoplanes
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const VERTEX_INDEX cube_index,
	const GRADIENT_COORD_TYPE max_small_mag,
	const SCALAR_TYPE isovalue,
	std::vector<COORD_TYPE> & point_coord,
	std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
	std::vector<SCALAR_TYPE> & scalar,
	NUM_TYPE & num_gradients,
	const OFFSET_VOXEL & voxel)
{
	const GRADIENT_COORD_TYPE max_small_mag_squared =
		max_small_mag * max_small_mag;
	IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
	IJK::PROCEDURE_ERROR error("get_large_cube_gradients");

	// Initialize num_gradients
	num_gradients = 0;

	scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

	for (NUM_TYPE k = 0; k < scalar_grid.NumCubeVertices(); k++) {
		VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
		add_selected_gradient
			(scalar_grid, gradient_grid, iv, cube_coord.PtrConst(),
			max_small_mag_squared, isovalue, voxel,
			point_coord, gradient_coord, scalar, num_gradients);
	}

}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   with small gradient magnitude.
/// @pre Array vertex_flag[] is preallocated with size 
///      at least vertex_list.size().
void SHARPISO::deselect_vertices_with_small_gradients
	(const GRADIENT_GRID_BASE & gradient_grid, 
	const VERTEX_INDEX * vertex_list, const NUM_TYPE num_vertices,
	const GRADIENT_COORD_TYPE max_small_mag_squared,
	bool vertex_flag[])
{
	for (NUM_TYPE i = 0; i < num_vertices; i++)
		if (vertex_flag[i]) {

			VERTEX_INDEX iv = vertex_list[i];
			GRADIENT_COORD_TYPE mag_squared = 
				gradient_grid.ComputeMagnitudeSquared(iv);

			if (mag_squared <= max_small_mag_squared) 
			{ vertex_flag[i] = false; }
		}
}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   with small gradient magnitude.
/// std::vector variation.
void SHARPISO::deselect_vertices_with_small_gradients
	(const GRADIENT_GRID_BASE & gradient_grid, 
	const std::vector<VERTEX_INDEX> & vertex_list,
	const GRADIENT_COORD_TYPE max_small_mag_squared,
	bool vertex_flag[])
{
	const NUM_TYPE num_vertices = vertex_list.size();

	if (num_vertices > 0) {

		deselect_vertices_with_small_gradients
			(gradient_grid,  &(vertex_list.front()), num_vertices, 
			max_small_mag_squared, vertex_flag);
	}
}


/// Set to false vertex_flag[i] for any vertex_list[i] 
///   determining an isoplane which does not intersect the cube.
/// @pre Array vertex_flag[] is preallocated with size 
///      at least vertex_list.size().
void SHARPISO::deselect_vertices_based_on_isoplanes
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid, 
	const GRID_COORD_TYPE * cube_coord, const OFFSET_VOXEL & voxel,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX * vertex_list, const NUM_TYPE num_vertices,
	bool vertex_flag[])
{
	typedef SHARPISO_SCALAR_GRID::DIMENSION_TYPE DTYPE;

	// static so not reallocated at each call
	static COORD_TYPE vertex_coord[DIM3];
	static COORD_TYPE coord[DIM3];

	for (NUM_TYPE i = 0; i < num_vertices; i++) {

		if (vertex_flag[i]) {

			VERTEX_INDEX iv = vertex_list[i];

			gradient_grid.ComputeCoord(iv, vertex_coord);
			for (DTYPE d = 0; d < DIM3; d++) {
				coord[d] = vertex_coord[d] - cube_coord[d] + voxel.OffsetFactor();
				coord[d] *= scalar_grid.Spacing(d); 
			}
			const GRADIENT_COORD_TYPE * vertex_gradient_coord =
				gradient_grid.VectorPtrConst(iv);
			SCALAR_TYPE s = scalar_grid.Scalar(iv);

			if (!iso_intersects_cube
				(voxel, coord, vertex_gradient_coord, s, isovalue))
			{ vertex_flag[i] = false; }
		}
	}

}

/// Set to false vertex_flag[i] for any vertex_list[i] 
///   determining an isoplane which does not intersect the cube.
/// std::vector variation.
void SHARPISO::deselect_vertices_based_on_isoplanes
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid, 
	const GRID_COORD_TYPE * cube_coord, const OFFSET_VOXEL & voxel,
	const SCALAR_TYPE isovalue,
	const std::vector<VERTEX_INDEX> & vertex_list,
	bool vertex_flag[])
{
	const NUM_TYPE num_vertices = vertex_list.size();

	if (num_vertices > 0) {

		deselect_vertices_based_on_isoplanes
			(scalar_grid, gradient_grid, cube_coord, voxel, isovalue,
			&(vertex_list.front()), num_vertices, vertex_flag);
	}
}

// **************************************************
// VOXEL
// **************************************************

SHARPISO::VOXEL::VOXEL()
{
	COORD_TYPE min_coord[DIM3], max_coord[DIM3];

	IJK::CUBE_INFO<int,NUM_TYPE>::SetDimension(DIM3);
	InitLocal();
	SetMaxVertexIndex();

	for (int d = 0; d < DIM3; d++) {
		min_coord[d] = 0;
		max_coord[d] = 1;
	}

	SetVertexCoord(min_coord, max_coord);
}

SHARPISO::VOXEL::~VOXEL()
{
	FreeLocal();
}

void SHARPISO::VOXEL::InitLocal()
{
	vertex_coord = new COORD_TYPE[DIM3*this->NumVertices()];
}

void SHARPISO::VOXEL::FreeLocal()
{
	if (vertex_coord != NULL) 
	{ delete [] vertex_coord; }

	vertex_coord = NULL;
}

void SHARPISO::VOXEL::SetMaxVertexIndex()
{
	if (this->NumVertices() > 0) 
	{ max_vertex_index = this->NumVertices()-1; }
	else
	{ max_vertex_index = 0; }
}

// Set voxel vertex coordinates.
void SHARPISO::VOXEL::SetVertexCoord             
	(const COORD_TYPE min_coord[DIM3], const COORD_TYPE max_coord[DIM3])
{
	for (long j = 0; j < this->NumVertices(); j++) {
		long j0 = j;
		for (int d = 0; d < DIM3; d++) {
			int iend = j0 % 2;
			if (iend == 0) 
			{ vertex_coord[j*dimension+d] = min_coord[d]; }
			else
			{ vertex_coord[j*dimension+d] = max_coord[d]; }
			j0 = j0/2;
		}
	}
}

/// Set coordinates of OFFSET_VOXEL
void SHARPISO::OFFSET_VOXEL::SetVertexCoord
	(const COORD_TYPE spacing[DIM3], const COORD_TYPE offset_factor)
{
	COORD_TYPE min_coord[DIM3];
	COORD_TYPE max_coord[DIM3];

	this->offset_factor = offset_factor;

	for (int d = 0; d < DIM3; d++) {
		min_coord[d] = 0;
		max_coord[d] = spacing[d]*(1 + 2*offset_factor);
	}

	VOXEL::SetVertexCoord(min_coord, max_coord);
}


// **************************************************
// GET_GRADIENTS_PARAM
// **************************************************

/// Initialize GET_GRADIENTS_PARAM
void SHARPISO::GET_GRADIENTS_PARAM::Init()
{
	grad_selection_method = UNKNOWN_GRAD_SELECTION_METHOD;
	use_only_cube_gradients = false;
	use_large_neighborhood = false;
	use_zero_grad_boundary = false;
	use_diagonal_neighbors = false;
	use_selected_gradients = true;
	select_based_on_grad_dir = false;
	use_intersected_edge_endpoint_gradients = false;
	use_gradients_determining_edge_intersections = false;
	allow_duplicates = false;
	flag_sort_gradients = false;
	grad_selection_cube_offset = 1;
	max_small_magnitude = 0.001;
	zero_tolerance = 0.0000001;
	max_grad_dist = 1;
}


void SHARPISO::GET_GRADIENTS_PARAM::SetGradSelectionMethod
	(const GRAD_SELECTION_METHOD grad_selection_method)
{
	IJK::PROCEDURE_ERROR error("GET_GRADIENTS_PARAM::SetGradSelectionMethod");

	this->grad_selection_method = grad_selection_method;

	// set defaults
	use_only_cube_gradients = false;
	use_selected_gradients = true;
	use_diagonal_neighbors = false;
	use_intersected_edge_endpoint_gradients = false;
	use_gradients_determining_edge_intersections = false;
	select_based_on_grad_dir = false;
	allow_duplicates = false;

	switch(grad_selection_method) {

	case GRAD_C:
		use_only_cube_gradients = true;
		use_selected_gradients = false;
		break;

	case GRAD_N:
		use_only_cube_gradients = false;
		use_selected_gradients = false;
		break;

	case GRAD_CS:
		use_only_cube_gradients = true;
		use_selected_gradients = true;
		break;

	case GRAD_NS:
		use_only_cube_gradients = false;
		use_selected_gradients = true;
		break;

	case GRAD_XS:
		use_only_cube_gradients = false;
		use_diagonal_neighbors = true;
		use_selected_gradients = true;
		break;

	case GRAD_IE:
		use_only_cube_gradients = true;
		use_selected_gradients = false;
		use_intersected_edge_endpoint_gradients = true;
		break;

	case GRAD_IES:
		use_only_cube_gradients = true;
		use_selected_gradients = true;
		use_intersected_edge_endpoint_gradients = true;
		break;

	case GRAD_IE_DIR:
		use_only_cube_gradients = true;
		use_selected_gradients = true;
		select_based_on_grad_dir = true;
		use_intersected_edge_endpoint_gradients = true;
		break;

	case GRAD_CD:
		use_only_cube_gradients = true;
		use_selected_gradients = false;
		use_intersected_edge_endpoint_gradients = true;
		use_intersected_edge_endpoint_gradients = false;
		use_gradients_determining_edge_intersections = true;
		break;

	case GRAD_CD_DUP:
		use_only_cube_gradients = true;
		use_selected_gradients = false;
		use_intersected_edge_endpoint_gradients = true;
		use_intersected_edge_endpoint_gradients = false;
		use_gradients_determining_edge_intersections = true;
		allow_duplicates = true;
		break;

	case GRAD_NIE:
		use_only_cube_gradients = false;
		use_selected_gradients = false;
		use_intersected_edge_endpoint_gradients = true;
		break;

  case GRAD_3x3x3:
	case GRAD_NIES:
		use_only_cube_gradients = false;
		use_selected_gradients = true;
		use_intersected_edge_endpoint_gradients = true;
		break;

  case GRAD_5x5x5:
  case GRAD_7x7x7:
  case GRAD_9x9x9:
	case GRAD_BIES:
		use_only_cube_gradients = false;
		use_selected_gradients = true;
		use_intersected_edge_endpoint_gradients = true;
		use_zero_grad_boundary = true;
		use_large_neighborhood = true;
		use_new_version = true;
		break;

	case GRAD_BIES_OLD:
		use_only_cube_gradients = false;
		use_selected_gradients = true;
		use_intersected_edge_endpoint_gradients = true;
		use_zero_grad_boundary = true;
		use_large_neighborhood = true;
		use_new_version = false;
		break;

	case GRAD_EDGEI_INTERPOLATE:
		break;

	case GRAD_EDGEI_GRAD:
		break;

	default:
		error("Programming error.  Unknown gradient selection type.");
		throw error;
	};

  // Set max_grad_dist, if necessary.
  switch(grad_selection_method) {

  case GRAD_5x5x5:
    max_grad_dist = 2;
    break;

  case GRAD_7x7x7:
    max_grad_dist = 3;
    break;

  case GRAD_9x9x9:
    max_grad_dist = 4;
    break;
  }
}

// **************************************************
// MAP GRAD_SELECTION_METHOD TO/FROM C++ string
// **************************************************

GRAD_SELECTION_METHOD SHARPISO::get_grad_selection_method
	(const std::string & s)
{
  if (s == "grad3") { return(GRAD_3x3x3); }
  else if (s == "grad5") { return(GRAD_5x5x5); }
  else if (s == "grad7") { return(GRAD_7x7x7); }
  else if (s == "grad9") { return(GRAD_9x9x9); }
	else if (s == "gradC") { return(GRAD_C); }
	else if (s == "gradN") { return(GRAD_N); }
	else if (s == "gradCS") { return(GRAD_CS); }
	else if (s == "gradNS") { return(GRAD_NS); }
	else if (s == "gradXS") { return(GRAD_XS); }
	else if (s == "gradIE") { return(GRAD_IE); }
	else if (s == "gradIES") { return(GRAD_IES); }
	else if (s == "gradIEDir") { return(GRAD_IE_DIR); }
	else if (s == "gradCD") { return(GRAD_CD); }
	else if (s == "gradCDdup") { return(GRAD_CD_DUP); }
	else if (s == "gradNIE") { return(GRAD_NIE); }
	else if (s == "gradNIES") { return(GRAD_NIES); }
	else if (s == "gradBIES") { return(GRAD_BIES); }
	else if (s == "gradBIES_old") { return(GRAD_BIES_OLD); }
	else if (s == "gradEIinterp") { return(GRAD_EDGEI_INTERPOLATE); }
	else if (s == "gradEIgrad") { return(GRAD_EDGEI_GRAD); }
	else if (s == "gradES") { return(GRAD_EDGEI_INTERPOLATE); }
	else if (s == "gradEC") { return(GRAD_EDGEI_GRAD); }
	else { return(UNKNOWN_GRAD_SELECTION_METHOD); }
}

void SHARPISO::get_grad_selection_string
	(const GRAD_SELECTION_METHOD grad_sel, std::string & s)
{
	if (grad_sel == GRAD_3x3x3) { s = "grad3"; }
	else if (grad_sel == GRAD_5x5x5) { s = "grad5"; }
	else if (grad_sel == GRAD_7x7x7) { s = "grad7"; }
	else if (grad_sel == GRAD_9x9x9) { s = "grad9"; }
  else if (grad_sel == GRAD_C) { s = "gradCD"; }
	else if (grad_sel == GRAD_CS) { s = "gradCS"; }
	else if (grad_sel == GRAD_NS) { s = "gradNS"; }
	else if (grad_sel == GRAD_XS) { s = "gradXS"; }
	else if (grad_sel == GRAD_IE) { s = "gradIE"; }
	else if (grad_sel == GRAD_IES) { s = "gradIES"; }
	else if (grad_sel == GRAD_IE_DIR) { s = "gradIEDir"; }
	else if (grad_sel == GRAD_CD) { s = "gradCD"; }
	else if (grad_sel == GRAD_CD_DUP) { s = "gradCDdup"; }
	else if (grad_sel == GRAD_NIE) { s = "gradNIE"; }
	else if (grad_sel == GRAD_NIES) { s = "gradNIES"; }
	else if (grad_sel == GRAD_BIES) { s = "gradBIES"; }
	else if (grad_sel == GRAD_BIES_OLD) { s = "gradBIES_old"; }
	else if (grad_sel == GRAD_EDGEI_INTERPOLATE) { s = "gradEIinterp"; }
	else if (grad_sel == GRAD_EDGEI_GRAD) { s = "gradEIgrad"; }
	else { s = "gradUnknown"; };
}

/// \file sharpiso_feature.cxx
/// Compute sharp isosurface vertices and edges.
/// Version 0.1.0

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


#include "sharpiso_feature.h"
#include "sharpiso_svd.h"
#include "sharpiso_intersect.h"
#include "sharpiso_findIntersect.h"
#include "sh_point_find.h"

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_scalar.txx"

// Local functions.
bool is_dist_to_cube_le
(const COORD_TYPE coord[DIM3], const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE max_dist);
void snap_to_cube
(const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE snap_dist, COORD_TYPE coord[DIM3]);
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE point_coord[DIM3]);
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE point_coord[DIM3],
		VERTEX_INDEX & conflicting_cube);
void diff_coord
(const GRID_COORD_TYPE coordA[DIM3], const GRID_COORD_TYPE coordB[DIM3],
		NUM_TYPE & num_diff, VERTEX_INDEX & icoord);


// **************************************************
// SVD ROUTINES TO COMPUTE SHARP VERTEX/EDGE
// **************************************************

/// Compute sharp isosurface vertex using singular valued decomposition.
void SHARPISO::svd_compute_sharp_vertex_for_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const OFFSET_CUBE_111 & cube_111,
		COORD_TYPE sharp_coord[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{

	if (sharpiso_param.use_lindstrom) {

		svd_compute_sharp_vertex_for_cube_lindstrom
		(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
				cube_111, sharp_coord, eigenvalues, num_large_eigenvalues,
				svd_info);
	}
	else {
		svd_compute_sharp_vertex_for_cube_lc_intersection
		(scalar_grid, gradient_grid, cube_index, isovalue, sharpiso_param,
				cube_111, sharp_coord, eigenvalues, num_large_eigenvalues,
				svd_info);
	}
}


/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use Lindstrom's formula.
void SHARPISO::svd_compute_sharp_vertex_for_cube_lindstrom
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const OFFSET_CUBE_111 & cube_111,
		COORD_TYPE sharp_coord[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{
	const EIGENVALUE_TYPE max_small_eigenvalue =
			sharpiso_param.max_small_eigenvalue;

	NUM_TYPE num_gradients = 0;
	std::vector<COORD_TYPE> point_coord;
	std::vector<GRADIENT_COORD_TYPE> gradient_coord;
	std::vector<SCALAR_TYPE> scalar;

	// Compute coord of the cube.
	GRID_COORD_TYPE cube_coord[DIM3];
	scalar_grid.ComputeCoord(cube_index, cube_coord);

	get_gradients
	(scalar_grid, gradient_grid, cube_index, isovalue,
			sharpiso_param, cube_111, sharpiso_param.flag_sort_gradients,
			point_coord, gradient_coord, scalar, num_gradients);

	svd_info.location = LOC_SVD;
	svd_info.flag_conflict = false;
	svd_info.flag_Linf_iso_vertex_location = false;

	// svd_calculate_sharpiso vertex using lindstrom
	COORD_TYPE central_point[DIM3];
	if (sharpiso_param.flag_dist2centroid) {
		compute_edgeI_centroid
		(scalar_grid, gradient_grid, isovalue, cube_index,
				sharpiso_param.use_sharp_edgeI, central_point);
	}
	else {
		scalar_grid.ComputeCubeCenterCoord(cube_index, central_point);
	}
	IJK::copy_coord(DIM3, central_point, svd_info.central_point);

	/// use the sharp version with the garlnd heckbert way of storing normals
	if (sharpiso_param.use_lindstrom_fast){

		svd_calculate_sharpiso_vertex_using_lindstrom_fast
      (num_gradients, max_small_eigenvalue,isovalue, &(scalar[0]), 
       &(point_coord[0]),	&(gradient_coord[0]), 
       num_large_eigenvalues, eigenvalues, central_point, sharp_coord);
	}
	else{

    svd_calculate_sharpiso_vertex_using_lindstrom
      (sharpiso_param.use_lindstrom2, &(point_coord[0]),
       &(gradient_coord[0]), &(scalar[0]),
       num_gradients, isovalue, max_small_eigenvalue,
       num_large_eigenvalues, eigenvalues, central_point, sharp_coord);
	}

	postprocess_isovert_location
	(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
			sharpiso_param, sharp_coord, svd_info);
}

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
/// when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void local_svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const std::vector<COORD_TYPE> & point_coord,
		const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		const std::vector<SCALAR_TYPE> & scalar,
		const NUM_TYPE num_gradients,
		COORD_TYPE sharp_coord[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info);

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
///    when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void SHARPISO::svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const OFFSET_CUBE_111 & cube_111,
		COORD_TYPE sharp_coord[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{
	const EIGENVALUE_TYPE max_small_eigenvalue =
			sharpiso_param.max_small_eigenvalue;
	const COORD_TYPE max_dist = sharpiso_param.max_dist;

	NUM_TYPE num_gradients = 0;
	std::vector<COORD_TYPE> point_coord;
	std::vector<GRADIENT_COORD_TYPE> gradient_coord;
	std::vector<SCALAR_TYPE> scalar;
	GRADIENT_COORD_TYPE line_direction[DIM3];

	// Compute coord of the cube.
	GRID_COORD_TYPE cube_coord[DIM3];
	scalar_grid.ComputeCoord(cube_index, cube_coord);

	get_gradients
	(scalar_grid, gradient_grid, cube_index, isovalue,
			sharpiso_param, cube_111, sharpiso_param.flag_sort_gradients,
			point_coord, gradient_coord, scalar, num_gradients);

	svd_info.location = LOC_SVD;
	svd_info.flag_conflict = false;
	svd_info.flag_Linf_iso_vertex_location = false;

	local_svd_compute_sharp_vertex_for_cube_lc_intersection
	(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
			sharpiso_param, point_coord, gradient_coord, scalar, num_gradients,
			sharp_coord, eigenvalues, num_large_eigenvalues, svd_info);

	postprocess_isovert_location
	(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
			sharpiso_param, sharp_coord, svd_info);
}

/// Compute sharp isosurface vertex using singular valued decomposition.
/// Use line-cube intersection to compute sharp isosurface vertex
///    when number of eigenvalues is 2.
/// Use centroid to compute isosurface vertex when number of eigenvalues is 1.
void local_svd_compute_sharp_vertex_for_cube_lc_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const std::vector<COORD_TYPE> & point_coord,
		const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		const std::vector<SCALAR_TYPE> & scalar,
		const NUM_TYPE num_gradients,
		COORD_TYPE sharp_coord[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{
	const EIGENVALUE_TYPE max_eigen = sharpiso_param.max_small_eigenvalue;
	const COORD_TYPE max_dist = sharpiso_param.max_dist;
	GRADIENT_COORD_TYPE line_direction[DIM3];

	svd_calculate_sharpiso_vertex_unit_normals
	(&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
			num_gradients, isovalue, max_eigen, num_large_eigenvalues,
			eigenvalues, sharp_coord, line_direction);

	if (num_large_eigenvalues == 2) {
		bool flag_conflict;
		COORD_TYPE line_origin[DIM3];
		// if there are 2 sing vals then the coord acts as the line origin
		IJK::copy_coord_3D(sharp_coord, line_origin);

		compute_vertex_on_line
		(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
				sharpiso_param, line_origin, line_direction, sharp_coord, svd_info);
		svd_info.SetRayInfo(line_origin, line_direction, sharp_coord);
	}
	else if (num_large_eigenvalues  == 1) {
		compute_edgeI_centroid
		(scalar_grid, gradient_grid, isovalue, cube_index,
				sharpiso_param.use_sharp_edgeI, sharp_coord);
		svd_info.location = CENTROID;
	}
	else if (num_large_eigenvalues == 0 ) {
		COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};
		IJK::add_coord_3D(cube_coord, cube_center, sharp_coord);
		svd_info.location = CUBE_CENTER;
	}

}


/// Compute sharp isosurface vertex on the line
/// @param flag_conflict True if sharp_coord conflicts with other cube.
void SHARPISO::compute_vertex_on_line
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		const COORD_TYPE line_origin[DIM3],
		const COORD_TYPE line_direction[DIM3],
		COORD_TYPE sharp_coord[DIM3],
		SVD_INFO & svd_info)
{
	const SIGNED_COORD_TYPE max_dist = sharpiso_param.max_dist;
	const GRADIENT_COORD_TYPE zero_tolerance = sharpiso_param.zero_tolerance;
	COORD_TYPE central_point[DIM3];
	VERTEX_INDEX conflicting_cube;
	static GRID_COORD_TYPE conflicting_cube_coord[DIM3];
	static COORD_TYPE Linf_coord[DIM3];
	VERTEX_INDEX icoord;
	NUM_TYPE num_diff;

	if (sharpiso_param.flag_dist2centroid) {
		compute_edgeI_centroid
		(scalar_grid, gradient_grid, isovalue, cube_index,
				sharpiso_param.use_sharp_edgeI, central_point);
	}
	else {
		for (int d = 0; d < DIM3; d++)
		{ central_point[d] = cube_coord[d] + 0.5; }
	}

	// Compute the closest point (L2 distance) on the line to central_point.
	compute_closest_point_on_line
	(central_point, line_origin, line_direction, zero_tolerance, sharp_coord);
	IJK::copy_coord(DIM3, central_point, svd_info.central_point);

	if (is_dist_to_cube_le(sharp_coord, cube_coord, max_dist)) {

		snap_to_cube(cube_coord, sharpiso_param.snap_dist, sharp_coord);
		if (check_conflict(scalar_grid, isovalue, cube_coord,
				sharp_coord, conflicting_cube)) {

			if (sharpiso_param.use_Linf_dist) {
				scalar_grid.ComputeCoord(conflicting_cube, conflicting_cube_coord);

				diff_coord(cube_coord, conflicting_cube_coord, num_diff, icoord);
				if (num_diff == 1 &&
						abs(line_direction[icoord]) >= sharpiso_param.max_small_grad_coord_Linf) {

					// Use Linf distance instead of L1 distance.
					compute_closest_point_on_line_linf
					(central_point, line_origin, line_direction, zero_tolerance, Linf_coord);

					snap_to_cube(cube_coord, sharpiso_param.snap_dist, Linf_coord);
					if (!check_conflict(scalar_grid, isovalue, cube_coord, Linf_coord)) {
						// No conflict.  Use Linf coord.
						IJK::copy_coord_3D(Linf_coord, sharp_coord);
						svd_info.flag_Linf_iso_vertex_location = true;
					}
				}
			}
		}
	}

}


/// Compute sharp isosurface vertex using edge-isosurface intersections.
/// Approximate gradients using linear interpolation on the edges.
void SHARPISO::svd_compute_sharp_vertex_edgeI_interpolate_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{
	const EIGENVALUE_TYPE max_small_eigenvalue =
			sharpiso_param.max_small_eigenvalue;
	GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3];
	SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D];
	NUM_TYPE num_gradients(0);

	get_cube_gradients
	(scalar_grid, gradient_grid, cube_index,
			gradient_coord, scalar);

	// Ray Direction to calculate intersection if there are 2 singular values.
	GRADIENT_COORD_TYPE ray_direction[3]={0.0};

	//tobe added as a parameters
	bool use_cmplx_interp = false;

	bool cube_create = shFindPoint
			(&(gradient_coord[0]), &(scalar[0]), isovalue, use_cmplx_interp,
					max_small_eigenvalue, eigenvalues, num_large_eigenvalues,
					svd_info, coord);

	COORD_TYPE cube_coord[DIM3];
	COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};

	scalar_grid.ComputeCoord(cube_index, cube_coord);
	//check if cube creation failed.
	if(cube_create){
		IJK::add_coord(DIM3, cube_coord, coord, coord);
		svd_info.location = LOC_SVD;
	}
	else{
		IJK::add_coord(DIM3, cube_coord, cube_center, coord);
		svd_info.location = CUBE_CENTER;
	}

	if (sharpiso_param.flag_round)
	{ IJK::round_coord
		(sharpiso_param.round_denominator, DIM3, coord, coord);
	}

}

/// Compute sharp isosurface vertex using edge-isosurface intersections.
/// Use sharp formula for computing gradient at intersection.
void SHARPISO::svd_compute_sharp_vertex_edgeI_sharp_gradient
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		COORD_TYPE coord[DIM3], EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues,
		SVD_INFO & svd_info)
{
	const EIGENVALUE_TYPE max_small_eigenvalue =
			sharpiso_param.max_small_eigenvalue;
	NUM_TYPE num_gradients(0);
	GRADIENT_COORD_TYPE gradient_coord[NUM_CUBE_VERTICES3D*DIM3];
	SCALAR_TYPE scalar[NUM_CUBE_VERTICES3D];

	svd_info.flag_conflict = false;

	get_cube_gradients
	(scalar_grid, gradient_grid, cube_index, gradient_coord, scalar);

	// Ray Direction to calculate intersection if there are 2 singular values.
	GRADIENT_COORD_TYPE ray_direction[3]={0.0};

	//tobe added as a parameters
	bool use_cmplx_interp = true;
	bool cube_create = shFindPoint
			(&(gradient_coord[0]), &(scalar[0]), isovalue, use_cmplx_interp,
					max_small_eigenvalue, eigenvalues, num_large_eigenvalues,
					svd_info, coord);


	COORD_TYPE cube_center[DIM3] = {0.5,0.5,0.5};
	GRID_COORD_TYPE cube_coord[DIM3];
	scalar_grid.ComputeCoord(cube_index, cube_coord);

	//check if cube creation failed.
	if(cube_create) {
		IJK::add_coord(DIM3, cube_coord, coord, coord);
		svd_info.location = LOC_SVD;
	}
	else {
		IJK::add_coord(DIM3, cube_coord, cube_center, coord);
		svd_info.location = CUBE_CENTER;
	}

	if (is_dist_to_cube_le(coord, cube_coord, sharpiso_param.max_dist)) {

		if (!sharpiso_param.flag_allow_conflict) {
			VERTEX_INDEX conflicting_cube;

			snap_to_cube(cube_coord, sharpiso_param.snap_dist, coord);
			bool flag_conflict =
					check_conflict(scalar_grid, isovalue, cube_coord, coord,
							conflicting_cube);

			if (flag_conflict) {
				svd_info.flag_conflict = true;
				svd_info.cube_containing_coord = conflicting_cube;
				process_conflict
				(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
						sharpiso_param, coord, svd_info);
			}
		}
	}
	else {
		process_far_point(scalar_grid, cube_index, cube_coord, isovalue,
				sharpiso_param, coord, svd_info);
	}


	// Clamp point to cube + max_dist.
	clamp_point(sharpiso_param.max_dist, cube_coord, coord);

	if (sharpiso_param.flag_round)
	{ IJK::round_coord
		(sharpiso_param.round_denominator, DIM3, coord, coord);
	}
}


/// Compute sharp isosurface vertex near facet
///    using singular valued decomposition.
/// @param facet_v0 Index of primary (lowest-left) facet vertex.
/// @param orth_dir Direction orthogonal to facet.
/// @pre Facet is in interior of grid.
/// @param[out] sharp_vertex_location  -1, 0, or 1.
///    -1: Sharp vertex is below/left of facet.
///    0: Vertex not sharp or relative location undetermined.
void SHARPISO::svd_compute_sharp_vertex_near_facet
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX facet_v0,
		const NUM_TYPE facet_orth_dir,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		NUM_TYPE & sharp_vertex_location,
		COORD_TYPE sharp_coord[DIM3],
		GRADIENT_COORD_TYPE line_direction[DIM3],
		EIGENVALUE_TYPE eigenvalues[DIM3],
		NUM_TYPE & num_large_eigenvalues)
{
	const GRADIENT_COORD_TYPE max_small_mag =
			sharpiso_param.max_small_magnitude;
	const EIGENVALUE_TYPE max_eigen = sharpiso_param.max_small_eigenvalue;
	const COORD_TYPE snap_dist = sharpiso_param.snap_dist;
	NUM_TYPE num_gradients = 0;
	std::vector<COORD_TYPE> point_coord;
	std::vector<GRADIENT_COORD_TYPE> gradient_coord;
	std::vector<SCALAR_TYPE> scalar;
	static COORD_TYPE v0_coord[DIM3];
	static COORD_TYPE end[2][DIM3];
	static COORD_TYPE endc[2];

	// Initialize
	sharp_vertex_location = 0;

	get_two_cube_gradients
	(scalar_grid, gradient_grid, facet_v0, facet_orth_dir,
			isovalue, sharpiso_param,
			point_coord, gradient_coord, scalar, num_gradients);

	svd_calculate_sharpiso_vertex_unit_normals
	(&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
			num_gradients, isovalue, max_eigen, num_large_eigenvalues,
			eigenvalues, sharp_coord, line_direction);

	if (num_large_eigenvalues > 1) {

		scalar_grid.ComputeCoord(facet_v0, v0_coord);

		if (num_large_eigenvalues == 3) {

			if (sharp_coord[facet_orth_dir] > v0_coord[facet_orth_dir]+snap_dist)
			{ sharp_vertex_location = 1; }
			else if (sharp_coord[facet_orth_dir]+snap_dist <
					v0_coord[facet_orth_dir]+snap_dist)
			{ sharp_vertex_location = -1; }
		}
		else {
			// case: num_large_eigenvalues == 2

					// Large number for long/narrow cylinder.
			const COORD_TYPE CYLINDER_LENGTH = 100;
			bool flag_intersect;

			intersect_line_square_cylinder
			(v0_coord, facet_orth_dir, 1, CYLINDER_LENGTH,
					sharp_coord, line_direction, max_small_mag, end, flag_intersect);

			if (flag_intersect) {
				endc[0] = end[0][facet_orth_dir];
				endc[1] = end[1][facet_orth_dir];

				if (endc[0] > endc[1]) { std::swap(endc[0], endc[1]); }

				if (endc[0] > v0_coord[facet_orth_dir])
				{ sharp_vertex_location = 1; }
				else if (endc[1] < v0_coord[facet_orth_dir])
				{ sharp_vertex_location = -1; }
				// else default (sharp_vertex_location = 0).
			}
			// else default (sharp_vertex_location = 0).
		}
	}
	// else default (sharp_vertex_location = 0).

}

/// Compute sharp isosurface vertex near facet
///    using singular valued decomposition.
/// Short parameter list.
void SHARPISO::svd_compute_sharp_vertex_near_facet
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX facet_v0,
		const NUM_TYPE facet_orth_dir,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		NUM_TYPE & sharp_vertex_location)
{
	NUM_TYPE num_large_eigenvalues;
	COORD_TYPE sharp_coord[DIM3];
	GRADIENT_COORD_TYPE line_direction[DIM3];
	EIGENVALUE_TYPE eigenvalues[DIM3];

	svd_compute_sharp_vertex_near_facet
	(scalar_grid, gradient_grid, facet_v0, facet_orth_dir,
			isovalue, sharpiso_param, sharp_vertex_location,
			sharp_coord, line_direction, eigenvalues, num_large_eigenvalues);
}


/// Compute sharp isosurface vertex using singular valued decomposition.
void SHARPISO::svd_compute_sharp_vertex_for_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
 const SHARPISO_EDGE_INDEX_GRID & edge_index,
 const VERTEX_INDEX cube_index,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & sharpiso_param,
 const OFFSET_CUBE_111 & cube_111,
 COORD_TYPE sharp_coord[DIM3],
 EIGENVALUE_TYPE eigenvalues[DIM3],
 NUM_TYPE & num_large_eigenvalues,
 SVD_INFO & svd_info)
{
 	const EIGENVALUE_TYPE max_small_eigenvalue =
			sharpiso_param.max_small_eigenvalue;
	const COORD_TYPE max_dist = sharpiso_param.max_dist;
  IJK::PROCEDURE_ERROR error("svd_compute_sharp_vertex_for_cube");

	NUM_TYPE num_gradients = 0;
	std::vector<COORD_TYPE> point_coord;
	std::vector<GRADIENT_COORD_TYPE> gradient_coord;
	std::vector<SCALAR_TYPE> scalar;

  if (!sharpiso_param.use_lindstrom) {
    error.AddMessage("Programming error.  Normal based computation only implemented with Lindstrom.");
    throw error;
  }

	// Compute coord of the cube.
	GRID_COORD_TYPE cube_coord[DIM3];
	scalar_grid.ComputeCoord(cube_index, cube_coord);

	get_gradients_from_list
    (scalar_grid, edgeI_coord, edgeI_normal_coord, edge_index, 
     cube_index, isovalue, sharpiso_param, 
     point_coord, gradient_coord, scalar, num_gradients);

	svd_info.location = LOC_SVD;
	svd_info.flag_conflict = false;
	svd_info.flag_Linf_iso_vertex_location = false;

	// svd_calculate_sharpiso vertex using lindstrom
	COORD_TYPE central_point[DIM3];
	if (sharpiso_param.flag_dist2centroid) {
		compute_edgeI_centroid
      (scalar_grid, edgeI_coord, edge_index, isovalue, cube_index, 
       central_point);
	}
	else {
		scalar_grid.ComputeCubeCenterCoord(cube_index, central_point);
	}
	IJK::copy_coord(DIM3, central_point, svd_info.central_point);

	/// use the sharp version with the garlnd heckbert way of storing normals
	if (sharpiso_param.use_lindstrom_fast){

		svd_calculate_sharpiso_vertex_using_lindstrom_fast
      (num_gradients, max_small_eigenvalue, isovalue, &(scalar[0]), 
       &(point_coord[0]),	&(gradient_coord[0]), 
       num_large_eigenvalues, eigenvalues, central_point, sharp_coord);
	}
	else{

    svd_calculate_sharpiso_vertex_using_lindstrom
      (sharpiso_param.use_lindstrom2, &(point_coord[0]),
       &(gradient_coord[0]), &(scalar[0]),
       num_gradients, isovalue, max_small_eigenvalue,
       num_large_eigenvalues, eigenvalues, central_point, sharp_coord);
	}

	svd_info.cube_containing_coord = cube_index;

	if (!is_dist_to_cube_le(sharp_coord, cube_coord, max_dist)) {
		process_far_point(scalar_grid, cube_index, cube_coord, isovalue,
                      sharpiso_param, sharp_coord, svd_info);
  }

}

// **************************************************
// SUBGRID ROUTINES TO COMPUTE SHARP VERTEX/EDGE
// **************************************************


/// Compute sharp vertex.
/// Use subgrid sampling to locate isosurface vertex on sharp edge/corner.
void SHARPISO::subgrid_compute_sharp_vertex_in_cube
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const SCALAR_TYPE isovalue,
		const GET_GRADIENTS_PARAM & get_gradients_param,
		const OFFSET_CUBE_111 & cube_111,
		const NUM_TYPE subgrid_axis_size,
		COORD_TYPE sharp_coord[DIM3],
		SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
	NUM_TYPE num_gradients = 0;
	std::vector<COORD_TYPE> point_coord;
	std::vector<GRADIENT_COORD_TYPE> gradient_coord;
	std::vector<SCALAR_TYPE> scalar;

	get_gradients
	(scalar_grid, gradient_grid, cube_index, isovalue,
			get_gradients_param, cube_111, false,
			point_coord, gradient_coord, scalar, num_gradients);

	IJK::ARRAY<GRID_COORD_TYPE> cube_coord(DIM3);
	scalar_grid.ComputeCoord(cube_index, cube_coord.Ptr());

	subgrid_calculate_iso_vertex_in_cube
	(point_coord, gradient_coord, scalar,
			num_gradients, cube_coord.PtrConst(), isovalue, subgrid_axis_size,
			sharp_coord, scalar_stdev, max_abs_scalar_error);
}

/// Calculate isosurface vertex using regular subgrid of the cube.
void SHARPISO::subgrid_calculate_iso_vertex_in_cube
(const COORD_TYPE * point_coord, const GRADIENT_COORD_TYPE * gradient_coord,
		const SCALAR_TYPE * scalar, const NUM_TYPE num_points,
		const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
		const NUM_TYPE subgrid_axis_size,
		COORD_TYPE sharp_coord[DIM3],
		SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
	COORD_TYPE coord[DIM3];
	COORD_TYPE center_coord[DIM3];
	SCALAR_TYPE sharp_stdev_squared(0);
	SCALAR_TYPE sharp_max_abs_error(0);
	COORD_TYPE sharp_dist2center_squared(0);
	IJK::PROCEDURE_ERROR error("subgrid_calculate_iso_vertex_in_cube");

	if (subgrid_axis_size < 1) {
		error.AddMessage
		("Programming error. Subgrid axis size must be at least 1.");
		error.AddMessage("  Subgrid axis size = ", subgrid_axis_size, ".");
		throw error;
	}

	// Compute center coordinate
	for (NUM_TYPE d = 0; d < DIM3; d++)
	{ center_coord[d] = cube_coord[d] + 0.5; }

	const COORD_TYPE h = 1.0/(subgrid_axis_size+1);

	bool flag_set_sharp(false);
	for (NUM_TYPE ix = 0; ix < subgrid_axis_size; ix++) {
		coord[0] = cube_coord[0] + (ix+1)*h;
		for (NUM_TYPE iy = 0; iy < subgrid_axis_size; iy++) {
			coord[1] = cube_coord[1] + (iy+1)*h;
			for (NUM_TYPE iz = 0; iz < subgrid_axis_size; iz++) {

				coord[2] = cube_coord[2] + (iz+1)*h;
				SCALAR_TYPE s, stdev_squared, max_abs_error;

				compute_gradient_based_scalar_diff
				(coord, isovalue, point_coord, gradient_coord, scalar, num_points,
						stdev_squared, max_abs_error);

				if (!flag_set_sharp ||
						stdev_squared < sharp_stdev_squared) {

					IJK::copy_coord(DIM3, coord, sharp_coord);
					sharp_stdev_squared = stdev_squared;
					sharp_max_abs_error = max_abs_error;
					IJK::compute_distance_squared
					(DIM3, coord, center_coord, sharp_dist2center_squared);
					flag_set_sharp = true;
				}
				else if (stdev_squared == sharp_stdev_squared) {
					COORD_TYPE dist2center_squared;
					IJK::compute_distance_squared
					(DIM3, coord, center_coord, dist2center_squared);
					if (dist2center_squared < sharp_dist2center_squared) {
						IJK::copy_coord(DIM3, coord, sharp_coord);
						sharp_stdev_squared = stdev_squared;
						sharp_max_abs_error = max_abs_error;
						sharp_dist2center_squared = dist2center_squared;
						flag_set_sharp = true;
					}
				}

			}
		}
	}

	scalar_stdev = std::sqrt(sharp_stdev_squared);
	max_abs_scalar_error = sharp_max_abs_error;
}


// **************************************************
// COMPUTE ISO VERTEX AT CENTROID
// **************************************************

/// Compute centroid of intersections of isosurface and grid edges.
/// @param use_sharp_edgeI If true, use sharp formula for
///          intersection of isosurface and grid edges.
void SHARPISO::compute_edgeI_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
		const bool use_sharp_edgeI, COORD_TYPE * coord)
{
	if (use_sharp_edgeI) {
		compute_edgeI_sharp_centroid
      (scalar_grid, gradient_grid, isovalue, iv, coord);
	}
	else {
		compute_edgeI_centroid(scalar_grid, isovalue, iv, coord);
	}

}

/// Compute centroid of intersections of isosurface and grid edges
void SHARPISO::compute_edgeI_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
		COORD_TYPE * coord)
{
	const int dimension = scalar_grid.Dimension();
	COORD_TYPE vcoord[dimension];
	COORD_TYPE coord0[dimension];
	COORD_TYPE coord1[dimension];
	COORD_TYPE coord2[dimension];

	int num_intersected_edges = 0;
	IJK::set_coord(dimension, 0.0, vcoord);

	for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
		for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
			VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
			VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

			if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {
				SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
				SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

				scalar_grid.ComputeCoord(iend0, coord0);
				scalar_grid.ComputeCoord(iend1, coord1);

				IJK::linear_interpolate_coord
				(dimension, s0, coord0, s1, coord1, isovalue, coord2);

				IJK::add_coord(dimension, vcoord, coord2, vcoord);

				num_intersected_edges++;
			}
		}

	if (num_intersected_edges > 0) {
		IJK::multiply_coord
		(dimension, 1.0/num_intersected_edges, vcoord, vcoord);
	}
	else {
		scalar_grid.ComputeCoord(iv, vcoord);
		for (int d = 0; d < dimension; d++)
		{ vcoord[d] += 0.5; };
	}

	IJK::copy_coord(dimension, vcoord, coord);
}

// Compute centroid of intersections of isosurface and grid edges.
//   Edge-isosurface intersections are given in an input list.
void SHARPISO::compute_edgeI_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const SHARPISO_EDGE_INDEX_GRID & edge_index,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX cube_index,
 COORD_TYPE * coord)
{
	const int dimension = scalar_grid.Dimension();
	COORD_TYPE vcoord[dimension];
  IJK::PROCEDURE_ERROR error("compute_edgeI_centroid");

	int num_intersected_edges = 0;
	IJK::set_coord(dimension, 0.0, vcoord);

	for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
		for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
			VERTEX_INDEX iend0 = scalar_grid.FacetVertex(cube_index, edge_dir, k);
			VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

			if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

        INDEX_DIFF_TYPE j = edge_index.Vector(iend0, edge_dir);

        if (j < 0) {
          error.AddMessage
            ("Error.  Missing edge-isosurface intersection for edge (",
             iend0, ",", iend1, ").");
          throw error;
        }

				IJK::add_coord(dimension, vcoord, &(edgeI_coord[j*DIM3]), vcoord);

				num_intersected_edges++;
			}
		}

	if (num_intersected_edges > 0) {
		IJK::multiply_coord
		(dimension, 1.0/num_intersected_edges, vcoord, vcoord);
	}
	else {
		scalar_grid.ComputeCoord(cube_index, vcoord);
		for (int d = 0; d < dimension; d++)
		{ vcoord[d] += 0.5; };
	}

	IJK::copy_coord(dimension, vcoord, coord);
}

/// Compute centroid of intersections of isosurface and grid edges.
/// Use sharp formula for computing intersection of isosurface and grid edge.
void SHARPISO::compute_edgeI_sharp_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
		COORD_TYPE * coord)
{
	const int dimension = scalar_grid.Dimension();
	GRID_COORD_TYPE grid_coord[dimension];
	COORD_TYPE vcoord[dimension];
	COORD_TYPE coord2[dimension];
	int num_intersected_edges = 0;

	IJK::set_coord(dimension, 0.0, vcoord);

	for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
		for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
			VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
			VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

			if (is_gt_min_le_max(scalar_grid, iend0, iend1, isovalue)) {

				compute_isosurface_grid_edge_intersection
          (scalar_grid, gradient_grid, isovalue,
           iend0, iend1, edge_dir, coord2);

				IJK::add_coord(dimension, vcoord, coord2, vcoord);

				num_intersected_edges++;
			}
		}

	if (num_intersected_edges > 0) {
		IJK::multiply_coord
		(dimension, 1.0/num_intersected_edges, vcoord, vcoord);
	}
	else {
		scalar_grid.ComputeCoord(iv, vcoord);
		for (int d = 0; d < dimension; d++)
		{ vcoord[d] += 0.5; };
	}

	IJK::copy_coord(dimension, vcoord, coord);
}



/// Calculate isosurface vertex using regular subgrid of the cube.
void SHARPISO::subgrid_calculate_iso_vertex_in_cube
(const std::vector<COORD_TYPE> & point_coord,
		const std::vector<GRADIENT_COORD_TYPE> & gradient_coord,
		const std::vector<SCALAR_TYPE> & scalar,
		const NUM_TYPE num_points,
		const GRID_COORD_TYPE cube_coord[DIM3], const SCALAR_TYPE isovalue,
		const NUM_TYPE subgrid_axis_size,
		COORD_TYPE sharp_coord[DIM3],
		SCALAR_TYPE & scalar_stdev, SCALAR_TYPE & max_abs_scalar_error)
{
	subgrid_calculate_iso_vertex_in_cube
	(&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]),
			num_points, cube_coord, isovalue, subgrid_axis_size,
			sharp_coord, scalar_stdev, max_abs_scalar_error);
}


// **************************************************
// ROUTINES TO MOVE POINTS
// **************************************************

// Clamp to cube cube_coord[]
void SHARPISO::clamp_point
(const float cube_offset,
		const GRID_COORD_TYPE cube_coord[DIM3],
		COORD_TYPE point[DIM3])
{
	for (int i=0; i<DIM3; i++) {
		float p = point[i] - cube_coord[i];
		if (p < (-cube_offset))
		{ point[i] = cube_coord[i] - cube_offset; }
		if (p > 1+cube_offset)
		{ point[i]=  cube_coord[i] + 1.0 + cube_offset; }
	}
}

// Clamp to unit cube (0,0,0) to (1,1,1).
void SHARPISO::clamp_point
(const float cube_offset,
		COORD_TYPE point[DIM3])
{
	for (int i=0; i<DIM3; i++) {
		if (point[i]< (-cube_offset))
		{ point[i] = -cube_offset; }
		if (point[i] > 1.0 + cube_offset)
		{ point[i]=  1.0 + cube_offset; }
	}
}

void SHARPISO::postprocess_isovert_location
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		COORD_TYPE iso_coord[DIM3],
		SVD_INFO & svd_info)
{
	const COORD_TYPE max_dist = sharpiso_param.max_dist;

	svd_info.cube_containing_coord = cube_index;
	if (is_dist_to_cube_le(iso_coord, cube_coord, max_dist)) {

		if (!sharpiso_param.flag_allow_conflict) {
			bool flag_conflict;
			VERTEX_INDEX conflicting_cube;
			snap_to_cube(cube_coord, sharpiso_param.snap_dist, iso_coord);
			flag_conflict =
					check_conflict(scalar_grid, isovalue, cube_coord, iso_coord,
							conflicting_cube);

			if (flag_conflict) {
				process_conflict
				(scalar_grid, gradient_grid, cube_index, cube_coord, isovalue,
						sharpiso_param, iso_coord, svd_info);
				svd_info.flag_conflict = true;
				svd_info.cube_containing_coord = conflicting_cube;
			}
		}
	}
	else {
		process_far_point(scalar_grid, cube_index, cube_coord, isovalue,
				sharpiso_param, iso_coord, svd_info);
	}

	if (sharpiso_param.flag_round)
	{ IJK::round_coord
		(sharpiso_param.round_denominator, DIM3, iso_coord, iso_coord); }
}

void SHARPISO::process_conflict(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		COORD_TYPE iso_coord[DIM3],
		SVD_INFO & svd_info)
{
	if (!sharpiso_param.flag_allow_conflict) {

		if (sharpiso_param.flag_clamp_conflict) {
			// Clamp point to the cube.
			clamp_point(0, cube_coord, iso_coord);
		}
		else {
			// Use the centroid.
			compute_edgeI_centroid
			(scalar_grid, isovalue, cube_index, iso_coord);
			svd_info.location = CENTROID;
		}
	}

}

void SHARPISO::process_far_point
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX cube_index,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & sharpiso_param,
		COORD_TYPE iso_coord[DIM3],
		SVD_INFO & svd_info)
{

	if (sharpiso_param.flag_clamp_far) {
		// Clamp point to cube + max_dist.
		clamp_point(sharpiso_param.max_dist, cube_coord, iso_coord);
	}
	else {
		// Use the centroid.
		compute_edgeI_centroid
		(scalar_grid, isovalue, cube_index, iso_coord);
		svd_info.location = CENTROID;
	}
}

// **************************************************
// MISC ROUTINES
// **************************************************

/// Return true if (Linf) distance from coord to cube is at most
bool is_dist_to_cube_le
(const COORD_TYPE coord[DIM3], const GRID_COORD_TYPE cube_coord[DIM3],
		const SCALAR_TYPE max_dist)
{
	for (int d=0; d<DIM3; d++) {
		if (coord[d]+max_dist < cube_coord[d]) { return(false); }
		if (coord[d] > cube_coord[d]+1+max_dist) { return(false); }
	}

	return(true);
}

/// Return index of cube containing point.
/// Set flag_boundary to true if point is on cube boundary.
/// @pre Point is contained in grid.
/// @pre grid.AxisSize(d) > 0 for every axis d.
void get_cube_containing_point
(const SHARPISO_GRID & grid, const COORD_TYPE * coord,
		VERTEX_INDEX & cube_index, bool & flag_boundary)
{
	static COORD_TYPE coord2[DIM3];

	flag_boundary = false;
	for (int d = 0; d < DIM3; d++) {
		coord2[d] = floor(coord[d]);

		if (coord2[d] == coord[d])
		{ flag_boundary = true; }

		if (coord2[d] >= grid.AxisSize(d))
		{ coord2[d] = grid.AxisSize(d)-1; }
	}

	cube_index = grid.ComputeVertexIndex(coord2);
}

/// Return list of cubes containing point.
/// Set flag_boundary to true if point is on cube boundary.
/// @pre Point is contained in grid.
/// @pre grid.AxisSize(d) > 0 for every axis d.
void get_all_cubes_containing_point
(const SHARPISO_GRID & grid, const COORD_TYPE * coord,
		std::vector<VERTEX_INDEX> & cube_list)
{
	static COORD_TYPE coord2[DIM3];

	cube_list.clear();


	for (int index = 0; index < 8; index++) {

		bool flag_skip(false);
		for (int d = 0; d < DIM3; d++) {
			coord2[d] = floor(coord[d]);

			if (coord2[d] == coord[d]) {
				int mask = (1L << d);
				if ((mask & index) == 0) {
					if (coord2[d] > 0)
					{ coord2[d]--; }
					else
					{ flag_skip = true; }
				}
			}

			if (coord2[d] >= grid.AxisSize(d))
			{ flag_skip = true; }
		}

		if (!flag_skip) {
			VERTEX_INDEX cube_index = grid.ComputeVertexIndex(coord2);
			cube_list.push_back(cube_index);
		}
	}

	if (cube_list.size() == 0) {
		// Point is not on cube boundary.
		VERTEX_INDEX cube_index;
		bool flag_boundary;
		get_cube_containing_point(grid, coord, cube_index, flag_boundary);
		cube_list.push_back(cube_index);
	}

}

// Snap coordinate to cube.
void snap_to_cube
(const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE snap_dist, COORD_TYPE coord[DIM3])
{
	for (int d = 0; d < DIM3; d++) {
		if (coord[d] < cube_coord[d]) {
			if (coord[d] + snap_dist >= cube_coord[d])
			{ coord[d] = cube_coord[d]; }
		}
		else {
			GRID_COORD_TYPE right_coord = cube_coord[d]+1;
			if (coord[d] > right_coord) {
				if (coord[d] <= right_coord+snap_dist)
				{ coord[d] = right_coord; }
			}
		}
	}
}

// Count number of different coordinates
void diff_coord
(const GRID_COORD_TYPE coordA[DIM3], const GRID_COORD_TYPE coordB[DIM3],
		NUM_TYPE & num_diff, VERTEX_INDEX & icoord)
{
	icoord = 0;
	num_diff = 0;
	for (NUM_TYPE d = 0; d < DIM3; d++) {
		if (coordA[d] != coordB[d]) {
			num_diff++;
			icoord = d;
		}
	}
}

// **************************************************
// Check for conflict
// **************************************************

/// Return true if cube contains point.
/// *** MOVE TO IJKGRID. ***
bool cube_contains_point
(const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE point_coord[DIM3])
{
	for (int d = 0; d < DIM3; d++) {
		if (point_coord[d] < cube_coord[d])
		{ return(false); }
		if (point_coord[d] > cube_coord[d]+1)
		{ return(false); }
	}

	return(true);
}

/// Return true if point lies in an occupied cube
///   other than the one given by cube_coord[].
/// Returns conflicting cube.
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE point_coord[DIM3],
		VERTEX_INDEX & conflicting_cube)
{
	if (cube_contains_point(cube_coord, point_coord)) {
		// Point is in cube.  No conflict.
		return(false);
	}

	if (scalar_grid.ContainsPoint(point_coord)) {
		// Point is inside grid.
		VERTEX_INDEX icube2;
		bool flag_boundary;
		get_cube_containing_point(scalar_grid, point_coord, icube2, flag_boundary);

		if (IJK::is_gt_cube_min_le_cube_max(scalar_grid, icube2, isovalue)) {
			// Isosurface intersects icube2. Conflict.
			conflicting_cube = icube2;
			return(true);
		}

		if (flag_boundary) {
			std::vector<VERTEX_INDEX> cube_list;
			get_all_cubes_containing_point(scalar_grid, point_coord, cube_list);

			for (int i = 0; i < cube_list.size(); i++) {
				if (IJK::is_gt_cube_min_le_cube_max
						(scalar_grid, cube_list[i], isovalue)) {
					// Isosurface intersects cube_list[i]. Conflict.
					conflicting_cube = cube_list[i];
					return(true);
				}
			}
		}
	}

	// No conflict.
	return(false);
}

/// Return true if point lies in an occupied cube
///   other than the one given by cube_coord[].
/// Returns conflicting cube.
bool check_conflict
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const GRID_COORD_TYPE cube_coord[DIM3],
		const COORD_TYPE point_coord[DIM3])
{
	VERTEX_INDEX conflicting_cube;

	return(check_conflict(scalar_grid, isovalue, cube_coord, point_coord,
			conflicting_cube));
}

// **************************************************
// SHARP_ISOVERT_PARAM
// **************************************************

/// Initialize SHARP_ISOVERT_PARAM
void SHARPISO::SHARP_ISOVERT_PARAM::Init()
{
	use_lindstrom = true;
	use_lindstrom2 = false;
	use_lindstrom_fast = true;
	flag_allow_conflict = false;
	flag_clamp_conflict = true;
	flag_clamp_far = false;
	flag_round = false;
	flag_dist2centroid = true;
	flag_check_disk = false;
	use_sharp_edgeI = false;
	use_Linf_dist = true;
	max_dist = 1.0;
	snap_dist = 1.0/16.0;
	max_small_eigenvalue = 0.1;
	round_denominator = 16;
	max_small_grad_coord_Linf = 0.2;
	linf_dist_thresh_merge_sharp = 1.0;
	bin_width = 5;
}

// **************************************************
// Class SHARP_ISOVERT_PARAM member functions
// **************************************************

/// Set
void SHARP_ISOVERT_PARAM::Set(const SHARP_ISOVERT_PARAM & param)
{
	*this = param;
}


// **************************************************
// SVD_INFO
// **************************************************

/// Initialize.
void SVD_INFO::Init()
{
	flag_conflict = false;
	flag_Linf_iso_vertex_location = false;
	location = LOC_SVD;
}

/// Set ray information.
void SVD_INFO::SetRayInfo
(const COORD_TYPE origin[DIM3], const COORD_TYPE direction[DIM3],
		const COORD_TYPE intersect[DIM3])

{
	IJK::copy_coord_3D(origin, ray_initial_point);
	IJK::copy_coord_3D(direction, ray_direction);
	IJK::copy_coord_3D(intersect, ray_cube_intersection);
}

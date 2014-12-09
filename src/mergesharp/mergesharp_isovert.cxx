/// \file mergesharp_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.

/*
Copyright (C) 2012-2014 Arindam Bhattacharya

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

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>  
#include <stdio.h>
#include <stdio.h>

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkgrid_macros.h"
#include "ijkscalar_grid.txx"

#include "mergesharp_types.h"
#include "mergesharp_isovert.h"
#include "mergesharp_position.h"

using namespace std;
using namespace SHARPISO;
using namespace MERGESHARP;
using namespace IJK;


// *** DEBUG ***
#include "ijkprint.txx"
bool flag_debug(false);


namespace {


	void compute_linf_dist
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE isovert_coord[DIM3], SCALAR_TYPE & linf_dist);

	void compute_linf_dist
  (const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert);

  bool find_3x3x3_overlap
  (const GRID_COORD_TYPE cubeA_coord[DIM3],
   const GRID_COORD_TYPE cubeB_coord[DIM3],
   GRID_COORD_TYPE rmin[DIM3], GRID_COORD_TYPE rmax[DIM3],
   int & overlap_dim);

  template <typename GRID_TYPE>
  bool find_3x3x3_overlap
  (const GRID_TYPE & grid,
   const VERTEX_INDEX cubeA_index, const VERTEX_INDEX cubeB_index,
   GRID_COORD_TYPE rmin[DIM3], GRID_COORD_TYPE rmax[DIM3],
   int & overlap_dim);

	bool are_connected
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2, const SCALAR_TYPE isovalue);

	bool is_angle_large
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const ISOVERT & isovertData, const VERTEX_INDEX iv,
		const SCALAR_TYPE threshold, const VERTEX_INDEX v1,
		const VERTEX_INDEX v2);

  bool is_facet_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, VERTEX_INDEX & cube_index1);

  bool is_facet_or_edge_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   VERTEX_INDEX & cube_index2);

  bool is_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   const VERTEX_INDEX & cube_index2);

  bool is_cube_or_neighbor_covered_by
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1);

  void store_svd_info
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX cube_index,
   const NUM_TYPE gcube_index,
   const NUM_TYPE num_large_eigenvalues,
   const SVD_INFO svd_info,
   ISOVERT & isovert, ISOVERT_INFO & isovert_info);

  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, ISOVERT & isovert);

  void output_gcube_list
  (const SHARPISO_GRID & grid, const ISOVERT & isovert, 
   const std::vector<NUM_TYPE> & gcube_index_list);

}

// **************************************************
// ISOVERT ROUTINES
// **************************************************

/**
* **Compute Isovertex Positions from gradients. 
*/
void compute_isovert_positions 
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert,
	ISOVERT_INFO & isovert_info)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
		isovert_param.grad_selection_cube_offset;
	OFFSET_VOXEL voxel;
	SVD_INFO svd_info;

	voxel.SetVertexCoord
		(scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);

	if (isovert_param.use_lindstrom) {

		IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
			NUM_TYPE index = isovert.sharp_ind_grid.Scalar(iv);

			if (index != ISOVERT::NO_INDEX) {
				// this is an active cube
				isovert.gcube_list[index].flag_centroid_location = false;

				// compute the sharp vertex for this cube
				EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
				NUM_TYPE num_large_eigenvalues;

				svd_compute_sharp_vertex_for_cube_lindstrom
					(scalar_grid, gradient_grid, iv, isovalue, isovert_param, voxel,
					isovert.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);
				store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
					svd_info, isovert, isovert_info);
			}
		}
	}
	else {

		IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
			NUM_TYPE index = isovert.sharp_ind_grid.Scalar(iv);

			if (index != ISOVERT::NO_INDEX) {
				// this is an active cube
				isovert.gcube_list[index].flag_centroid_location = false;

				// compute the sharp vertex for this cube
				EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
				NUM_TYPE num_large_eigenvalues;

				svd_compute_sharp_vertex_for_cube_lc_intersection
					(scalar_grid, gradient_grid, iv, isovalue, isovert_param, voxel,
					isovert.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);

				store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
					svd_info, isovert, isovert_info);
			}
		}

	}

	store_boundary_bits(scalar_grid, isovert.gcube_list);
  set_cube_containing_isovert(scalar_grid, isovalue, isovert);
}

// Compute isosurface vertex positions using isosurface-edge intersections.
void compute_isovert_positions_edgeI
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const VERTEX_POSITION_METHOD vertex_position_method,
	ISOVERT & isovert,
	ISOVERT_INFO & isovert_info)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
		isovert_param.grad_selection_cube_offset;
	SVD_INFO svd_info;

	if (vertex_position_method == EDGEI_GRADIENT) {

		IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
			NUM_TYPE index = isovert.sharp_ind_grid.Scalar(iv);

			if (index != ISOVERT::NO_INDEX) {
				// this is an active cube
				isovert.gcube_list[index].flag_centroid_location = false;

				// compute the sharp vertex for this cube
				EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
				NUM_TYPE num_large_eigenvalues;

				svd_compute_sharp_vertex_edgeI_sharp_gradient
					(scalar_grid, gradient_grid, iv, isovalue, isovert_param,
					isovert.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);


				store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
					svd_info, isovert, isovert_info);
			}
		}

	}
	else {
		// vertex_position_method == EDGEI_INTERPOLATE

		IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
			NUM_TYPE index = isovert.sharp_ind_grid.Scalar(iv);

			if (index != ISOVERT::NO_INDEX) {
				// this is an active cube
				isovert.gcube_list[index].flag_centroid_location = false;

				// compute the sharp vertex for this cube
				EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
				NUM_TYPE num_large_eigenvalues;

				svd_compute_sharp_vertex_edgeI_interpolate_gradients
					(scalar_grid, gradient_grid, iv, isovalue, isovert_param,
					isovert.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);

				store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
					svd_info, isovert, isovert_info);

			}
		}

	}

	store_boundary_bits(scalar_grid, isovert.gcube_list);
}

/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered.
/// takes isovert_info as parameter.
// NOT USED
void MERGESHARP::recompute_isovert_positions 
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert,
	ISOVERT_INFO & isovert_info)
{
	for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
		GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

		if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE)) 
		{
			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

			compute_edgeI_centroid
				(scalar_grid, gradient_grid, isovalue, cube_index,
				isovert_param.use_sharp_edgeI, isovert.gcube_list[i].isovert_coord);
			isovert.gcube_list[i].flag_centroid_location = true;
		}

		// recompute the covered point
		if (cube_flag == COVERED_POINT)
		{
			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

			const SIGNED_COORD_TYPE grad_selection_cube_offset =
				isovert_param.grad_selection_cube_offset;
			OFFSET_VOXEL voxel;
			SVD_INFO svd_info;

			voxel.SetVertexCoord
				(scalar_grid.SpacingPtrConst(), grad_selection_cube_offset);
			//DEBUG
			// compute the sharp vertex for this cube
			EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
			NUM_TYPE num_large_eigenvalues;

			//cout <<"isovert param "<<isovert_param.use_intersected_edge_endpoint_gradients <<" "
			//	<<isovert_param.use_zero_grad_boundary <<endl;

			SHARP_ISOVERT_PARAM  temp_isovert_param = isovert_param;

			temp_isovert_param.max_small_eigenvalue = temp_isovert_param.max_small_eigenvalue*4;
			//temp_isovert_param.max_small_eigenvalue = 1.0;
			//temp_isovert_param.max_grad_dist = 2;

			//cout <<"temp isovert param "<<temp_isovert_param.use_intersected_edge_endpoint_gradients<<" "
			//	<<temp_isovert_param.use_zero_grad_boundary << endl;
			////shift to gradNIES 
			//temp_isovert_param.use_only_cube_gradients = false;
			//temp_isovert_param.use_selected_gradients = true;

			//temp_isovert_param.use_intersected_edge_endpoint_gradients = true;
			//temp_isovert_param.use_zero_grad_boundary = false;
			//temp_isovert_param.use_large_neighborhood = false;

			/*
			//find the local gradients. 
			if(cube_index == 22185){
			cout <<"test"<<endl;
			COORD_TYPE coord1[DIM3];

			for (NUM_TYPE k = 0; k < NUM_CUBE_VERTICES3D; k++) 
			{
			int cc =  scalar_grid.CubeVertex(cube_index, k);

			scalar_grid.ComputeCoord(cc, coord1);
			cout <<"ci "<<cc<<" ["<<coord1[0]<<","<<coord1[1]<<","<<coord1[2]<<"] ";
			scalar_grid.ComputeScaledCoord(cc, coord1);
			cout <<"scaled "<<cc<<" ["<<coord1[0]<<","<<coord1[1]<<","<<coord1[2]<<"]";
			cout<<"cube center "<< coord1[0]+scalar_grid.Spacing(0)/2 <<" "
			<< coord1[1]+scalar_grid.Spacing(1)/2<<" "
			<<coord1[2]+scalar_grid.Spacing(2)/2<<endl;
			//red

			}


			vector<VERTEX_INDEX> vertex_list(NUM_CUBE_VERTICES3D,0);
			get_cube_vertices(scalar_grid, cube_index, &vertex_list[0]);
			vector<COORD_TYPE> point_coord(NUM_CUBE_VERTICES3D*DIM3,0);
			vector<COORD_TYPE> gradient_coord(NUM_CUBE_VERTICES3D*DIM3,-1);
			vector<COORD_TYPE> scalar(NUM_CUBE_VERTICES3D,0);

			get_vertex_gradients
			(scalar_grid, gradient_grid, &(vertex_list[0]), NUM_CUBE_VERTICES3D,
			&(point_coord[0]), &(gradient_coord[0]), &(scalar[0]));

			//DEBUG
			for (int c = 0; c < NUM_CUBE_VERTICES3D; c++)
			{
			int v = scalar_grid.ComputeVertexIndex(&(point_coord[DIM3*c]));
			cout <<vertex_list[c]<<" ["<<point_coord[DIM3*c+0]<<","<<point_coord[DIM3*c+1]<<","<<point_coord[DIM3*c+2]<<"]ci: "<<v;
			cout <<" gradients ["<<gradient_coord[DIM3*c+0]<<","<<gradient_coord[DIM3*c+1]
			<<","<<gradient_coord[DIM3*c+2]<<"] scalar : "<< scalar[c]<<endl;

			}

			cout <<"test end "<<endl;
			}

			svd_compute_sharp_vertex_for_cube_lindstrom
			(scalar_grid, gradient_grid, cube_index, isovalue, temp_isovert_param, voxel,
			isovert.gcube_list[i].isovert_coord,
			eigenvalues, num_large_eigenvalues, svd_info);
			store_svd_info(scalar_grid, cube_index, i, num_large_eigenvalues,
			svd_info, isovert, isovert_info);
			//DEBUG
			cout <<"new coveredpoint "<< isovert.gcube_list[i].isovert_coord[0]
			<<" "<<isovert.gcube_list[i].isovert_coord[1]
			<<" "<<isovert.gcube_list[i].isovert_coord[2]<<endl;
			cout <<"cube index "<< cube_index <<" ";


			cout <<"num large eigen "<< num_large_eigenvalues<<" [";
			for (int y = 0; y < num_large_eigenvalues; y++)
			{
			cout <<eigenvalues[y]<<" ";
			}
			cout <<"] \n\n";
			*/


			compute_edgeI_centroid
				(scalar_grid, gradient_grid, isovalue, cube_index,
				isovert_param.use_sharp_edgeI, isovert.gcube_list[i].isovert_coord);

			isovert.gcube_list[i].flag_centroid_location = true;
		}
	}
}

/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered.
void MERGESHARP::recompute_isovert_positions 
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert)
{

	for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
		GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

		if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE) 
			|| (cube_flag == COVERED_POINT)) 
		{
			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
			//DEBUG
			//if (cube_flag == COVERED_POINT)
			//{
			// cout <<"Remaining Covered point ["<<isovert.gcube_list[i].cube_index<<"] "<< isovert.gcube_list[i].isovert_coord[0]
			// <<" "<<isovert.gcube_list[i].isovert_coord[1]
			// <<" "<<isovert.gcube_list[i].isovert_coord[2]<<endl;

			//}

			compute_edgeI_centroid
				(scalar_grid, gradient_grid, isovalue, cube_index,
				isovert_param.use_sharp_edgeI, isovert.gcube_list[i].isovert_coord);
			isovert.gcube_list[i].flag_centroid_location = true;
			//DEBUG
			//if (cube_flag == COVERED_POINT)
			//{
			// cout <<"to sharp centroid "<< isovert.gcube_list[i].isovert_coord[0]
			// <<" "<<isovert.gcube_list[i].isovert_coord[1]
			// <<" "<<isovert.gcube_list[i].isovert_coord[2]<<endl;

			//}

		}
	}

}


/// Recompute isosurface vertex positions for cubes 
/// which are not selected or covered.
/// This is not called by the sharp computations
void MERGESHARP::recompute_isovert_positions 
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert)
{
	SHARPISO_EDGE_INDEX_GRID edge_index
		(DIM3, scalar_grid.AxisSize(), DIM3);

	set_edge_index(edgeI_coord, edge_index);

	for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
		GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

		if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE)) {

			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
			compute_edgeI_centroid
				(scalar_grid, edgeI_coord, edge_index, isovalue, cube_index,
				isovert.gcube_list[i].isovert_coord);



			isovert.gcube_list[i].flag_centroid_location = true;
		}
	}

}


void round_down(const COORD_TYPE * coord, GRID_COORD_TYPE min_coord[DIM3])
{
	for (NUM_TYPE d = 0; d < DIM3; d++)
	{ min_coord[d] = int(std::floor(coord[d])); }
}

void MERGESHARP::set_edge_index
	(const std::vector<COORD_TYPE> & edgeI_coord,
	SHARPISO_EDGE_INDEX_GRID & edge_index)
{
	GRID_COORD_TYPE min_coord[DIM3];

	edge_index.SetAllCoord(ISOVERT::NO_INDEX);

	for (NUM_TYPE i = 0; i < edgeI_coord.size()/DIM3; i++) {
		round_down(&(edgeI_coord[i*DIM3]), min_coord);

		int edge_dir = 0;
		for (int d = 1; d < DIM3; d++) {
			if (edgeI_coord[i*DIM3+d] != min_coord[d]) { edge_dir = d; }
		}

		VERTEX_INDEX iv0 = edge_index.ComputeVertexIndex(min_coord);

		edge_index.Set(iv0, edge_dir, i);
	}

	// Set edge index from edgeI_coord[] which are on grid vertices.
	for (NUM_TYPE i = 0; i < edgeI_coord.size()/DIM3; i++) {
		round_down(&(edgeI_coord[i*DIM3]), min_coord);

		if (is_coord_equal_3D(&(edgeI_coord[i*DIM3]), min_coord)) {

			VERTEX_INDEX iv0 = edge_index.ComputeVertexIndex(min_coord);
			for (int edge_dir = 0; edge_dir < DIM3; edge_dir++) {
				if (edge_index.Vector(iv0, edge_dir) == ISOVERT::NO_INDEX)
				{ edge_index.Set(iv0, edge_dir, i); }

				if (min_coord[edge_dir] > 0) {
					VERTEX_INDEX iv1 = edge_index.PrevVertex(iv0, edge_dir);
					if (edge_index.Vector(iv1, edge_dir) == ISOVERT::NO_INDEX)
					{ edge_index.Set(iv1, edge_dir, i); }
				}
			}
		}


	}

}


// Compute isosurface vertex positions using hermite data.
void compute_isovert_positions 
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const std::vector<COORD_TYPE> & edgeI_normal_coord,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
		isovert_param.grad_selection_cube_offset;
	SHARPISO_EDGE_INDEX_GRID edge_index
		(DIM3, scalar_grid.AxisSize(), DIM3);
	SVD_INFO svd_info;

	set_edge_index(edgeI_coord, edge_index);

	IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX) {
		NUM_TYPE index = isovert.sharp_ind_grid.Scalar(iv);
		if (index!=ISOVERT::NO_INDEX) {
			// this is an active cube
			// compute the sharp vertex for this cube
			EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
			VERTEX_INDEX num_large_eigenvalues;

			svd_compute_sharp_vertex_for_cube_hermite
				(scalar_grid, edgeI_coord, edgeI_normal_coord, edge_index,
				iv, isovalue, isovert_param,
				isovert.gcube_list[index].isovert_coord,
				eigenvalues, num_large_eigenvalues, svd_info);

			//set num eigen
			isovert.gcube_list[index].num_eigenvalues =  
				(unsigned char)num_large_eigenvalues;
			//set the sharp vertex type to be *AVAILABLE*
			if(num_large_eigenvalues > 1)
			{ isovert.gcube_list[index].flag=AVAILABLE_GCUBE; }
			else
			{ isovert.gcube_list[index].flag=SMOOTH_GCUBE; }

			// store distance
			compute_linf_dist(scalar_grid, iv,
				isovert.gcube_list[index].isovert_coord,
				isovert.gcube_list[index].linf_dist);
		}
	}

	store_boundary_bits(scalar_grid, isovert.gcube_list);
}


// Copy isovert position from from_gcube to to_gcube.
void copy_isovert_position
(const SHARPISO_GRID & grid, 
 const NUM_TYPE from_gcube, const NUM_TYPE to_gcube,
 ISOVERT & isovert)
{
  // *** DEBUG ***
  /*
    using namespace std;
    VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
    cerr << "Setting " << conflicting_cube
    << " from " << cubeA_index << endl;
  */

  IJK::copy_coord_3D(isovert.gcube_list[from_gcube].isovert_coord,
                     isovert.gcube_list[to_gcube].isovert_coord);
  isovert.gcube_list[to_gcube].cube_containing_isovert = 
    isovert.gcube_list[to_gcube].cube_index;
  isovert.gcube_list[to_gcube].num_eigenvalues =
    isovert.gcube_list[from_gcube].num_eigenvalues;
  if (isovert.gcube_list[to_gcube].flag == COVERED_POINT)
    { isovert.gcube_list[to_gcube].flag = AVAILABLE_GCUBE; }
  isovert.gcube_list[to_gcube].flag_conflict = false;
  isovert.gcube_list[to_gcube].flag_coord_from_other = true;

  // store new Linf distance
  compute_linf_dist(grid, to_gcube, isovert);
}


// Change isovert positions in cubes with conflicts.
void swap_isovert_positions
(const SHARPISO_GRID & grid, ISOVERT & isovert)
{
  std::vector<NUM_TYPE> gcube_sharp_list;

  get_corner_or_edge_cubes(isovert.gcube_list, gcube_sharp_list);

  for (NUM_TYPE i = 0; i < gcube_sharp_list.size(); i++) {

    NUM_TYPE gcubeA_index = gcube_sharp_list[i];
    if (isovert.gcube_list[gcubeA_index].flag_conflict) {

      VERTEX_INDEX conflicting_cube = 
        isovert.gcube_list[gcubeA_index].cube_containing_isovert;
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(conflicting_cube);

      if (gcubeB_index != ISOVERT::NO_INDEX) {

        if (isovert.gcube_list[gcubeB_index].flag_conflict) {

          VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
          if (isovert.gcube_list[gcubeB_index].cube_containing_isovert == 
              cubeA_index) {
            // swap isovert_coord for gcubeA_index and gcubeB_index

            // *** DEBUG ***
            /*
            using namespace std;
            VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
            cerr << "Swapping " << conflicting_cube
                 << " and " << cubeA_index << endl;
            */

            swap(isovert.gcube_list[gcubeA_index].isovert_coord[0],
                 isovert.gcube_list[gcubeB_index].isovert_coord[0]);
            swap(isovert.gcube_list[gcubeA_index].isovert_coord[1],
                 isovert.gcube_list[gcubeB_index].isovert_coord[1]);
            swap(isovert.gcube_list[gcubeA_index].isovert_coord[2],
                 isovert.gcube_list[gcubeB_index].isovert_coord[2]);
            swap(isovert.gcube_list[gcubeA_index].num_eigenvalues,
                 isovert.gcube_list[gcubeB_index].num_eigenvalues);

            isovert.gcube_list[gcubeA_index].flag_conflict = false;
            isovert.gcube_list[gcubeB_index].flag_conflict = false;
            isovert.gcube_list[gcubeA_index].flag_coord_from_other = true;
            isovert.gcube_list[gcubeB_index].flag_coord_from_other = true;
            if (isovert.gcube_list[gcubeA_index].flag == COVERED_POINT)
              { isovert.gcube_list[gcubeA_index].flag = AVAILABLE_GCUBE; }
            if (isovert.gcube_list[gcubeB_index].flag == COVERED_POINT)
              { isovert.gcube_list[gcubeB_index].flag = AVAILABLE_GCUBE; }

            // store new Linf distance
            compute_linf_dist(grid, gcubeA_index, isovert);
            compute_linf_dist(grid, gcubeB_index, isovert);
          }
          else {
            copy_isovert_position
              (grid, gcubeA_index, gcubeB_index, isovert);
          }
        }
      }
      else {
        // Should never happen, but...
        IJK::PROCEDURE_ERROR error("swap_isovert_positions");
        error.AddMessage("Programming error.  Illegal conflicting cube index: ",
                         conflicting_cube, " for cube ",
                         isovert.CubeIndex(gcubeA_index), ".");
        throw error;

      }
    }
  }

}

// Change isovert positions which lie in covered cubes.
void reset_covered_isovert_positions
(const	SHARPISO_GRID_NEIGHBORS & grid, 
 const SHARPISO_BOOL_GRID & covered_grid,
 ISOVERT & isovert)
{
  std::vector<NUM_TYPE> gcube_sharp_list;

  get_corner_or_edge_cubes(isovert.gcube_list, gcube_sharp_list);

  for (NUM_TYPE i = 0; i < gcube_sharp_list.size(); i++) {

    NUM_TYPE gcubeA_index = gcube_sharp_list[i];

    if (isovert.gcube_list[gcubeA_index].flag_conflict) {

      VERTEX_INDEX conflicting_cube = 
        isovert.gcube_list[gcubeA_index].cube_containing_isovert;
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(conflicting_cube);

      if (gcubeB_index != ISOVERT::NO_INDEX) {

        if (isovert.gcube_list[gcubeB_index].flag != SELECTED_GCUBE) {

          VERTEX_INDEX cubeB_index = isovert.CubeIndex(gcubeB_index);
          VERTEX_INDEX cubeC_index =
            isovert.gcube_list[gcubeB_index].cube_containing_isovert;

          if (cubeC_index != cubeB_index) {


            if (covered_grid.Scalar(cubeC_index)) {

              // *** DEBUG ***
              /*
              using namespace std;
              VERTEX_INDEX cubeA_index = isovert.CubeIndex(gcubeA_index);
              cerr << "Copying position from " << cubeA_index << " ";
              ijkgrid_output_vertex_coord(cerr, grid, cubeA_index);
              cerr << " to " << cubeB_index << " ";
              ijkgrid_output_vertex_coord(cerr, grid, cubeB_index);
              cerr << endl;
              cerr << "  Old position: ";
              IJK::print_coord3D(cerr, isovert.gcube_list[gcubeB_index].isovert_coord);
              cerr << ".  New position: ";
              IJK::print_coord3D(cerr, isovert.gcube_list[gcubeA_index].isovert_coord);
              cerr << endl;
              */

              copy_isovert_position(grid, gcubeA_index, gcubeB_index, isovert);
            }
          }
        }
      }
      else {
        // Should never happen, but...
        IJK::PROCEDURE_ERROR error("reset_covered_isovert_positions");
        error.AddMessage("Programming error.  Illegal conflicting cube index: ",
                         conflicting_cube, " for cube ",
                         isovert.CubeIndex(gcubeA_index), ".");
        throw error;

      }
    }
  }

}

class GCUBE_COMPARE {

public:
	const std::vector<GRID_CUBE> * gcube_list;

	GCUBE_COMPARE(const std::vector<GRID_CUBE> & gcube_list)
	{ this->gcube_list = &gcube_list; };

	bool operator () (int i,int j)
	{
		if (gcube_list->at(i).num_eigenvalues == 
			gcube_list->at(j).num_eigenvalues) {
				return ((gcube_list->at(i).linf_dist) < (gcube_list->at(j).linf_dist)); 
		}
		else {
			return ((gcube_list->at(i).num_eigenvalues) > 
				(gcube_list->at(j).num_eigenvalues)); 
		}
	}
};

/// Get cubes with more than one eigenvalue.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to sharp cubes
///    sorted by number of large eigenvalues and by distance 
///    of sharp coord from cube center.
void MERGESHARP::get_corner_or_edge_cubes
(const std::vector<GRID_CUBE> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list)
{
  GCUBE_COMPARE gcube_compare(gcube_list);

  for (int i=0;i<gcube_list.size();i++)
    {
      // *** SHOULD CHECK THAT POINT IS NOT COMPUTED USING CENTROID ***
      if (gcube_list[i].num_eigenvalues > 1)
        { gcube_index_list.push_back(i); }
    }

  sort(gcube_index_list.begin(), gcube_index_list.end(), gcube_compare);
}

void MERGESHARP::store_boundary_bits
	(const SHARPISO_GRID & grid, GRID_CUBE_ARRAY & gcube_list)
{
	for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {
		grid.ComputeBoundaryCubeBits
			(gcube_list[i].cube_index, gcube_list[i].boundary_bits);
	}
}

/// Compute the cube index from the gc index
VERTEX_INDEX cube_ind_frm_gc_ind
	(const ISOVERT &isoData, const  NUM_TYPE &gc_ind){
		return isoData.gcube_list[gc_ind].cube_index;
}

/// get the *connected_list* of vertices which are selected and
/// connected to vertex *vert_ind*
void get_connected(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX vert_index,
	const vector<VERTEX_INDEX> &selected_list,
	vector<VERTEX_INDEX> &connected_list)
{
	for (int iv=0; iv < selected_list.size();iv++){
		if (are_connected(scalar_grid, vert_index, selected_list[iv], isovalue)){
			connected_list.push_back(selected_list[iv]);
		}
	}
}

inline void divide_coord_3D
	(const AXIS_SIZE_TYPE bin_width, GRID_COORD_TYPE coord[DIM3])
{
	for (int d = 0; d < DIM3; d++)
		coord[d] = IJK::integer_divide(coord[d], bin_width);
}

/// Get the selected vertices around iv.
void get_selected
	(const SHARPISO_GRID & grid,
	const VERTEX_INDEX iv,
	const BIN_GRID<VERTEX_INDEX> & bin_grid,
	const AXIS_SIZE_TYPE bin_width,
	std::vector<VERTEX_INDEX> & selected_list)
{
	const int dimension = grid.Dimension();
	static GRID_COORD_TYPE coord[DIM3];
	static GRID_COORD_TYPE min_coord[DIM3];
	static GRID_COORD_TYPE max_coord[DIM3];

	long boundary_bits;

	grid.ComputeCoord(iv, coord);
	divide_coord_3D(bin_width, coord);
	VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
	bin_grid.ComputeBoundaryBits(ibin, boundary_bits);

	selected_list.resize(bin_grid.ListLength(ibin));
	std::copy(bin_grid(ibin).list.begin(), bin_grid(ibin).list.end(),
            selected_list.begin());

	if (boundary_bits == 0) {

		for (NUM_TYPE k = 0; k < bin_grid.NumVertexNeighborsC(); k++) {
			VERTEX_INDEX jbin = bin_grid.VertexNeighborC(ibin, k);
			for (int i = 0; i < bin_grid.ListLength(jbin); i++) {
        selected_list.push_back(bin_grid.List(jbin, i)); 
      }
		}
	}
	else {

		for (int d = 0; d < dimension; d++) {
			if (coord[d] > 0)
			{ min_coord[d] = coord[d] - 1; }
			else
			{ min_coord[d] = 0; }

			if (coord[d]+1 < bin_grid.AxisSize(d))
			{ max_coord[d] = coord[d] + 1; }
			else
			{ max_coord[d] = coord[d]; }
		}

		for (coord[0] = min_coord[0]; coord[0] <= max_coord[0];
			coord[0]++)
			for (coord[1] = min_coord[1]; coord[1] <= max_coord[1];
				coord[1]++)
				for (coord[2] = min_coord[2]; coord[2] <= max_coord[2];
					coord[2]++) {

						VERTEX_INDEX jbin = bin_grid.ComputeVertexIndex(coord);
						if (ibin != jbin) {
							for (int i = 0; i < bin_grid.ListLength(jbin); i++)
							{ selected_list.push_back(bin_grid.List(jbin, i)); }
						}
				}
	}
}


// Check if selecting this vertex creates a triangle with a large angle.
// @param check_triangle_angle If true, check it triangle has large angles.
// @param bin_grid Contains the already selected vertices.
// @param[out] v1,v2 vertex indices which form a triangle with iv
bool MERGESHARP::creates_triangle 
	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const bool check_triangle_angle,
	const ISOVERT & isovert,
	const VERTEX_INDEX iv,
	const SCALAR_TYPE isovalue,
	const BIN_GRID<VERTEX_INDEX> & bin_grid,
	const AXIS_SIZE_TYPE bin_width,
	VERTEX_INDEX & v1,
	VERTEX_INDEX & v2)
{
	const SCALAR_TYPE threshold = cos(140*M_PI/180);
	vector <VERTEX_INDEX> selected_list;
	vector <VERTEX_INDEX> connected_list;

	// get the selected vertices around iv
	get_selected(scalar_grid, iv, bin_grid, bin_width, selected_list);

	// get the list of vertices connected to the vertex iv
	get_connected(scalar_grid, isovalue, iv, selected_list, connected_list);

	int limit = connected_list.size();
	// for each pair jv1 jv2 in the connected list
	for (int i=0; i< limit-1; ++i){
		for(int j=i+1; j <= (limit-1); ++j)
		{
			if (are_connected(scalar_grid, connected_list[i],
				connected_list[j], isovalue))
			{
				v1 = connected_list[i];
				v2 = connected_list[j];

				if (check_triangle_angle){
					// checking angle
					bool flag_large_angle = 
						is_angle_large(scalar_grid, isovert, iv, threshold, v1 , v2);

					if (flag_large_angle)
					{ 
            // *** DEBUG ***
            if (flag_debug) {
              using namespace std;
              cerr << "   Cube " << iv << " ";
              ijkgrid_output_vertex_coord(cerr, scalar_grid, iv);
              cerr << " creates large angle triangle with:" << endl;
              cerr << "      " << v1 << " ";
              ijkgrid_output_vertex_coord(cerr, scalar_grid, v1);
              cerr << " and " << v2 << " ";
              ijkgrid_output_vertex_coord(cerr, scalar_grid, v2);
              cerr << endl;
            }

						return true; 
					}
				}
				else 
				{
					// returning true without checking angles angles
					return true;
				}
			}
		}
	}


	return false;
}

/// Initialize bin_grid.
/// @param bin_width = number of cubes along each axis.
void MERGESHARP::init_bin_grid
	(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
	BIN_GRID<VERTEX_INDEX> & bin_grid)
{
	const int dimension = grid.Dimension();
	IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);

	for (int d = 0; d < dimension; d++) {
		axis_size[d] = IJK::compute_subsample_size(grid.AxisSize(d), bin_width);
	}

	//bin_grid.SetSize(dimension, axis_size);
	bin_grid.SetSize(dimension, axis_size.PtrConst());

}

void MERGESHARP::bin_grid_insert
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid)
{
	static GRID_COORD_TYPE coord[DIM3];

	grid.ComputeCoord(cube_index, coord);
	divide_coord_3D(bin_width, coord);
	VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
	bin_grid.Insert(ibin, cube_index);
}

void MERGESHARP::bin_grid_remove
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid)
{
	static GRID_COORD_TYPE coord[DIM3];

	grid.ComputeCoord(cube_index, coord);
	divide_coord_3D(bin_width, coord);
	VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
	bin_grid.Remove(ibin, cube_index);
}


/*
Check if the sharp vertex is in a cube which is already covered.
if yes then return true else return false.
*/
bool check_covered_point(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SHARPISO_BOOL_GRID &covered_grid,
	ISOVERT &isovert,
	const VERTEX_INDEX &gcube_index)
{
  VERTEX_INDEX cube_index;

  cube_index = isovert.gcube_list[gcube_index].cube_containing_isovert;
  
	if (covered_grid.Scalar(cube_index)) { return true; }

  // *** IS THIS NECESSARY? ***
  // *** SHOULD USE IsCoveredOrSelected() ***
	// if the sharp vertex is not in the same cube. 
	if (isovert.sharp_ind_grid.Scalar(cube_index) != ISOVERT::NO_INDEX){
		if(isovert.gcube_list[gcube_index].cube_index != cube_index){
			if (isovert.isFlag(cube_index, COVERED_A_GCUBE))
        { return true; }
		}
	}

	return false;
}



void select_cube
	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
  const NUM_TYPE gcube_index,
	const GRID_CUBE_FLAG flag,
	ISOVERT & isovert
	)
{
	const int bin_width = isovert_param.bin_width;
  const VERTEX_INDEX cube_index = 
    isovert.gcube_list[gcube_index].cube_index;

  isovert.gcube_list[gcube_index].flag = SELECTED_GCUBE;

  covered_grid.Set(cube_index, true);

  // *** DEBUG ***
  if (flag_debug) {
    using namespace std;
    cerr << "Selecting cube " << cube_index << "  ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << "  Linf_dist: " << isovert.gcube_list[gcube_index].linf_dist;
    cerr << endl;
  }

  bin_grid_insert(scalar_grid, bin_width, cube_index, bin_grid);

  // mark all the neighbors as covered
  for (int i=0;i < gridn.NumVertexNeighborsC(); i++) {
    VERTEX_INDEX cube_index2 = gridn.VertexNeighborC(cube_index, i);

    covered_grid.Set(cube_index2, true);

    NUM_TYPE gcube_index2 = isovert.sharp_ind_grid.Scalar(cube_index2);
    if(gcube_index2 != ISOVERT::NO_INDEX) {
      isovert.gcube_list[gcube_index2].flag = flag;

      if (isovert.gcube_list[gcube_index2].covered_by ==
          isovert.gcube_list[gcube_index2].cube_index) {
        isovert.gcube_list[gcube_index2].covered_by = cube_index;
      }
    }
	}
}

void check_and_select_cube
	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
  const NUM_TYPE gcube_index,
	const GRID_CUBE_FLAG flag,
	ISOVERT & isovert)
{
	const int bin_width = isovert_param.bin_width;
  VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;
	VERTEX_INDEX v1, v2;

	//Check if the sharp vertex is a covered point.
	//Covered point: the point is inside a covered cube. 
	bool flag_covered_point = check_covered_point
    (scalar_grid, covered_grid, isovert, gcube_index);

	if(flag_covered_point)
	{   
		isovert.gcube_list[gcube_index].flag = COVERED_POINT;
		return;
	}

	bool flag_check_angle = isovert_param.flag_check_triangle_angle;
	bool triangle_flag =
		creates_triangle(scalar_grid, flag_check_angle, isovert, cube_index,
                     isovalue, bin_grid, bin_width, v1, v2);

	if (!triangle_flag) {
    select_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
       isovert_param, gcube_index, flag, isovert);
	}
	else
	{
		isovert.gcube_list[gcube_index].flag = UNAVAILABLE_GCUBE;
	}

}


/// Unselect cube grid_cube.  Uncover all neighbors.
/// @pre Neigboring cubes are covered only grid_cube.
void unselect_cube
(const SHARPISO_GRID_NEIGHBORS & grid,
 const AXIS_SIZE_TYPE bin_width,
 const NUM_TYPE gcube_index,
 SHARPISO_BOOL_GRID & covered_grid,
 BIN_GRID<int> & bin_grid,
 ISOVERT & isovert)
{
  VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE;
  bin_grid_remove(grid, bin_width, cube_index, bin_grid);

  if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < grid.NumVertexNeighborsC(); k++) {

      VERTEX_INDEX cubeB_index = grid.VertexNeighborC(cube_index, k);
      INDEX_DIFF_TYPE gcubeB_index = isovert.GCubeIndex(cubeB_index);

      covered_grid.Set(cubeB_index, false);
      if (gcubeB_index != ISOVERT::NO_INDEX) {

        if (isovert.gcube_list[gcubeB_index].num_eigenvalues > 1 &&
            !isovert.gcube_list[gcubeB_index].flag_centroid_location) {
          isovert.gcube_list[gcubeB_index].flag = AVAILABLE_GCUBE;
        }
        else {
          isovert.gcube_list[gcubeB_index].flag = SMOOTH_GCUBE;
        }
        isovert.gcube_list[gcubeB_index].covered_by = cubeB_index;
      }
    }
  }
  else {
    // *** NOT YET IMPLEMENTED ***
  }
}

/*
* Select corner cubes (eigen value 3)
*/
void select_corner_cubes	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const vector<NUM_TYPE> & sortd_ind2gcube_list,
	ISOVERT &isovert)
{
	const COORD_TYPE linf_dist_threshold = 
		isovert_param.linf_dist_thresh_merge_sharp;
	const int bin_width = isovert_param.bin_width;

	for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    NUM_TYPE gcube_index = sortd_ind2gcube_list[ind];

		GRID_CUBE c;
		c = isovert.gcube_list[gcube_index];

		// check boundary
		if (c.boundary_bits == 0)

			// select corners first
				if (isovert.isFlag
            (cube_ind_frm_gc_ind(isovert, sortd_ind2gcube_list[ind]), AVAILABLE_GCUBE)
            && c.linf_dist < linf_dist_threshold && c.num_eigenvalues > 2) 
				{
					check_and_select_cube
						(scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_CORNER_GCUBE, isovert);
				}
	}

}

/*
* Select edge cubes (eigenvalue 2)
*/

/// Select edge cubes (eigenvalue 2)
void select_edge_cubes	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const vector<NUM_TYPE> & sortd_ind2gcube_list,
	ISOVERT & isovert)
{
	const int bin_width = isovert_param.bin_width;
  COORD_TYPE linf_dist = isovert_param.linf_dist_thresh_merge_sharp;

	for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    NUM_TYPE gcube_index = sortd_ind2gcube_list[ind];

		GRID_CUBE c;
		c = isovert.gcube_list[gcube_index];

		// check boundary
		if (c.boundary_bits == 0)
      if (c.flag == AVAILABLE_GCUBE &&
          c.linf_dist < linf_dist  && c.num_eigenvalues == 2 &&
          !c.flag_conflict)
				{
					check_and_select_cube
						(scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_A_GCUBE, isovert);
				}
	}
}

/// Select one edge cube (eigenvalue 2).
/// @param from_list List of gcube indices sorted by increasing
///    distance from isovert_coord to cube centers.
///    Select cube from list from_list.
bool select_one_edge_cube	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const vector<NUM_TYPE> & from_list,
  INDEX_DIFF_TYPE & selected_gcube_index,
	ISOVERT & isovert)
{
	const int bin_width = isovert_param.bin_width;
  COORD_TYPE linf_dist = isovert_param.linf_dist_thresh_merge_sharp;

  selected_gcube_index = ISOVERT::NO_INDEX;

	for (int ind=0; ind < from_list.size(); ind++) {

    NUM_TYPE gcube_index = from_list[ind];

		GRID_CUBE c;
		c = isovert.gcube_list[gcube_index];

		if (c.boundary_bits == 0)
				if (c.flag == AVAILABLE_GCUBE &&
            c.linf_dist < linf_dist && c.num_eigenvalues == 2 &&
            !c.flag_conflict)
				{
					check_and_select_cube
						(scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_A_GCUBE, isovert);

          if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {

            // *** DEBUG ***
            /*
            cerr << "  Selected new cube: " << c.cube_index << "  ";
            ijkgrid_output_vertex_coord(cerr, scalar_grid, c.cube_index);
            cerr << " Linf dist: " << c.linf_dist << endl;
            cerr << "  isovert_position: ";
            IJK::print_coord3D
              (cerr, isovert.gcube_list[gcube_index].isovert_coord);
            cerr << endl;
            */

            selected_gcube_index = gcube_index;
            return(true); 
          }
				}
	}

  return(false);
}

/// Get corner or edge cubes around given cube.
void get_corner_or_edge_cubes_around_cube
(const SHARPISO_GRID_NEIGHBORS & grid,
 const ISOVERT & isovert,
 const NUM_TYPE gcube_index,
 std::vector<NUM_TYPE> & neighbor_list)
{
  const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);

  neighbor_list.clear();
  if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

    for (NUM_TYPE k = 0; k < grid.NumVertexNeighborsC(); k++) {

      VERTEX_INDEX adjacent_cube_index = 
        grid.VertexNeighborC(cube_index, k);

      INDEX_DIFF_TYPE adjacent_gcube_index = 
        isovert.GCubeIndex(adjacent_cube_index);

      if (adjacent_gcube_index != ISOVERT::NO_INDEX) {
        if (isovert.gcube_list[adjacent_gcube_index].num_eigenvalues > 1) {
          neighbor_list.push_back(adjacent_gcube_index);
        }
      }
    }
  }
}

/// Reselect two edge cubes around gcube_index.
/// Select one cube from list from_list.
/// @param from_list List of gcube indices sorted by increasing
///    distance from isovert_coord to cube centers.
void reselect_two_edge_cubes (
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
  const VERTEX_INDEX gcube_index,
  const std::vector<NUM_TYPE> from_list,
	ISOVERT & isovert)
{
  const int bin_width = isovert_param.bin_width;
  vector<NUM_TYPE> neighbor_list;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);

  get_corner_or_edge_cubes_around_cube
    (gridn, isovert, gcube_index, neighbor_list);
  sort(neighbor_list.begin(), neighbor_list.end(), gcube_compare);

  unselect_cube
    (gridn, bin_width, gcube_index, covered_grid, bin_grid, isovert);

  INDEX_DIFF_TYPE gcubeB_index, gcubeC_index;
  if (select_one_edge_cube(scalar_grid, covered_grid, bin_grid, gridn, isovalue,
                           isovert_param, from_list, gcubeB_index, isovert)) {

    bool flag = select_one_edge_cube
      (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
       isovert_param, neighbor_list, gcubeC_index, isovert);

    // *** DEBUG ***
    /*
    using namespace std;
    if (!flag) {
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
      cerr << "*** Reselection of second cube failed. ***" << endl;
      cerr << "  Original cube: " << cube_index << "  ";
      ijkgrid_output_vertex_coord(cerr, gridn, cube_index);
      cerr << endl;
      cerr << "  Neighbor list: " << endl;
      for (int i = 0; i < neighbor_list.size(); i++) {
        NUM_TYPE gcubeB_index = neighbor_list[i];
        VERTEX_INDEX cubeB_index = isovert.CubeIndex(gcubeB_index);
        cerr << " " << cubeB_index << " ";
        ijkgrid_output_vertex_coord(cerr, scalar_grid, cubeB_index);
        cerr << " flag centroid: " 
             << int(isovert.gcube_list[gcubeB_index].flag_centroid_location);
        cerr << " num_eigen: " 
             << int(isovert.gcube_list[gcubeB_index].num_eigenvalues);
        cerr << " flag: "
             << int(isovert.gcube_list[gcubeB_index].flag);
        cerr << " Linf dist: " << isovert.gcube_list[gcubeB_index].linf_dist;
        cerr << endl;
        cerr << "    flag_conflict: "
             << int(isovert.gcube_list[gcubeB_index].flag_conflict);
        cerr << "  isovert_coord: ";
        IJK::print_coord3D
          (cerr, isovert.gcube_list[gcubeB_index].isovert_coord);
        cerr << endl;
      }
      cerr << endl;
    }
    */
  }
  else {

    // Reselect cube gcube_index
    select_cube(scalar_grid, covered_grid, bin_grid, gridn, isovalue,
                isovert_param, gcube_index, COVERED_A_GCUBE, isovert);

    // *** DEBUG ***
    /*
    VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    cerr << "*** Reselection around cube " << cube_index << " ";
    ijkgrid_output_vertex_coord(cerr, scalar_grid, cube_index);
    cerr << " failed." << endl;
    cerr << "  From list: ";
    for (int i = 0; i < from_list.size(); i++) {
        cerr << " " << isovert.CubeIndex(from_list[i]) << " ";
        ijkgrid_output_vertex_coord
          (cerr, scalar_grid, isovert.CubeIndex(from_list[i]));
    }
    cerr << endl;
    */
  }

}

/// Reselect edge cubes around gcube_index.
/// @param nearby_sharp_list List of nearby sharp cubes.
void reselect_edge_cubes (
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
  const VERTEX_INDEX gcube_index,
  const std::vector<VERTEX_INDEX> & nearby_selected_list,
	ISOVERT & isovert
)
{
  const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
  GRID_COORD_TYPE cube_coord[DIM3], cubeB_coord[DIM3], cubeC_coord[DIM3];
  GRID_COORD_TYPE rmin[DIM3], rmax[DIM3];
  std::vector<NUM_TYPE> from_list;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);

  scalar_grid.ComputeCoord(cube_index, cube_coord);

  int max_overlap_dim = -1;
  VERTEX_INDEX max_overlap_cube = 0;
  for (NUM_TYPE i = 0; i < nearby_selected_list.size(); i++) {
    VERTEX_INDEX cubeB_index = nearby_selected_list[i];

    if (cube_index != cubeB_index) {

      scalar_grid.ComputeCoord(cubeB_index, cubeB_coord);

      int overlap_dim;
      find_3x3x3_overlap(cube_coord, cubeB_coord, rmin, rmax, overlap_dim);
      if (overlap_dim > max_overlap_dim) { 
        max_overlap_dim = overlap_dim; 
        max_overlap_cube = cubeB_index;
      }
    }
  }

  if (max_overlap_dim >= 3) { 
    // Some cube is covered by both gcube_index and another selected cube.
    return; 
  }
  else {

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "Reselecting around cube " << cube_index << ".  ";
    IJK::print_coord3D(cerr, cube_coord);
    cerr << "  linf_dist: " << isovert.gcube_list[gcube_index].linf_dist
         << endl;
    cerr << "  Overlaps cube: " << max_overlap_cube << "  ";
    ijkgrid_output_vertex_coord(cerr, gridn, max_overlap_cube);
    cerr << "  overlap_dim: " << max_overlap_dim << endl;
    */

    if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

      scalar_grid.ComputeCoord(max_overlap_cube, cubeB_coord);

      for (NUM_TYPE k = 0; k < gridn.NumVertexNeighborsC(); k++) {

        VERTEX_INDEX cubeC_index = gridn.VertexNeighborC(cube_index, k);
        INDEX_DIFF_TYPE gcubeC_index = isovert.GCubeIndex(cubeC_index);

        if (gcubeC_index != ISOVERT::NO_INDEX) {

          scalar_grid.ComputeCoord(cubeC_index, cubeC_coord);

          int overlap_dim;
          find_3x3x3_overlap(cubeB_coord, cubeC_coord, rmin, rmax, overlap_dim);
          if (overlap_dim >= 3) 
            { from_list.push_back(gcubeC_index); }
        }
      }

      if (from_list.size() > 0) {
        sort(from_list.begin(), from_list.end(), gcube_compare);

        reselect_two_edge_cubes
          (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
           isovert_param, gcube_index, from_list, isovert);
      }
    }
    else {
      // NOT YET IMPLEMENTED.
    }
  }

}

/// Reselect edge cubes around gcube_index.
/// @param nearby_selected_list List of nearby sharp cubes.
void reselect_edge_cubes (
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID & covered_grid,
	BIN_GRID<VERTEX_INDEX> & bin_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const vector<NUM_TYPE> & sortd_ind2gcube_list,
	ISOVERT & isovert
)
{
	const int bin_width = isovert_param.bin_width;
  std::vector<VERTEX_INDEX> nearby_selected_list;

  for (NUM_TYPE i = 0; i < sortd_ind2gcube_list.size(); i++) {

    NUM_TYPE gcube_index = sortd_ind2gcube_list[i];
    if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE &&
        isovert.gcube_list[gcube_index].num_eigenvalues == 2 &&
        !isovert.gcube_list[gcube_index].flag_near_corner) {

      VERTEX_INDEX cube_index = isovert.gcube_list[gcube_index].cube_index;

      get_selected(scalar_grid, cube_index, bin_grid, bin_width, 
                   nearby_selected_list);

      reselect_edge_cubes
        (scalar_grid, covered_grid, bin_grid, gridn, isovalue,
         isovert_param, gcube_index, nearby_selected_list, isovert);
    }
  }
}


/*
* is th neighbor of type flag
*/
bool is_neighbor(	
	GRID_CUBE c,
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_GRID_NEIGHBORS & gridn,
	ISOVERT &isovert,
	GRID_CUBE_FLAG flag
	)
{
	VERTEX_INDEX neighbor_cube_index;
	for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
	{
		neighbor_cube_index = gridn.CubeNeighborF(c.cube_index, j);
		//neighbor_gcube_index is an entry into the gcube list 
		INDEX_DIFF_TYPE neighbor_gcube_index
			= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

		if(neighbor_gcube_index == ISOVERT::NO_INDEX)
		{
			continue;
		}
		if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
		{
			return true;
		}
	}

	//edge
	for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
	{
		neighbor_cube_index = gridn.CubeNeighborE(c.cube_index, j);
		//neighbor_gcube_index is an entry into the gcube list 
		INDEX_DIFF_TYPE neighbor_gcube_index
			= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

		if(neighbor_gcube_index == ISOVERT::NO_INDEX)
		{
			continue;
		}
		if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
		{
			return true;
		}
	}

	//edge
	for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
	{
		neighbor_cube_index = gridn.CubeNeighborV(c.cube_index, j);
		//neighbor_gcube_index is an entry into the gcube list 
		INDEX_DIFF_TYPE neighbor_gcube_index
			= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);

		if(neighbor_gcube_index == ISOVERT::NO_INDEX)
		{
			continue;
		}
		if(isovert.gcube_list[neighbor_gcube_index].flag == flag)
		{
			return true;
		}
	}

	return false;
}


/*
* Select near corners. neighbor of covered corner
*/
void select_cubes_near_corners	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	SHARPISO_BOOL_GRID &covered_grid,
	BIN_GRID<VERTEX_INDEX> &bin_grid,
	SHARPISO_GRID_NEIGHBORS &gridn,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	vector<NUM_TYPE> sortd_ind2gcube_list,
	ISOVERT &isovert)
{
	const COORD_TYPE linf_dist_threshold = 
		isovert_param.linf_dist_thresh_merge_sharp;
	const int bin_width = isovert_param.bin_width;


	for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    NUM_TYPE gcube_index = sortd_ind2gcube_list[ind];

		GRID_CUBE c;
		c = isovert.gcube_list[gcube_index];

    if (isovert.isFlag(cube_ind_frm_gc_ind(isovert, sortd_ind2gcube_list[ind]), AVAILABLE_GCUBE)) {

      // check boundary
      if(c.boundary_bits == 0 && c.linf_dist <= 0.5 ) {
        bool flag = is_neighbor
          (c, scalar_grid, gridn,  isovert,  COVERED_CORNER_GCUBE );

        if (flag)
          check_and_select_cube
            (scalar_grid, covered_grid, bin_grid, gridn, isovalue, 
             isovert_param, gcube_index, COVERED_A_GCUBE, isovert );

        if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) 
          { isovert.gcube_list[gcube_index].flag_near_corner = true; }
      }
    }
	}
}

// Set gcube_list[i].covered_by to gcube_list[i].cube_index for each i.
void initialize_covered_by(ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) 
    { isovert.gcube_list[i].covered_by = isovert.gcube_list[i].cube_index; }
}


void MERGESHARP::select_sharp_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert)
{
	const int dimension = scalar_grid.Dimension();
	const int bin_width = isovert_param.bin_width;
  GCUBE_COMPARE gcube_compare(isovert.gcube_list);
	std::vector<NUM_TYPE> sharp_gcube_list;

	// get corner or edge cubes
	get_corner_or_edge_cubes(isovert.gcube_list, sharp_gcube_list);

	BIN_GRID<VERTEX_INDEX> bin_grid;
	init_bin_grid(scalar_grid, bin_width, bin_grid);

  initialize_covered_by(isovert);

	SHARPISO_GRID_NEIGHBORS gridn;
	gridn.SetSize(scalar_grid);

	//Keeps track of all cubes which are covered.
	SHARPISO_BOOL_GRID covered_grid;
	covered_grid.SetSize(scalar_grid);
	covered_grid.SetAll(false);

  // *** DEBUG ***
  /*
  using namespace std;
  flag_debug = true;
  */

  // *** DEBUG ***
  if (flag_debug) {
    output_gcube_list(scalar_grid, isovert, sharp_gcube_list);
    cerr << endl << "*** Selecting corner cubes." << endl;
  }

	// select corner cubes
	select_corner_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  reset_covered_isovert_positions(gridn, covered_grid, isovert);

  // Resort sharp gcube_list
  sort(sharp_gcube_list.begin(), sharp_gcube_list.end(), gcube_compare);

  // *** DEBUG ***
  if (flag_debug) { 
    output_gcube_list(scalar_grid, isovert, sharp_gcube_list);
    cerr << endl << "*** Selecting near corner cubes." << endl; 
  }

	// select cubes near corners
	select_cubes_near_corners
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);

  // *** DEBUG ***
  if (flag_debug) { cerr << endl << "*** Selecting edge cubes." << endl; }

	// select edge cubes.
	select_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param,
     sharp_gcube_list, isovert);

  // reselect edge cubes
	reselect_edge_cubes
    (scalar_grid, covered_grid, bin_grid, gridn, isovalue, isovert_param, 
     sharp_gcube_list, isovert);
}


/**
Sets the sharp_ind_grid and gcube_list

Traverse the scalar grid,
set everything to be ISOVERT::NO_INDEX
if isosurface intersects the cube then set *index*
push the corresponding *grid_cube* into the *gcube_list*
*/
void create_active_cubes (
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	ISOVERT &isovert)
{
	NUM_TYPE index = 0;
	//set the size of sharp index grid
	isovert.sharp_ind_grid.SetSize(scalar_grid);
	IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX)
	{
		if (is_gt_cube_min_le_cube_max(scalar_grid, iv, isovalue))
		{
			index = isovert.gcube_list.size();
			isovert.sharp_ind_grid.Set(iv,index);
			GRID_CUBE gc;
			gc.cube_index = iv;
			isovert.gcube_list.push_back(gc);
		}
		else
		{
			// if it is not intersected then mark as ISOVERT::NO_INDEX
			isovert.sharp_ind_grid.Set(iv,ISOVERT::NO_INDEX);
		}
	}
}

/// process edge called from are connected
void process_edge(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
				  const VERTEX_INDEX v0,
				  const VERTEX_INDEX v1,
				  const SCALAR_TYPE isovalue,
				  bool &is_intersect)
{
	is_intersect = is_gt_min_le_max(scalar_grid, v0, v1, isovalue);
}

/**
* Compute dual isovert.
*/
void MERGESHARP::compute_dual_isovert
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	const VERTEX_POSITION_METHOD vertex_position_method,
	ISOVERT & isovert,
	ISOVERT_INFO & isovert_info)
{
	IJK::PROCEDURE_ERROR error("compute_dual_isovert");

	if (!gradient_grid.Check
		(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	create_active_cubes(scalar_grid, isovalue, isovert);

	if (vertex_position_method == GRADIENT_POSITIONING) {

		compute_isovert_positions 
			(scalar_grid, gradient_grid, isovalue, isovert_param, 
			isovert, isovert_info);

    swap_isovert_positions(scalar_grid, isovert);
	}
	else {
		compute_isovert_positions_edgeI
			(scalar_grid, gradient_grid, isovalue, isovert_param, 
			vertex_position_method, isovert, isovert_info);
	}
}

void MERGESHARP::compute_dual_isovert
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const std::vector<COORD_TYPE> & edgeI_coord,
	const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT &isovert)
{
	create_active_cubes(scalar_grid, isovalue, isovert);

	compute_isovert_positions 
		(scalar_grid, edgeI_coord, edgeI_normal_coord,
		isovalue, isovert_param, isovert);
}


// Select all grid cubes which are not smooth.
void MERGESHARP::select_non_smooth(ISOVERT & isovert)
{
	for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
		if (isovert.gcube_list[i].flag != SMOOTH_GCUBE) 
		{ isovert.gcube_list[i].flag = SELECTED_GCUBE; };
	}
}

// Get list of grid cubes from isovert.
void MERGESHARP::get_cube_list
	(const ISOVERT & isovert, std::vector<VERTEX_INDEX> & cube_list)
{
	const NUM_TYPE num_gcube = isovert.gcube_list.size();

	cube_list.resize(num_gcube);
	for (NUM_TYPE i = 0; i < num_gcube; i++) {
		cube_list[i] = isovert.gcube_list[i].cube_index;
	}
}

// Store isosurface lookup table index in gcube_list.
void MERGESHARP::store_table_index
	(const std::vector<IJKDUALTABLE::TABLE_INDEX> & table_index,
	GRID_CUBE_ARRAY & gcube_list)
{
	IJK::PROCEDURE_ERROR error("store_table_index");

	if (table_index.size() != gcube_list.size()) {
		error.AddMessage("Programming error.  Numbers of elements in table_index and gcube_list differ.");
		error.AddMessage("  table_index.size() = ", table_index.size(), ".");
		error.AddMessage("  gcube_list.size() = ", gcube_list.size(), ".");
		throw error;
	}

	for (NUM_TYPE i = 0; i < table_index.size(); i++) 
	{ gcube_list[i].table_index = table_index[i]; }
}

// **************************************************
// GRID_CUBE member functions
// **************************************************

void GRID_CUBE::Init()
{
  num_eigenvalues = 0;
  flag = AVAILABLE_GCUBE;
  boundary_bits = 0;
  cube_index = 0;
  linf_dist = 0;
  flag_centroid_location = false;
  flag_conflict = false;
  flag_near_corner = false;
  cube_containing_isovert = 0;
  flag_coord_from_other = false;
  table_index = 0;
  covered_by = 0;
}

bool GRID_CUBE::IsCoveredOrSelected() const
{
  if (flag == COVERED_A_GCUBE || flag == COVERED_B_GCUBE ||
      flag == COVERED_CORNER_GCUBE || flag == SELECTED_GCUBE) 
    { return true; }
  else
    { return false; }
}

// **************************************************
// ISOVERT member functions
// **************************************************

bool ISOVERT::isFlag(const int cube_index,  GRID_CUBE_FLAG _flag){
	if (isActive(cube_index)){
		if (gcube_list[sharp_ind_grid.Scalar(cube_index)].flag == _flag)
			return true;
		else
			return false;
	}
	else
		return false;
}

bool ISOVERT::isActive(const int cube_index){
	if (sharp_ind_grid.Scalar(cube_index) != NO_INDEX)
		return true;
	else
		return false;
}

// **************************************************
// ISOVERT_INFO member functions
// **************************************************

void ISOVERT_INFO::Clear()
{
	num_sharp_corners = 0;
	num_sharp_edges = 0;
	num_smooth_vertices = 0;
	num_merged_iso_vertices = 0;
	num_conflicts = 0;
	num_Linf_iso_vertex_locations = 0;
}

// **************************************************
// Set ISOVERT_INFO
// **************************************************

/// Count number of vertices on sharp corners or sharp edges.
/// Count number of smooth vertices.
void MERGESHARP::count_vertices
	(const ISOVERT & isovert, ISOVERT_INFO & isovert_info)
{
	isovert_info.num_sharp_corners = 0;
	isovert_info.num_sharp_edges = 0;
	isovert_info.num_smooth_vertices = 0;
	for (VERTEX_INDEX i = 0; i < isovert.gcube_list.size(); i++) {
		if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
			if (isovert.gcube_list[i].num_eigenvalues == 2)
			{ isovert_info.num_sharp_edges++; }
			else if (isovert.gcube_list[i].num_eigenvalues == 3)
			{ isovert_info.num_sharp_corners++; }
		}
		else if (isovert.gcube_list[i].flag == SMOOTH_GCUBE) {
			isovert_info.num_smooth_vertices++;
		}
	}
}

// **************************************************
// Local routines
// **************************************************

namespace {


	// Find the linf distance between the sharp vertex for the cube
	// and the center of the cube.
	void compute_linf_dist
  (const SHARPISO_GRID & grid, const VERTEX_INDEX iv,
   const COORD_TYPE isovert_coord[DIM3], SCALAR_TYPE & linf_dist)
	{
		SCALAR_TYPE squareDist=0.0;
		COORD_TYPE cc[DIM3];

		grid.ComputeCubeCenterScaledCoord(iv,cc);
		SCALAR_TYPE temp_d   = 0.0;
		SCALAR_TYPE max_dist = -1.0;

		for (int d=0; d<DIM3; d++) {
			temp_d = abs(isovert_coord[d]-cc[d])/grid.Spacing(d);
			if (temp_d > max_dist)
			{ max_dist=temp_d; }
		}
		linf_dist = max_dist;
	}

	// Find the linf distance between the sharp vertex for the cube
	// and the center of the cube.
	void compute_linf_dist
		(const SHARPISO_GRID & grid, const NUM_TYPE gcube_index, ISOVERT & isovert)
  {
    const VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
    compute_linf_dist(grid, cube_index, 
                      isovert.gcube_list[gcube_index].isovert_coord,
                      isovert.gcube_list[gcube_index].linf_dist);
  }

  /// Compute the overlap region between two 3x3x3 regions around cubes.
  /// @param[out] Overlap dimension.  No overlap has overlap dimension -1.
  bool find_3x3x3_overlap
  (const GRID_COORD_TYPE cubeA_coord[DIM3],
   const GRID_COORD_TYPE cubeB_coord[DIM3],
   GRID_COORD_TYPE rmin[DIM3],
   GRID_COORD_TYPE rmax[DIM3],
   int & overlap_dim)
  {
    for (int d=0;d<DIM3;d++){
      rmin[d]=max(cubeA_coord[d]-1, cubeB_coord[d]-1);
      rmax[d]=min(cubeA_coord[d]+2, cubeB_coord[d]+2);
    }

    overlap_dim = 0;

    for (int d=0; d<DIM3; d++) {

      if(rmin[d] > rmax[d]) {
        overlap_dim = -1;
        return(false);
      }

      if (rmin[d] < rmax[d]) 
        { overlap_dim++; }
    }

    return(true);
  }

  /// Compute the overlap region between two 3x3x3 regions around cubes.
  /// @param[out] Overlap dimension.  No overlap has overlap dimension -1.
  template <typename GRID_TYPE>
  bool find_3x3x3_overlap
  (const GRID_TYPE & grid,
   const VERTEX_INDEX cubeA_index,
   const VERTEX_INDEX cubeB_index,
   GRID_COORD_TYPE rmin[DIM3],
   GRID_COORD_TYPE rmax[DIM3],
   int & overlap_dim)
  {
    GRID_COORD_TYPE cubeA_coord[DIM3], cubeB_coord[DIM3];
    bool is_overlap;

    grid.ComputeCoord(cubeA_index, cubeA_coord);
    grid.ComputeCoord(cubeB_index, cubeB_coord);

    is_overlap = 
      find_3x3x3_overlap(cubeA_coord, cubeB_coord, rmin, rmax, overlap_dim);

    return(is_overlap);
  }

	// Return true if two  cube-indices are connected
	bool are_connected 
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2,
		const SCALAR_TYPE isovalue )
	{
		// find the overlap region
		GRID_COORD_TYPE rmin[DIM3], rmax[DIM3];

    int overlap_dim;
    find_3x3x3_overlap
      (scalar_grid, cube_index1, cube_index2, rmin, rmax, overlap_dim);

		if (overlap_dim < 2) { return false; }

		COORD_TYPE vbase = scalar_grid.ComputeVertexIndex(rmin);
		COORD_TYPE base_coord[DIM3];
		VERTEX_INDEX base_index;
		VERTEX_INDEX vnext;
		COORD_TYPE diff[DIM3]={0.0,0.0,0.0};
		//visit each edge along the way
		VERTEX_INDEX v0,v1,v2;

		int d0=0,d1=0,d2=0;
		bool is_intersect=false;

		int num=0;
		// if there is an overlap
		for(int d0=0;d0<DIM3;d0++){
			d1 = (d0+1)%3;
			d2 = (d0+2)%3;
			v1 = vbase;
			for (int i1=rmin[d1]; i1 <=rmax[d1];i1++){
				v2=v1;
				for (int i2=rmin[d2];i2<=rmax[d2];i2++){
					v0=v2;
					for(int i0=rmin[d0];i0<rmax[d0];i0++){
						vnext = scalar_grid.NextVertex(v0,d0);
						process_edge(scalar_grid, v0,vnext,isovalue, is_intersect);
						num++;
						if (is_intersect)
							return true;
						v0=vnext;
					}
					v2=scalar_grid.NextVertex(v2,d2);
				}
				v1=scalar_grid.NextVertex(v1,d1);
			}
		}

		return false;
	}

	// Compute the cosine of the angle between (v2,v1) and (v2,v3)
	void compute_cos_angle
		(const ISOVERT & isovert,
		const VERTEX_INDEX gcube_list_index_v1,
		const VERTEX_INDEX gcube_list_index_v2,
		const VERTEX_INDEX gcube_list_index_v3,
		SCALAR_TYPE & cos_angle)
	{
		COORD_TYPE coord0[DIM3], coord1[DIM3],
			coord2[DIM3],vec12[DIM3], vec32[DIM3];

		subtract_coord(DIM3, isovert.gcube_list[gcube_list_index_v2].isovert_coord,
			isovert.gcube_list[gcube_list_index_v1].isovert_coord, vec12);
		normalize_vector(DIM3, vec12, 0.001, vec12);

		subtract_coord(DIM3, isovert.gcube_list[gcube_list_index_v2].isovert_coord,
			isovert.gcube_list[gcube_list_index_v3].isovert_coord, vec32);
		normalize_vector(DIM3, vec32, 0.001, vec32);

		compute_inner_product(DIM3, vec12, vec32, cos_angle);
	}


	// Return true if triangle(iv,v1,v2) has a large angle
	// iv, v1 and v2 are cube_indices.

	bool is_angle_large
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const ISOVERT &isovert,
		const VERTEX_INDEX iv,
		const SCALAR_TYPE threshold,
		const VERTEX_INDEX v1,
		const VERTEX_INDEX v2)
	{

		SCALAR_TYPE cos_angle_iv, cos_angle_v1, cos_angle_v2;
		VERTEX_INDEX gcube_list_index_v1;
		VERTEX_INDEX gcube_list_index_v2;
		VERTEX_INDEX gcube_list_index_v3;

		gcube_list_index_v1=isovert.sharp_ind_grid.Scalar(v1);
		gcube_list_index_v2=isovert.sharp_ind_grid.Scalar(v2);
		gcube_list_index_v3=isovert.sharp_ind_grid.Scalar(iv);


		compute_cos_angle(isovert, gcube_list_index_v1,
			gcube_list_index_v3,  gcube_list_index_v2, cos_angle_iv);

		compute_cos_angle(isovert, gcube_list_index_v3,
			gcube_list_index_v1,  gcube_list_index_v2, cos_angle_v1);

		compute_cos_angle(isovert, gcube_list_index_v3,
			gcube_list_index_v2,  gcube_list_index_v1, cos_angle_v2);

		if ((cos_angle_iv < threshold) || (cos_angle_v1 < threshold)
			|| (cos_angle_v2 < threshold))
		{
			return true;
		}
		else
			return false;
	}

  // Return true if some facet neighbor is covered by a sharp cube.
  bool is_facet_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, VERTEX_INDEX & cube_index1)
  {
    cube_index1 = 0;
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborF(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.sharp_ind_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        if (isovert.gcube_list[index_gcube].covered_by != icube) {
          cube_index1 = isovert.gcube_list[index_gcube].covered_by;
          return(true);
        }
      }
    }
    return(false);
  }

  // Return true if some facet or edge neighbor is covered by a sharp cube
  //   other than cube_index1.
  // @pre cube_index0 is not on grid boundary.
  bool is_facet_or_edge_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   VERTEX_INDEX & cube_index2)
  {
    cube_index2 = 0;

    // Check facet neighbors
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborF(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.sharp_ind_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1) {
          cube_index2 = covered_by;
          return(true);
        }
      }
    }

    // Check edge neighbors
    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
      VERTEX_INDEX icube = grid.CubeNeighborE(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.sharp_ind_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1) {
          cube_index2 = covered_by;
          return(true);
        }
      }
    }

    return(false);
  }

  // Return true if some facet, edge or vertex neighbor is covered 
  //   by a sharp cube other than cube_index1 or cube_index2
  // @pre cube_index0 is not on grid boundary.
  bool is_neighbor_covered
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1,
   const VERTEX_INDEX & cube_index2)
  {
    for (NUM_TYPE j = 0; j < grid.NumVertexNeighborsC(); j++) {
      VERTEX_INDEX icube = grid.VertexNeighborC(cube_index0, j);

      VERTEX_INDEX index_gcube = isovert.sharp_ind_grid.Scalar(icube);

      if (index_gcube != ISOVERT::NO_INDEX) {
        VERTEX_INDEX covered_by = isovert.gcube_list[index_gcube].covered_by;
        if (covered_by != icube && covered_by != cube_index1 && 
            covered_by != cube_index2) {
          return(true);
        }
      }
    }
    return(false);
  }

  // Return true if cube or some cube neighbor of cube_index0 
  //   is covered by cube_index1
  bool is_cube_or_neighbor_covered_by
  (const SHARPISO_GRID_NEIGHBORS & grid, const ISOVERT & isovert,
   const VERTEX_INDEX cube_index0, const VERTEX_INDEX cube_index1)
  {
    VERTEX_INDEX gcube_index;

    gcube_index = isovert.sharp_ind_grid.Scalar(cube_index0);
    if (gcube_index != ISOVERT::NO_INDEX) {
      VERTEX_INDEX covered_by = isovert.gcube_list[gcube_index].covered_by;
      if (covered_by == cube_index1) { return(true); }
    }

    if (isovert.gcube_list[gcube_index].boundary_bits == 0) {

      for (NUM_TYPE j = 0; j < grid.NumVertexNeighborsC(); j++) {
        VERTEX_INDEX icube = grid.VertexNeighborC(cube_index0, j);

        gcube_index = isovert.sharp_ind_grid.Scalar(icube);

        if (gcube_index != ISOVERT::NO_INDEX) {
          VERTEX_INDEX covered_by = isovert.gcube_list[gcube_index].covered_by;
          if (covered_by == cube_index1) { return(true); }
        }
      }
    }
    else {
      // *** NOT YET IMPLEMENTED. ***
    }

    return(false);
  }

	// Store information from singular valued decomposition.
  void store_svd_info
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX cube_index,
		const NUM_TYPE gcube_index,
		const NUM_TYPE num_large_eigenvalues,
		const SVD_INFO svd_info,
		ISOVERT & isovert, ISOVERT_INFO & isovert_info)
	{
		if (svd_info.flag_conflict) 
		{ isovert_info.num_conflicts++; }
		else if (svd_info.flag_Linf_iso_vertex_location)
		{ isovert_info.num_Linf_iso_vertex_locations++; }

		// set number of eigenvalues.
		isovert.gcube_list[gcube_index].num_eigenvalues =
			(unsigned char) num_large_eigenvalues;

		// set the sharp vertex type to be *AVAILABLE*.
		if (num_large_eigenvalues > 1 && svd_info.location == LOC_SVD)
      { isovert.gcube_list[gcube_index].flag = AVAILABLE_GCUBE; }
		else 
      { isovert.gcube_list[gcube_index].flag = SMOOTH_GCUBE; }

    if (svd_info.location == CENTROID) 
      { isovert.gcube_list[gcube_index].flag_centroid_location = true; }
    else
      { isovert.gcube_list[gcube_index].flag_centroid_location = false; }

    isovert.gcube_list[gcube_index].flag_conflict = svd_info.flag_conflict;

		// store distance.
		compute_linf_dist(scalar_grid, cube_index, 
			isovert.gcube_list[gcube_index].isovert_coord,
			isovert.gcube_list[gcube_index].linf_dist);
	}

	// Set cube_containing_isovert
  void set_cube_containing_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue, ISOVERT & isovert)
  {
    COORD_TYPE cube_coord[DIM3];
    VERTEX_INDEX cube_containing_point;

    for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {

      const COORD_TYPE * point_coord = isovert.gcube_list[i].isovert_coord;
      VERTEX_INDEX cube_index = isovert.CubeIndex(i);
      scalar_grid.ComputeCoord(cube_index, cube_coord);

      if (cube_contains_point
          (cube_coord, point_coord, scalar_grid.SpacingPtrConst()) ||
          !scalar_grid.ContainsPoint(point_coord)) {

        isovert.gcube_list[i].cube_containing_isovert = cube_index;
        isovert.gcube_list[i].flag_conflict = false;
      }
      else {

        bool flag_boundary;
        bool flag_active;

        get_cube_containing_point
          (scalar_grid, isovalue, point_coord, cube_containing_point,
           flag_boundary, flag_active);

        isovert.gcube_list[i].cube_containing_isovert = 
          cube_containing_point;

        if (flag_active)
          { isovert.gcube_list[i].flag_conflict = true; }
        else if (!flag_boundary) 
          { isovert.gcube_list[i].flag_conflict = false; }
        else {

          VERTEX_INDEX conflicting_cube;
          if (check_conflict(scalar_grid, isovalue, cube_coord, 
                             isovert.gcube_list[i].isovert_coord, 
                             conflicting_cube)) {
            
            isovert.gcube_list[i].cube_containing_isovert = 
              conflicting_cube;
            isovert.gcube_list[i].flag_conflict = true;
          }
          else {
            isovert.gcube_list[i].flag_conflict = false;
          }
        }
      }
    }
  }

  // DEBUG routine
  void output_gcube_list
  (const SHARPISO_GRID & grid, const ISOVERT & isovert, 
   const std::vector<NUM_TYPE> & gcube_index_list)
  {
    for (int i = 0; i < gcube_index_list.size(); i++) {
      NUM_TYPE gcube_index = gcube_index_list[i];
      VERTEX_INDEX cube_index = isovert.CubeIndex(gcube_index);
      cerr << "Cube : " << cube_index << " ";
      ijkgrid_output_vertex_coord(cerr, grid, cube_index);
      cerr << "  isovert_coord: ";
      IJK::print_coord3D(cerr, isovert.gcube_list[gcube_index].isovert_coord);
      cerr << endl;
      cerr << "    linf dist: " << isovert.gcube_list[gcube_index].linf_dist;
      cerr << "  flag conflict: " 
           << int(isovert.gcube_list[gcube_index].flag_conflict);
      cerr << "  flag_coord_from_other: "
           << int(isovert.gcube_list[gcube_index].flag_coord_from_other);
      cerr << "  num_eigen: "
           << int(isovert.gcube_list[gcube_index].num_eigenvalues)
           << endl;
    }
  }

}

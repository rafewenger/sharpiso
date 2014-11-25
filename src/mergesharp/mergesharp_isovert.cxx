/// \file mergesharp_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.

/*
Copyright (C) 2012-2013 Arindam Bhattacharya

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


namespace {


	void compute_linf_dist
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],SCALAR_TYPE & linf_dist);

	bool are_connected
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2, const SCALAR_TYPE isovalue);

	bool is_angle_large
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const ISOVERT & isovertData, const VERTEX_INDEX iv,
		const SCALAR_TYPE threshold, const VERTEX_INDEX v1,
		const VERTEX_INDEX v2);

	// Store information from singular valued decomposition.
	inline void store_svd_info
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

		// store distance.
		compute_linf_dist(scalar_grid, cube_index, 
			isovert.gcube_list[gcube_index].isovert_coord,
			isovert.gcube_list[gcube_index].linf_dist);
	}

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
				//DEBUG
				/*using namespace std;
				cout <<"\n---------------------------\ncube index eigen "<< iv<<","; 
				*/
				svd_compute_sharp_vertex_for_cube_lindstrom
					(scalar_grid, gradient_grid, iv, isovalue, isovert_param, voxel,
					isovert.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);
				store_svd_info(scalar_grid, iv, index, num_large_eigenvalues,
					svd_info, isovert, isovert_info);

				////DEBUG
				//cout <<" eigen values ";

				//if(svd_info.location!=CENTROID)
				//{
				//	for (int j = 0; j < num_large_eigenvalues; j++)
				//	{
				//		cout << eigenvalues[j]<<" ";
				//	}
				//	cout << endl;
				//}
				//else
				//{
				//	cout <<"location computed using centroid.\n";
				//}
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
		//DEBUG
		// recompute the covered point
		if (cube_flag == COVERED_POINT)
		{
			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

			/*cout <<"coveredpoint "<< isovert.gcube_list[i].isovert_coord[0]
			<<" "<<isovert.gcube_list[i].isovert_coord[1]
			<<" "<<isovert.gcube_list[i].isovert_coord[2]<<endl;*/

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


	//DEBUG
	//new computations for COVERED_POINT
	/*for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
	GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

	if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE)) 
	{
	VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
	//DEBUG
	std::cerr <<" RED "<< int (isovert_param.use_sharp_edgeI)<<std::endl;
	compute_edgeI_centroid
	(scalar_grid, gradient_grid, isovalue, cube_index,
	isovert_param.use_sharp_edgeI, isovert.gcube_list[i].isovert_coord);
	isovert.gcube_list[i].flag_centroid_location = true;
	}
	if (cube_flag == COVERED_POINT)
	{
	cout <<"coveredpoint "<< isovert.gcube_list[i].isovert_coord[0]
	<<" "<<isovert.gcube_list[i].isovert_coord[1]
	<<" "<<isovert.gcube_list[i].isovert_coord[2]<<endl;

	//recompute using gradNS 
	isovert.gcube_list[i].flag_centroid_location = true;
	}
	}*/

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

// Create and set the sortd_ind2gcube_list.
void MERGESHARP::sort_gcube_list
	(const std::vector<GRID_CUBE> & gcube_list,
	std::vector<NUM_TYPE> & sortd_ind2gcube_list)
{
  GCUBE_COMPARE gcube_compare(gcube_list);

  for (int i=0;i<gcube_list.size();i++)
    {
      // *** SHOULD CHECK THAT POINT IS NOT COMPUTED USING CENTROID ***
      if (gcube_list[i].num_eigenvalues > 1)
        sortd_ind2gcube_list.push_back(i);
    }

  sort (sortd_ind2gcube_list.begin(),sortd_ind2gcube_list.end(), 
        gcube_compare);
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
			for (int i = 0; i < bin_grid.ListLength(jbin); i++)
			{ selected_list.push_back(bin_grid.List(jbin, i)); }
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
				//Debug
				//cout <<"create triangle v1 "<< v1 <<",create triangle v2 "<< v2 << endl;

				if (check_triangle_angle){
					// checking angle
					bool flag_large_angle = 
						is_angle_large(scalar_grid, isovert, iv, threshold, v1 , v2);

					if (flag_large_angle)
					{ 
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


/**
Check if the sharp vertex is in a cube which is already covered.
if yes then return true else return false.
*/
bool check_covered_point(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SHARPISO_BOOL_GRID &covered_grid,
	ISOVERT &isovert,
	const VERTEX_INDEX &gcube_index)
{
	//stor the coordinated of the cube which has point A
	COORD_TYPE cube_coord[DIM3];
	for (int d = 0; d < DIM3; d++)
	{
		COORD_TYPE c = isovert.gcube_list[gcube_index].isovert_coord[d];
		cube_coord[d] = int(c/scalar_grid.Spacing(d));
	}

	// Cube index of the actual cube in which the vertex is present.
	VERTEX_INDEX cube_index = scalar_grid.ComputeVertexIndex(cube_coord);
	//DEBUG
	//cout <<"[ "<<isovert.gcube_list[gcube_index].isovert_coord[0] <<" "
	//	<< isovert.gcube_list[gcube_index].isovert_coord[1]<<" "
	//	<< isovert.gcube_list[gcube_index].isovert_coord[2]<<"]\n";
	COORD_TYPE coord1[DIM3];
	//scalar_grid.ComputeScaledCoord(cube_index, coord1);
	//cout <<"cube_index1 scaled "<< cube_index <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
	//scalar_grid.ComputeCoord(cube_index, coord1);
	//cout <<"cube_index1 unscal "<< cube_index <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;

	if(covered_grid.Scalar(cube_index))
	{
		//cout <<"covered\n";
		return true;
	}
	else
	{
		//cout <<"not covered "<<endl;
	}

	//if the sharp vertex is not in the same cube. 
	if (isovert.sharp_ind_grid.Scalar(cube_index) != ISOVERT::NO_INDEX){
		if(isovert.gcube_list[gcube_index].cube_index != cube_index){
			//DEBUG
			/*VERTEX_INDEX cube_index2 = isovert.gcube_list[gcube_index].cube_index;
			scalar_grid.ComputeScaledCoord(cube_index2, coord1);
			cout <<"cube_index original scaled "<<cube_index2 <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;

			scalar_grid.ComputeCoord(cube_index2, coord1);
			cout <<"cube_index original unscal "<<cube_index2 <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;*/

			if (isovert.isFlag(cube_index, COVERED_A_GCUBE)){
				//debug
				//cout <<"return true;"<<endl;
				return true;
			}
			else{
				//debug
				//cout << "Flag "<<isovert.gcube_list[isovert.sharp_ind_grid.Scalar(cube_index2)].flag ;
				//cout <<" return false;"<<endl;
				return false;
			}
		}
	}

	return false;
}
/**
* Select the 3x3x3 regions
*/

void select_3x3x3_regions
	(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	vector<NUM_TYPE> sortd_ind2gcube_list,
	ISOVERT &isovert)
{
	const int dimension = scalar_grid.Dimension();
	const int bin_width = isovert_param.bin_width;

	const COORD_TYPE linf_dist_threshold = 
		isovert_param.linf_dist_thresh_merge_sharp;

	BIN_GRID<VERTEX_INDEX> bin_grid;
	init_bin_grid(scalar_grid, bin_width, bin_grid);

	// list of selected vertices
	vector<VERTEX_INDEX> selected_list;

	SHARPISO_GRID_NEIGHBORS gridn;
	gridn.SetSize(scalar_grid);

	//Keeps track of all cubes which are covered.
	SHARPISO_BOOL_GRID covered_grid;
	covered_grid.SetSize(scalar_grid);
	covered_grid.SetAll(false);

	for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

		GRID_CUBE c;
		c = isovert.gcube_list[sortd_ind2gcube_list[ind]];
		// check boundary
		if(c.boundary_bits == 0)
			if (isovert.isFlag(cube_ind_frm_gc_ind(isovert, sortd_ind2gcube_list[ind]), AVAILABLE_GCUBE)
				&& c.linf_dist < linf_dist_threshold) 
			{

				VERTEX_INDEX v1, v2;
				//Check if the sharp vertex is a covered point.
				//Covered point: the point is inside a covered cube. 

				bool flag_covered_point 
					= check_covered_point(scalar_grid, covered_grid, isovert, sortd_ind2gcube_list[ind]);
				if(flag_covered_point)
				{   
					isovert.gcube_list[sortd_ind2gcube_list[ind]].flag = COVERED_POINT;
					//debug
					/*COORD_TYPE coord1[DIM3];
					scalar_grid.ComputeCoord(isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index, coord1);
					cout <<"covered point cube "<< isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index
					<<" coord "<< coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<" , ";
					scalar_grid.ComputeScaledCoord(isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index, coord1);
					cout <<" scaledcoord "<< coord1[0] + scalar_grid.Spacing(0)/2<<" "<<coord1[1]
					+scalar_grid.Spacing(1)/2<<" "<<coord1[2] + scalar_grid.Spacing(2)/2<<". \n ";*/

					continue;
				}


				bool flag_check_angle = isovert_param.flag_check_triangle_angle;
				bool triangle_flag =
					creates_triangle(scalar_grid, flag_check_angle, isovert, c.cube_index,
					isovalue, bin_grid, bin_width, v1, v2);

				if (!triangle_flag) {
					//The vertex is selcted as sharp.
					isovert.gcube_list[sortd_ind2gcube_list[ind]].flag =
						SELECTED_GCUBE;

					//DEBUG 4
					COORD_TYPE coord1[DIM3];
					scalar_grid.ComputeCoord(
						isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index, coord1);

					//cout << "*** "<<isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[0] <<" "
					//  << isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[1]<<" "
					//  <<isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[2]<<" *** "
					//  <<" linf "<<isovert.gcube_list[sortd_ind2gcube_list[ind]].linf_dist<<" ind "<< ind<<
					//  " cubeindex  "<< isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index<<endl;

					// add to selected list
					VERTEX_INDEX cube_index =
						isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index;

					selected_list.push_back(cube_index);
					covered_grid.Set(cube_index, true);

					bin_grid_insert(scalar_grid, bin_width, cube_index, bin_grid);

					// mark all the neighbors as covered
					for (int i=0;i < gridn.NumVertexNeighborsC(); i++)
					{
						VERTEX_INDEX n = gridn.VertexNeighborC
							(isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index, i);

						//DEBUG
						covered_grid.Set(n, true);

						if(isovert.sharp_ind_grid.Scalar(n)!=ISOVERT::NO_INDEX)
						{
							VERTEX_INDEX neighbor_index_2_gclist = isovert.sharp_ind_grid.Scalar(n);
							isovert.gcube_list[neighbor_index_2_gclist].flag = COVERED_A_GCUBE;
							//DEBUG
							/*VERTEX_INDEX ci = isovert.gcube_list[neighbor_index_2_gclist].cube_index;
							COORD_TYPE coord1[DIM3];
							scalar_grid.ComputeScaledCoord(ci, coord1);
							cout <<"neighbor"<< ci << " scaled "<< ci <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];
							scalar_grid.ComputeCoord(ci, coord1);
							cout <<"neighbor unscal "<< ci <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
							cout <<" flag " << isovert.gcube_list[neighbor_index_2_gclist].flag << endl;*/
						}
						//DEBUG
						else
						{

							/* COORD_TYPE coord1[DIM3];
							scalar_grid.ComputeScaledCoord(n, coord1);
							cout <<"neighbor no index "<< n << " scaled "<< n <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];
							scalar_grid.ComputeCoord(n, coord1);
							cout <<"neighbor unscal "<< n <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;*/

						}
					}
				}
				else
				{
					isovert.gcube_list[sortd_ind2gcube_list[ind]].flag = UNAVAILABLE_GCUBE;
				}
			}
	}
/*
  const int dimension = scalar_grid.Dimension();
  const int bin_width = isovert_param.bin_width;

  const COORD_TYPE linf_dist_threshold = 
    isovert_param.linf_dist_thresh_merge_sharp;

  BIN_GRID<VERTEX_INDEX> bin_grid;
  init_bin_grid(scalar_grid, bin_width, bin_grid);

  // list of selected vertices
  vector<VERTEX_INDEX> selected_list;

  SHARPISO_GRID_NEIGHBORS gridn;
  gridn.SetSize(scalar_grid);

  //Keeps track of all cubes which are covered.
  SHARPISO_BOOL_GRID covered_grid;
  covered_grid.SetSize(scalar_grid);
  covered_grid.SetAll(false);

  for (int ind=0; ind < sortd_ind2gcube_list.size(); ind++) {

    GRID_CUBE c;
    c = isovert.gcube_list[sortd_ind2gcube_list[ind]];
    // check boundary
    if(c.boundary_bits == 0)
      if (isovert.isFlag
          (cube_ind_frm_gc_ind(isovert, sortd_ind2gcube_list[ind]), 
           AVAILABLE_GCUBE) && c.linf_dist < linf_dist_threshold) {

        VERTEX_INDEX v1, v2;
		//Check if the sharp vertex is a covered point.
		//Covered point: the point is inside a covered cube. 
		
		bool flag_covered_point 
			= check_covered_point(scalar_grid, covered_grid, isovert, sortd_ind2gcube_list[ind]);
		if(flag_covered_point)
		{   
			isovert.gcube_list[sortd_ind2gcube_list[ind]].flag = COVERED_POINT;
			continue;
		}
		

        bool flag_check_angle = isovert_param.flag_check_triangle_angle;
        bool triangle_flag =
          creates_triangle(scalar_grid, flag_check_angle, isovert, c.cube_index,
                           isovalue, bin_grid, bin_width, v1, v2);

        if (!triangle_flag) {
			//The vertex is selcted as sharp.
          isovert.gcube_list[sortd_ind2gcube_list[ind]].flag =
            SELECTED_GCUBE;
		
		  //DEBUG 4
		  COORD_TYPE coord1[DIM3];
		  scalar_grid.ComputeCoord(
			  isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index,
			  coord1);
		 
		 cout << "*** "<<isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[0] <<" "
			  << isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[1]<<" "
			  <<isovert.gcube_list[sortd_ind2gcube_list[ind]].isovert_coord[2]<<" *** "
			  <<" linf "<<isovert.gcube_list[sortd_ind2gcube_list[ind]].linf_dist<<" ind "<< ind<<
			  " cubeindex  "<< isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index<<endl;

		  // add to selected list
          VERTEX_INDEX cube_index =
            isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index;

          selected_list.push_back(cube_index);
		  covered_grid.Set(cube_index, true);

          bin_grid_insert(scalar_grid, bin_width, cube_index, bin_grid);

          // mark all the neighbors as covered
          for (int i=0;i < gridn.NumVertexNeighborsC(); i++)
            {
              VERTEX_INDEX n = gridn.VertexNeighborC
                (isovert.gcube_list[sortd_ind2gcube_list[ind]].cube_index, i);

			  //DEBUG
			  covered_grid.Set(n, true);

              if(isovert.sharp_ind_grid.Scalar(n)!=ISOVERT::NO_INDEX)
                {
                  VERTEX_INDEX neighbor_index_2_gclist = isovert.sharp_ind_grid.Scalar(n);
				  isovert.gcube_list[neighbor_index_2_gclist].flag = COVERED_A_GCUBE;
				  //DEBUG
				  VERTEX_INDEX ci = isovert.gcube_list[neighbor_index_2_gclist].cube_index;
				  COORD_TYPE coord1[DIM3];
				  scalar_grid.ComputeScaledCoord(ci, coord1);
				  cout <<"neighbor"<< ci << " scaled "<< ci <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];
				  scalar_grid.ComputeCoord(ci, coord1);
				  cout <<"neighbor unscal "<< ci <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
				  cout <<" flag " << isovert.gcube_list[neighbor_index_2_gclist].flag << endl;
                }
			  //DEBUG
			  else
			  {
				 
				  COORD_TYPE coord1[DIM3];
				  scalar_grid.ComputeScaledCoord(n, coord1);
				  cout <<"neighbor no index "<< n << " scaled "<< n <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];
				  scalar_grid.ComputeCoord(n, coord1);
				  cout <<"neighbor unscal "<< n <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
				  
			  }
            }
        }
        else
          {
            isovert.gcube_list[sortd_ind2gcube_list[ind]].flag = UNAVAILABLE_GCUBE;
          }
      }
  }
*/

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
/// Compute the overlap region between two cube indices
bool find_overlap(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX &cube_index1,
	const VERTEX_INDEX &cube_index2,
	COORD_TYPE *rmin,
	COORD_TYPE *rmax)
{
	COORD_TYPE coord1[DIM3], coord2[DIM3];
	scalar_grid.ComputeCoord(cube_index1, coord1);
	scalar_grid.ComputeCoord(cube_index2, coord2);
	for (int d=0;d<DIM3;d++){
		rmin[d]=max(coord1[d]-1, coord2[d]-1);
		rmax[d]=min(coord1[d]+2, coord2[d]+2);
	}

	// track the dimension of the tracked,
	// if the tracked regions has at least 2 dimension
	int dim_of_overlap=0;
	for (int d=0;d<DIM3;d++)
	{
		if(rmin[d] > rmax[d])
		{ return false; }
		if(rmin[d] < rmax[d])
			dim_of_overlap++;
	}
	if (dim_of_overlap>=2)
		return true;
	else
		return false;
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
	}
	else {
		compute_isovert_positions_edgeI
			(scalar_grid, gradient_grid, isovalue, isovert_param, 
			vertex_position_method, isovert, isovert_info);
	}
}

void MERGESHARP::select_sharp_isovert
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovert)
{
	// keep track of the sorted indices
	std::vector<NUM_TYPE> sortd_ind2gcube_list;
	sort_gcube_list(isovert.gcube_list, sortd_ind2gcube_list);
	//DEBUG 
	//for (int i = 0; i < sortd_ind2gcube_list.size(); i++)
	//{
	// //DEBUG
	// cout <<"isovert "<< isovert.gcube_list[sortd_ind2gcube_list[i]].isovert_coord[0]
	// <<" "<< isovert.gcube_list[sortd_ind2gcube_list[i]].isovert_coord[1]
	// <<" "<< isovert.gcube_list[sortd_ind2gcube_list[i]].isovert_coord[2]<<endl;


	// COORD_TYPE coord1[DIM3];
	// VERTEX_INDEX cube_index = isovert.gcube_list[sortd_ind2gcube_list[i]].cube_index;

	// scalar_grid.ComputeScaledCoord(cube_index, coord1);
	// cout <<"cube_index scaled   ["<< cube_index <<"] " 
	//  << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
	// cout <<"cube center coord "<< coord1[0]-scalar_grid.Spacing(0)/2.0<<" "
	//  <<coord1[1]-scalar_grid.Spacing(1)/2.0 <<" "
	//  <<coord1[2]-scalar_grid.Spacing(2)/2.0 <<"\n";
	// scalar_grid.ComputeCoord(cube_index, coord1);
	// cout <<"cube_index unscaled ["<< cube_index <<"] " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl<<endl;

	//}
	select_3x3x3_regions (scalar_grid, isovalue, isovert_param, 
		sortd_ind2gcube_list, isovert);
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
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],
		SCALAR_TYPE & linf_dist)
	{
		SCALAR_TYPE squareDist=0.0;
		COORD_TYPE cc[DIM3];
		scalar_grid.ComputeCubeCenterScaledCoord(iv,cc);
		SCALAR_TYPE temp_d   = 0.0;
		SCALAR_TYPE max_dist = -1.0;

		for (int d=0; d<DIM3; d++) {
			temp_d = abs(isovert_coord[d]-cc[d])/scalar_grid.Spacing(d);
			if (temp_d > max_dist)
			{ max_dist=temp_d; }
		}
		linf_dist = max_dist;
	}

	// Return true if two  cube-indices are connected
	bool are_connected 
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX & cube_index1,
		const VERTEX_INDEX & cube_index2,
		const SCALAR_TYPE isovalue )
	{
		// find the overlap region
		COORD_TYPE rmin[DIM3], rmax[DIM3];

		bool is_overlap =
			find_overlap(scalar_grid, cube_index1, cube_index2, rmin, rmax);

		if (!is_overlap) { return false; }

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

}

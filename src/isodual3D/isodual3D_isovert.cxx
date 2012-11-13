/// \file isodual3D_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.
/// Created on: Oct 29, 2012
/// Author: arindam


#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdio.h>
#include "isodual3D_isovert.h"
#include "isodual3D_position.h"
#include "isodual3D_datastruct.h"
#include "ijkgrid_macros.h"
#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkscalar_grid.txx"
#include "sharpiso_grids.h"



using namespace std;
using namespace SHARPISO;
using namespace ISODUAL3D;
using namespace IJK;



void compute_linf_dist
( const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],SCALAR_TYPE &linf_dist);
/// Decide if two  cube-indices are connected
bool are_connected (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,const VERTEX_INDEX &cube_index1,
		const VERTEX_INDEX &cube_index2, const SCALAR_TYPE isovalue );

/*
 * Traverse the *isovertData.sharp_ind_grid*
 * if index is not ISOVERT::NO_INDEX then compute the sharp vertex for the cube
 * set grid_cube_flag to be avaliable
 */
void compute_isovert_positions (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT &isovertData)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
			isovert_param.grad_selection_cube_offset;
	OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

	SVD_INFO svd_info;
	IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX)
	{
		NUM_TYPE index = isovertData.sharp_ind_grid.Scalar(iv);
		if (index!=ISOVERT::NO_INDEX)
		{
			// this is an active cube
			// compute the sharp vertex for this cube
			EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
			VERTEX_INDEX num_large_eigenvalues;

			svd_compute_sharp_vertex_for_cube
			(scalar_grid, gradient_grid, iv, isovalue, isovert_param, cube_111,
					isovertData.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);

			//set num eigen
			isovertData.gcube_list[index].num_eigen =  (unsigned char)num_large_eigenvalues;
			//set the sharp vertex type to be *AVAILABLE*
			if(num_large_eigenvalues > 1)
				isovertData.gcube_list[index].flag=AVAILABLE_GCUBE;
			else
				isovertData.gcube_list[index].flag=SMOOTH_GCUBE;

			// store distance
			compute_linf_dist( scalar_grid, iv,
					isovertData.gcube_list[index].isovert_coord,
					isovertData.gcube_list[index].linf_dist);
			// store boundary bits
			scalar_grid.ComputeBoundaryCubeBits
			(iv,isovertData.gcube_list[index].boundary_bits);
		}
	}
}


class gcube_compare{
public:
	std::vector<GRID_CUBE> gcube_list;
	gcube_compare(vector<GRID_CUBE> & gcube_list_ ){gcube_list= gcube_list_;};
	bool operator () (int i,int j)
	{
		return gcube_list[i].linf_dist < gcube_list[j].linf_dist;
	}

};
void sort_gcube_list(vector<NUM_TYPE> &sortd_ind2gcube_list, vector<GRID_CUBE> &gcube_list)
{
	//set up the *sortd_ind2gcube_list*
	sortd_ind2gcube_list.resize(gcube_list.size(), 0);

	for (int i=0;i<sortd_ind2gcube_list.size();i++)
	{
		sortd_ind2gcube_list[i]=i;
	}

	sort (sortd_ind2gcube_list.begin(),sortd_ind2gcube_list.end(), gcube_compare(gcube_list));

	//debug
	/*
	int num3=0;
	cout <<gcube_list.size()<<" "<<sortd_ind2gcube_list.size()<<endl;

	for (int i=0;i<sortd_ind2gcube_list.size();i++)
	{
		if((int)gcube_list[sortd_ind2gcube_list[i]].num_eigen==3)
		cout <<"["<<(int)gcube_list[sortd_ind2gcube_list[i]].num_eigen<<"] "<<
				gcube_list[sortd_ind2gcube_list[i]].linf_dist<<" ["<<
				gcube_list[sortd_ind2gcube_list[i]].isovert_coord[0]<<" "<<
				gcube_list[sortd_ind2gcube_list[i]].isovert_coord[1]<<" "<<
				gcube_list[sortd_ind2gcube_list[i]].isovert_coord[2]<<"] "
				<<endl;

		if((int)gcube_list[sortd_ind2gcube_list[i]].num_eigen==3)
			num3++;
	}
	cout <<"num3 "<<num3<<endl;
	 */
}

/// Compute the cube index from the gc index
VERTEX_INDEX cube_ind_frm_gc_ind
(const ISOVERT &isoData,const  NUM_TYPE &gc_ind){
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
/// Check if selecting this vertex creates a triangle
/// *selected_list* is a list of already selected vertices.
bool creates_triangle (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv,
		const SCALAR_TYPE isovalue,
		vector<VERTEX_INDEX> &selected_list
){
	vector <VERTEX_INDEX> connected_list;
	// get the list of vertices connected to the vertex iv
	get_connected(scalar_grid, isovalue,
			iv, selected_list, connected_list);
	int limit = connected_list.size();
	// for each pair jv1 jv2 in the connected list
	for (int i=0; i< limit-1;++i){
		for(int j=i+1;j<= (limit-1);++j)
		{
			if (are_connected(scalar_grid, connected_list[i],
					connected_list[j], isovalue)){
				return true;
			}
		}
	}
	return false;
}


/*
 * select the 3x3 regions
 */

void select_3x3_regions
(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		vector<NUM_TYPE> sortd_ind2gcube_list,
		ISOVERT &isovertData)

{
	const int dimension = scalar_grid.Dimension();
	int num3=0;
	// list of selected vertices
	vector<VERTEX_INDEX> selected_list;
	SHARPISO_GRID_NEIGHBORS gridn;
	gridn.SetSize(scalar_grid);
	for (int neigen=3;neigen>0;neigen--)
		for (int ind=0;ind<sortd_ind2gcube_list.size();ind++)
		{
			GRID_CUBE c;
			c = isovertData.gcube_list[sortd_ind2gcube_list[ind]];
			// check boundary
			if(c.boundary_bits == 0)
				if (isovertData.isFlag
						(cube_ind_frm_gc_ind(isovertData, sortd_ind2gcube_list[ind]),AVAILABLE_GCUBE)
						&& c.linf_dist < isovertData.linf_dist_threshold
						&& (int)c.num_eigen == neigen)
				{
					/*
					// DEBUG
					cout <<"cube index  "<< c.cube_index;
					COORD_TYPE test[DIM3];
					scalar_grid.ComputeCoord(c.cube_index, test);
					cout <<" ["<<test[0]<<" "<<test[1]<<" "<<test[2]<<"] ";
					cout <<" ne ["<<(int)c.num_eigen<<"] ";
					cout <<" ["<<c.isovert_coord[0]<<" "<<c.isovert_coord[1]<<" "<<c.isovert_coord[2]<<"]"<<endl;
					 */

					if (creates_triangle(scalar_grid, c.cube_index, isovalue, selected_list) == false)
					{
						isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag= SELECTED_GCUBE;
						// add to selected list
						selected_list.push_back(isovertData.gcube_list[sortd_ind2gcube_list[ind]].cube_index);

						for (int i=0;i<gridn.NumVertexNeighborsC();i++)
						{
							VERTEX_INDEX n = gridn.VertexNeighborC
									(isovertData.gcube_list[sortd_ind2gcube_list[ind]].cube_index,i);

							if(isovertData.sharp_ind_grid.Scalar(n)!=ISOVERT::NO_INDEX)
							{
								VERTEX_INDEX neighbor_index_2_gclist
								= isovertData.sharp_ind_grid.Scalar(n);
								//if covered and not boundary
								if (isovertData.gcube_list[neighbor_index_2_gclist].flag == COVERED_GCUBE &&
										(isovertData.gcube_list[neighbor_index_2_gclist].boundary_bits == 0))
								{
									for(int j=0;j<gridn.NumVertexNeighborsC();j++)
									{
										VERTEX_INDEX k=
												gridn.VertexNeighborC(isovertData.gcube_list[neighbor_index_2_gclist].cube_index,j);
										if(isovertData.sharp_ind_grid.Scalar(k)!=ISOVERT::NO_INDEX){
											if(isovertData.gcube_list[isovertData.sharp_ind_grid.Scalar(k)].flag == AVAILABLE_GCUBE)
												isovertData.gcube_list[isovertData.sharp_ind_grid.Scalar(k)].flag = UNAVAILABLE_GCUBE;
										}
									}

								}
								else
								{
									isovertData.gcube_list[neighbor_index_2_gclist].flag = COVERED_GCUBE;
								}
							}
						}
					}
					else
					{
						// what to do if creates triangle
						isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag= UNAVAILABLE_GCUBE;
					}
				}
		}
}



/*
 * Traverse the scalar grid,
 * set everything to be ISOVERT::NO_INDEX
 * if isosurface intersects the cube then set *index*
 * push the corresponding *grid_cube* into the *gcube_list*
 */
void create_active_cubes (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		ISOVERT &isovertData)
{
	NUM_TYPE index=0;
	//set the size of sharp index grid
	isovertData.sharp_ind_grid.SetSize(scalar_grid);
	IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX)
	{
		if (is_gt_cube_min_le_cube_max(scalar_grid, iv, isovalue))
		{
			index=isovertData.gcube_list.size();
			isovertData.sharp_ind_grid.Set(iv,index);
			GRID_CUBE gc;
			gc.cube_index = iv;
			isovertData.gcube_list.push_back(gc);
		}
		else
		{
			// if it is not intersected then mark as ISOVERT::NO_INDEX
			isovertData.sharp_ind_grid.Set(iv,ISOVERT::NO_INDEX);
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

	for (int d=0;d<DIM3;d++)
	{
		if(rmin[d] > rmax[d])
			{ return false; }
	}
	return true;
}
/// process edge called from are connected
void process_edge(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX v0,
		const VERTEX_INDEX v1,
		const SCALAR_TYPE isovalue,
		bool &is_intersect)
{
	COORD_TYPE coord1[DIM3];
	scalar_grid.ComputeCoord(v0,coord1);
	//cout <<"["<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<"]";
	scalar_grid.ComputeCoord(v1,coord1);
	//cout <<"["<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<"]";
	is_intersect = is_gt_min_le_max(scalar_grid, v0, v1, isovalue);
	//cout <<"  -- "<<(int)is_intersect<<" "<<endl;
}

/// Decide if two  cube-indices are connected
bool are_connected (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX &cube_index1,
		const VERTEX_INDEX &cube_index2,
		const SCALAR_TYPE isovalue ){

	// find the overlap region
	COORD_TYPE rmin[DIM3], rmax[DIM3];
	//debug

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
		for (int i1=rmin[d1]; i1 <=rmax[d1]; i1++){
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
	/*	DEBUG
		if(num >0){
			cout <<" rmin "<<rmin[0]<<" "<<rmin[1]<<" "<<rmin[2]<<endl;
			cout <<" rmax "<<rmax[0]<<" "<<rmax[1]<<" "<<rmax[2]<<endl;
			cout <<" num of edges "<<num<<endl;
		}
	 */
	return false;
}

void ISODUAL3D::compute_dual_isovert(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT &isovertData)
{
	create_active_cubes(scalar_grid, gradient_grid, isovalue, isovertData);

	compute_isovert_positions (scalar_grid, gradient_grid, isovalue, isovert_param,
			isovertData);
	// keep track of the sorted indices
	std::vector<NUM_TYPE> sortd_ind2gcube_list;
	sort_gcube_list(sortd_ind2gcube_list,isovertData.gcube_list);
	select_3x3_regions (scalar_grid, gradient_grid, isovalue,
			isovert_param, sortd_ind2gcube_list, isovertData);




	/*  DEBUG
	VERTEX_INDEX v1,v2;
//	COORD_TYPE coord1[DIM3]={16.0,20.0,14.0};
//	COORD_TYPE coord2[DIM3]={19.0,18.0,13.0};

	COORD_TYPE coord1[DIM3]={16.0,20.0,14.0};
	COORD_TYPE coord2[DIM3]={17.0,19.0,12.0};

	v1=scalar_grid.ComputeVertexIndex(coord1);
	v2 =scalar_grid.ComputeVertexIndex(coord2);
	cout <<"---- v1 "<<v1<<"v2 "<<v2<<endl;
	are_connected(scalar_grid,v1,v2,isovalue);
	 */
}


/*
 * Find the linf distance between the sharp vertex for the cube and the center of the cube
 */
void compute_linf_dist( const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],
		SCALAR_TYPE &linf_dist)
{
	SCALAR_TYPE squareDist=0.0;
	COORD_TYPE cc[DIM3];
	scalar_grid.ComputeCubeCenterCoord(iv,cc);
	SCALAR_TYPE temp_d   = 0.0;
	SCALAR_TYPE max_dist = -1.0;
	for (int d=0;d<DIM3;d++){
		temp_d = abs(isovert_coord[d]-cc[d]);
		if (temp_d > max_dist)
			max_dist=temp_d;
	}
	linf_dist = max_dist;
};

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
// Set ISOVERT_INFO
// **************************************************

/// Count number of vertices on sharp corners or sharp edges.
/// Count number of smooth vertices.
void ISODUAL3D::count_vertices
(const ISOVERT & isovert, ISOVERT_INFO & isovert_info)
{
	isovert_info.num_sharp_corners = 0;
	isovert_info.num_sharp_edges = 0;
	isovert_info.num_smooth_vertices = 0;
	for (VERTEX_INDEX i = 0; i < isovert.gcube_list.size(); i++) {
		if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
			if (isovert.gcube_list[i].num_eigen == 2)
			{ isovert_info.num_sharp_edges++; }
			else if (isovert.gcube_list[i].num_eigen == 3) {
				{ isovert_info.num_sharp_corners++; }
			}
		}
		else if (isovert.gcube_list[i].flag == SMOOTH_GCUBE) {
			isovert_info.num_smooth_vertices++;
		}
	}
}

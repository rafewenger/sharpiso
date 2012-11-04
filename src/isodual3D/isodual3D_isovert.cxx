/// \file isodual3D_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.
/// Created on: Oct 29, 2012
/// Author: arindam


#include <iostream>
#include <cmath>
#include <algorithm>
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



void compute_l2_dist
( const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],SCALAR_TYPE &l2dist);

/*
 * Traverse the *isovertData.sharp_ind_grid*
 * if index is not -1 then compute the sharp vertex for the cube
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
		if (index!=-1)
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
			isovertData.gcube_list[index].num_eigen = num_large_eigenvalues;
			//set the sharp vertex type to be *AVAILABLE*
			if(num_large_eigenvalues > 1)
				isovertData.gcube_list[index].flag=AVAILABLE_GCUBE;
			else
				isovertData.gcube_list[index].flag=SMOOTH_GCUBE;

			// STORE DISTANCE
			compute_l2_dist( scalar_grid, iv,
					isovertData.gcube_list[index].isovert_coord,
					isovertData.gcube_list[index].l2dist);
		}
	}
}


class gcube_compare{
public:
	std::vector<GRID_CUBE> gcube_list;
	gcube_compare(vector<GRID_CUBE> & gcube_list_ ){gcube_list= gcube_list_;};
	bool operator () (int i,int j)
	{
		return gcube_list[i].l2dist < gcube_list[j].l2dist;
	}

};
void sort_gcube_list(vector<NUM_TYPE> &sortd_ind2gcube_list, vector<GRID_CUBE> &gcube_list)
{
	//set up the *sortd_ind2gcube_list*
	sortd_ind2gcube_list.resize(gcube_list.size(), 0);
	vector<NUM_TYPE>::iterator it;
	for (int i=0;i<sortd_ind2gcube_list.size();i++)
	{
		sortd_ind2gcube_list[i]=i;
	}

	sort (sortd_ind2gcube_list.begin(),sortd_ind2gcube_list.end(), gcube_compare(gcube_list));
}


/*
 * Helper functions to select_3x3_regions
 */

bool is_cube(const GRID_CUBE_FLAG flag,
		const VERTEX_INDEX iv,
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		ISOVERT &isovertData)
{
	//check if intersected
	if(isovertData.sharp_ind_grid.Scalar(iv)!=-1)
	{
		NUM_TYPE index = isovertData.sharp_ind_grid.Scalar(iv);
		if (isovertData.gcube_list[index].flag == flag)
			return true;
	}
	return false;
}

/*
 * select the 3x3 regions
 */
/// does not check for boundary regions

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

	SHARPISO_GRID_NEIGHBORS gridn;
	gridn.SetSize(scalar_grid);

	for (int ind=0;ind<sortd_ind2gcube_list.size();ind++)
	{
		GRID_CUBE c;
		c = isovertData.gcube_list[sortd_ind2gcube_list[ind]];

		if ( c.flag == AVAILABLE_GCUBE && c.l2dist < 0.8)
		{
			isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag= SELECTED_GCUBE;
			for (int i=0;i<gridn.NumVertexNeighborsC();i++)
			{
				VERTEX_INDEX n = gridn.VertexNeighborC
						(isovertData.gcube_list[sortd_ind2gcube_list[ind]].index2sg,i);
				if(isovertData.sharp_ind_grid.Scalar(n)!=-1)
				{
					VERTEX_INDEX neighbor_index_2_gclist
					= isovertData.sharp_ind_grid.Scalar(n);
					if (isovertData.gcube_list[neighbor_index_2_gclist].flag == COVERED_GCUBE)
					{
						// set unavailable : not implemented yet
					}
					else
					{
						isovertData.gcube_list[neighbor_index_2_gclist].flag = COVERED_GCUBE;
					}
				}
			}
		}
	}

}



/*
 * Traverse the scalar grid,
 * set everything to be -1
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
			gc.index2sg = iv;
			isovertData.gcube_list.push_back(gc);
		}
		else
		{
			// if it is not intersected then mark as -1
			isovertData.sharp_ind_grid.Set(iv,-1);
		}
	}
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


}


/*
 * Find the distance between the sharp vertex for the cube and the center of the cube
 */
void compute_l2_dist( const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],
		SCALAR_TYPE &l2dist)
{
	SCALAR_TYPE squareDist=0.0;
	COORD_TYPE cc[DIM3];
	scalar_grid.ComputeCubeCenterCoord(iv,cc);
	for (int d=0;d<DIM3;d++	)
		squareDist = squareDist + (isovert_coord[d]-cc[d])*(isovert_coord[d]-cc[d]);
	l2dist=sqrt(squareDist);
};

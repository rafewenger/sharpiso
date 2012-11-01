/*
 * isodual3D_isovert.cxx
 *
 *  Created on: Oct 29, 2012
 *      Author: arindam
 */

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



using namespace std;
using namespace SHARPISO;
using namespace ISODUAL3D;
using namespace IJK;



void compute_l2_dist
( const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],SCALAR_TYPE &l2dist);

/*
 * Traverse the *isovertData.sharp_ind_grid*
 * if index is not -1 then compute the sharp vertex for the cube
 * set grid_cube_flag to be avaliable
 */
void compute_isovert_positions (
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const ISODUAL_PARAM & isodual_param,
		ISOVERT &isovertData)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
			isodual_param.grad_selection_cube_offset;
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
			(scalar_grid, gradient_grid, iv, isovalue, isodual_param, cube_111,
					isovertData.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);
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
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
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
void select_3x3_regions
(
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const ISODUAL_PARAM & isodual_param,
		ISOVERT &isovertData)

{
	const int dimension = scalar_grid.Dimension();
	for  (VERTEX_INDEX iv=0;iv<scalar_grid.NumVertices();iv++)
	{
		/*
		//check if the vertex is being intersected
		if(isovertData.sharp_ind_grid.Scalar(iv)!=-1)
		{
			NUM_TYPE index = isovertData.sharp_ind_grid.Scalar(iv);
			if( ! is_in_boundary (scalar_grid,iv))
			{
				//not in boudary
				if(isovertData.gcube_list[index].flag == AVAILABLE_GCUBE)
				{
					for (int d = 0; d < dimension; d++) {
						VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv, d);
						if(is_cube(AVAILABLE_GCUBE, iv0, scalar_grid, isovertData))
						{

						}
						VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv, d);

					}
				}
			}
		}

		 */
	}
}



/*
 * Traverse the scalar grid,
 * set everything to be -1
 * if isosurface intersects the cube then set *index*
 * push the corresponding *grid_cube* into the *gcube_list*
 */
void create_active_cubes (
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
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
			isovertData.gcube_list.push_back(gc);
		}
		else
		{
			// if it is not intersected then mark as -1
			isovertData.sharp_ind_grid.Set(iv,-1);
		}
	}
}


void ISODUAL3D::computeDualIsovert(
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const ISODUAL_PARAM & isodual_param,
		ISOVERT &isovertData)
{
	create_active_cubes(scalar_grid, gradient_grid, isovalue, isovertData);

	compute_isovert_positions (scalar_grid, gradient_grid, isovalue, isodual_param,
			isovertData);

	std::vector<NUM_TYPE> sortd_ind2gcube_list;

	sort_gcube_list(sortd_ind2gcube_list,isovertData.gcube_list);

	select_3x3_regions (scalar_grid, gradient_grid, isovalue,
			isodual_param, isovertData);
}


/*
 * Find the distance between the sharp vertex for the cube and the center of the cube
 */
void compute_l2_dist( const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
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

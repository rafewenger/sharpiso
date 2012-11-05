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



void compute_linf_dist
( const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],SCALAR_TYPE &linf_dist);

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
			scalar_grid.ComputeBoundaryBits
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
    /*
		if( gcube_list[i].num_eigen == 3 && gcube_list[j].num_eigen < 3)
		{
		  cout <<" "<<(int)gcube_list[i].num_eigen<<endl;
		  return true;
		}
		else
		{
		  return gcube_list[i].linf_dist < gcube_list[j].linf_dist;
		}
		*/
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

	SHARPISO_GRID_NEIGHBORS gridn;
	gridn.SetSize(scalar_grid);
  for (int neigen=3;neigen>0;neigen--)
	for (int ind=0;ind<sortd_ind2gcube_list.size();ind++)
	{
		GRID_CUBE c;
		c = isovertData.gcube_list[sortd_ind2gcube_list[ind]];

    // *** USE PARAMETER, NOT 0.8.  (PARAMETER SHOULD BE 1.4.)
		// check boundary
		if(c.boundary_bits == 0)
		if ( c.flag == AVAILABLE_GCUBE && c.linf_dist < 1.4 && c.num_eigen == neigen)
		{
			isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag= SELECTED_GCUBE;
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
						  VERTEX_INDEX k=gridn.VertexNeighborC(isovertData.gcube_list[neighbor_index_2_gclist].cube_index,j);
						  if(isovertData.sharp_ind_grid.Scalar(k)!=ISOVERT::NO_INDEX){
						    if(isovertData.gcube_list[isovertData.sharp_ind_grid.Scalar(k)].flag == AVAILABLE_GCUBE)
						    isovertData.gcube_list[isovertData.sharp_ind_grid.Scalar(k)].flag == UNAVAILABLE_GCUBE;
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

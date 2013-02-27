/// \file isodual3D_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.
/// Created on: Oct 29, 2012
/// Author: arindam


#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdio.h>

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkgrid_macros.h"
#include "ijkscalar_grid.txx"

#include "isodual3D_isovert.h"
#include "isodual3D_position.h"


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
bool is_angle_large(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const ISOVERT &isovertData, const VERTEX_INDEX iv,
		const SCALAR_TYPE threshold, const VERTEX_INDEX v1,
		const VERTEX_INDEX v2);

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
    isovertData.gcube_list[index].flag_centroid_location = false;

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

			// set number of eigenvalues.
			isovertData.gcube_list[index].num_eigenvalues =
        (unsigned char) num_large_eigenvalues;
			// set the sharp vertex type to be *AVAILABLE*.
			if(num_large_eigenvalues > 1 && svd_info.location == LOC_SVD)
        { isovertData.gcube_list[index].flag = AVAILABLE_GCUBE; }
			else 
        { isovertData.gcube_list[index].flag = SMOOTH_GCUBE; }

			// store distance.
			compute_linf_dist( scalar_grid, iv,
					isovertData.gcube_list[index].isovert_coord,
					isovertData.gcube_list[index].linf_dist);
			// store boundary bits
			scalar_grid.ComputeBoundaryCubeBits
			(iv,isovertData.gcube_list[index].boundary_bits);
		}
	}
}

/// Recompute isosurface vertex positions for cubes 
///   which are not selected or covered.
void ISODUAL3D::recompute_isovert_positions (
    const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT & isovert)
{
	for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
		GRID_CUBE_FLAG cube_flag = isovert.gcube_list[i].flag;

		if ((cube_flag == AVAILABLE_GCUBE) || (cube_flag == UNAVAILABLE_GCUBE)) {

			VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

			compute_edgeI_centroid
        (scalar_grid, gradient_grid, isovalue, cube_index,
         isovert_param.use_sharp_edgeI, isovert.gcube_list[i].isovert_coord);

      isovert.gcube_list[i].flag_centroid_location = true;
		}
	}

}


/// Recompute isosurface vertex positions for cubes 
///   which are not selected or covered.
void ISODUAL3D::recompute_isovert_positions (
    const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
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

void ISODUAL3D::set_edge_index
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


/*
 * Traverse the *isovertData.sharp_ind_grid*
 * if index is not ISOVERT::NO_INDEX then compute the sharp vertex for the cube
 * set grid_cube_flag to be avaliable
 */
void compute_isovert_positions (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const std::vector<COORD_TYPE> & edgeI_coord,
		const std::vector<COORD_TYPE> & edgeI_normal_coord,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT & isovertData)
{
	const SIGNED_COORD_TYPE grad_selection_cube_offset =
			isovert_param.grad_selection_cube_offset;
	OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);
	SHARPISO_EDGE_INDEX_GRID edge_index
	(DIM3, scalar_grid.AxisSize(), DIM3);
	SVD_INFO svd_info;

	set_edge_index(edgeI_coord, edge_index);

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
			(scalar_grid, edgeI_coord, edgeI_normal_coord, edge_index,
					iv, isovalue, isovert_param, cube_111,
					isovertData.gcube_list[index].isovert_coord,
					eigenvalues, num_large_eigenvalues, svd_info);

			//set num eigen
			isovertData.gcube_list[index].num_eigenvalues =  
        (unsigned char)num_large_eigenvalues;
			//set the sharp vertex type to be *AVAILABLE*
			if(num_large_eigenvalues > 1)
				{
					isovertData.gcube_list[index].flag=AVAILABLE_GCUBE;
				}
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


void ISODUAL3D::sort_gcube_list
(const std::vector<GRID_CUBE> & gcube_list,
 std::vector<NUM_TYPE> & sortd_ind2gcube_list)
{
	GCUBE_COMPARE gcube_compare(gcube_list);

	for (int i=0;i<gcube_list.size();i++)
	{
		if (gcube_list[i].num_eigenvalues > 1)
			sortd_ind2gcube_list.push_back(i);
	}

	sort (sortd_ind2gcube_list.begin(),sortd_ind2gcube_list.end(), 
        gcube_compare);
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
/// @param check_triangl_angle If true, check it triangle has large angles.
// @param bin_grid Contains the already selected vertices.
// @param[out] v1,v2 vertex indices which form a triangle with iv
bool ISODUAL3D::creates_triangle (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const bool check_triangle_angle,
		const ISOVERT & isovertData,
		const VERTEX_INDEX iv,
		const SCALAR_TYPE isovalue,
		const BIN_GRID<VERTEX_INDEX> & bin_grid,
		const AXIS_SIZE_TYPE bin_width,
		VERTEX_INDEX & v1,
		VERTEX_INDEX & v2
)
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
	for (int i=0; i< limit-1;++i){
		for(int j=i+1;j<= (limit-1);++j)
		{
			if (are_connected(scalar_grid, connected_list[i],
					connected_list[j], isovalue)){
				v1 = connected_list[i];
				v2 = connected_list[j];
				if (check_triangle_angle){
					// checking angle
					bool flag_large_angle = is_angle_large(scalar_grid, isovertData, iv, threshold, v1 , v2);

					if (flag_large_angle)
					{
						return true;
					}
				}
				else{
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
void ISODUAL3D::init_bin_grid
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
		BIN_GRID<VERTEX_INDEX> & bin_grid)
{
	const int dimension = grid.Dimension();
	AXIS_SIZE_TYPE axis_size[dimension];

	for (int d = 0; d < dimension; d++) {
		axis_size[d] = IJK::compute_subsample_size(grid.AxisSize(d), bin_width);
	}

	bin_grid.SetSize(dimension, axis_size);
}

void ISODUAL3D::bin_grid_insert
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid)
{
	static GRID_COORD_TYPE coord[DIM3];

	grid.ComputeCoord(cube_index, coord);
	divide_coord_3D(bin_width, coord);
	VERTEX_INDEX ibin = bin_grid.ComputeVertexIndex(coord);
	bin_grid.Insert(ibin, cube_index);
}

/*
 * Compute the cosine of the angle formed by the
 * vector v1v2 and v3v2
 */
void compute_cos_angle(
		const ISOVERT &isovertData,
		const VERTEX_INDEX gcube_list_index_v1,
		const VERTEX_INDEX gcube_list_index_v2,
		const VERTEX_INDEX gcube_list_index_v3,
		SCALAR_TYPE & cos_angle
)
{
	COORD_TYPE coord0[DIM3], coord1[DIM3],
	coord2[DIM3],vec12[DIM3], vec32[DIM3];

	subtract_coord(DIM3, isovertData.gcube_list[gcube_list_index_v2].isovert_coord,
			isovertData.gcube_list[gcube_list_index_v1].isovert_coord, vec12);
	normalize_vector(DIM3, vec12, 0.001, vec12);

	subtract_coord(DIM3, isovertData.gcube_list[gcube_list_index_v2].isovert_coord,
			isovertData.gcube_list[gcube_list_index_v3].isovert_coord, vec32);
	normalize_vector(DIM3, vec32, 0.001, vec32);

	compute_inner_product(DIM3, vec12, vec32, cos_angle);
}

/*
 * is angle large
 * is the angle between the triangle formed by iv,v1,v2 too large(or cosine angle less)
 */
bool is_angle_large(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const ISOVERT &isovertData,
		const VERTEX_INDEX iv,
		const SCALAR_TYPE threshold,
		const VERTEX_INDEX v1,
		const VERTEX_INDEX v2)
{

	SCALAR_TYPE cos_angle_iv, cos_angle_v1, cos_angle_v2;
	VERTEX_INDEX gcube_list_index_v1;
	VERTEX_INDEX gcube_list_index_v2;
	VERTEX_INDEX gcube_list_index_v3;

	gcube_list_index_v1=isovertData.sharp_ind_grid.Scalar(v1);
	gcube_list_index_v2=isovertData.sharp_ind_grid.Scalar(v2);
	gcube_list_index_v3=isovertData.sharp_ind_grid.Scalar(iv);


	compute_cos_angle(isovertData, gcube_list_index_v1,
			gcube_list_index_v3,  gcube_list_index_v2, cos_angle_iv);

	compute_cos_angle(isovertData, gcube_list_index_v3,
			gcube_list_index_v1,  gcube_list_index_v2, cos_angle_v1);

	compute_cos_angle(isovertData, gcube_list_index_v3,
			gcube_list_index_v2,  gcube_list_index_v1, cos_angle_v2);

	if ((cos_angle_iv < threshold) || (cos_angle_v1 < threshold)
			|| (cos_angle_v2 < threshold))
	{
		return true;
	}
	else
		return false;
}
/*
 * select the 3x3x3 regions
 */

void select_3x3x3_regions
(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		vector<NUM_TYPE> sortd_ind2gcube_list,
		ISOVERT &isovertData)
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
	for (int ind=0;ind<sortd_ind2gcube_list.size();ind++) {

		GRID_CUBE c;
		c = isovertData.gcube_list[sortd_ind2gcube_list[ind]];
		// check boundary
		if(c.boundary_bits == 0)
			if (isovertData.isFlag
					(cube_ind_frm_gc_ind(isovertData, sortd_ind2gcube_list[ind]),AVAILABLE_GCUBE)
					&& c.linf_dist < linf_dist_threshold) {

				VERTEX_INDEX v1, v2;

				bool flag_check_angle = isovert_param.flag_check_triangle_angle;
				bool triangle_flag =
						creates_triangle(scalar_grid, flag_check_angle, isovertData, c.cube_index,
								isovalue, bin_grid, bin_width, v1, v2);

				if (!triangle_flag) {

					isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag =
							SELECTED_GCUBE;

					// add to selected list
					VERTEX_INDEX cube_index =
							isovertData.gcube_list[sortd_ind2gcube_list[ind]].cube_index;

					selected_list.push_back(cube_index);

					bin_grid_insert(scalar_grid, bin_width, cube_index, bin_grid);

					// mark all the neighbors as covered
					for (int i=0;i<gridn.NumVertexNeighborsC();i++)
					{

						VERTEX_INDEX n = gridn.VertexNeighborC
								(isovertData.gcube_list[sortd_ind2gcube_list[ind]].cube_index,i);

						if(isovertData.sharp_ind_grid.Scalar(n)!=ISOVERT::NO_INDEX)
						{

							VERTEX_INDEX neighbor_index_2_gclist = isovertData.sharp_ind_grid.Scalar(n);
							isovertData.gcube_list[neighbor_index_2_gclist].flag = COVERED_GCUBE;
						}
					}
				}
				else
				{
					isovertData.gcube_list[sortd_ind2gcube_list[ind]].flag = UNAVAILABLE_GCUBE;
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
	COORD_TYPE coord1[DIM3];
	scalar_grid.ComputeCoord(v0,coord1);
	scalar_grid.ComputeCoord(v1,coord1);
	is_intersect = is_gt_min_le_max(scalar_grid, v0, v1, isovalue);
}

/// Decide if two  cube-indices are connected
bool are_connected (
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX &cube_index1,
		const VERTEX_INDEX &cube_index2,
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

	return false;
}

void ISODUAL3D::compute_dual_isovert(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT & isovertData)
{
	IJK::PROCEDURE_ERROR error("compute_dual_isovert");

	if (!gradient_grid.Check
			(scalar_grid, "gradient grid", "scalar grid", error))
	{ throw error; }

	create_active_cubes(scalar_grid, isovalue, isovertData);

	compute_isovert_positions 
    (scalar_grid, gradient_grid, isovalue, isovert_param, isovertData);
}

void ISODUAL3D::select_sharp_isovert(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT & isovertData)
{
	// keep track of the sorted indices
	std::vector<NUM_TYPE> sortd_ind2gcube_list;
	sort_gcube_list(isovertData.gcube_list, sortd_ind2gcube_list);
	select_3x3x3_regions (scalar_grid, isovalue, isovert_param, 
                        sortd_ind2gcube_list, isovertData);
}

void ISODUAL3D::compute_dual_isovert(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const std::vector<COORD_TYPE> & edgeI_coord,
		const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
		ISOVERT &isovertData)
{
	create_active_cubes(scalar_grid, isovalue, isovertData);

	compute_isovert_positions 
	(scalar_grid, edgeI_coord, edgeI_normal_coord,
   isovalue, isovert_param, isovertData);
}


/*
 * Find the linf distance between the sharp vertex for the cube and the center of the cube
 */
void compute_linf_dist(
		const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv,
		COORD_TYPE isovert_coord[DIM3],
		SCALAR_TYPE &linf_dist)
{
	SCALAR_TYPE squareDist=0.0;
	COORD_TYPE cc[DIM3];
	scalar_grid.ComputeCubeCenterCoord(iv,cc);
	SCALAR_TYPE temp_d   = 0.0;
	SCALAR_TYPE max_dist = -1.0;
	for (int d=0; d<DIM3; d++){
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

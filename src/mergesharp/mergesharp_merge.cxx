/// \file mergesharp_merge.cxx
/// Merge cubes containing sharp vertices.

/*
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


#include <vector>
#include <algorithm>
#include <set>

#include "ijkcoord.txx"
#include "ijkgraph.txx"
#include "ijkmesh.txx"
#include "ijkgrid_macros.h"

#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"
#include "mergesharp_merge.h"
#include "mergesharp_extract.h"



// forward declarations
namespace {

	void determine_gcube_map
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		const SHARP_ISOVERT_PARAM & sharp_isovert_param,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info);

	void determine_gcube_map
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		const SHARP_ISOVERT_PARAM & sharp_isovert_param,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info);

	void map_gcube_indices(const std::vector<VERTEX_INDEX> & gcube_map,
		std::vector<VERTEX_INDEX> & gcube_index);

	void set_first_gcube_isov
		(const MERGESHARP::ISOVERT & isovert, 
		const std::vector<MERGESHARP::DUAL_ISOVERT> & iso_vlist,
		SHARPISO::NUM_TYPE * first_gcube_isov);

	void insert_vertex_list
		(const std::vector<VERTEX_INDEX> & vlist,
		MERGESHARP::VERTEX_HASH_TABLE & vertex_hash);

	void remap_vertex_list
		(const MERGESHARP::VERTEX_HASH_TABLE & vertex_hash,
		const std::vector<VERTEX_INDEX> & vlist,
		std::vector<VERTEX_INDEX> & new_vlist);
}

// **************************************************
// Merge some isosurface vertices
// **************************************************

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void MERGESHARP::merge_sharp_iso_vertices
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
	const SCALAR_TYPE isovalue,
	ISOVERT & isovert,
	const SHARP_ISOVERT_PARAM & sharp_isovert_param,
	std::vector<VERTEX_INDEX> & isoquad_cube,
	std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info)
{
	const NUM_TYPE num_gcube = isovert.gcube_list.size();

	determine_gcube_map
		(scalar_grid, isovalue, isovert, sharp_isovert_param, 
		gcube_map, sharpiso_info);

	// Count number merged isosurface vertices.
	NUM_TYPE num_merged = 0;
	for (int i = 0; i < num_gcube; i++) {
		if (gcube_map[i] != i) { num_merged++; }
	}
	sharpiso_info.num_merged_iso_vertices = num_merged;

	map_gcube_indices(gcube_map, isoquad_cube);
}

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void MERGESHARP::merge_sharp_iso_vertices
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
	ISOVERT & isovert, 
	const SHARP_ISOVERT_PARAM & sharp_isovert_param,
	std::vector<VERTEX_INDEX> & isoquad_cube,
	SHARPISO_INFO & sharpiso_info)
{
	const NUM_TYPE num_gcube = isovert.gcube_list.size();
	std::vector<VERTEX_INDEX> gcube_map(num_gcube);

	merge_sharp_iso_vertices
		(scalar_grid, isovalue, isovert, sharp_isovert_param,
		isoquad_cube, gcube_map, sharpiso_info);
}

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
// Allows multiple isosurface vertices per cube.
void MERGESHARP::merge_sharp_iso_vertices_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	const SCALAR_TYPE isovalue,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	ISOVERT & isovert,
	const SHARP_ISOVERT_PARAM & sharp_isovert_param,
	std::vector<VERTEX_INDEX> & poly_vert,
	std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info)
{
	const NUM_TYPE num_gcube = isovert.gcube_list.size();
	IJK::ARRAY<NUM_TYPE> first_gcube_isov(num_gcube);

	determine_gcube_map
		(scalar_grid, isodual_table, isovalue, isovert, sharp_isovert_param, 
		gcube_map, sharpiso_info);

	// Count number merged isosurface vertices.
	NUM_TYPE num_merged = 0;
	for (int i = 0; i < num_gcube; i++) {
		if (gcube_map[i] != i) { num_merged++; }
	}
	sharpiso_info.num_merged_iso_vertices = num_merged;

	set_first_gcube_isov(isovert, iso_vlist, first_gcube_isov.Ptr());

	for (NUM_TYPE i = 0; i < poly_vert.size(); i++) {
		NUM_TYPE k = poly_vert[i];
		VERTEX_INDEX cube_index0 = iso_vlist[k].cube_index;
		NUM_TYPE gcube_index0 = isovert.sharp_ind_grid.Scalar(cube_index0);
		VERTEX_INDEX gcube_index1 = gcube_map[gcube_index0];
		if (isovert.gcube_list[gcube_index1].flag == SELECTED_GCUBE) {
			// Reset poly_vert[i] to first isosurface vertex for cube_index1.
			poly_vert[i] = first_gcube_isov[gcube_index1];
		}
	}
}

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
// Allows multiple isosurface vertices per cube.
void MERGESHARP::merge_sharp_iso_vertices_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	const SCALAR_TYPE isovalue,
	const std::vector<DUAL_ISOVERT> & iso_vlist,
	ISOVERT & isovert,
	const SHARP_ISOVERT_PARAM & sharp_isovert_param,
	std::vector<VERTEX_INDEX> & poly_vert,
	SHARPISO_INFO & sharpiso_info)
{
	const NUM_TYPE num_gcube = isovert.gcube_list.size();
	std::vector<VERTEX_INDEX> gcube_map(num_gcube);

	merge_sharp_iso_vertices_multi
		(scalar_grid, isodual_table, isovalue, iso_vlist, isovert, 
		sharp_isovert_param, poly_vert, gcube_map, sharpiso_info);
}

// **************************************************
// Map isosurface vertices
// **************************************************

namespace {

	using namespace MERGESHARP;

	void map_iso_vertex(const std::vector<GRID_CUBE> & gcube_list,
		const INDEX_DIFF_TYPE from_cube,
		const INDEX_DIFF_TYPE to_cube,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		if (from_cube != ISOVERT::NO_INDEX) {
			if (gcube_list[from_cube].flag != SELECTED_GCUBE && 
				gcube_map[from_cube] == from_cube) {
					if (gcube_list[from_cube].boundary_bits == 0) {
						//DEBUG
						using namespace std;
						//cout <<"successfully mapped"<< endl;
						// Map gcube_list[from_cube] to isosurface vertex in cube to_cube.
						gcube_map[from_cube] = to_cube;
					}
					else
					{//DEBUG
						using namespace std;
						//cout <<"\n";
					}
			}
		}
	}

	// Find the list of gcube_index of vertices which are "connected"
	// to a vertex v with gcube_index_v (the index to the gcube_list of vertex v.)
	//[out] connected_sharp, list of gcube_indices
	void find_connected_sharp(
		const VERTEX_INDEX gcube_index_v,// entry into the gcube_list.
		const MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
		std::vector<VERTEX_INDEX> & connected_sharp // list of gcube_indices.
		)
	{
		//DEBUG
		using namespace std;
		SHARPISO_GRID_NEIGHBORS gridn;
		// Set size of grid neighbors grid.
		gridn.SetSize(isovert.sharp_ind_grid);


		if (isovert.gcube_list[gcube_index_v].boundary_bits == 0) {

			VERTEX_INDEX cube_index_v = isovert.gcube_list[gcube_index_v].cube_index;
			//DEBUG
			//cout <<" cube_index_v " << cube_index_v << endl;  
			for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
			{  
				//cube_index of neighbor of v.
				VERTEX_INDEX neighbor_cube_index_v;
				neighbor_cube_index_v = gridn.CubeNeighborF(cube_index_v, j);
				//DEBUG
				//std::cout <<"neighbor cube index " << neighbor_cube_index_v << std::endl;
				INDEX_DIFF_TYPE neighbor_gcube_index_v = isovert.sharp_ind_grid.Scalar(neighbor_cube_index_v);
				//NOT_ACTIVE
				if(neighbor_gcube_index_v == ISOVERT::NO_INDEX)
					continue;
				//DEBUG
				//std::cout <<"gcube_map size "<<
				//	  gcube_map.size() <<" neighbor_gcube_index_v neighbor "<< neighbor_gcube_index_v <<std::endl;
				// the neighbor maps to x;
				VERTEX_INDEX x = gcube_map[neighbor_gcube_index_v]; 

				bool is_x_selected = false;
				if(isovert.gcube_list[x].flag == SELECTED_GCUBE)
				{
					is_x_selected = true;
				}
				if (is_x_selected)
				{
					std::vector<VERTEX_INDEX>::iterator it;
					it = find (connected_sharp.begin(), connected_sharp.end(), x);
					// if x not in connected_sharp, add it to connected_sharp
					if (it == connected_sharp.end())
					{
						connected_sharp.push_back(x);
					}
				}
			}//for_end

			for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
			{
				//cube_index of neighbor of v.
				VERTEX_INDEX neighbor_cube_index_v;
				neighbor_cube_index_v = gridn.CubeNeighborE(cube_index_v, j);
				INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_cube_index_v);
				//NOT_ACTIVE
				if(k == ISOVERT::NO_INDEX)
					continue;
				// the neighbor maps to x;
				VERTEX_INDEX x = gcube_map[k];

				//initlize is_x_selected to false; 
				bool is_x_selected = false;
				if(isovert.gcube_list[x].flag == SELECTED_GCUBE)
				{
					is_x_selected = true;
				}
				if (is_x_selected)
				{
					std::vector<VERTEX_INDEX>::iterator it;
					it = find (connected_sharp.begin(), connected_sharp.end(), x);
					// if x not in connected_sharp, add it to connected_sharp
					if (it == connected_sharp.end())
					{
						connected_sharp.push_back(x);
					}
				}
			}//for_end
			//DEBUG
			//find connected sharp. 
			/*
			//For edge neighbors
			for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
			{
			//cube_index of neighbor of v.
			VERTEX_INDEX neighbor_cube_index_v;
			neighbor_cube_index_v = gridn.CubeNeighborE(cube_index_v, j);
			INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_cube_index_v);
			//NOT_ACTIVE
			if(k == ISOVERT::NO_INDEX)
			continue;
			// the neighbor maps to x;
			VERTEX_INDEX x = gcube_map[k];

			//initlize is_x_selected to false; 
			bool is_x_selected = false;
			if(isovert.gcube_list[x].flag == SELECTED_GCUBE)
			{
			is_x_selected = true;
			}
			if (is_x_selected)
			{
			std::vector<VERTEX_INDEX>::iterator it;
			it = find (connected_sharp.begin(), connected_sharp.end(), x);
			// if x not in connected_sharp, add it to connected_sharp
			if (it == connected_sharp.end())
			{
			connected_sharp.push_back(x);
			}
			}
			}//for_end

			//For vertex neighbors
			for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
			{
			//cube_index of neighbor of v.
			VERTEX_INDEX neighbor_cube_index_v;
			neighbor_cube_index_v = gridn.CubeNeighborV(cube_index_v, j);
			INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_cube_index_v);
			//NOT_ACTIVE
			if(k == ISOVERT::NO_INDEX)
			continue;
			// the neighbor maps to x;
			VERTEX_INDEX x = gcube_map[k];

			//initlize is_x_selected to false; 
			bool is_x_selected = false;

			if(isovert.gcube_list[x].flag == SELECTED_GCUBE)
			{
			is_x_selected = true;
			}
			if (is_x_selected)
			{
			std::vector<VERTEX_INDEX>::iterator it;
			it = find (connected_sharp.begin(), connected_sharp.end(), x);
			// if x not in connected_sharp, add it to connected_sharp
			if (it == connected_sharp.end())
			{
			connected_sharp.push_back(x);
			}
			}
			}//for_end
			*/

		}
	}
	//****************************************
	/// are connected and required functions 
	//****************************************
	/// process edge called from are connected
	void process_edge
		(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX v0,
		const VERTEX_INDEX v1,
		const SCALAR_TYPE isovalue,
		bool &is_intersect)
	{
		is_intersect = is_gt_min_le_max(scalar_grid, v0, v1, isovalue);
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
			rmin[d]=std::max(coord1[d]-1, coord2[d]-1);
			rmax[d]=std::min(coord1[d]+2, coord2[d]+2);
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
	//
	// Find a good mapping for a vertex. 
	bool find_good_mapping(
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		const MERGESHARP::ISOVERT & isovert, 
		const VERTEX_INDEX gcube_index_from_v, // gcube_index of the from_vertex
		const std::vector<VERTEX_INDEX> & connected_sharp, // list of gcube_indices.
		VERTEX_INDEX & gcube_index_to_vertex // gcube_index of the vertex which is a good mapping
		)
	{
		//debug
		/*using namespace std;
		cout <<"connected sharp to "<< isovert.gcube_list[gcube_index_from_v].cube_index<<" this is the neighbor being mapped"<<endl;
		for (int i = 0; i < connected_sharp.size(); i++)
		{
			cout <<isovert.gcube_list[connected_sharp[i]].cube_index<<endl;
		}*/


		if(connected_sharp.size() < 1)
		{
			return false;
		}
		else if(connected_sharp.size() <= 2)
		{


			for (int i = 0; i < connected_sharp.size(); i++)
			{
				if(connected_sharp[i] == gcube_index_to_vertex)
				{
					gcube_index_to_vertex = connected_sharp[i];

					return true;
				}
			}
			return false;
		}
		else if(connected_sharp.size() > 3)
		{  
			return false;
		}
		else
		{
			//check for degree 3
			short num_of_deg3 = 0;
			int deg3_index = 0;
			for (int i = 0; i < connected_sharp.size(); i++)
			{
				if( isovert.gcube_list[connected_sharp[i]].num_eigenvalues == 3)
				{
					num_of_deg3++;
					deg3_index = i;
				}
			}
			if(num_of_deg3 == 1)
			{
				gcube_index_to_vertex = connected_sharp[deg3_index];
				return true;
			}
			else if(num_of_deg3 > 1)
			{
				using namespace std;
				//COORD_TYPE coord1[DIM3]= {0.0,0.0,0.0};
				//scalar_grid.ComputeCoord(
				//	isovert.gcube_list[gcube_index_from_v].cube_index, coord1);
				//cout <<"ERROR "<< isovert.gcube_list[gcube_index_from_v].cube_index
				//	<<" "<< coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]
				//<<" connected to more than one degree 3 vertex.\n";
				//cerr <<"ERROR "<< isovert.gcube_list[gcube_index_from_v].cube_index
				//	<<" "<< coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]
				//<<" connected to more than one degree 3 vertex.\n";
				return false;
			}
			else{

				for (int k = 0; k < 3; k++)
				{
					if (connected_sharp[k] == gcube_index_to_vertex)
					{

						VERTEX_INDEX k1 = (k+1)%3;
						VERTEX_INDEX k2 = (k+2)%3;

						VERTEX_INDEX cube_index1 = 
							isovert.gcube_list[connected_sharp[k1]].cube_index;
						VERTEX_INDEX cube_index2 = 
							isovert.gcube_list[connected_sharp[k2]].cube_index;
						bool flag = are_connected(scalar_grid, cube_index1,
							cube_index2, isovalue);

						//debug
						/*VERTEX_INDEX cube_index_temp = 
						isovert.gcube_list[connected_sharp[k]].cube_index;
						bool flag2 = are_connected(scalar_grid, cube_index1,
						cube_index_temp, isovalue);
						cout << cube_index1 <<" and " << cube_index_temp <<" are among them connected and "
						<<flag2 <<" ";

						bool flag3 = are_connected(scalar_grid, cube_index2,
						cube_index_temp, isovalue);
						cout << cube_index2 <<" and " << cube_index_temp <<" are among them connected and "
						<<flag2 <<" \n";*/
						//debug_end

						//debug
						//if(flag == false)// The other two vertices are not connected.
						if(flag == false)
						{
							gcube_index_to_vertex = connected_sharp[k];
							return true;
						}
						/*else if (flag == true)
						{
						bool temp = flag2 && flag3;
						cout <<" temp "<< temp <<" ";
						if (temp == false)
						{
						gcube_index_to_vertex = connected_sharp[k];
						return true;
						}
						else
						return false;
						}*/
						else
						{
							//cout << cube_index1 <<"  and " << cube_index2 <<" are conencted."<<endl;
						}
					}
				}
				//cout <<"all mappping failed"<<endl;
				return false;
			}
		}
	}

	// Map isosurface vertices adjacent to sharp and convered cubes to 
	// Selected cubes. 
	void map_adjacent_cubes_extended
		(
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		std::vector<NUM_TYPE> sorted_gcube_list;
		using namespace MERGESHARP;
		const NUM_TYPE num_gcube = isovert.gcube_list.size();
		VERTEX_INDEX cube_index_covered, neighbor_index;

		SHARPISO_GRID_NEIGHBORS gridn;
		// Set size of grid neighbors grid.
		gridn.SetSize(isovert.sharp_ind_grid);

		//setup the  sorted_gcube_list
		sort_gcube_list(isovert.gcube_list, sorted_gcube_list);


		// Set cubes which share facets with selected cubes.
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) 
		{
			//index to sharp cube in the sorted gcube list
			NUM_TYPE gcube_index = sorted_gcube_list[i];

			//debug
			using namespace std;


			//covered and not selected.
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {
					//cube index of the sharp cube
					cube_index_covered = isovert.gcube_list[gcube_index].cube_index;

					//debug
					COORD_TYPE coord1[DIM3];
					scalar_grid.ComputeScaledCoord(cube_index_covered, coord1);
					cout <<"******covered sharpF  "<< coord1[0]<<" "<< coord1[1]<<" "<<  coord1[2] << endl;


					if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
						// Cube cube_index is an interior cube.



						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
						{

							neighbor_index = gridn.CubeNeighborF(cube_index_covered, j);
							//k is an entry into the gcube list 
							INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
							if(k == ISOVERT::NO_INDEX)
							{
								continue;
							}


							//map to cube where cube index covered is mapped.
							VERTEX_INDEX to_cube = gcube_map[gcube_index];

							//DEBUG

							VERTEX_INDEX v = isovert.gcube_list[to_cube].cube_index;
							COORD_TYPE coord1[DIM3];
							using namespace std;
							scalar_grid.ComputeScaledCoord(v, coord1);
							cout <<"\n---------------------------------------------\n";
							cout <<"to_cube scaled "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
							cout <<"linf "<< isovert.gcube_list[gcube_index].linf_dist<<endl;
							cout <<"neighbor gcubeindex F "<<k <<" neighbor cubeindex " << neighbor_index <<endl;
							cout <<" neighbor from "<< isovert.gcube_list[k].isovert_coord[0]<<" "
								<< isovert.gcube_list[k].isovert_coord[1]<<" "<<
								isovert.gcube_list[k].isovert_coord[2];

							scalar_grid.ComputeCoord(v, coord1);
							cout <<" to_cube unscal  {"<< v <<"} with coord {" 
								<< coord1[0]<<" "<<coord1[1]<<" "<< coord1[2]<<"}";

							//debug end

							std::vector<VERTEX_INDEX> connected_sharp;
							//NOT_ACTIVE

							find_connected_sharp(k, isovert, gcube_map, connected_sharp);

							bool flag = find_good_mapping(scalar_grid, isovalue, 
								isovert, k, connected_sharp, to_cube);
							if(flag)
							{
								//debug
								cout <<"map Vertex: "<< endl;
								map_iso_vertex(isovert.gcube_list, k, to_cube, gcube_map);
								//isovert.gcube_list[gcube_index].flag = COVERED_B_GCUBE;
							}
							//debug
							else
							{
								//cout <<"no good mapping" << endl;
							}
						}


						/*
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
						{

						neighbor_index = gridn.CubeNeighborE(cube_index_covered, j);

						INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
						//NOT_ACTIVE
						if(k == -1)
						continue;
						//map to cube where cube index covered is mapped.
						VERTEX_INDEX to_cube = gcube_map[gcube_index];
						VERTEX_INDEX v = isovert.gcube_list[to_cube].cube_index;
						COORD_TYPE coord1[DIM3];
						using namespace std;
						scalar_grid.ComputeScaledCoord(v, coord1);
						cout <<"to_cube  scaled "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;


						scalar_grid.ComputeCoord(v, coord1);
						cout <<"neighbor gcubeindex E "<<k <<" neighbor cubeindex " << neighbor_index <<endl;
						cout <<"from "<< isovert.gcube_list[k].isovert_coord[0]<<" "
						<< isovert.gcube_list[k].isovert_coord[1]<<" "<<
						isovert.gcube_list[k].isovert_coord[2];
						cout <<" to_cube unscal  "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];

						std::vector<VERTEX_INDEX> connected_sharp;

						find_connected_sharp(k, isovert, gcube_map, connected_sharp);
						bool flag = find_good_mapping(scalar_grid, isovalue, 
						isovert, k, connected_sharp, to_cube);
						if(flag)
						{
						//cout <<"successfully mapped"<< endl;
						map_iso_vertex(isovert.gcube_list, k, to_cube, gcube_map);
						//isovert.gcube_list[gcube_index].flag = COVERED_B_GCUBE;
						}
						else
						{
						cout <<"no good mapping" << endl;
						}
						}//

						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
						{

						neighbor_index = gridn.CubeNeighborV(cube_index_covered, j);

						INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
						//NOT_ACTIVE
						if(k == -1)
						continue;
						//map to cube where cube index covered is mapped.
						VERTEX_INDEX to_cube = gcube_map[gcube_index];
						VERTEX_INDEX v = isovert.gcube_list[to_cube].cube_index;
						COORD_TYPE coord1[DIM3];
						using namespace std;
						scalar_grid.ComputeScaledCoord(v, coord1);
						cout <<"to_cube  scaled "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
						scalar_grid.ComputeCoord(v, coord1);
						cout <<"neighbor gcubeindex V "<<k <<" neighbor cubeindex " << neighbor_index <<endl;
						cout <<"from "<< isovert.gcube_list[k].isovert_coord[0]<<" "
						<< isovert.gcube_list[k].isovert_coord[1]<<" "<<
						isovert.gcube_list[k].isovert_coord[2];
						cout <<" to_cube unscal  "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2];
						std::vector<VERTEX_INDEX> connected_sharp;
						find_connected_sharp(k, isovert, gcube_map, connected_sharp);
						bool flag = find_good_mapping(scalar_grid, isovalue, 
						isovert, k, connected_sharp, to_cube);
						if(flag)
						{
						//cout <<"successfully mapped"<< endl;
						map_iso_vertex(isovert.gcube_list, k, to_cube, gcube_map);
						//isovert.gcube_list[gcube_index].flag = COVERED_B_GCUBE;
						}
						else
						{
						cout <<"no good mapping" << endl;
						}
						}*/


					}
					else {
						// *** Handle boundary case. ***
					}
			}
		}

		// Set cubes which share edges with selected cubes.
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			//covered and not selected.
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {

					cube_index_covered = isovert.gcube_list[gcube_index].cube_index;

					//debug
					COORD_TYPE coord1[DIM3];
					using namespace std;
					scalar_grid.ComputeScaledCoord(cube_index_covered, coord1);
					cout <<"******covered sharpE  "<< coord1[0]<<" "<< coord1[1]<<" "<<  coord1[2] << endl;
					cout <<"linf "<< isovert.gcube_list[gcube_index].linf_dist<<endl;
					if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
						// Cube cube_index is an interior cube.

						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) {

							neighbor_index = gridn.CubeNeighborE(cube_index_covered, j);

							INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
							//NOT_ACTIVE
							if(k == -1)
								continue;
							//map to cube where cube index covered is mapped.
							VERTEX_INDEX to_cube = gcube_map[gcube_index];
							//debug

							VERTEX_INDEX v = isovert.gcube_list[to_cube].cube_index;
							COORD_TYPE coord1[DIM3];
							using namespace std;
							scalar_grid.ComputeScaledCoord(v, coord1);
							cout <<"\n---------------------------------------------\n";
							cout <<"to_cube  scaled "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;


							scalar_grid.ComputeCoord(v, coord1);
							cout <<"neighbor gcubeindex E "<<k <<" neighbor cubeindex " << neighbor_index <<endl;
							cout <<" neighbor from "<< isovert.gcube_list[k].isovert_coord[0]<<" "
								<< isovert.gcube_list[k].isovert_coord[1]<<" "<<
								isovert.gcube_list[k].isovert_coord[2];


							cout <<" to_cube unscal  {"<< v <<"} with coord {" 
								<< coord1[0]<<" "<<coord1[1]<<" "<< coord1[2]<<"}";

							//debug_end

							std::vector<VERTEX_INDEX> connected_sharp;

							find_connected_sharp(k, isovert, gcube_map, connected_sharp);
							bool flag = find_good_mapping(scalar_grid, isovalue, 
								isovert, k, connected_sharp, to_cube);
							if(flag)
							{
								//debug
								cout <<"map Vertex: "<< endl;

								map_iso_vertex(isovert.gcube_list, k, to_cube, gcube_map);
							}
							//debug
							else
							{
								cout <<"no good mapping" << endl;
							}
						}//
					}
					else {
						// *** Handle boundary case. ***
					}
			}
		}
		//covered and not selected.

		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {

					cube_index_covered = isovert.gcube_list[gcube_index].cube_index;
					//debug
					using namespace std;
					COORD_TYPE coord1[DIM3];
					scalar_grid.ComputeScaledCoord(cube_index_covered, coord1);
					cout <<"******covered sharpV  "<< coord1[0]<<" "<< coord1[1]<<" "<<  coord1[2] << endl;
					cout <<"linf "<< isovert.gcube_list[gcube_index].linf_dist<<endl;
					//debug_end

					if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
						// Cube cube_index is an interior cube.

						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) {

							neighbor_index = gridn.CubeNeighborV(cube_index_covered, j);

							INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
							//NOT_ACTIVE
							if(k == ISOVERT::NO_INDEX)
								continue;
							//map to cube where cube index covered is mapped.
							VERTEX_INDEX to_cube = gcube_map[gcube_index];
							//debug

							VERTEX_INDEX v = isovert.gcube_list[to_cube].cube_index;
							COORD_TYPE coord1[DIM3];
							using namespace std;
							scalar_grid.ComputeScaledCoord(v, coord1);
							cout <<"\n---------------------------------------------\n";
							cout <<"to_cube  scaled "<< v <<" " << coord1[0]<<" "<<coord1[1]<<" "<< coord1[2] << endl;
							scalar_grid.ComputeCoord(v, coord1);
							cout <<"neighbor gcubeindex V ["<<k <<"], neighbor cubeindex [" << neighbor_index<<"]"<<endl;
							cout <<" neighbor from "<< isovert.gcube_list[k].isovert_coord[0]<<" "
								<< isovert.gcube_list[k].isovert_coord[1]<<" "<<
								isovert.gcube_list[k].isovert_coord[2];
							cout <<" to_cube unscal  {"<< v <<"} with coord {" 
								<< coord1[0]<<" "<<coord1[1]<<" "<< coord1[2]<<"}";

							//debug_end

							std::vector<VERTEX_INDEX> connected_sharp;
							find_connected_sharp(k, isovert, gcube_map, connected_sharp);
							bool flag = find_good_mapping(scalar_grid, isovalue, 
								isovert, k, connected_sharp, to_cube);
							if(flag)
							{
								//debug
								cout <<"map Vertex: "<< endl;
								map_iso_vertex(isovert.gcube_list, k, to_cube, gcube_map);
							}
							else
							{
								cout <<"no good mapping" << endl;
							}
						}
					}
					else {
						// *** Fill in. ***
					}
			}
		}
	}
	/*
	* Compute mapping original
	*/
	void compute_mapping_temp(
		const VERTEX_INDEX neighbor_cube_index,
		const VERTEX_INDEX neighbor_gcube_index,
		VERTEX_INDEX & tocube_gcube_index,
		const SCALAR_TYPE isovalue,
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		MERGESHARP::ISOVERT & isovert,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map
		)
	{
		std::vector<VERTEX_INDEX> connected_sharp;
		find_connected_sharp(neighbor_gcube_index, isovert, gcube_map, connected_sharp);

		bool found_mapping = find_good_mapping(scalar_grid, isovalue, 
			isovert, neighbor_gcube_index, connected_sharp, tocube_gcube_index);
		if(found_mapping)
		{
			map_iso_vertex(isovert.gcube_list, neighbor_gcube_index, tocube_gcube_index, gcube_map);
		}
	}
	/*
	* Compute mapping version 2 for new extended merge 
	*/
	void compute_mapping(
		const VERTEX_INDEX neighbor_cube_index,
		const VERTEX_INDEX neighbor_gcube_index,
		VERTEX_INDEX & tocube_gcube_index,
		const SCALAR_TYPE isovalue,
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		MERGESHARP::ISOVERT & isovert,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map
		)
	{
		COORD_TYPE c1coord[DIM3]={0.0,0.0,0.0};
		COORD_TYPE c2coord[DIM3]={0.0,0.0,0.0};

		VERTEX_INDEX tocube_cube_index = isovert.gcube_list[tocube_gcube_index].cube_index;
		scalar_grid.ComputeScaledCoord( tocube_cube_index, c1coord);

		scalar_grid.ComputeScaledCoord( neighbor_cube_index, c2coord);


		COORD_TYPE vdiff[DIM3]={0.0,0.0,0.0};
		IJK::subtract_coord_3D(c1coord,c2coord,vdiff);
		int m=0;
		for (int i = 0; i < DIM3; i++)
		{
			if(abs(int(vdiff[i])) == 2)
			{
				m++;
			}
		}
		if (m<2)
		{

			std::vector<VERTEX_INDEX> connected_sharp;
			find_connected_sharp(neighbor_gcube_index, isovert, gcube_map, connected_sharp);

			bool found_mapping = find_good_mapping(scalar_grid, isovalue, 
				isovert, neighbor_gcube_index, connected_sharp, tocube_gcube_index);
			if(found_mapping)
			{
				//debug
				/*using namespace std;
				COORD_TYPE cc[DIM3] = {0.0,0.0,0.0};
				scalar_grid.ComputeCoord(neighbor_cube_index, cc);
				cout <<"\nmapping neighbor "<<neighbor_cube_index<<" {"<<cc[0]<<","<<cc[1]<<","<<cc[2]<<"}"
					<<" to cube " << tocube_cube_index<<" ";*/

				map_iso_vertex(isovert.gcube_list, neighbor_gcube_index, tocube_gcube_index, gcube_map);

			}
		}
		else{
			using namespace std;
			//cout <<"coord differ by 2 vdiff "<<vdiff[0]<<" "<<vdiff[1]<<" "<<vdiff[2]<<endl;
		}
	}
	/*
	* Check FACE neighbors
	* @param covered_gcube_index, gcube_index of the sharp covered cube
	* @param covered_cube_index, cube index of the sharp covered cube
	*/
	void check_face_neighbors(
		const VERTEX_INDEX covered_gcube_index,
		const SCALAR_TYPE isovalue,
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SHARPISO_GRID_NEIGHBORS &gridn,
		MERGESHARP::ISOVERT & isovert,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		VERTEX_INDEX covered_cube_index= 
			isovert.gcube_list[covered_gcube_index].cube_index;
		////debug
		using namespace std;
		//cout <<"\n**********covered sharp Face "<<covered_cube_index<<"\n";

		VERTEX_INDEX neighbor_cube_index;
		for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
		{
			neighbor_cube_index = gridn.CubeNeighborF(covered_cube_index, j);
			//neighbor_gcube_index is an entry into the gcube list 
			INDEX_DIFF_TYPE neighbor_gcube_index
				= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
			//debug
			//cout <<"j "<< j <<" cube index "<< neighbor_cube_index<<endl;
			
			if(neighbor_gcube_index == ISOVERT::NO_INDEX)
			{
				continue;
			}
			//cout <<"intersects surface "<< endl;
			//map to cube where cube index covered is mapped.
			VERTEX_INDEX tocube_gcube_index = gcube_map[covered_gcube_index];
			compute_mapping (neighbor_cube_index, neighbor_gcube_index, 
				tocube_gcube_index, isovalue, scalar_grid, isovert, gcube_map);
		}
	}
	/*
	* Check EDGE neighbors
	* @param covered_gcube_index, gcube_index of the sharp covered cube
	* @param covered_cube_index, cube index of the sharp covered cube
	*/
	void check_edge_neighbors(
		const VERTEX_INDEX covered_gcube_index,
		const SCALAR_TYPE isovalue,
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SHARPISO_GRID_NEIGHBORS &gridn,
		MERGESHARP::ISOVERT & isovert,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		VERTEX_INDEX covered_cube_index= 
			isovert.gcube_list[covered_gcube_index].cube_index;
		////debug
		using namespace std;
		//cout <<"\n**********covered sharp Edge "<<covered_cube_index<<"\n";
		VERTEX_INDEX neighbor_cube_index;
		for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
		{
			neighbor_cube_index = gridn.CubeNeighborE(covered_cube_index, j);
			//neighbor_gcube_index is an entry into the gcube list 
			INDEX_DIFF_TYPE neighbor_gcube_index
				= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
			//debug
			COORD_TYPE  cc[DIM3]={0.0,0.0,0.0};
			scalar_grid.ComputeCoord(neighbor_cube_index,cc);
			//cout <<"j "<< j <<" cube index "<< neighbor_cube_index
			//	<<" "<<cc[0]<<" "<<cc[1]<<" "<<cc[2]
			//	<<endl;
			if(neighbor_gcube_index == ISOVERT::NO_INDEX)
			{
				continue;
			}
			//map to cube where cube index covered is mapped.
			VERTEX_INDEX tocube_gcube_index = gcube_map[covered_gcube_index];
			compute_mapping (neighbor_cube_index, neighbor_gcube_index, 
				tocube_gcube_index, isovalue, scalar_grid, isovert, gcube_map);
		}
	}
	/*
	* Check VERTEX neighbors
	* @param covered_gcube_index, gcube_index of the sharp covered cube
	* @param covered_cube_index, cube index of the sharp covered cube
	*/
	void check_vertex_neighbors(
		const VERTEX_INDEX covered_gcube_index,
		const SCALAR_TYPE isovalue,
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const SHARPISO_GRID_NEIGHBORS &gridn,
		MERGESHARP::ISOVERT & isovert,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		VERTEX_INDEX covered_cube_index= 
			isovert.gcube_list[covered_gcube_index].cube_index;
		using namespace std;
		//cout <<"\n**********covered sharp Vertex "<<covered_cube_index<<"\n";
		VERTEX_INDEX neighbor_cube_index;
		for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
		{
			neighbor_cube_index = gridn.CubeNeighborV(covered_cube_index, j);
			//neighbor_gcube_index is an entry into the gcube list 
			INDEX_DIFF_TYPE neighbor_gcube_index
				= isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
			//debug
			COORD_TYPE  cc[DIM3]={0.0,0.0,0.0};
			scalar_grid.ComputeCoord(neighbor_cube_index,cc);
			/*cout <<"j "<< j <<" cube index "<< neighbor_cube_index
				<<" "<<cc[0]<<" "<<cc[1]<<" "<<cc[2]
				<<endl;*/
			if(neighbor_gcube_index == ISOVERT::NO_INDEX)
			{
				//cout <<"not intersected"<<endl;
				continue;
			}
			//map to cube where cube index covered is mapped.
			VERTEX_INDEX tocube_gcube_index = gcube_map[covered_gcube_index];
			compute_mapping (neighbor_cube_index, neighbor_gcube_index, 
				tocube_gcube_index, isovalue, scalar_grid, isovert, gcube_map);
		}
	}
	/*
	* Map extended version 3
	*/
	// Map isosurface vertices adjacent to sharp and convered cubes to 
	// Selected cubes. 
	void map_adjacent_cubes_extended_version3
		(
		const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		std::vector<NUM_TYPE> sorted_gcube_list;
		using namespace MERGESHARP;
		const NUM_TYPE num_gcube = isovert.gcube_list.size();
		VERTEX_INDEX cube_index_covered, neighbor_cube_index, neighbor_index;

		SHARPISO_GRID_NEIGHBORS gridn;
		// Set size of grid neighbors grid.
		gridn.SetSize(isovert.sharp_ind_grid);

		//setup the  sorted_gcube_list
		sort_gcube_list(isovert.gcube_list, sorted_gcube_list);

		//FACE
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) 
		{
			//index to sharp cube in the sorted gcube list
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			//covered and not selected.
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {
					//check boundary 
					if (isovert.gcube_list[gcube_index].boundary_bits == 0) 
					{
						check_face_neighbors(gcube_index, isovalue,
							scalar_grid, gridn, isovert, gcube_map);
					}
					else {
						// *** Handle boundary case. ***
					}
			}
		}

		//EDGE
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) 
		{
			//index to sharp cube in the sorted gcube list
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			//covered and not selected.
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {
					//check boundary 
					if (isovert.gcube_list[gcube_index].boundary_bits == 0) 
					{
						check_edge_neighbors(gcube_index, isovalue, 
							scalar_grid, gridn, isovert, gcube_map);
					}
					else {
						// *** Handle boundary case. ***
					}
			}
		}

		//VERTICES
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) 
		{
			//index to sharp cube in the sorted gcube list
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			//covered and not selected.
			if (isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {
					//check boundary 
					if (isovert.gcube_list[gcube_index].boundary_bits == 0) 
					{
						check_vertex_neighbors(gcube_index, isovalue, 
							scalar_grid, gridn, isovert, gcube_map);
					}
					else {
						// *** Handle boundary case. ***
					}
			}
		}
	}
	//map_extended version 3 end.
	// 2nd version of the map cube extension 

	void map_adjacent_cubes_extended_version2
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		const MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map,
		const SHARP_ISOVERT_PARAM & isovert_param)
	{
		std::vector<NUM_TYPE> sorted_gcube_list;
		using namespace MERGESHARP;
		const NUM_TYPE num_gcube = isovert.gcube_list.size();
		//Cube indices
		VERTEX_INDEX covered_cube_index, neighbor_cube_index;

		SHARPISO_GRID_NEIGHBORS gridn;
		// Set size of grid neighbors grid.
		gridn.SetSize(isovert.sharp_ind_grid);

		//setup the  sorted_gcube_list
		sort_gcube_list(isovert.gcube_list, sorted_gcube_list);

		// Set cubes which share facets with selected cubes.
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) 
		{
			//index to covered sharp cube in the sorted gcube list
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			//covered and not selected.
			if ( isovert.gcube_list[gcube_index].flag == COVERED_A_GCUBE
				&& isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) 
			{
				//cube index of the covered sharp cube
				covered_cube_index = isovert.gcube_list[gcube_index].cube_index;

				if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
					// Cube cube_index is an interior cube.
					COORD_TYPE real_cube_coord[DIM3];
					for (int d = 0; d < DIM3; d++)
					{
						COORD_TYPE c = isovert.gcube_list[gcube_index].isovert_coord[d];
						real_cube_coord[d] = int(c/scalar_grid.Spacing(d));
					}

					// Cube index of the actual cube in which the vertex is present.
					VERTEX_INDEX real_cube_index =  scalar_grid.ComputeVertexIndex(real_cube_coord);
					//Where the current covered sharp cube is going
					VERTEX_INDEX to_cube_gcube_index = gcube_map[gcube_index];
					VERTEX_INDEX to_cube_cube_index = isovert.gcube_list[to_cube_gcube_index].cube_index;
					// The second check ensures that the covered cube puts a point somewhere other than itself.
					if(( to_cube_cube_index != real_cube_index) && (real_cube_index != covered_cube_index))
					{ 
						using namespace std;
						COORD_TYPE coord1[DIM3];
						//Find the adjacent cubes of real_cube_index 
						set<VERTEX_INDEX> real_cube_adjacencies_cube_index;
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
						{
							VERTEX_INDEX v = gridn.CubeNeighborF(real_cube_index, j);
							VERTEX_INDEX v_gcube_index = isovert.sharp_ind_grid.Scalar(v);
							if(v_gcube_index != ISOVERT::NO_INDEX)
							{
								real_cube_adjacencies_cube_index.insert(v);
							}
						}
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
						{
							VERTEX_INDEX v = gridn.CubeNeighborE(real_cube_index, j);
							VERTEX_INDEX v_gcube_index = isovert.sharp_ind_grid.Scalar(v);
							if(v_gcube_index != ISOVERT::NO_INDEX)
							{
								real_cube_adjacencies_cube_index.insert(v);
							}
						}
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
						{
							VERTEX_INDEX v = gridn.CubeNeighborV(real_cube_index, j);
							VERTEX_INDEX v_gcube_index = isovert.sharp_ind_grid.Scalar(v);
							if(v_gcube_index != ISOVERT::NO_INDEX)
							{
								real_cube_adjacencies_cube_index.insert(v);
							}
						}

						const int bin_width = isovert_param.bin_width;
						BIN_GRID<VERTEX_INDEX> bin_grid;

						init_bin_grid(scalar_grid, bin_width, bin_grid);
						bool check_angle = true;
						VERTEX_INDEX v1,v2; // not used


						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
						{
							neighbor_cube_index = gridn.CubeNeighborF(covered_cube_index, j);
							//k is an entry into the gcube list 
							INDEX_DIFF_TYPE neighbor_gcube_index = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
							if(neighbor_gcube_index == ISOVERT::NO_INDEX)
							{
								//Not intersected by the surface
								continue;
							}
							else
							{
								//
								if(isovert.gcube_list[neighbor_gcube_index].flag != COVERED_A_GCUBE && 
									isovert.gcube_list[neighbor_gcube_index].flag != SELECTED_GCUBE)
								{

									//c3 is adjacent to c2
									set<VERTEX_INDEX>::iterator it;
									it = real_cube_adjacencies_cube_index.find(neighbor_cube_index);
									if(it != real_cube_adjacencies_cube_index.end())
									{
										std::vector<VERTEX_INDEX> connected_sharp;
										find_connected_sharp(neighbor_gcube_index, isovert, gcube_map, connected_sharp);
										bool flag = find_good_mapping(scalar_grid, isovalue, 
											isovert, neighbor_gcube_index, connected_sharp, 
											to_cube_gcube_index);
										if(flag)
										{
											map_iso_vertex(isovert.gcube_list, neighbor_gcube_index, to_cube_gcube_index, gcube_map);

										}
									}
								}
							}
						}//endfor
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) 
						{
							neighbor_cube_index = gridn.CubeNeighborE(covered_cube_index, j);
							//k is an entry into the gcube list 
							INDEX_DIFF_TYPE neighbor_gcube_index = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
							if(neighbor_gcube_index == ISOVERT::NO_INDEX)
							{
								//Not intersected by the surface
								continue;
							}
							else
							{
								//
								if(isovert.gcube_list[neighbor_gcube_index].flag != COVERED_A_GCUBE && 
									isovert.gcube_list[neighbor_gcube_index].flag != SELECTED_GCUBE)
								{

									//c3 is adjacent to c2
									set<VERTEX_INDEX>::iterator it;
									it = real_cube_adjacencies_cube_index.find(neighbor_cube_index);
									if(it != real_cube_adjacencies_cube_index.end())
									{
										std::vector<VERTEX_INDEX> connected_sharp;
										find_connected_sharp(neighbor_gcube_index, isovert, gcube_map, connected_sharp);
										bool flag = find_good_mapping(scalar_grid, isovalue, 
											isovert, neighbor_gcube_index, connected_sharp, 
											to_cube_gcube_index);
										if(flag)
										{
											map_iso_vertex(isovert.gcube_list, neighbor_gcube_index, to_cube_gcube_index, gcube_map);
										}

									}

								}

							}
						}//endfor
						for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) 
						{
							neighbor_cube_index = gridn.CubeNeighborV(covered_cube_index, j);
							//k is an entry into the gcube list 
							INDEX_DIFF_TYPE neighbor_gcube_index = isovert.sharp_ind_grid.Scalar(neighbor_cube_index);
							if(neighbor_gcube_index == ISOVERT::NO_INDEX)
							{
								//Not intersected by the surface
								continue;
							}
							else
							{
								//
								if(isovert.gcube_list[neighbor_gcube_index].flag != COVERED_A_GCUBE && 
									isovert.gcube_list[neighbor_gcube_index].flag != SELECTED_GCUBE)
								{
									//c3 is adjacent to c2
									set<VERTEX_INDEX>::iterator it;
									it = real_cube_adjacencies_cube_index.find(neighbor_cube_index);
									if(it != real_cube_adjacencies_cube_index.end())
									{
										std::vector<VERTEX_INDEX> connected_sharp;
										find_connected_sharp(neighbor_gcube_index, isovert, gcube_map, connected_sharp);
										bool flag = find_good_mapping(scalar_grid, isovalue, 
											isovert, neighbor_gcube_index, connected_sharp, 
											to_cube_gcube_index);
										if(flag)
										{
											map_iso_vertex(isovert.gcube_list, neighbor_gcube_index, to_cube_gcube_index, gcube_map);
										}
									}
								}
							}
						}//endfor

					}//endif
				}

			}
		}
	}


	// Map isosurface vertices adjacent to selected cubes
	//  to the isosurface vertex in the selected cube.
	void map_adjacent_cubes
		(const MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		std::vector<NUM_TYPE> sorted_gcube_list;
		using namespace MERGESHARP;

		const NUM_TYPE num_gcube = isovert.gcube_list.size();
		VERTEX_INDEX cube_index, neighbor_index;
		SHARPISO_GRID_NEIGHBORS gridn;

		// Set size of grid neighbors grid.
		gridn.SetSize(isovert.sharp_ind_grid);

		sort_gcube_list(isovert.gcube_list, sorted_gcube_list);

		for (NUM_TYPE i = 0; i < num_gcube; i++)
		{ gcube_map[i] = i; }

		// Set cubes which share facets with selected cubes.
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
			//index to sharp cube in sorted gcube_list
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {
				// cube index of the sharp cube.
				cube_index = isovert.gcube_list[gcube_index].cube_index;

				if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
					// Cube cube_index is an interior cube.
					//find the neighbors
					for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) 
					{
						//neighbor to the sharp cube. (cube_index)
						neighbor_index = gridn.CubeNeighborF(cube_index, j);

						INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
						// k is the 'from-cube', and gcube_index is the 'to-cube'.
						map_iso_vertex(isovert.gcube_list, k, gcube_index, gcube_map);
					}
				}
				else {
					// *** Handle boundary case. ***
				}
			}
		}

		// Set cubes which share edges with selected cubes.
		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
			NUM_TYPE gcube_index = sorted_gcube_list[i];

			if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {

				cube_index = isovert.gcube_list[gcube_index].cube_index;

				if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
					// Cube cube_index is an interior cube.

					for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) {

						neighbor_index = gridn.CubeNeighborE(cube_index, j);

						INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
						map_iso_vertex(isovert.gcube_list, k, gcube_index, gcube_map);
					}
				}
				else {
					// *** Handle boundary case. ***
				}
			}
		}

		for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
			NUM_TYPE gcube_index = sorted_gcube_list[i];
			if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {

				cube_index = isovert.gcube_list[gcube_index].cube_index;

				if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
					// Cube cube_index is an interior cube.

					for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) {

						neighbor_index = gridn.CubeNeighborV(cube_index, j);

						INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
						map_iso_vertex(isovert.gcube_list, k, gcube_index, gcube_map);
					}

				}
				else {
					// *** Fill in. ***
				}

			}
		}
	}

	void set_first_gcube_isov
		(const ISOVERT & isovert, const std::vector<DUAL_ISOVERT> & iso_vlist,
		NUM_TYPE * first_gcube_isov)
	{
		// Scan iso_vlist from back to front.
		for (int i = 0; i < iso_vlist.size(); i++) {
			int j = iso_vlist.size()-i-1;
			VERTEX_INDEX cube_index = iso_vlist[j].cube_index;
			VERTEX_INDEX gcube_index = isovert.sharp_ind_grid.Scalar(cube_index);
			first_gcube_isov[gcube_index] = j;
		}
	}

	void map_gcube_indices(const std::vector<VERTEX_INDEX> & gcube_map,
		std::vector<VERTEX_INDEX> & gcube_index)
	{
		for (VERTEX_INDEX i = 0; i < gcube_index.size(); i++) {
			VERTEX_INDEX gcube_index_i = gcube_index[i];
			gcube_index[i] = gcube_map[gcube_index_i];
		}
	}

	// Map cube_list[i] to target cube determined by gcube_map.
	void map_cube_list(const ISOVERT & isovert,
		const std::vector<VERTEX_INDEX> & gcube_map,
		std::vector<VERTEX_INDEX> & cube_list)
	{
		for (VERTEX_INDEX i = 0; i < cube_list.size(); i++) {
			VERTEX_INDEX cube_index0 = cube_list[i];
			VERTEX_INDEX gcube_index0 = isovert.sharp_ind_grid.Scalar(cube_index0);
			VERTEX_INDEX gcube_index1 = gcube_map[gcube_index0];
			cube_list[i] = isovert.gcube_list[gcube_index1].cube_index;
		}
	}

	// If cube0 containing isosurface vertex isov_index[i] is merged 
	//   with some other cube cube1, then map isov_index[i] 
	//   to the isosurface vertex in cube1.
	void map_isov_indices(const ISOVERT & isovert,
		const std::vector<VERTEX_INDEX> & gcube_map,
		const std::vector<DUAL_ISOVERT> & iso_vlist,
		std::vector<ISO_VERTEX_INDEX> & isov_index)
	{
		VERTEX_HASH_TABLE cube_hash;

		for (NUM_TYPE i = 0; i < iso_vlist.size(); i++) {
			VERTEX_INDEX cube_index = iso_vlist[i].cube_index;

			VERTEX_HASH_TABLE::iterator cube_iter = cube_hash.find(cube_index);
			if (cube_iter == cube_hash.end()) {
				cube_hash.insert(VERTEX_HASH_TABLE::value_type(cube_index, i));
			}
		}

		for (NUM_TYPE i = 0; i < isov_index.size(); i++) {
			ISO_VERTEX_INDEX j0 = isov_index[i];
			VERTEX_INDEX cube_index0 = iso_vlist[j0].cube_index;
			VERTEX_INDEX gcube_index0 = isovert.sharp_ind_grid.Scalar(cube_index0);
			VERTEX_INDEX gcube_index1 = gcube_map[gcube_index0];
			if (gcube_index1 != gcube_index0 ||
				isovert.gcube_list[gcube_index1].flag == SELECTED_GCUBE) {
					VERTEX_INDEX cube_index1 = isovert.gcube_list[gcube_index1].cube_index;
					VERTEX_HASH_TABLE::iterator cube1_iter = cube_hash.find(cube_index1);
					isov_index[i] = cube1_iter->second;
			}
		}
	}


	// Forward declaration
	void unmap_non_disk_isopatches
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info);

	void unmap_non_disk_isopatches
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info);

	void determine_gcube_map
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		const SHARP_ISOVERT_PARAM & sharp_isovert_param,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info)
	{
		map_adjacent_cubes(isovert, gcube_map);
		//DEBUG
		if(sharp_isovert_param.flag_map_extended){

			
			//map_adjacent_cubes_extended(scalar_grid, isovalue, isovert, gcube_map);
			map_adjacent_cubes_extended_version3(scalar_grid, isovalue, isovert, gcube_map);
		}

		if (sharp_isovert_param.flag_check_disk) {
			unmap_non_disk_isopatches
				(scalar_grid, isovalue, isovert, gcube_map, sharpiso_info);
		}
	}

	void determine_gcube_map
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		const SHARP_ISOVERT_PARAM & sharp_isovert_param,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info)
	{
		map_adjacent_cubes(isovert, gcube_map);
		if(sharp_isovert_param.flag_map_extended){

			//DEBUG going back to the old 
			//map_adjacent_cubes_extended(scalar_grid, isovalue, isovert, gcube_map);
			map_adjacent_cubes_extended_version3(scalar_grid, isovalue, isovert, gcube_map);
		}
		if (sharp_isovert_param.flag_check_disk) {
			unmap_non_disk_isopatches
				(scalar_grid, isodual_table, isovalue, isovert, gcube_map, 
				sharpiso_info);
		}
	}

}

// **************************************************
// Unmap non-disk isopatches
// **************************************************

namespace {

	void unmap_merged_cubes
		(MERGESHARP::ISOVERT & isovert, const VERTEX_INDEX cube_index0,
		const AXIS_SIZE_TYPE dist2cube,
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
	{
		const NUM_TYPE gcube_index0 = isovert.sharp_ind_grid.Scalar(cube_index0);
		std::vector<VERTEX_INDEX> merged_cube_list;

		isovert.gcube_list[gcube_index0].flag = NON_DISK_GCUBE;

		get_merged_cubes(isovert.sharp_ind_grid, isovert, cube_index0,
			gcube_map, dist2cube,  merged_cube_list);
		for (NUM_TYPE i = 0; i < merged_cube_list.size(); i++) {
			VERTEX_INDEX cube_index1 = merged_cube_list[i];
			NUM_TYPE gcube_index1 = isovert.sharp_ind_grid.Scalar(cube_index1);
			if (gcube_map[gcube_index1] == gcube_index0) {
				gcube_map[gcube_index1] = gcube_index1; 
				//DEBUG : COVERED_GCUBE is being replaced by COVERED_A_GCUBE and COVERED_B_GCUBE
				if (isovert.gcube_list[gcube_index1].flag == COVERED_A_GCUBE) {
					isovert.gcube_list[gcube_index1].flag = SMOOTH_GCUBE;
				}
			}
		}
	}

	// Reverse merges which create isopatches which are not disks.
	void unmap_non_disk_isopatches
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info)
	{
		const NUM_TYPE num_gcube = isovert.gcube_list.size();

		const int dist2cube = 1;
		std::vector<ISO_VERTEX_INDEX> tri_vert;
		std::vector<ISO_VERTEX_INDEX> quad_vert;

		bool passed_all_disk_checks;
		do {
			passed_all_disk_checks = true;

			for (NUM_TYPE i = 0; i < num_gcube; i++) {
				if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
					VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

					extract_dual_isopatch_incident_on
						(scalar_grid, isovalue, isovert, cube_index,
						gcube_map, dist2cube, tri_vert, quad_vert);

					IJK::reorder_quad_vertices(quad_vert);

					if (!is_isopatch_disk3D(tri_vert, quad_vert)) {
						unmap_merged_cubes(isovert, cube_index, dist2cube, gcube_map);
						sharpiso_info.num_non_disk_isopatches++;
						passed_all_disk_checks = false;
					}
				}
			}
		}
		while (!passed_all_disk_checks);
	}

	// Renumber vertices in vlist so that they are numbered starting at 0.
	// Construct cube list of distinct vertices in vlist.
	template <typename VTYPE0, typename VTYPE1>
	void construct_cube_list
		(std::vector<VTYPE0> & vlist, std::vector<VTYPE1> & cube_list)
	{
		VERTEX_HASH_TABLE vertex_hash;

		cube_list.clear();
		insert_vertex_list(vlist, vertex_hash);

		remap_vertex_list(vertex_hash, vlist, vlist);

		cube_list.resize(vertex_hash.size());
		for (VERTEX_HASH_TABLE::const_iterator vertex_iter = vertex_hash.begin();
			vertex_iter != vertex_hash.end(); vertex_iter++) {
				VERTEX_INDEX iv = vertex_iter->first;
				NUM_TYPE n = vertex_iter->second;
				cube_list[n] = iv;
		}
	}

	// Reverse merges which create isopatches which are not disks.
	void unmap_non_disk_isopatches
		(const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
		const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
		const SCALAR_TYPE isovalue,
		MERGESHARP::ISOVERT & isovert, 
		std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
		MERGESHARP::SHARPISO_INFO & sharpiso_info)
	{
		const NUM_TYPE num_gcube = isovert.gcube_list.size();

		const int dist2cube = 1;
		std::vector<ISO_VERTEX_INDEX> tri_vert;
		std::vector<ISO_VERTEX_INDEX> quad_vert;

		bool passed_all_disk_checks;
		do {
			passed_all_disk_checks = true;

			for (NUM_TYPE i = 0; i < num_gcube; i++) {
				if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {
					VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

					extract_dual_isopatch_incident_on_multi
						(scalar_grid, isodual_table, isovalue, isovert, 
						cube_index, gcube_map, dist2cube, tri_vert, quad_vert);
					IJK::reorder_quad_vertices(quad_vert);

					if (!is_isopatch_disk3D(tri_vert, quad_vert)) {
						unmap_merged_cubes(isovert, cube_index, dist2cube, gcube_map);
						isovert.gcube_list[i].flag = NON_DISK_GCUBE;
						sharpiso_info.num_non_disk_isopatches++;
						passed_all_disk_checks = false;
					}
				}
			}
		}
		while (!passed_all_disk_checks);
	}

	// Set table_index[i] of cube cube_list[i].
	// Table index is stored in isovert.
	template <typename VTYPE, typename ITYPE>
	void set_table_index
		(const MERGESHARP::ISOVERT & isovert, const std::vector<VTYPE> & cube_list,
		std::vector<ITYPE> & table_index)
	{
		table_index.resize(cube_list.size());

		for (NUM_TYPE i = 0; i < cube_list.size(); i++) {
			VTYPE cube_index = cube_list[i];
			VTYPE gcube_index = isovert.sharp_ind_grid.Scalar(cube_index);
			table_index[i] = isovert.gcube_list[gcube_index].table_index;
		}
	}


}



// **************************************************
// ROUTINE: is_isopatch_disk3D
// **************************************************

// Forward declarations:
template <typename VTYPE>
void renumber_tri_quad_vertices
	(const std::vector<VTYPE> & tri_vert,
	const std::vector<VTYPE> & quad_vert,
	std::vector<VTYPE> & new_tri_vert,
	std::vector<VTYPE> & new_quad_vert,
	NUM_TYPE & num_vert);

// Search cycle starting at iv0.
// @pre All vertices in cycle containing iv0 have is_visited set to false.
// @pre All vertices have degree two.
void search_cycle
	(const VERTEX_INDEX iv0, std::vector<CYCLE_VERTEX> & cycle_vertex);

// Return true if isopatch incident on vertex is a disk.
// @param tri_vert Triangle vertices.
// @param quad_vert Quadrilateral vertices in order around quadrilateral.
// @pre Assumes the boundary of the isopatch is the link of some vertex.
bool MERGESHARP::is_isopatch_disk3D
	(const std::vector<ISO_VERTEX_INDEX> & tri_vert,
	const std::vector<ISO_VERTEX_INDEX> & quad_vert)
{
	const NUM_TYPE num_tri = tri_vert.size()/NUM_VERT_PER_TRI;
	const NUM_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
	std::vector<ISO_VERTEX_INDEX> tri_vert2;
	std::vector<ISO_VERTEX_INDEX> quad_vert2;
	NUM_TYPE num_vert;

	VERTEX_HASH_TABLE vertex_hash;
	EDGE_HASH_TABLE edge_hash;

	// Renumber tri and quad vertices
	renumber_tri_quad_vertices
		(tri_vert, quad_vert, tri_vert2, quad_vert2, num_vert);

	// Check for edges in more than two isosurface polygons.
	insert_tri_quad_edges(tri_vert2, quad_vert2, edge_hash);

	for (EDGE_HASH_TABLE::const_iterator edge_iter = edge_hash.begin();
		edge_iter != edge_hash.end(); edge_iter++) {
			if (edge_iter->second > 2) {
				return(false); 
			}
	}

	// Check that boundary is a cycle or edge.
	std::vector<CYCLE_VERTEX> cycle_vertex(num_vert);

	construct_boundary_cycle(edge_hash, cycle_vertex);
	NUM_TYPE num_boundary_vertices = 0;
	NUM_TYPE first_adjacent = 0;
	for (NUM_TYPE i = 0; i < cycle_vertex.size(); i++) {
		NUM_TYPE num_adjacent = cycle_vertex[i].num_adjacent;
		if (num_adjacent == 2) {
			first_adjacent = i;
			num_boundary_vertices++;
		}
		else if (num_adjacent != 0) {
			return(false); 
		}
	}

	if (num_boundary_vertices < 3) { 
		// Disk must have at least three boundary cycle vertices.
		return(false); 
	}

	search_cycle(first_adjacent, cycle_vertex);

	for (NUM_TYPE i = 0; i < cycle_vertex.size(); i++) { 
		if (cycle_vertex[i].num_adjacent == 2) {
			if (!cycle_vertex[i].is_visited) { return(false); }
		}
	}

	return(true);
}


/// Get list of cubes merged with icube.
void MERGESHARP::get_merged_cubes
	(const SHARPISO_GRID & grid,
	const ISOVERT & isovert,
	const VERTEX_INDEX cube_index0,
	const std::vector<VERTEX_INDEX> & gcube_map,
	const AXIS_SIZE_TYPE dist2cube,
	std::vector<VERTEX_INDEX> & merged_cube_list)
{
	const int dimension = grid.Dimension();
	const VERTEX_INDEX gcube_index0 = 
		isovert.sharp_ind_grid.Scalar(cube_index0);
	VERTEX_INDEX region_iv0;
	IJK::ARRAY<AXIS_SIZE_TYPE> region_axis_size(dimension);

	IJK::compute_region_around_cube
		(cube_index0, dimension, grid.AxisSize(), dist2cube, 
		region_iv0, region_axis_size.Ptr());

	NUM_TYPE num_region_cubes;
	IJK::compute_num_grid_cubes
		(dimension, region_axis_size.PtrConst(), num_region_cubes);

	IJK::ARRAY<VERTEX_INDEX> region_cube_list(num_region_cubes);
	IJK::get_subgrid_cubes
		(dimension, grid.AxisSize(), region_iv0, region_axis_size.PtrConst(),
		region_cube_list.Ptr());

	for (NUM_TYPE i = 0; i < num_region_cubes; i++) {
		VERTEX_INDEX cube_index1 = region_cube_list[i];
		VERTEX_INDEX gcube_index1 = isovert.sharp_ind_grid.Scalar(cube_index1);
		if (gcube_map[gcube_index1] == gcube_index0) 
		{ merged_cube_list.push_back(cube_index1); }
	}
}

/// Get edges on boundary of merged cubes.
void MERGESHARP::get_merged_boundary_edges
	(const SHARPISO_GRID & grid,
	const std::vector<VERTEX_INDEX> & merged_cube_list,
	std::vector<EDGE_INDEX> & boundary_edge_list)
{
	typedef SHARPISO_GRID::DIMENSION_TYPE DTYPE;
	typedef std::unordered_map<EDGE_INDEX, NUM_TYPE> HASH_TABLE;

	HASH_TABLE edge_hash;

	for (NUM_TYPE i = 0; i < merged_cube_list.size(); i++) {
		VERTEX_INDEX cube_index = merged_cube_list[i];
		for (DTYPE edge_dir = 0; edge_dir < DIM3; edge_dir++) {
			for (NUM_TYPE k = 0; k < NUM_CUBE_FACET_VERTICES3D; k++) {
				VERTEX_INDEX iv = grid.FacetVertex(cube_index, edge_dir, k);
				EDGE_INDEX iedge = iv*DIM3+edge_dir;

				HASH_TABLE::iterator edge_iter = edge_hash.find(iedge);
				if (edge_iter == edge_hash.end()) {
					edge_hash.insert(HASH_TABLE::value_type(iedge, 1));
				}
				else {
					edge_iter->second++;
				}
			}
		}
	}

	for (HASH_TABLE::const_iterator edge_iter = edge_hash.begin();
		edge_iter != edge_hash.end(); edge_iter++) {
			if (edge_iter->second != NUM_QUAD_VERTICES) 
			{ boundary_edge_list.push_back(edge_iter->first); }
	}

}

/// Select edges from edge list which are in the grid interior.
void select_interior_grid_edges
	(const SHARPISO_GRID & grid,
	const std::vector<EDGE_INDEX> & edge_list,
	std::vector<EDGE_INDEX> & interior_edge_list)
{
	const int dimension = grid.Dimension();

	for (NUM_TYPE i = 0; i < edge_list.size(); i++) {
		EDGE_INDEX edge_index = edge_list[i];
		EDGE_INDEX iend0 = edge_index/dimension;
		int edge_dir = edge_index%dimension;
		if (!IJK::is_edge_on_grid_boundary
			(iend0, edge_dir, dimension, grid.AxisSize())) {
				interior_edge_list.push_back(edge_index);
		}
	}
}

// Extract dual isosurface patch with vertex in merged cube.
void MERGESHARP::extract_dual_quad_isopatch_incident_on
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const ISOVERT & isovert,
	const VERTEX_INDEX cube_index0,
	const std::vector<VERTEX_INDEX> & gcube_map,
	const AXIS_SIZE_TYPE dist2cube,
	std::vector<ISO_VERTEX_INDEX> & isoquad_cube,
	std::vector<FACET_VERTEX_INDEX> & facet_vertex)
{
	std::vector<VERTEX_INDEX> merged_cube_list;
	std::vector<EDGE_INDEX> boundary_edge_list;
	std::vector<EDGE_INDEX> edge_list;

	get_merged_cubes(scalar_grid, isovert, cube_index0, gcube_map, dist2cube, 
		merged_cube_list);
	get_merged_boundary_edges(scalar_grid, merged_cube_list, boundary_edge_list);
	select_interior_grid_edges(scalar_grid, boundary_edge_list, edge_list);
	extract_dual_isopoly_from_list(scalar_grid, isovalue, edge_list, 
		isoquad_cube, facet_vertex);
}

// Extract dual isosurface patch with vertex in merged cube.
void MERGESHARP::extract_dual_isopatch_incident_on
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const SCALAR_TYPE isovalue,
	const ISOVERT & isovert,
	const VERTEX_INDEX cube_index0,
	const std::vector<VERTEX_INDEX> & gcube_map,
	const AXIS_SIZE_TYPE dist2cube,
	std::vector<ISO_VERTEX_INDEX> & tri_vert,
	std::vector<ISO_VERTEX_INDEX> & quad_vert)
{
	std::vector<ISO_VERTEX_INDEX> isoquad_cube;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;

	tri_vert.clear();
	quad_vert.clear();

	extract_dual_quad_isopatch_incident_on
		(scalar_grid, isovalue, isovert, cube_index0, gcube_map, dist2cube,
		isoquad_cube, facet_vertex);

	map_cube_list(isovert, gcube_map, isoquad_cube);
	IJK::get_non_degenerate_quad_btlr(isoquad_cube, tri_vert, quad_vert);
}

// Extract dual isosurface patch with vertex in merged cube.
// Allow multiple isosurface vertices in each cube.
// Version returning iso_vlist.
void MERGESHARP::extract_dual_isopatch_incident_on_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	const SCALAR_TYPE isovalue,
	const ISOVERT & isovert,
	const VERTEX_INDEX cube_index,
	const std::vector<VERTEX_INDEX> & gcube_map,
	const AXIS_SIZE_TYPE dist2cube,
	std::vector<ISO_VERTEX_INDEX> & tri_vert,
	std::vector<ISO_VERTEX_INDEX> & quad_vert,
	std::vector<DUAL_ISOVERT> & iso_vlist)
{
	std::vector<ISO_VERTEX_INDEX> isoquad_cube;
	std::vector<FACET_VERTEX_INDEX> facet_vertex;
	std::vector<VERTEX_INDEX> cube_list;
	std::vector<ISO_VERTEX_INDEX> quad_vert2;
	std::vector<IJKDUALTABLE::TABLE_INDEX> table_index;

	extract_dual_quad_isopatch_incident_on
		(scalar_grid, isovalue, isovert, cube_index, gcube_map, dist2cube,
		isoquad_cube, facet_vertex);

	map_cube_list(isovert, gcube_map, isoquad_cube);

	construct_cube_list(isoquad_cube, cube_list);

	NUM_TYPE num_split;

	set_table_index(isovert, cube_list, table_index);

	IJK::split_dual_isovert
		(isodual_table, cube_list, table_index,
		isoquad_cube, facet_vertex, iso_vlist, quad_vert2, num_split);

	map_isov_indices(isovert, gcube_map, iso_vlist, quad_vert2);

	tri_vert.clear();
	quad_vert.clear();
	IJK::get_non_degenerate_quad_btlr(quad_vert2, tri_vert, quad_vert);
}

// Extract dual isosurface patch with vertex in merged cube.
// Allow multiple isosurface vertices in each cube.
void MERGESHARP::extract_dual_isopatch_incident_on_multi
	(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
	const SCALAR_TYPE isovalue,
	const ISOVERT & isovert,
	const VERTEX_INDEX cube_index,
	const std::vector<VERTEX_INDEX> & gcube_map,
	const AXIS_SIZE_TYPE dist2cube,
	std::vector<ISO_VERTEX_INDEX> & tri_vert,
	std::vector<ISO_VERTEX_INDEX> & quad_vert)
{
	std::vector<DUAL_ISOVERT> iso_vlist;

	extract_dual_isopatch_incident_on_multi
		(scalar_grid, isodual_table, isovalue, isovert, cube_index, gcube_map,
		dist2cube, tri_vert, quad_vert, iso_vlist);
}

// Insert polygon edges in edge hash table.
void insert_poly_edges
	(const std::vector<ISO_VERTEX_INDEX> & poly_vert, 
	const NUM_TYPE num_vert_per_poly,
	EDGE_HASH_TABLE & edge_hash)
{
	NUM_TYPE num_poly = poly_vert.size()/num_vert_per_poly;

	for (NUM_TYPE i = 0; i < num_poly; i++) {
		for (int k0 = 0; k0 < num_vert_per_poly; k0++) {
			VERTEX_INDEX iv0 = poly_vert[i*num_vert_per_poly+k0];
			int k1 = (k0+1)%num_vert_per_poly;
			VERTEX_INDEX iv1 = poly_vert[i*num_vert_per_poly+k1];
			if (iv0 > iv1) { std::swap(iv0, iv1); }
			VERTEX_PAIR key = std::make_pair(iv0, iv1);

			EDGE_HASH_TABLE::iterator edge_iter = edge_hash.find(key);
			if (edge_iter == edge_hash.end()) {
				edge_hash.insert(EDGE_HASH_TABLE::value_type(key, 1));
			}
			else {
				edge_iter->second++;
			}
		}
	}
}

// Insert triangle and quadrilateral edges into edge hash table.
void MERGESHARP::insert_tri_quad_edges
	(const std::vector<ISO_VERTEX_INDEX> & tri_vert,
	const std::vector<ISO_VERTEX_INDEX> & quad_vert,
	EDGE_HASH_TABLE & edge_hash)
{
	insert_poly_edges(tri_vert, NUM_VERT_PER_TRI, edge_hash);
	insert_poly_edges(quad_vert, NUM_VERT_PER_QUAD, edge_hash);
}

namespace {

	// Insert or increment vertex count in vertex hash table.
	void insert_vertex(const VERTEX_INDEX iv,
		VERTEX_HASH_TABLE & vertex_hash)
	{
		VERTEX_HASH_TABLE::iterator vertex_iter = vertex_hash.find(iv);
		if (vertex_iter == vertex_hash.end()) {
			vertex_hash.insert(std::pair<VERTEX_INDEX, int>(iv, 1));
		}
		else {
			vertex_iter->second++;
		}
	}

	// Insert vertices from vlist into vertex hash table.
	// @param vertex_hash Maps vertex to unique number 
	//        in range [0-(vertex_hash.size()-1)].
	void insert_vertex_list
		(const std::vector<VERTEX_INDEX> & vlist,
		VERTEX_HASH_TABLE & vertex_hash)
	{
		for (NUM_TYPE i = 0; i < vlist.size(); i++) {
			VERTEX_INDEX iv = vlist[i];
			VERTEX_HASH_TABLE::iterator vertex_iter = vertex_hash.find(iv);
			if (vertex_iter == vertex_hash.end()) {
				NUM_TYPE n = vertex_hash.size();
				vertex_hash.insert(std::pair<VERTEX_INDEX, int >(iv, n));
			}
		}
	}

	// Remap vertices from vlist to values vertex hash table.
	// @pre Every vertex in vlist is in the vertex hash table.
	void remap_vertex_list
		(const VERTEX_HASH_TABLE & vertex_hash,
		const std::vector<VERTEX_INDEX> & vlist,
		std::vector<VERTEX_INDEX> & new_vlist)
	{
		new_vlist.resize(vlist.size());

		for (NUM_TYPE i = 0; i < vlist.size(); i++) {
			VERTEX_INDEX iv = vlist[i];
			VERTEX_HASH_TABLE::const_iterator vertex_iter = vertex_hash.find(iv);
			new_vlist[i] = vertex_iter->second;
		}
	}

	// Insert vertex in cycle.
	// @param vertex_loc Location of vertex in list cycle_vertex.
	void insert_cycle_vertex
		(const VERTEX_INDEX iv0, const VERTEX_INDEX iv1,
		std::vector<CYCLE_VERTEX> & cycle_vertex)
	{
		NUM_TYPE num_adjacent = cycle_vertex[iv0].num_adjacent;
		if (num_adjacent < 2) 
		{ cycle_vertex[iv0].adjacent[num_adjacent] = iv1; }
		cycle_vertex[iv0].num_adjacent++;
	}
}

// Construct list of boundary cycle vertices.
void MERGESHARP::construct_boundary_cycle
	(const EDGE_HASH_TABLE & edge_hash,
	std::vector<CYCLE_VERTEX> & cycle_vertex)
{
	for (EDGE_HASH_TABLE::const_iterator edge_iter = edge_hash.begin();
		edge_iter != edge_hash.end(); edge_iter++) {
			if (edge_iter->second == 1) {
				VERTEX_INDEX iv0 = (edge_iter->first).first;
				VERTEX_INDEX iv1 = (edge_iter->first).second;
				insert_cycle_vertex(iv0, iv1, cycle_vertex);
				insert_cycle_vertex(iv1, iv0, cycle_vertex);
			}
	}
}

// Renumber tri and quad vertices so that they are from 0 to num_vert-1.
template <typename VTYPE>
void renumber_tri_quad_vertices
	(const std::vector<VTYPE> & tri_vert,
	const std::vector<VTYPE> & quad_vert,
	std::vector<VTYPE> & new_tri_vert,
	std::vector<VTYPE> & new_quad_vert,
	NUM_TYPE & num_vert)
{
	VERTEX_HASH_TABLE vertex_hash;

	insert_vertex_list(tri_vert, vertex_hash);
	insert_vertex_list(quad_vert, vertex_hash);

	num_vert = vertex_hash.size();
	remap_vertex_list(vertex_hash, tri_vert, new_tri_vert);
	remap_vertex_list(vertex_hash, quad_vert, new_quad_vert);
}

// Search cycle starting at iv0.
// @pre All vertices in cycle containing iv0 have is_visited set to false.
// @pre All vertices have degree two.
void search_cycle
	(const VERTEX_INDEX iv0, std::vector<CYCLE_VERTEX> & cycle_vertex)
{
	VERTEX_INDEX iv = iv0;
	VERTEX_INDEX ivprev = cycle_vertex[iv0].adjacent[0];
	while (!cycle_vertex[iv].is_visited) {
		cycle_vertex[iv].is_visited = true;
		if (cycle_vertex[iv].adjacent[0] == ivprev) {
			ivprev = iv;
			iv = cycle_vertex[iv].adjacent[1];
		}
		else {
			ivprev = iv;
			iv = cycle_vertex[iv].adjacent[0];
		}
	}
}

/// Test merge_sharp_iso_vertices

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012 Rephael Wenger

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


#include "ijkgrid_macros.h"

#include "isodual3D_isovert.h"
#include "isodual3D_decimate.h"

using namespace IJK;
using namespace ISODUAL3D;
using namespace std;

// global variables
bool flag_output_entire_grid(false);

// forward declarations
void set_isovert(const SHARPISO_GRID & grid, ISOVERT & isovert);
void init_gcube_map(const ISOVERT & isovert,
                    std::vector<VERTEX_INDEX> & gcube_map);
void set_scalarA
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const NUM_TYPE num_zero);
void set_scalarE
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_gcube_mapE
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & gcube_map);
void map_neighbors_to_cube
(const ISOVERT & isovert, const VERTEX_INDEX cube_index, 
 std::vector<VERTEX_INDEX> & gcube_map);
void output_scalar_grid(const SHARPISO_SCALAR_GRID_BASE & scalar_grid);
void output_scalar_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index);
void output_gcube_map
(const ISOVERT & isovert, const std::vector<VERTEX_INDEX> & gcube_map);
void output_gcube_map
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 const std::vector<VERTEX_INDEX> & gcube_map);
void output_is_isopatch_disk
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map);


int main()
{
  SHARPISO_SCALAR_GRID scalar_grid;
  ISOVERT isovert;
  SHARP_ISOVERT_PARAM isovert_param;
  AXIS_SIZE_TYPE axis_sizeA[DIM3] = { 5, 6, 7 };
  GRID_COORD_TYPE cubeA_coord[DIM3] = { 2, 3, 4 };
  VERTEX_INDEX cube_index;
  std::vector<VERTEX_INDEX> gcube_map;
  const SCALAR_TYPE isovalue = 1;

  scalar_grid.SetSize(DIM3, axis_sizeA);

  set_isovert(scalar_grid, isovert);
  init_gcube_map(isovert, gcube_map);
  cube_index = scalar_grid.ComputeVertexIndex(cubeA_coord);
  map_neighbors_to_cube(isovert, cube_index, gcube_map);

  if (flag_output_entire_grid) {
    output_gcube_map(isovert, gcube_map);
    cout << endl;
  }
  else {
    /* DEBUG
    output_gcube_map(isovert, cube_index, gcube_map);
    cout << endl;
    */

    set_scalarE(scalar_grid, cube_index);
    set_gcube_mapE(isovert, cube_index, gcube_map);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_gcube_map(isovert, cube_index, gcube_map);
    output_is_isopatch_disk
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);
  }
}

void set_isovert(const SHARPISO_GRID & grid, ISOVERT & isovert)
{
  COORD_TYPE coord[DIM3];

  isovert.sharp_ind_grid.SetSize(grid);
  isovert.sharp_ind_grid.SetAll(ISOVERT::NO_INDEX);

  IJK_FOR_EACH_GRID_CUBE(icube, grid, VERTEX_INDEX) {
    GRID_CUBE grid_cube;
    NUM_TYPE i = isovert.gcube_list.size();

    grid_cube.cube_index = icube;
    grid.ComputeCubeCenterCoord(icube, grid_cube.isovert_coord);
    grid.ComputeBoundaryCubeBits(icube, grid_cube.boundary_bits);
    grid_cube.flag = SMOOTH_GCUBE;

    isovert.gcube_list.push_back(grid_cube);
    isovert.sharp_ind_grid.Set(icube, i);
  }
}

void init_gcube_map(const ISOVERT & isovert,
                    std::vector<VERTEX_INDEX> & gcube_map)
{
  gcube_map.resize(isovert.gcube_list.size());
  
  for (int i = 0; i < gcube_map.size(); i++) 
    { gcube_map[i] = i; }
}

void map_neighbors_to_cube
(const ISOVERT & isovert,
 const VERTEX_INDEX cube_index, std::vector<VERTEX_INDEX> & gcube_map)
{
  const int dimension = isovert.sharp_ind_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = isovert.sharp_ind_grid.AxisSize();
  SHARPISO_GRID_NEIGHBORS gridn;
  long boundary_bits;
  IJK::PROCEDURE_ERROR error("map_neighbors_to_cube");

  gridn.SetSize(dimension, axis_size);

  NUM_TYPE gcube_index = isovert.sharp_ind_grid.Scalar(cube_index);
  boundary_bits = isovert.gcube_list[gcube_index].boundary_bits;

  if (boundary_bits != 0) {
    error.AddMessage
      ("Error. Map neighbors to cube not implemented for boundary cubes.");
    error.AddMessage
      ("Cube ", cube_index, " is a boundary cube.");
    throw error;
  }

  for (int i = 0; i < gridn.NumVertexNeighborsC(); i++) {
    VERTEX_INDEX cube_neighbor = gridn.VertexNeighborC(cube_index, i);
    NUM_TYPE gcube_neighbor = isovert.sharp_ind_grid.Scalar(cube_neighbor);
    gcube_map[gcube_neighbor] = gcube_index;
  }

}

// Set first num_zero vertices in region to 0.
void set_scalarA
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const NUM_TYPE num_zero)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const int REGION_EDGE_LENGTH = 3;

  // Initialize all scalar values to 0.
  scalar_grid.SetAll(0);

  VERTEX_INDEX iv0 = 
    cube_index - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);
  VERTEX_INDEX num_vertices;
  compute_num_grid_vertices_in_region
    (dimension, REGION_EDGE_LENGTH, num_vertices);

  IJK::ARRAY<VERTEX_INDEX> vlist(num_vertices);
  get_grid_vertices_in_region
    (dimension, axis_size, iv0, REGION_EDGE_LENGTH, vlist.Ptr());

  for (NUM_TYPE i = 0; i < num_vertices; i++) {
    if (i >= num_zero) 
      { scalar_grid.Set(vlist[i], 2); }
  }

}

// Set first num_zero primary cube vertices of cubes in region to 0.
void set_scalarB
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const NUM_TYPE num_zero)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const int REGION_EDGE_LENGTH = 3;

  // Initialize all scalar values to 2.
  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = 
    cube_index - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);
  VERTEX_INDEX num_cubes;
  compute_num_grid_cubes_in_region
    (dimension, REGION_EDGE_LENGTH, num_cubes);

  IJK::ARRAY<VERTEX_INDEX> vlist(num_cubes);
  get_grid_cubes_in_region
    (dimension, axis_size, iv0, REGION_EDGE_LENGTH, vlist.Ptr(), num_cubes);

  for (NUM_TYPE i = 0; i < num_cubes; i++) {
    if (i < num_zero) 
      { scalar_grid.Set(vlist[i], 0); }
  }

}


// Set all vertices of cube cube_index to 2.
void set_scalarC
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  // Initialize all scalar values to 0.
  scalar_grid.SetAll(0);

  for (int k = 0; k < NUM_CUBE_VERTICES3D; k++) {
    VERTEX_INDEX iv = scalar_grid.CubeVertex(cube_index, k);
    scalar_grid.Set(iv, 2);
  }
}

// Set all vertices in box to 2.
void set_scalarD
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index,
 const IJK::BOX<GRID_COORD_TYPE> & box)
{
  // Initialize all scalar values to 0.
  scalar_grid.SetAll(0);
  scalar_grid.SetRegion(box, 1, 2);
}


template <typename GTYPE, typename STYPE>
void set_slice
(GTYPE & scalar_grid, const VERTEX_INDEX iv0,
 STYPE slice[], AXIS_SIZE_TYPE slice_axis_size[2])
{
  GRID_COORD_TYPE coord0[DIM3];
  GRID_COORD_TYPE coord[DIM3];

  scalar_grid.ComputeCoord(iv0, coord0);

  coord[2] = coord0[2];
  for (int y = 0; y < slice_axis_size[1]; y++) {
    coord[1] = coord0[1]+y;
    for (int x = 0; x < slice_axis_size[0]; x++) {
      coord[0] = coord0[0]+x;
      VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(coord);
      scalar_grid.Set(iv, slice[y*slice_axis_size[0] + x]);
    }
  }
}

void set_scalarE
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 0, 2 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 0, 0, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      0, 0, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_gcube_mapE
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
  std::vector<VERTEX_INDEX> & gcube_map)
{
  SHARPISO_INDEX_GRID gcube_map_grid;
  NUM_TYPE gcube_index = isovert.sharp_ind_grid.Scalar(cube_index);
  AXIS_SIZE_TYPE slice_axis_size[2] = { 3, 3 };
  SCALAR_TYPE slice1[9] =
    { 0, 0, 0,
      0, 1, 1,
      1, 1, 1 };
  SCALAR_TYPE slice2[9] =
    { 0, 0, 0,
      1, 1, 1,
      1, 1, 1 };
  SCALAR_TYPE slice3[9] =
    { 0, 0, 0,
      1, 1, 0,
      1, 1, 1 };

  gcube_map_grid.SetSize(isovert.sharp_ind_grid);

  for (int i = 0; i < gcube_map.size(); i++)
    { gcube_map[i] = i; };

  gcube_map_grid.SetAll(0);

  VERTEX_INDEX iv0 = gcube_map_grid.PrevVertex(cube_index, 0);
  iv0 = gcube_map_grid.PrevVertex(iv0, 1);
  iv0 = gcube_map_grid.PrevVertex(iv0, 2);
  set_slice(gcube_map_grid, iv0, slice1, slice_axis_size);

  iv0 = gcube_map_grid.NextVertex(iv0, 2);
  set_slice(gcube_map_grid, iv0, slice2, slice_axis_size);

  iv0 = gcube_map_grid.NextVertex(iv0, 2);
  set_slice(gcube_map_grid, iv0, slice3, slice_axis_size);

  for (int i = 0; i < isovert.gcube_list.size(); i++) {
    VERTEX_INDEX iv = isovert.gcube_list[i].cube_index;
    if (gcube_map_grid.Scalar(iv) == 1) {
      gcube_map[i] = gcube_index;
    }
  }

}


// **************************************************
// Check routines
// **************************************************


// **************************************************
// Output routines
// **************************************************

void output_scalar_grid(const SHARPISO_SCALAR_GRID_BASE & scalar_grid)
{
  cout << "Scalar grid:" << endl;
  IJK::output_scalar_grid(cout, scalar_grid);
}


void output_scalar_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE axis_size[DIM3] = { 4, 4, 4 };
  SHARPISO_SCALAR_GRID scalar_subgrid(dimension, axis_size);

  VERTEX_INDEX iv0 = 
    cube_index - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

  scalar_subgrid.CopyRegion(scalar_grid, iv0, axis_size, 0);

  cout << "Scalar values around cube " << cube_index << ":" << endl;
  IJK::output_scalar_grid(cout, scalar_subgrid);
}

void output_gcube_map
(const ISOVERT & isovert, const std::vector<VERTEX_INDEX> & gcube_map)
{
  const int dimension = isovert.sharp_ind_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = isovert.sharp_ind_grid.AxisSize();
  SHARPISO_INDEX_GRID gcube_map_grid;

  gcube_map_grid.SetSize(dimension, axis_size);

  for (VERTEX_INDEX iv = 0; iv < isovert.sharp_ind_grid.NumVertices(); iv++) {
    INDEX_DIFF_TYPE gcube_index = isovert.sharp_ind_grid.Scalar(iv);
    if (0 <= gcube_index && gcube_index < gcube_map.size()) 
      { gcube_map_grid.Set(iv, gcube_map[gcube_index]); }
    else
      { gcube_map_grid.Set(iv, gcube_index); }
  }

  cout << "gcube map: " << endl;
  output_scalar_grid(cout, gcube_map_grid);
}

// Output gcube map around cube.
void output_gcube_map
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 const std::vector<VERTEX_INDEX> & gcube_map)
{
  const int dimension = isovert.sharp_ind_grid.Dimension();
  const AXIS_SIZE_TYPE axis_size[DIM3] = { 3, 3, 3 };
  SHARPISO_INDEX_GRID gcube_map_grid(dimension, axis_size);

  VERTEX_INDEX iv0 = 
    cube_index - 
    isovert.sharp_ind_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

  gcube_map_grid.CopyRegion(isovert.sharp_ind_grid, iv0, axis_size, 0);

  for (VERTEX_INDEX iv = 0; iv < gcube_map_grid.NumVertices(); iv++) {
    INDEX_DIFF_TYPE gcube_index = gcube_map_grid.Scalar(iv);
    if (0 <= gcube_index && gcube_index < gcube_map.size())
      { gcube_map_grid.Set(iv, gcube_map[gcube_index]); }
  }

  cout << "gcube map around cube " << cube_index << ": " << endl;
  output_scalar_grid(cout, gcube_map_grid);
}

void output_is_isopatch_disk
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map)
{

}

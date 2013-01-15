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


#include <iomanip>

#include "ijkgrid_macros.h"
#include "ijkmesh.txx"

#include "ijkdualtable.h"

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
void set_scalarF1
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_scalarF2
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_scalarF3
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_scalarF4
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_scalarF5
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_scalarF6
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index);
void set_gcube_map
(const ISOVERT & isovert, const VERTEX_INDEX cube_index0,
 const bool A3x3x3[DIM3][DIM3][DIM3],
 std::vector<VERTEX_INDEX> & gcube_map);
void set_gcube_mapE
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & gcube_map);
void set_gcube_mapA1
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & gcube_map);
void set_gcube_mapA2
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & gcube_map);
void set_gcube_mapF6
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 std::vector<VERTEX_INDEX> & gcube_map);
void map_neighbors_to_cube
(const ISOVERT & isovert, const VERTEX_INDEX cube_index, 
 std::vector<VERTEX_INDEX> & gcube_map);
template <typename SGRID>
void output_scalar_grid(const SGRID & scalar_grid);
template <typename SGRID, typename VTYPE>
void output_scalar_subgrid(const SGRID & scalar_grid, const VTYPE cube_index);
void output_index_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index);
void output_gcube_map
(const ISOVERT & isovert, const std::vector<VERTEX_INDEX> & gcube_map);
void output_gcube_map
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 const std::vector<VERTEX_INDEX> & gcube_map);
void output_merged_cubes
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 const std::vector<VERTEX_INDEX> & gcube_map);
void output_is_isopatch_disk2
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map);
void check_scalarA
(const VERTEX_INDEX cube_index,
 const ISOVERT & isovert,
 std::vector<VERTEX_INDEX> & gcube_map);


int main()
{
  SHARPISO_SCALAR_GRID scalar_grid;
  ISOVERT isovert;
  SHARP_ISOVERT_PARAM isovert_param;
  AXIS_SIZE_TYPE axis_sizeA[DIM3] = { 6, 6, 6};
  GRID_COORD_TYPE cubeA_coord[DIM3] = { 2, 2, 2 };
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

    /*
    set_gcube_mapA2(isovert, cube_index, gcube_map);
    output_gcube_map(isovert, cube_index, gcube_map);
    output_index_subgrid(scalar_grid, cube_index);
    output_merged_cubes(isovert, cube_index, gcube_map);

    set_scalarF1(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);

    set_scalarF2(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);

    set_scalarF3(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);

    set_scalarF4(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);

    set_scalarF5(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
      (scalar_grid, isovalue, cube_index, isovert, gcube_map);
    */

    set_gcube_mapF6(isovert, cube_index, gcube_map);
    output_gcube_map(isovert, cube_index, gcube_map);
    output_index_subgrid(scalar_grid, cube_index);
    output_merged_cubes(isovert, cube_index, gcube_map);

    set_scalarF6(scalar_grid, cube_index);
    output_scalar_subgrid(scalar_grid, cube_index);
    output_is_isopatch_disk2
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

void set_scalarF1
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 0, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);
}

void set_scalarF2
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 2, 2,
      2, 0, 2, 2,
      2, 0, 2, 2,
      2, 2, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_scalarF3
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 2,
      0, 2, 2, 2,
      0, 0, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 2, 2,
      0, 0, 2, 2,
      0, 2, 2, 2,
      2, 2, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_scalarF4
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 0, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 2, 2,
      2, 0, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_scalarF5
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice2[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 2,
      2, 0, 2, 2,
      2, 2, 0, 2,
      2, 2, 2, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);
  set_slice(scalar_grid, iv0, slice2, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_scalarF6
(SHARPISO_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX cube_index)
{
  AXIS_SIZE_TYPE slice_axis_size[2] = { 4, 4 };
  SCALAR_TYPE slice3[16] =
    { 2, 2, 2, 0,
      2, 2, 0, 0,
      2, 0, 0, 2,
      0, 0, 2, 2 };
  SCALAR_TYPE slice4[16] =
    { 2, 2, 0, 0,
      2, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0 };

  scalar_grid.SetAll(2);

  VERTEX_INDEX iv0 = scalar_grid.PrevVertex(cube_index, 0);
  iv0 = scalar_grid.PrevVertex(iv0, 1);


  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice3, slice_axis_size);

  iv0 = scalar_grid.NextVertex(iv0, 2);
  set_slice(scalar_grid, iv0, slice4, slice_axis_size);
}

void set_gcube_map
(const ISOVERT & isovert, const VERTEX_INDEX cube_index0,
 const bool A3x3x3[DIM3][DIM3][DIM3],
 std::vector<VERTEX_INDEX> & gcube_map)
{
  const NUM_TYPE gcube_index0 = isovert.sharp_ind_grid.Scalar(cube_index0);
  GRID_COORD_TYPE cube_coord0[DIM3];
  GRID_COORD_TYPE cube_coord1[DIM3];

  isovert.sharp_ind_grid.ComputeCoord(cube_index0, cube_coord0);

  for (int i = 0; i < gcube_map.size(); i++)
    { gcube_map[i] = i; };

  for (int z = 0; z < DIM3; z++) {
    for (int y = 0; y < DIM3; y++) {
      for (int x = 0; x < DIM3; x++) {
        cube_coord1[0] = cube_coord0[0]+x-1;
        cube_coord1[1] = cube_coord0[1]+y-1;
        cube_coord1[2] = cube_coord0[2]+z-1;
        VERTEX_INDEX cube_index1 = 
          isovert.sharp_ind_grid.ComputeVertexIndex(cube_coord1);
        VERTEX_INDEX gcube_index1 =
          isovert.sharp_ind_grid.Scalar(cube_index1);
        if (A3x3x3[z][y][x]) {
          gcube_map[gcube_index1] = gcube_index0;
        }
      }
    }
  }
}

void clear_A3x3x3(bool A3x3x3[DIM3][DIM3][DIM3])
{
  for (int x = 0; x < DIM3; x++)
    for (int y = 0; y < DIM3; y++)
      for (int z = 0; z < DIM3; z++)
        { A3x3x3[x][y][z] = false; }
}

void set_slice(const int z, const bool slice[DIM3][DIM3],
               bool A3x3x3[DIM3][DIM3][DIM3])
{
  for (int x = 0; x < DIM3; x++)
    for (int y = 0; y < DIM3; y++)
      { A3x3x3[z][y][x] = slice[y][x]; }
}

void set_gcube_mapA1
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
  std::vector<VERTEX_INDEX> & gcube_map)
{
  bool A3x3x3[DIM3][DIM3][DIM3];

  clear_A3x3x3(A3x3x3);
  A3x3x3[1][1][1] = 1;
  set_gcube_map(isovert, cube_index, A3x3x3, gcube_map);
}

void set_gcube_mapA2
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
  std::vector<VERTEX_INDEX> & gcube_map)
{
  bool A3x3x3[DIM3][DIM3][DIM3];

  clear_A3x3x3(A3x3x3);
  A3x3x3[1][1][1] = 1;
  A3x3x3[2][1][1] = 1;
  set_gcube_map(isovert, cube_index, A3x3x3, gcube_map);
}

void set_gcube_mapE
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
  std::vector<VERTEX_INDEX> & gcube_map)
{
  bool A3x3x3[DIM3][DIM3][DIM3];
  bool slice0[DIM3][DIM3] = { { 0, 0, 0 },
                              { 0, 1, 1 },
                              { 1, 1, 1 } };
  bool slice1[DIM3][DIM3] = { { 0, 0, 0 },
                              { 1, 1, 1 },
                              { 1, 1, 1 } };
  bool slice2[DIM3][DIM3] = { { 0, 0, 0 },
                              { 1, 1, 0 },
                              { 1, 1, 1 } };

  clear_A3x3x3(A3x3x3);
  set_slice(0, slice0, A3x3x3);
  set_slice(1, slice1, A3x3x3);
  set_slice(2, slice2, A3x3x3);
  set_gcube_map(isovert, cube_index, A3x3x3, gcube_map);
}

void set_gcube_mapF6
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
  std::vector<VERTEX_INDEX> & gcube_map)
{
  bool A3x3x3[DIM3][DIM3][DIM3];
  bool slice0[DIM3][DIM3] = { { 0, 0, 0 },
                              { 0, 0, 0 },
                              { 0, 0, 0 } };
  bool slice1[DIM3][DIM3] = { { 0, 1, 0 },
                              { 1, 1, 1 },
                              { 0, 1, 1 } };
  bool slice2[DIM3][DIM3] = { { 1, 1, 0 },
                              { 1, 1, 1 },
                              { 0, 1, 1 } };

  clear_A3x3x3(A3x3x3);
  set_slice(0, slice0, A3x3x3);
  set_slice(1, slice1, A3x3x3);
  set_slice(2, slice2, A3x3x3);
  set_gcube_map(isovert, cube_index, A3x3x3, gcube_map);

  VERTEX_INDEX cube_index2 = isovert.sharp_ind_grid.PrevVertex(cube_index, 0);
  cube_index2 = isovert.sharp_ind_grid.NextVertex(cube_index2, 1);
  VERTEX_INDEX cube_index3 = isovert.sharp_ind_grid.NextVertex(cube_index2, 2);
  VERTEX_INDEX gcube_index2 = isovert.sharp_ind_grid.Scalar(cube_index2);
  VERTEX_INDEX gcube_index3 = isovert.sharp_ind_grid.Scalar(cube_index3);
  gcube_map[gcube_index3] = gcube_index2;

  cube_index2 = isovert.sharp_ind_grid.NextVertex(cube_index, 0);
  cube_index2 = isovert.sharp_ind_grid.PrevVertex(cube_index2, 1);
  cube_index3 = isovert.sharp_ind_grid.NextVertex(cube_index2, 2);
  gcube_index2 = isovert.sharp_ind_grid.Scalar(cube_index2);
  gcube_index3 = isovert.sharp_ind_grid.Scalar(cube_index3);
  gcube_map[gcube_index3] = gcube_index2;
}

// **************************************************
// Check routines
// **************************************************

void check_scalarA
(const VERTEX_INDEX cube_index,
 const ISOVERT & isovert,
 std::vector<VERTEX_INDEX> & gcube_map)
{
  const int dimension = isovert.sharp_ind_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = isovert.sharp_ind_grid.AxisSize();
  const int REGION_EDGE_LENGTH = 3;
  const int isovalue = 1;
  SHARPISO_SCALAR_GRID scalar_grid(dimension, axis_size);
  VERTEX_INDEX num_vertices;
  int num_non_disk(0);

  compute_num_grid_vertices_in_region
    (dimension, REGION_EDGE_LENGTH, num_vertices);

  for (int i = 0; i <= num_vertices; i++) {
    set_scalarA(scalar_grid, cube_index, i);

    // *** CHECK IF ISOPATCH IS A DISK ***
    /*
    {
      cout << "Isopatch formed by cube " << cube_index << " is not a disk." 
           << endl;

      output_scalar_subgrid(scalar_grid, cube_index);
      num_non_disk++;
    }
    */
  }

  if (num_non_disk == 0) {
    cout << "Passed check_scalarA." << endl;
  }
}


// **************************************************
// Output routines
// **************************************************

template <typename SGRID>
void output_scalar_grid(const SGRID & scalar_grid)
{
  cout << "Scalar grid:" << endl;
  IJK::output_scalar_grid(cout, scalar_grid);
}

template <typename SGRID>
void output_grid_horizontal
(const int num_width, const SGRID & scalar_grid)
{
  GRID_COORD_TYPE coord[DIM3];
  IJK::PROCEDURE_ERROR error("output_grid_horizontal");

  if (scalar_grid.Dimension() != 3) {
    error.AddMessage("Programming error.  Dimension should be 3.");
    throw error;
  }

  for (int y = 0; y < scalar_grid.AxisSize(1); y++) {
    for (int z = 0; z < scalar_grid.AxisSize(2); z++) {
      for (int x = 0; x < scalar_grid.AxisSize(0); x++) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
        VERTEX_INDEX iv = scalar_grid.ComputeVertexIndex(coord);
        cout << " " << setw(num_width) << scalar_grid.Scalar(iv);
      }
      cout << "  ";
    }
    cout << endl;
  }
}


template <typename SGRID, typename VTYPE>
void output_subgrid
(const char * label, const int num_width, 
 const SGRID & scalar_grid, const VTYPE cube_index)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE axis_size[DIM3] = { 4, 4, 4 };
  SGRID scalar_subgrid(dimension, axis_size);

  VERTEX_INDEX iv0 = 
    cube_index - scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

  scalar_subgrid.CopyRegion(scalar_grid, iv0, axis_size, 0);

  cout << label << " around cube " << cube_index << ":" << endl;

  output_grid_horizontal(num_width, scalar_subgrid);
  cout << endl;
}

template <typename SGRID, typename VTYPE>
void output_scalar_subgrid
(const SGRID & scalar_grid, const VTYPE cube_index)
{
  output_subgrid("Scalar values", 1, scalar_grid, cube_index);
}

void output_index_subgrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const VERTEX_INDEX cube_index)
{
  SHARPISO_INDEX_GRID index_grid;

  index_grid.SetSize(scalar_grid);

  for (VERTEX_INDEX iv = 0; iv < index_grid.NumVertices(); iv++)
    { index_grid.Set(iv, iv); }

  output_subgrid("Vertex indices", 3, index_grid, cube_index);
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
  output_grid_horizontal(3, gcube_map_grid);
  cout << endl;
}


// Output cubes merged with cube_index
void output_merged_cubes
(const ISOVERT & isovert, const VERTEX_INDEX cube_index,
 const std::vector<VERTEX_INDEX> & gcube_map)
{
  std::vector<VERTEX_INDEX> merged_cube_list;
  std::vector<EDGE_INDEX> boundary_edge_list;
  const int dist2cube = 1;

  get_merged_cubes
    (isovert.sharp_ind_grid, isovert, cube_index, gcube_map, dist2cube,
     merged_cube_list);

  cout << "Cubes merged with " << cube_index << ":";
  for (int i = 0; i < merged_cube_list.size(); i++) 
    { cout << " " << merged_cube_list[i]; }
  cout << endl;

  get_merged_boundary_edges
    (isovert.sharp_ind_grid, merged_cube_list, boundary_edge_list);

  /* DEBUG
  cout << "Merged boundary edges:" << endl;
  for (int i = 0; i < boundary_edge_list.size(); i++) {
    EDGE_INDEX edge_index = boundary_edge_list[i];
    VERTEX_INDEX iv0 = edge_index / DIM3;
    int edge_dir = edge_index%DIM3;
    VERTEX_INDEX iv1 = isovert.sharp_ind_grid.NextVertex(iv0, edge_dir);
    cout << "  ( " << iv0 << "," << iv1 << ") [Dir " << edge_dir << "]";
    cout << endl;
  }
  cout << endl;
  */

}


VERTEX_INDEX get_cube_index
(const ISOVERT & isovert, const NUM_TYPE gcube_index)
{
  return(isovert.gcube_list[gcube_index].cube_index);
}

template <typename T>
void output_poly_vert
(const char * label,
 const std::vector<T> & poly_vert, const NUM_TYPE num_vert_per_poly)
{
  NUM_TYPE num_poly = poly_vert.size()/num_vert_per_poly;
  cout << label << " (" << num_poly << "):" << endl;
  for (int i = 0; i < num_poly; i++) {
    cout << "  (";
    for (int j = 0; j < num_vert_per_poly; j++) {
      if (j > 0) { cout << " "; }
      cout <<  poly_vert[i*num_vert_per_poly+j];
    }
    cout << ")";
    if (i%2 == 1) { cout << endl; }
  }
  cout << endl;
}

void output_isovert(const ISODUAL3D::DUAL_ISOVERT & isov)
{
  cout << isov.cube_index << ":" << int(isov.patch_index);
}

template <typename T>
void output_poly_vert
(const char * label, const std::vector<ISODUAL3D::DUAL_ISOVERT> & iso_vlist,
 const std::vector<T> & poly_vert, const NUM_TYPE num_vert_per_poly)
{
  NUM_TYPE num_poly = poly_vert.size()/num_vert_per_poly;
  cout << label << " (" << num_poly << "):" << endl;
  for (int i = 0; i < num_poly; i++) {
    cout << "  (";
    for (int j = 0; j < num_vert_per_poly; j++) {
      if (j > 0) { cout << " "; }
      VERTEX_INDEX iv = poly_vert[i*num_vert_per_poly+j];
      output_isovert(iso_vlist[iv]);
    }
    cout << ")";
    if (i%2 == 1) { cout << endl; }
  }
  cout << endl;
}

// Renumber tri and quad vertices so that they are from 0 to num_vert-1.
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

// Insert vertices from vlist into vertex hash table.
// @param vertex_hash Maps vertex to unique number 
//        in range [0-(vertex_hash.size()-1)].
void insert_vertex_list
(const std::vector<VERTEX_INDEX> & vlist,
 VERTEX_HASH_TABLE & vertex_hash);

// Remap vertices from vlist to values vertex hash table.
// @pre Every vertex in vlist is in the vertex hash table.
void remap_vertex_list
(const VERTEX_HASH_TABLE & vertex_hash,
 const std::vector<VERTEX_INDEX> & vlist,
 std::vector<VERTEX_INDEX> & new_vlist);

// Construct cube list.
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

void output_is_isopatch_disk2
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX cube_index,
 const ISOVERT & isovert, 
 const std::vector<VERTEX_INDEX> & gcube_map)
{
  const int dimension = scalar_grid.Dimension();
  const int dist2cube = 1;
  std::vector<VERTEX_INDEX> isoquad_cube;
  std::vector<FACET_VERTEX_INDEX> facet_vertex;
  std::vector<ISO_VERTEX_INDEX> tri_vert;
  std::vector<ISO_VERTEX_INDEX> quad_vert;
  std::vector<ISO_VERTEX_INDEX> tri_vert2;
  std::vector<ISO_VERTEX_INDEX> quad_vert2;
  std::vector<ISO_VERTEX_INDEX> quad_vert3;
  std::vector<ISODUAL3D::DUAL_ISOVERT> iso_vlist;
  std::vector<VERTEX_INDEX> cube_list;

  bool flag_separate_opposite(true);
  bool flag_separate_neg(true);
  IJKDUALTABLE::ISODUAL_CUBE_TABLE 
    isodual_table(dimension, flag_separate_neg, flag_separate_opposite);

  extract_dual_quad_isopatch_incident_on
    (scalar_grid, isovalue, isovert, cube_index, gcube_map, dist2cube,
     isoquad_cube, facet_vertex);

  construct_cube_list(isoquad_cube, cube_list);

  NUM_TYPE num_split;

  // *** DEBUG ***
  using namespace std;
  cerr << "cube_list (" << cube_list.size() << "):" << endl;
  for (int j = 0; j < cube_list.size(); j++)
    { cout << " " << cube_list[j]; }
  cout << endl;

  /* DEBUG
  output_poly_vert("Quad cubes ", isoquad_cube, NUM_VERT_PER_QUAD);
  */

  IJK::split_dual_isovert
    (scalar_grid, isodual_table, isovalue, cube_list,
     isoquad_cube, facet_vertex, iso_vlist, quad_vert3, num_split);

  // *** DEBUG ***
  /*
  cerr << "iso_vlist: " << endl;
  for (int j = 0; j < iso_vlist.size(); j++) {
    cerr << "  cube: " << iso_vlist[j].cube_index
         << "  patch: " << int(iso_vlist[j].patch_index)
         << "  table: " << iso_vlist[j].table_index << endl;
  }
  */

  IJK::get_non_degenerate_quad_btlr(quad_vert3, tri_vert, quad_vert);
  reorder_quad_vertices(quad_vert);

  output_poly_vert
    ("Isosurface triangles ", iso_vlist, tri_vert, NUM_VERT_PER_TRI);
  output_poly_vert("Isosurface quadrilaterals ", iso_vlist,
                   quad_vert, NUM_VERT_PER_QUAD);

  if (is_isopatch_disk3D(tri_vert, quad_vert)) {
    cout << "Isopatch formed by cube " << cube_index << " is a disk." << endl;
  }
  else {
    EDGE_HASH_TABLE edge_hash;
    NUM_TYPE num_vert;

    // Renumber tri and quad vertices
    renumber_tri_quad_vertices
      (tri_vert, quad_vert, tri_vert2, quad_vert2, num_vert);

    std::vector<VERTEX_INDEX> original_vertex_id(num_vert);
    for (NUM_TYPE i = 0; i < tri_vert.size(); i++) 
      { original_vertex_id[tri_vert2[i]] = tri_vert[i]; }
    for (NUM_TYPE i = 0; i < quad_vert.size(); i++) 
      { original_vertex_id[quad_vert2[i]] = quad_vert[i]; }

    insert_tri_quad_edges(tri_vert2, quad_vert2, edge_hash);
    cout << "Isopatch formed by cube " << cube_index << " is not a disk." 
         << endl;

    bool flag_found_error(false);
    for (EDGE_HASH_TABLE::const_iterator edge_iter = edge_hash.begin();
         edge_iter != edge_hash.end(); edge_iter++) {
      if (edge_iter->second > 2) {
        VERTEX_INDEX i0 = original_vertex_id[(edge_iter->first).first];
        VERTEX_INDEX i1 = original_vertex_id[(edge_iter->first).second];
        cout << "  Edge: (";
        output_isovert(iso_vlist[i0]);
        cout << ",";
        output_isovert(iso_vlist[i1]);
        cout << ") is contained in "
             << edge_iter->second << " polygons." << endl;
        flag_found_error = true;
      }
    }

    if (flag_found_error) { 
      cout << endl;
      return; 
    };

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
        VERTEX_INDEX i0 = original_vertex_id[i];
        cout << "Vertex " << i0
             << " is incident on " << cycle_vertex[i].num_adjacent
             << " boundary edges." << endl;
        flag_found_error = true;
      }
    }

    if (flag_found_error) { 
      cout << endl;
      return; 
    };

    if (num_boundary_vertices < 3) { 
      // Disk must have at least three boundary cycle vertices.
      cout << "Patch has only " << num_boundary_vertices
           << " boundary vertices." << endl;
      cout << endl;
      return;
    }

    search_cycle(first_adjacent, cycle_vertex);

    int num_not_connected = 0;
    for (NUM_TYPE i = 0; i < cycle_vertex.size(); i++) { 
      if (cycle_vertex[i].num_adjacent == 2) {
        if (!cycle_vertex[i].is_visited) { 
          VERTEX_INDEX i0 = original_vertex_id[first_adjacent];
          VERTEX_INDEX i1 = original_vertex_id[i];
          if (num_not_connected == 0) {
            cout << "Vertices not in the same boundary cycles as ";
            output_isovert(iso_vlist[i0]);
            cout << " :" << endl;
          }
          cout << "  ";
          output_isovert(iso_vlist[i1]);

          num_not_connected++;
        }
      }
    }
    if (num_not_connected > 0) 
      { cout << endl; }
  }
  cout << endl;
}


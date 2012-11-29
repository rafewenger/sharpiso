/// \file isodual3D_decimate.cxx
/// Decimate dual isosurface.

/*
  IJK: Isosurface Jeneration Kode
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


#include <vector>

#include "ijkgraph.txx"
#include "ijkmesh.txx"
#include "ijkgrid_macros.h"

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"
#include "isodual3D_decimate.h"


// forward declarations
namespace {

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const ISODUAL3D::ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

}

// **************************************************
// Merge some isosurface vertices
// **************************************************

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::merge_sharp_iso_vertices
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const SHARP_ISOVERT_PARAM & sharp_isovert_param,
 std::vector<VERTEX_INDEX> & quad_vert,
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

  for (VERTEX_INDEX i = 0; i < quad_vert.size(); i++) {
    VERTEX_INDEX cube_index = quad_vert[i];
    quad_vert[i] = gcube_map[cube_index];
  }
}

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::merge_sharp_iso_vertices
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOVERT & isovert, 
 const SHARP_ISOVERT_PARAM & sharp_isovert_param,
 std::vector<VERTEX_INDEX> & quad_vert,
 SHARPISO_INFO & sharpiso_info)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();
  std::vector<VERTEX_INDEX> gcube_map(num_gcube);

  merge_sharp_iso_vertices
    (scalar_grid, isovalue, isovert, sharp_isovert_param,
     quad_vert, gcube_map, sharpiso_info);
}


// **************************************************
// Map isosurface vertices
// **************************************************

namespace {

  using namespace ISODUAL3D;

  void map_iso_vertex(const std::vector<GRID_CUBE> & gcube_list,
                      const INDEX_DIFF_TYPE from_cube,
                      const INDEX_DIFF_TYPE to_cube,
                      std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    if (from_cube != ISOVERT::NO_INDEX) {
      if (gcube_list[from_cube].flag != SELECTED_GCUBE && 
          gcube_map[from_cube] == from_cube) {
        if (gcube_list[from_cube].boundary_bits == 0) {
          // Map gcube_list[from_cube] to isosurface vertex in cube to_cube.
          gcube_map[from_cube] = to_cube;
        }
      }
    }
  }
    
  // Map isosurface vertices adjacent to selected cubes
  //  to the isosurface vertex in the selected cube.
  void map_adjacent_cubes
  (const ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    using namespace ISODUAL3D;

    const NUM_TYPE num_gcube = isovert.gcube_list.size();
    VERTEX_INDEX cube_index, neighbor_index;
    SHARPISO_GRID_NEIGHBORS gridn;
    GRID_COORD_TYPE cube_coord[DIM3];

    // Set size of grid neighbors grid.
    gridn.SetSize(isovert.sharp_ind_grid);

    for (NUM_TYPE i = 0; i < num_gcube; i++)
      { gcube_map[i] = i; }

    // Set cubes which share facets with selected cubes.
    for (NUM_TYPE i = 0; i < num_gcube; i++) {
      if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

        cube_index = isovert.gcube_list[i].cube_index;

        if (isovert.gcube_list[i].boundary_bits == 0) {
          // Cube cube_index is an interior cube.

          for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) {

            neighbor_index = gridn.CubeNeighborF(cube_index, j);

            INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
            map_iso_vertex(isovert.gcube_list, k, i, gcube_map);
          }
        }
        else {
          // *** Handle boundary case. ***
        }
      }
    }

    // Set cubes which share edges with selected cubes.
    for (NUM_TYPE i = 0; i < num_gcube; i++) {
      if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

        cube_index = isovert.gcube_list[i].cube_index;

        if (isovert.gcube_list[i].boundary_bits == 0) {
          // Cube cube_index is an interior cube.

          for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsE(); j++) {

            neighbor_index = gridn.CubeNeighborE(cube_index, j);

            INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
            map_iso_vertex(isovert.gcube_list, k, i, gcube_map);
          }
        }
        else {
          // *** Handle boundary case. ***
        }
      }
    }

    for (NUM_TYPE i = 0; i < num_gcube; i++) {
      if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

        cube_index = isovert.gcube_list[i].cube_index;

        if (isovert.gcube_list[i].boundary_bits == 0) {
          // Cube cube_index is an interior cube.

          for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsV(); j++) {

            neighbor_index = gridn.CubeNeighborV(cube_index, j);

            INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
            map_iso_vertex(isovert.gcube_list, k, i, gcube_map);
          }

        }
        else {
          // *** Fill in. ***
        }

      }
    }
  }

  // Forward declaration
  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const ISODUAL3D::ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
  {
    map_adjacent_cubes(isovert, gcube_map);

    if (sharp_isovert_param.flag_check_disk) {
      unmap_non_disk_isopatches
        (scalar_grid, isovalue, isovert, gcube_map, sharpiso_info);
    }
  }

}

// **************************************************
// Unmap non-disk isopatches
// **************************************************

namespace {

  // Reverse merges which create isopatches which are not disks.
  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   const ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
  {
    const NUM_TYPE num_gcube = isovert.gcube_list.size();

    sharpiso_info.num_non_disk_isopatches = 0;

    for (int d = 0; d < DIM3; d++) {
      if (scalar_grid.AxisSize(d) < 
          IS_ISOPATCH_DISK::num_vert_along_region_axis) 
        { return; }
    }

    IS_ISOPATCH_DISK is_isopatch_disk(scalar_grid);

    // Set cubes which share facets with selected cubes.
    for (NUM_TYPE i = 0; i < num_gcube; i++) {
      if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

        // *** Ignore boundary cubes for now. ***
        if (isovert.gcube_list[i].boundary_bits == 0) {

          VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;

          if (!is_isopatch_disk.IsIsopatchDisk
              (scalar_grid, isovalue, cube_index, isovert, gcube_map)) {
            is_isopatch_disk.UnmapAdjacent(cube_index, isovert, gcube_map); 

            sharpiso_info.num_non_disk_isopatches++;
          }
        }
      }
    }
  }

}

// **************************************************
// Function class IS_ISOPATCH_DISK
// **************************************************

namespace ISODUAL3D {

  /// Constructor for IS_ISOPATCH_DISK
  IS_ISOPATCH_DISK::IS_ISOPATCH_DISK(const SHARPISO_GRID & grid)
  {
    for (int d = 0; d < DIM3; d++) 
      { region_axis_size[d] = num_vert_along_region_axis; }

    cubeFlag.SetSize(DIM3, region_axis_size);
    selectedCubeBoundary.SetSize(DIM3, region_axis_size);
    regionScalar.SetSize(DIM3, region_axis_size);
    regionBoundary.SetSize(DIM3, region_axis_size);
    compute_boundary_grid(regionBoundary);

    region_boundary_cube.resize
      (selectedCubeBoundary.ComputeNumBoundaryCubes());
    IJK::get_boundary_grid_cubes
      (DIM3, region_axis_size, &(region_boundary_cube[0]));

    visited = new bool[regionScalar.NumVertices()];
    region_vertex_increment = new INDEX_DIFF_TYPE[regionScalar.NumVertices()];

    neighbor_grid.SetSize(grid);
    SetRegionVertexIncrement(grid);
  }

  /// Destructor for IS_ISOPATCH_DISK
  IS_ISOPATCH_DISK::~IS_ISOPATCH_DISK()
  {
    delete [] visited;
    visited = NULL;
    delete [] region_vertex_increment;
    region_vertex_increment = NULL;
  }

  bool IS_ISOPATCH_DISK::IsIsopatchDisk
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue,
     const VERTEX_INDEX cube_index,
     const ISOVERT & isovert,
     const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    std::vector<VERTEX_INDEX> vlist;
    std::vector<VERTEX_INDEX> elist;
    
    SetCubeFlag(cube_index, isovert, gcube_map);
    SetSelectedCubeBoundary();
    SetScalar(scalar_grid, cube_index);

    // Compute number of positive components.
    GetBoundaryPosVertices(isovalue, vlist);
    GetBoundaryPosEdges(isovalue, elist);

    int num_pos;
    IJK::compute_num_connected_components_using_elist
      (IJK::vector2pointer(vlist), vlist.size(), regionScalar.NumVertices(),
       IJK::vector2pointer(elist), elist.size()/2,
       num_pos, visited);

    // Compute number of negative components.
    GetBoundaryNegVertices(isovalue, vlist);
    GetBoundaryNegEdges(isovalue, elist);

    int num_neg;
    IJK::compute_num_connected_components_using_elist
      (IJK::vector2pointer(vlist), vlist.size(), regionScalar.NumVertices(),
       IJK::vector2pointer(elist), elist.size()/2,
       num_neg, visited);

    if (num_pos <= 1 && num_neg <= 1)
      { return(true); }
    else
      { return(false); }
  }

  // Reverse merges to isosurface vertex at cube_index.
  void IS_ISOPATCH_DISK::UnmapAdjacent
  (const NUM_TYPE cube_index, const ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map) const
  {
    for (NUM_TYPE i = 0; i < region_boundary_cube.size(); i++) {
      VERTEX_INDEX region_cube_index = region_boundary_cube[i];

      if (cubeFlag.Scalar(region_cube_index)) {
        INDEX_DIFF_TYPE iv = 
          cube_index + region_vertex_increment[region_cube_index];

        if (0 <= iv && iv < neighbor_grid.NumVertices()) {

          INDEX_DIFF_TYPE gcube_index = isovert.sharp_ind_grid.Scalar(iv);

          // Reset gcube_map[gcube_index] to gcube_index.
          gcube_map[gcube_index] = gcube_index;
        }
      }
    }
  }

  void IS_ISOPATCH_DISK::SetScalar
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const VERTEX_INDEX cube_index)
  {
    VERTEX_INDEX iv0 = 
      cube_index -  neighbor_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);
    regionScalar.CopyRegion(scalar_grid, iv0, region_axis_size, 0);
  }

  void IS_ISOPATCH_DISK::SetRegionVertexIncrement
  (const SHARPISO_GRID & grid)
  {
    const int dimension = grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = grid.AxisSize();
    IJK::PROCEDURE_ERROR error("SetRegionVertexIncrement");

    NUM_TYPE num_region_vertices;
    IJK::compute_num_grid_vertices_in_region
      (dimension, axis_size, 0, region_edge_length, num_region_vertices);

    if (num_region_vertices != cubeFlag.NumVertices()) {
      error.AddMessage
        ("Programming error. Scalar grid is smaller than region.");
      error.AddMessage
        ("  Grid size: ", grid.AxisSize(0), "x", grid.AxisSize(1),
         "x", grid.AxisSize(2), ".");
      error.AddMessage
        ("  Region size: ", cubeFlag.AxisSize(0), "x", cubeFlag.AxisSize(1),
         "x", cubeFlag.AxisSize(2), ".");
      throw error;
    }

    IJK::get_grid_vertices_in_region
      (dimension, axis_size, 0, region_edge_length, region_vertex_increment);

    VERTEX_INDEX offset 
      = neighbor_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);
    for (NUM_TYPE k = 0; k < num_region_vertices; k++)
      { region_vertex_increment[k] = region_vertex_increment[k] - offset; }
  }

  void IS_ISOPATCH_DISK::SetCubeFlag
  (const VERTEX_INDEX cube_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    const INDEX_DIFF_TYPE gcube_index = 
      isovert.sharp_ind_grid.Scalar(cube_index);

    cubeFlag.SetAll(false);

    for (NUM_TYPE i = 0; i < region_boundary_cube.size(); i++) {
      VERTEX_INDEX region_cube_index = region_boundary_cube[i];
      INDEX_DIFF_TYPE iv = 
        cube_index + region_vertex_increment[region_cube_index];

      if (0 <= iv && iv < neighbor_grid.NumVertices()) {

        INDEX_DIFF_TYPE neighbor_gcube_index = 
          isovert.sharp_ind_grid.Scalar(iv);

        if (gcube_index != ISOVERT::NO_INDEX) {
          if (gcube_map[neighbor_gcube_index] == gcube_index) {
            cubeFlag.Set(region_cube_index, true);
          }
        }
      }
    }

  }

  /// Set region boundary.
  /// @pre Region is a 4x4x4 patch and center cube is in region.
  void IS_ISOPATCH_DISK::SetSelectedCubeBoundary()
  {
    selectedCubeBoundary.SetAll(false);

    for (NUM_TYPE i = 0; i < region_boundary_cube.size(); i++) {
      VERTEX_INDEX icube = region_boundary_cube[i];

      if (cubeFlag.Scalar(icube)) {

        for (NUM_TYPE j = 0; j < NUM_CUBE_VERTICES3D; j++) {
          VERTEX_INDEX iv = selectedCubeBoundary.CubeVertex(icube, j);
          if (regionBoundary.Scalar(iv)) 
            { selectedCubeBoundary.Set(iv, true); }
        }
      }
      else {
        for (NUM_TYPE j = 0; j < NUM_CUBE_VERTICES3D; j++) {
          VERTEX_INDEX iv = selectedCubeBoundary.CubeVertex(icube, j);
          if (!regionBoundary.Scalar(iv)) 
            { selectedCubeBoundary.Set(iv, true); }
        }
      }
    }
  }

  void IS_ISOPATCH_DISK::GetBoundaryVertices
  (const SCALAR_TYPE isovalue, const bool flag_pos,
   std::vector<int> & vlist) const 
  {
    if (flag_pos) 
      { GetBoundaryPosVertices(isovalue, vlist); }
    else 
      { GetBoundaryNegVertices(isovalue, vlist); }
  }

  void IS_ISOPATCH_DISK::GetBoundaryPosVertices
  (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const 
  {
    vlist.clear();
    for (NUM_TYPE iv = 0; iv < selectedCubeBoundary.NumVertices(); iv++) {
      if (selectedCubeBoundary.Scalar(iv) && 
          regionScalar.Scalar(iv) >= isovalue)
        { vlist.push_back(iv); }
    }
  }

  void IS_ISOPATCH_DISK::GetBoundaryNegVertices
  (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const 
  {
    vlist.clear();
    for (NUM_TYPE iv = 0; iv < selectedCubeBoundary.NumVertices(); iv++) {
      if (selectedCubeBoundary.Scalar(iv) && 
          regionScalar.Scalar(iv) < isovalue)
        { vlist.push_back(iv); }
    }
  } 

  void IS_ISOPATCH_DISK::GetBoundaryEdges
  (const SCALAR_TYPE isovalue, const bool flag_pos,
   std::vector<int> & elist) const 
  {
    if (flag_pos) 
      { GetBoundaryPosEdges(isovalue, elist); }
    else 
      { GetBoundaryNegEdges(isovalue, elist); }
  }


  void IS_ISOPATCH_DISK::GetBoundaryPosEdges
  (const SCALAR_TYPE isovalue, std::vector<int> & elist) const 
  {
    elist.clear();
    IJK_FOR_EACH_GRID_EDGE(iv0, dir, selectedCubeBoundary, VERTEX_INDEX) 
      {
      VERTEX_INDEX iv1 = selectedCubeBoundary.NextVertex(iv0, dir);
      if (selectedCubeBoundary.Scalar(iv0) && 
          regionScalar.Scalar(iv0) >= isovalue &&
          selectedCubeBoundary.Scalar(iv1) &&
          regionScalar.Scalar(iv1) >= isovalue ) {
        elist.push_back(iv0);
        elist.push_back(iv1);
      }
    }
  }

  void IS_ISOPATCH_DISK::GetBoundaryNegEdges
  (const SCALAR_TYPE isovalue, std::vector<int> & elist) const 
  {
    elist.clear();
    IJK_FOR_EACH_GRID_EDGE(iv0, dir, selectedCubeBoundary, VERTEX_INDEX) 
      {
      VERTEX_INDEX iv1 = selectedCubeBoundary.NextVertex(iv0, dir);
      if (selectedCubeBoundary.Scalar(iv0) && 
          regionScalar.Scalar(iv0) < isovalue &&
          selectedCubeBoundary.Scalar(iv1) &&
          regionScalar.Scalar(iv1) < isovalue ) {
        elist.push_back(iv0);
        elist.push_back(iv1);
      }
    }
  }

}


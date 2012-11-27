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
 std::vector<VERTEX_INDEX> & quad_vert,
 std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();

  determine_gcube_map(scalar_grid, isovalue,
                      isovert, gcube_map, sharpiso_info);

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
 std::vector<VERTEX_INDEX> & quad_vert,
 SHARPISO_INFO & sharpiso_info)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();
  std::vector<VERTEX_INDEX> gcube_map(num_gcube);

  merge_sharp_iso_vertices
    (scalar_grid, isovalue, isovert, quad_vert, gcube_map, sharpiso_info);
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
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
  {
    map_adjacent_cubes(isovert, gcube_map);

    /* DEBUG
    unmap_non_disk_isopatches
      (scalar_grid, isovalue, isovert, gcube_map, sharpiso_info);
    */
  }

}

// **************************************************
// Unmap non-disk isopatches
// **************************************************

namespace {

  /// Function class for determining if isopatch is a disk.
  class IS_ISOPATCH_DISK {

  protected:
    static const AXIS_SIZE_TYPE num_vert_along_region_axis = 4;
    AXIS_SIZE_TYPE region_axis_size[DIM3];

    SHARPISO_SCALAR_GRID regionScalar;
    SHARPISO_BOOL_GRID cubeFlag;
    SHARPISO_BOOL_GRID gridBoundary;
    SHARPISO_BOOL_GRID selectedCubeBoundary;
    SHARPISO_GRID_NEIGHBORS neighbor_grid;
    std::vector<VERTEX_INDEX> selected_cube_boundary_cube_list;
    bool * visited;

    void SetCubeFlag
    (const VERTEX_INDEX cube_index,
     const ISOVERT & isovert,
     const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map);
    void SetScalar
    (const SHARPISO_SCALAR_GRID & scalar_grid, const VERTEX_INDEX cube_index);

    void SetRegionBoundary();
    void GetBoundaryVertices
    (const SCALAR_TYPE isovalue, const bool flag_pos,
     std::vector<int> & vlist) const;
    void GetBoundaryPosVertices
    (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const;
    void GetBoundaryNegVertices
    (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const;
    void GetBoundaryEdges
    (const SCALAR_TYPE isovalue, const bool flag_pos,
     std::vector<int> & elist) const;
    void GetBoundaryPosEdges
    (const SCALAR_TYPE isovalue, std::vector<int> & elist) const;
    void GetBoundaryNegEdges
    (const SCALAR_TYPE isovalue, std::vector<int> & elist) const;

  public:
    IS_ISOPATCH_DISK();
    ~IS_ISOPATCH_DISK();

    bool IsIsopatchDisk
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     SCALAR_TYPE isovalue,
     const VERTEX_INDEX cube_index,
     const ISOVERT & isovert,
     const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map);
  };

}

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
    IS_ISOPATCH_DISK is_isopatch_disk;

    // Set cubes which share facets with selected cubes.
    for (NUM_TYPE i = 0; i < num_gcube; i++) {
      if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

        VERTEX_INDEX cube_index = isovert.gcube_list[i].cube_index;
        if (!is_isopatch_disk.IsIsopatchDisk
            (scalar_grid, isovalue, cube_index, isovert, gcube_map))
          {
            // *** UNMAP ADJACENT VERTICES ***
          }
      }
    }

  }

}

// **************************************************
// Function class IS_ISOPATCH_DISK
// **************************************************

namespace {

  /// Constructor for IS_ISOPATCH_DISK
  IS_ISOPATCH_DISK::IS_ISOPATCH_DISK()
  {
    for (int d = 0; d < DIM3; d++) 
      { region_axis_size[d] = num_vert_along_region_axis; }

    cubeFlag.SetSize(DIM3, region_axis_size);
    selectedCubeBoundary.SetSize(DIM3, region_axis_size);
    regionScalar.SetSize(DIM3, region_axis_size);

    gridBoundary.SetSize(DIM3, region_axis_size);
    compute_boundary_grid(gridBoundary);

    selected_cube_boundary_cube_list.resize
      (selectedCubeBoundary.ComputeNumBoundaryCubes());
    IJK::get_boundary_grid_cubes
      (DIM3, region_axis_size, &(selected_cube_boundary_cube_list[0]));

    visited = new bool[regionScalar.NumVertices()];
  }

  /// Destructor for IS_ISOPATCH_DISK
  IS_ISOPATCH_DISK::~IS_ISOPATCH_DISK()
  {
    delete [] visited;
    visited = NULL;
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

    neighbor_grid.SetSize(scalar_grid);

    SetCubeFlag(cube_index, isovert, gcube_map);
    SetRegionBoundary();

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

  void IS_ISOPATCH_DISK::SetScalar
  (const SHARPISO_SCALAR_GRID & scalar_grid, const VERTEX_INDEX cube_index)
  {
    VERTEX_INDEX iv0 = 
      cube_index -  neighbor_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);
    regionScalar.CopyRegion(scalar_grid, iv0, region_axis_size, 0);

  }

  void IS_ISOPATCH_DISK::SetCubeFlag
  (const VERTEX_INDEX cube_index,
   const ISOVERT & isovert,
   const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map)
  {
    cubeFlag.SetAll(false);

    VERTEX_INDEX iv0 = 
      cube_index -  neighbor_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

    VERTEX_INDEX iv_z = iv0;
    VERTEX_INDEX icube_z = 0;
    for (int z = 0; z < cubeFlag.AxisSize(2); z++) {
      VERTEX_INDEX iv_y = iv_z;
      VERTEX_INDEX icube_y = icube_z;
      for (int y = 0; y < cubeFlag.AxisSize(1); y++) {
        VERTEX_INDEX iv_x = iv_y;
        VERTEX_INDEX icube_x = icube_y;
        for (int x = 0; x < cubeFlag.AxisSize(1); x++) {

          INDEX_DIFF_TYPE gcube_index = isovert.sharp_ind_grid.Scalar(iv_x);
          if (gcube_index != ISOVERT::NO_INDEX) {
            if (gcube_map[gcube_index] == cube_index) {
              cubeFlag.Set(icube_x, true);
            }
          }

          iv_x = neighbor_grid.NextVertex(0, iv_x);
          icube_x = cubeFlag.NextVertex(0, icube_x);
        }
        iv_y = neighbor_grid.NextVertex(1, iv_y);
        icube_y = cubeFlag.NextVertex(1, icube_y);
      }
      iv_z = neighbor_grid.NextVertex(2, iv_z);
      icube_z = cubeFlag.NextVertex(2, icube_z);
    }

  }


  /// Set region boundary.
  /// @pre Region is a 4x4x4 patch and center cube is in region.
  void IS_ISOPATCH_DISK::SetRegionBoundary()
  {
    selectedCubeBoundary.SetAll(false);

    for (NUM_TYPE i = 0; i < selected_cube_boundary_cube_list.size(); i++) {
      VERTEX_INDEX icube = selected_cube_boundary_cube_list[i];

      if (cubeFlag.Scalar(icube)) {

        for (NUM_TYPE j = 0; j < NUM_CUBE_VERTICES3D; j++) {
          VERTEX_INDEX iv = selectedCubeBoundary.CubeVertex(icube, j);
          if (gridBoundary.Scalar(iv)) 
            { selectedCubeBoundary.Set(iv, true); }
        }
      }
      else {
        for (NUM_TYPE j = 0; j < NUM_CUBE_VERTICES3D; j++) {
          VERTEX_INDEX iv = selectedCubeBoundary.CubeVertex(icube, j);
          if (!gridBoundary.Scalar(iv)) 
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


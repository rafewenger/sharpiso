/// \file isodual3D_decimate.cxx
/// Decimate dual isosurface.

/*
  IJK: Isosurface Jeneration Kode
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

#include "ijkgraph.txx"
#include "ijkmesh.txx"
#include "ijkgrid_macros.h"

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"
#include "isodual3D_decimate.h"
#include "isodual3D_extract.h"


// forward declarations
namespace {

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

  void map_gcube_indices(const std::vector<VERTEX_INDEX> & gcube_map,
                         std::vector<VERTEX_INDEX> & gcube_index);

  void set_first_gcube_isov
  (const ISODUAL3D::ISOVERT & isovert, 
   const std::vector<ISODUAL3D::DUAL_ISOVERT> & iso_vlist,
   SHARPISO::NUM_TYPE * first_gcube_isov);
}

// **************************************************
// Merge some isosurface vertices
// **************************************************

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::merge_sharp_iso_vertices
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
void ISODUAL3D::merge_sharp_iso_vertices
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
void ISODUAL3D::merge_sharp_iso_vertices_multi
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
void ISODUAL3D::merge_sharp_iso_vertices_multi
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
    std::vector<NUM_TYPE> sorted_gcube_list;
    using namespace ISODUAL3D;

    const NUM_TYPE num_gcube = isovert.gcube_list.size();
    VERTEX_INDEX cube_index, neighbor_index;
    SHARPISO_GRID_NEIGHBORS gridn;
    GRID_COORD_TYPE cube_coord[DIM3];

    // Set size of grid neighbors grid.
    gridn.SetSize(isovert.sharp_ind_grid);

    sort_gcube_list(isovert.gcube_list, sorted_gcube_list);

    for (NUM_TYPE i = 0; i < num_gcube; i++)
      { gcube_map[i] = i; }

    // Set cubes which share facets with selected cubes.
    for (NUM_TYPE i = 0; i < sorted_gcube_list.size(); i++) {
      NUM_TYPE gcube_index = sorted_gcube_list[i];

      if (isovert.gcube_list[gcube_index].flag == SELECTED_GCUBE) {

        cube_index = isovert.gcube_list[gcube_index].cube_index;

        if (isovert.gcube_list[gcube_index].boundary_bits == 0) {
          // Cube cube_index is an interior cube.

          for (NUM_TYPE j = 0; j < gridn.NumCubeNeighborsF(); j++) {

            neighbor_index = gridn.CubeNeighborF(cube_index, j);

            INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
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
      VERTEX_INDEX cube_index1 = isovert.gcube_list[gcube_index1].cube_index;
      VERTEX_HASH_TABLE::iterator cube1_iter = cube_hash.find(cube_index1);

      // *** DEBUG ***
      if (cube1_iter == cube_hash.end()) {
        using namespace std;
        cerr << "*** ERROR.  Cube index: " << cube_index1
             << " is not in the hash table." << endl;
        exit(1000);
      }
      isov_index[i] = cube1_iter->second;
    }
  }


  // Forward declaration
  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info);

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
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

  void determine_gcube_map
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
  {
    map_adjacent_cubes(isovert, gcube_map);

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
  (ISODUAL3D::ISOVERT & isovert, const VERTEX_INDEX cube_index0,
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
        if (isovert.gcube_list[gcube_index1].flag == COVERED_GCUBE) {
          isovert.gcube_list[gcube_index1].flag = SMOOTH_GCUBE;
        }
      }
    }
  }

  // Reverse merges which create isopatches which are not disks.
  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
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

  // Reverse merges which create isopatches which are not disks.
  void unmap_non_disk_isopatches
  (const SHARPISO::SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
   const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   ISODUAL3D::ISOVERT & isovert, 
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map, 
   ISODUAL3D::SHARPISO_INFO & sharpiso_info)
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
bool ISODUAL3D::is_isopatch_disk3D
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
void ISODUAL3D::get_merged_cubes
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
void ISODUAL3D::get_merged_boundary_edges
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
void ISODUAL3D::extract_dual_quad_isopatch_incident_on
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
void ISODUAL3D::extract_dual_isopatch_incident_on
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
void ISODUAL3D::extract_dual_isopatch_incident_on_multi
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

  extract_dual_quad_isopatch_incident_on
    (scalar_grid, isovalue, isovert, cube_index, gcube_map, dist2cube,
     isoquad_cube, facet_vertex);

  map_cube_list(isovert, gcube_map, isoquad_cube);

  construct_cube_list(isoquad_cube, cube_list);

  NUM_TYPE num_split;

  IJK::split_dual_isovert
    (scalar_grid, isodual_table, isovalue, cube_list,
     isoquad_cube, facet_vertex, iso_vlist, quad_vert2, num_split);

  map_isov_indices(isovert, gcube_map, iso_vlist, quad_vert2);

  tri_vert.clear();
  quad_vert.clear();
  IJK::get_non_degenerate_quad_btlr(quad_vert2, tri_vert, quad_vert);
}

// Extract dual isosurface patch with vertex in merged cube.
// Allow multiple isosurface vertices in each cube.
void ISODUAL3D::extract_dual_isopatch_incident_on_multi
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
void ISODUAL3D::insert_tri_quad_edges
(const std::vector<ISO_VERTEX_INDEX> & tri_vert,
 const std::vector<ISO_VERTEX_INDEX> & quad_vert,
 EDGE_HASH_TABLE & edge_hash)
{
  insert_poly_edges(tri_vert, NUM_VERT_PER_TRI, edge_hash);
  insert_poly_edges(quad_vert, NUM_VERT_PER_QUAD, edge_hash);
}

// Insert or increment vertex count in vertex hash table.
void insert_vertex(const VERTEX_INDEX iv,
                   VERTEX_HASH_TABLE & vertex_hash)
{
  VERTEX_HASH_TABLE::iterator vertex_iter = vertex_hash.find(iv);
  if (vertex_iter == vertex_hash.end()) {
    vertex_hash.insert({iv, 1});
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
      vertex_hash.insert({iv, n});
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

// Construct list of boundary cycle vertices.
void ISODUAL3D::construct_boundary_cycle
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


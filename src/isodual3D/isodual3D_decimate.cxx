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

#include "ijkmesh.txx"

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"
#include "isodual3D_decimate.h"


// forward declarations
namespace {

  void determine_gcube_map
  (const ISODUAL3D::ISOVERT & isovert, 
   ISODUAL3D::DUAL_ISOSURFACE & dual_isosurface,
   std::vector<SHARPISO::VERTEX_INDEX> & gcube_map);

}

// **************************************************
// Merge some isosurface vertices
// **************************************************

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::merge_sharp_iso_vertices
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const ISOVERT & isovert, 
 DUAL_ISOSURFACE & dual_isosurface, 
 std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();

  determine_gcube_map(isovert, dual_isosurface, gcube_map);

  // Count number merged isosurface vertices.
  NUM_TYPE num_merged = 0;
  for (int i = 0; i < num_gcube; i++) {
    if (gcube_map[i] != i) { num_merged++; }
  }
  sharpiso_info.num_merged_iso_vertices = num_merged;

  const NUM_TYPE quad_vert_size = dual_isosurface.quad_vert.size();
  std::vector<VERTEX_INDEX> quad_vert2(quad_vert_size);

  for (VERTEX_INDEX i = 0; i < quad_vert_size; i++) {
    VERTEX_INDEX cube_index = dual_isosurface.quad_vert[i];
    quad_vert2[i] = gcube_map[cube_index];
  }

  dual_isosurface.tri_vert.clear();
  dual_isosurface.quad_vert.clear();

  // Change lower-left, lower-right, upper-left, upper-right order
  //   to counter-clockwise order.
  IJK::reorder_quad_vertices(quad_vert2);

  IJK::get_non_degenerate_quad
    (quad_vert2, dual_isosurface.tri_vert, dual_isosurface.quad_vert);

  // Change counter-clockwise order to lower-left, lower-right,
  //   upper-left, upper-right order.
  IJK::reorder_quad_vertices(dual_isosurface.quad_vert);
}

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::merge_sharp_iso_vertices
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const ISOVERT & isovert, 
 DUAL_ISOSURFACE & dual_isosurface, SHARPISO_INFO & sharpiso_info)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();
  std::vector<VERTEX_INDEX> gcube_map(num_gcube);

  merge_sharp_iso_vertices
    (scalar_grid, isovert, dual_isosurface, gcube_map, sharpiso_info);
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
    
  void determine_gcube_map
  (const ISODUAL3D::ISOVERT & isovert, 
   ISODUAL3D::DUAL_ISOSURFACE & dual_isosurface,
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

}

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


// **************************************************
// Merge some isosurface vertices
// **************************************************

// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
void ISODUAL3D::decimate_dual_isopoly
(const ISOVERT & isovert, DUAL_ISOSURFACE & dual_isosurface)
{
  const NUM_TYPE num_gcube = isovert.gcube_list.size();
  std::vector<VERTEX_INDEX> gcube_map(num_gcube);
  VERTEX_INDEX cube_index, neighbor_index;
  SHARPISO_GRID_NEIGHBORS gridn;
  GRID_COORD_TYPE cube_coord[DIM3];
  int boundary_bits, boundary_bits2;

  // Set size of grid neighbors grid.
  gridn.SetSize(isovert.sharp_ind_grid);

  for (NUM_TYPE i = 0; i < num_gcube; i++)
    { gcube_map[i] = i; }

  for (NUM_TYPE i = 0; i < num_gcube; i++) {
    if (isovert.gcube_list[i].flag == SELECTED_GCUBE) {

      cube_index = isovert.gcube_list[i].index2sg;

      isovert.sharp_ind_grid.ComputeBoundaryBits
        (cube_index, boundary_bits);

      if (boundary_bits == 0) {
        // Cube cube_index is an interior cube.

        for (NUM_TYPE j = 0; j < gridn.NumVertexNeighborsC(); j++) {

          neighbor_index = gridn.VertexNeighborC(cube_index, j);

          INDEX_DIFF_TYPE k = isovert.sharp_ind_grid.Scalar(neighbor_index);
          if (k != ISOVERT::NO_INDEX) {

            if (isovert.gcube_list[k].flag != SELECTED_GCUBE) {

              isovert.sharp_ind_grid.ComputeBoundaryBits
                (neighbor_index, boundary_bits2);

              if (boundary_bits2 == 0) {
                // Map gcube_list[k] to isosurface vertex in cube i.
                gcube_map[k] = i;
              }
            }
          }
        }
      }
      else {
        // *** Fill in. ***
      }
      
    }
  }

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

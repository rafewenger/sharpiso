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
  int boundary_bits;

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
          }
        }
      }
      else {
        // *** Fill in. ***
      }
      
    }
  }
}

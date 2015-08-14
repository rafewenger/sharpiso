/// \file shrec_isovert.cxx
/// Data structures for creating and processing sharp isosurface vertices.

/*
Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iomanip>  
#include <string>  
#include <stdio.h>

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkgrid_macros.h"
#include "ijkscalar_grid.txx"


#include "shrec_types.h"
#include "shrec_isovert.h"
#include "shrec_position.h"

#include "sharpiso_closest.h"

#include "shrec_debug.h"

using namespace std;
using namespace SHARPISO;
using namespace SHREC;
using namespace IJK;


// *** DEBUG ***
#include "ijkprint.txx"
bool flag_debug(false);

inline bool is_uncovered_sharp(const GRID_CUBE_DATA & gcube_data)
{
  GRID_CUBE_FLAG gcube_flag = gcube_data.flag;
  if (gcube_flag == COVERED_POINT || gcube_flag == UNAVAILABLE_GCUBE) {

    if (gcube_data.cube_index == gcube_data.maps_to_cube)
      { return(true); }
  }

  return(false);
}



void recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, ISOVERT & isovert)
{
  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); isovert++) {

    if (is_uncovered_sharp(isovert.gcube_list[i])) {
      bool flag_loc_recomputed;

      recompute_using_adjacent
        (scalar_grid, cube_index, isovert, flag_loc_recomputed);
    }
  }

}

void recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const VERTEX_INDEX cube_index, ISOVERT & isovert,
 bool & flag_loc_recomputed)
{
  const gcube_index = isovert.GCubeIndex(cube_index);
  const BOUNDARY_BITS_TYPE boundary_bits = 
    isovert.gcube_list[gcube_index].boundary_bits;
  SHARPISO_FIXED_ARRAY<NUM_CUBE_NEIGHBORS3D, VERTEX_INDEX, NUM_TYPE>
    neighbor_list;
  COORD_TYPE new_coord[DIM3];

  flag_loc_recomputed = false;
  if (boundary_bits == 0) {

    for (int d = 0; d < DIM3; d++) {

      for (int iadj = -1; iadj < 2; iadj += 2) {

        cube2_index = scalar_grid.AdjacentVertex(cube_index, d, iadj);
        NUM_TYPE gcube2_index = isovert.GCubeIndex(cube2_index);

        if (gcube2_index == ISOVERT::NO_INDEX) { continue; }

        if (is_uncovered_sharp(isovert.gcube_list[gcube2_index])) {

          // SHOULD CHECK THAT SHARED FACET IS ACTIVE

          VERTEX_INDEX cube3_index = 
            isovert.gcube_list[gcube2_index].maps_to_cube;

          if (neighbor_list.Contains(cube3_index) == false) 
            { neighbor_list.PushBack(cube3_index); }
        }
      }
    }

    if (neighbor_list.NumElements() < 3) { return; }

    IJK::set_coord_3D(0, new_coord);
    for (NUM_TYPE i = 0; i < neighbor_list.NumElements(); i++) {
      VERTEX_INDEX cube3_index = neighbor_list[i];
      NUM_TYPE gcube3_index = isovert.GCubeIndex(cube3_index);
      IJK::add_coord_3D(new_coord, isovert.IsovertCoord(gcube3_index), new_coord);
    }

    IJK::divide_coord
      (DIM3, COORD_TYPE(neighbor_list.NumElements()), new_coord, new_coord);

    IJK::copy_coord_3D(new_coord, isovert.gcube_list[gcube_index].isovert_coord);

    // SET SOME FLAG IN GRID_CUBE_DATA
  }
  else {
    // HANDLE BOUNDARY_CUBE
  }

}



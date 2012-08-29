/// \file sharpiso_intersect.h
/// Use sharp formulat to compute intersection of isosurface and grid edge.
/// Version v0.1.1

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

#ifndef _SHARPISO_INTERSECT_
#define _SHARPISO_INTERSECT_

#include "sharpiso_grids.h"
#include "sharpiso_types.h"

namespace SHARPISO {

  void intersect_isosurface_grid_edge_sharp3D
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const GRADIENT_GRID_BASE & gradient_grid,
     const SCALAR_TYPE isovalue, 
     const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
     VERTEX_INDEX & iv2, COORD_TYPE * coord);

};

#endif


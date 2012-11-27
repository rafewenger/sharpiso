/// \file isodual3D_decimate.h
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


#ifndef _ISODUAL3D_DECIMATE_
#define _ISODUAL3D_DECIMATE_

#include <vector>

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"


namespace ISODUAL3D {


  // **************************************************
  // Merge some isosurface vertices
  // **************************************************

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & quad_vert,
   std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info);

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & quad_vert,
   SHARPISO_INFO & sharpiso_info);

};


#endif

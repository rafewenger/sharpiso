/// \file isodual3D_extract.h
/// Subroutines for extracting dual isosurface mesh

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

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

#ifndef _ISODUAL3D_EXTRACT_
#define _ISODUAL3D_EXTRACT_

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"

namespace ISODUAL3D {

  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  void extract_dual_isopoly
    (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
     ISODUAL_INFO & isodual_info);
}

#endif

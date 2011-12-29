/// \file sharpiso_cubes.h
/// Definitions of sharp isosurface cube types.
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
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

#ifndef _SHARPISO_CUBES_
#define _SHARPISO_CUBES_

#include "sharpiso_types.h"

#include "ijkcube.txx"

/// Definitions for sharp isosurface processing.
namespace SHARPISO {

  // **************************************************
  // CUBE DATA STRUCTURES
  // **************************************************

  typedef IJK::CUBE3D<NUM_TYPE, NUM_TYPE, COORD_TYPE, COORD_TYPE>
  SHARPISO_CUBE;                        ///< Cube

  typedef IJK::UNIT_CUBE3D<NUM_TYPE, NUM_TYPE, GRID_COORD_TYPE>
  SHARPISO_UNIT_CUBE;                   ///< Unit cube

};

#endif

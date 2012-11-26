/// \file sharpiso_grids.h
/// Definitions of sharp isosurface scalar and gradient grids.
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

#ifndef _SHARPISO_GRIDS_
#define _SHARPISO_GRIDS_

#include "sharpiso_types.h"

#include "ijkgrid.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"

/// Definitions for sharp isosurface processing.
namespace SHARPISO {

  // **************************************************
  // TYPES
  // **************************************************

  typedef int AXIS_SIZE_TYPE;
  typedef int GRADIENT_LENGTH_TYPE;

  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJK::GRID_PLUS<NUM_TYPE, AXIS_SIZE_TYPE, VERTEX_INDEX, NUM_TYPE>
    SHARPISO_GRID;                  ///< Regular grid.
  typedef IJK::GRID_NEIGHBORS
    <NUM_TYPE, AXIS_SIZE_TYPE, VERTEX_INDEX, INDEX_DIFF_TYPE, NUM_TYPE>
    SHARPISO_GRID_NEIGHBORS;        ///< Grid with neighbor data.

  typedef IJK::SCALAR_GRID_BASE<SHARPISO_GRID, SCALAR_TYPE>
    SHARPISO_SCALAR_GRID_BASE;              ///< sharpiso base scalar grid.
  typedef IJK::SCALAR_GRID_WRAPPER<SHARPISO_GRID, SCALAR_TYPE>
    SHARPISO_SCALAR_GRID_WRAPPER;   ///< sharpiso scalar grid wrapper.
  typedef IJK::SCALAR_GRID<SHARPISO_GRID, SCALAR_TYPE>
    SHARPISO_SCALAR_GRID;           ///< sharpiso scalar grid.
  typedef IJK::VECTOR_GRID_BASE
    <SHARPISO_GRID, GRADIENT_LENGTH_TYPE, GRADIENT_COORD_TYPE>
    GRADIENT_GRID_BASE;             ///<  sharpiso base gradient grid
  typedef IJK::VECTOR_GRID
    <SHARPISO_GRID, GRADIENT_LENGTH_TYPE, GRADIENT_COORD_TYPE>
    GRADIENT_GRID;                  ///< sharpiso gradient grid

  /// Index grid.  Signed to allow for -1.
  typedef IJK::SCALAR_GRID<SHARPISO_GRID, INDEX_DIFF_TYPE> SHARPISO_INDEX_GRID;

  /// Edge index grid.
  typedef IJK::VECTOR_GRID
    <SHARPISO_GRID, GRADIENT_LENGTH_TYPE, INDEX_DIFF_TYPE>
    SHARPISO_EDGE_INDEX_GRID;

  typedef IJK::BOOL_GRID_BASE<SHARPISO_GRID> 
    SHARPISO_BOOL_GRID_BASE;        ///< Boolean grid base.
  typedef IJK::BOOL_GRID<SHARPISO_GRID> 
    SHARPISO_BOOL_GRID;             ///< Boolean grid.
};

#endif

#ifndef _RELIGRADIENT_DATASTRUCT_
#define _RELIGRADIENT_DATASTRUCT_

#include <string>

#include "ijk.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"
#include "ijkmerge.txx"

#include "sharpiso_types.h"

using namespace SHARPISO;

namespace RELIGRADIENT {
  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJK::GRID_PLUS<int, AXIS_SIZE_TYPE, VERTEX_INDEX, VERTEX_INDEX> 
    RELIGRADIENT_GRID;                  ///< Marching Cubes grid.

  typedef IJK::SCALAR_GRID_BASE<RELIGRADIENT_GRID, SCALAR_TYPE> 
    RELIGRADIENT_SCALAR_GRID_BASE;      ///< isodual base scalar grid.
  typedef IJK::SCALAR_GRID_WRAPPER<RELIGRADIENT_GRID, SCALAR_TYPE>
    RELIGRADIENT_SCALAR_GRID_WRAPPER;   ///< isodual scalar grid wrapper.
  typedef IJK::SCALAR_GRID<RELIGRADIENT_GRID, SCALAR_TYPE> 
    RELIGRADIENT_SCALAR_GRID;           ///< Marching Cubes scalar grid.
  typedef IJK::VECTOR_GRID_BASE
  <RELIGRADIENT_GRID, GRADIENT_LENGTH_TYPE, GRADIENT_COORD_TYPE> 
    GRADIENT_GRID_BASE;            ///< isodual base gradient grid
  typedef IJK::VECTOR_GRID
  <RELIGRADIENT_GRID, GRADIENT_LENGTH_TYPE, GRADIENT_COORD_TYPE> 
    GRADIENT_GRID;                 ///< isodual gradient grid
}
#endif

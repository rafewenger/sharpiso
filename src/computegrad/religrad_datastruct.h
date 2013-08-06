#ifndef _RELIGRADIENT_DATASTRUCT_
#define _RELIGRADIENT_DATASTRUCT_

#include <string>

#include "ijk.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"
#include "ijkmerge.txx"

#include "sharpiso_types.h"
#include "sharpiso_grids.h"

namespace RELIGRADIENT {
  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef SHARPISO::SHARPISO_GRID RELIGRADIENT_GRID;
  typedef SHARPISO::SHARPISO_SCALAR_GRID_BASE RELIGRADIENT_SCALAR_GRID_BASE;
  typedef SHARPISO::SHARPISO_SCALAR_GRID RELIGRADIENT_SCALAR_GRID;
}
#endif

/// \file sharpiso_linear_alg.txx
/// sharpiso templates for linear algebra
/// Version 0.1.0

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

#ifndef _SHARPISO_LINEAR_ALG_
#define _SHARPISO_LINEAR_ALG_

#include <cmath>

/// Coordinate arithmetic functions.
namespace SHARPISO {

  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename ZTYPE>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 v0[], CTYPE1 v1[],
   bool & is_v0_zero, const ZTYPE zero_tolerance_squared)
  {
    CTYPE0 magnitude = 0;
    for (DTYPE d = 0; d < dimension; d++)
      { magnitude += v0[d] * v0[d]; }

    if (magnitude > zero_tolerance_squared) {
      magnitude = std::sqrt(magnitude);

      for (DTYPE d = 0; d < dimension; d++)
        { v1[d] = v0[d]/magnitude; }
      is_v0_zero = false;
    }
    else {
      for (DTYPE d = 0; d < dimension; d++)
        { v1[d] = 0; };
      is_v0_zero = true;
    }
  }

  template <typename CTYPE0, typename CTYPE1,
            typename ZTYPE>
  void normalize_vector_3D
  (const CTYPE0 v0[], CTYPE1 v1[],
   bool & is_v0_zero, const ZTYPE zero_tolerance_squared)
  {
    const int DIM3 = 3;

    normalize_vector(DIM3, v0, v1, is_v0_zero, zero_tolerance_squared);
  }

}

#endif

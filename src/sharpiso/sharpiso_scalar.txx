/// \file sharpiso_scalar.txx
/// Compute scalar values based on gradients.
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011,2012 Rephael Wenger

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

#ifndef _SHARPISO_SCALAR_TXX_
#define _SHARPISO_SCALAR_TXX_

#include <cmath>

#include "ijkcube.txx"
#include "sharpiso_types.h"

/// Definitions
namespace SHARPISO {

  // **************************************************
  // TEMPLATE ROUTINES TO COMPUTE SCALAR VALUES
  // **************************************************

  /// Compute scalar at coord0 based on gradient and scalar at coord1.
  template <typename CTYPE0, typename CTYPE1, typename GTYPE,
            typename STYPE>
  inline STYPE compute_gradient_based_scalar
  (const CTYPE0 coord0[DIM3], const CTYPE1 coord1[DIM3],
   const GTYPE gradient1[DIM3], const STYPE scalar1)
  {
    STYPE scalar0 =
      (coord0[0]-coord1[0])*gradient1[0] +
      (coord0[1]-coord1[1])*gradient1[1] +
      (coord0[2]-coord1[2])*gradient1[2] + scalar1;

    return(scalar0);
  }

  /// Compute maximum absolute value and average sum of squares
  ///   of difference between given and predicted scalar values.
  template <typename CTYPE0, typename CTYPE1, typename GTYPE,
            typename STYPE0, typename STYPE1, typename NTYPE,
            typename ETYPE>
  void compute_gradient_based_scalar_diff
  (const CTYPE0 coord0[DIM3], const STYPE0 scalar0,
   const CTYPE1 * coord1, const GTYPE * gradient1, const STYPE1 * scalar1,
   const NTYPE num_point1,
   ETYPE & avg_sum_of_diff_squared, ETYPE & max_abs_diff)
  {
    // Initialize output.
    max_abs_diff = 0;
    avg_sum_of_diff_squared = 0;

    for (NTYPE i = 0; i < num_point1; i++) {
      STYPE1 s = compute_gradient_based_scalar
        (coord0, coord1+i*DIM3, gradient1+i*DIM3, scalar1[i]);
      ETYPE diff = s - scalar0;
      if (diff < 0) { diff = -diff; };
      avg_sum_of_diff_squared += diff*diff;
      if (diff > max_abs_diff) { max_abs_diff = diff; }
    }

    if (num_point1 > 0)
      { avg_sum_of_diff_squared = avg_sum_of_diff_squared/num_point1; }
  }

  // **************************************************
  // INTERSECTION TEMPLATE ROUTINES
  // **************************************************

  /// Return true if isosurface intersects line segment
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename GTYPE, typename STYPE0, typename STYPE1>
  bool iso_intersects_line_segment
  (const CTYPE0 endpoint0[DIM3], const CTYPE1 endpoint1[DIM3],
   const CTYPE2 coord[DIM3], const GTYPE gradient[DIM3],
   const STYPE0 scalar, const STYPE1 isovalue)
  {
    STYPE0 s0 =
        compute_gradient_based_scalar(endpoint0, coord, gradient, scalar);
    STYPE1 s1 =
        compute_gradient_based_scalar(endpoint1, coord, gradient, scalar);

    if (s0 < isovalue) {
      if (s1 >= isovalue) { return(true); }
    }
    else {
      if (s1 <= isovalue) { return(true); }
      if (s0 == isovalue) { return(true); }
    }

    return(false);
  }

  /// Return true if isosurface intersects cube
  template <typename CUBE_TYPE, typename CTYPE0, typename GTYPE0,
            typename STYPE0, typename STYPE1>
  bool iso_intersects_cube
  (const CUBE_TYPE & cube,
   const CTYPE0 coord[DIM3], const GTYPE0 gradient[DIM3],
   const STYPE0 scalar, const STYPE1 isovalue)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;
    typedef typename CUBE_TYPE::COORD_TYPE CTYPE2;

    NTYPE icorner;
    IJK::compute_corner_nearest_direction(DIM3, gradient, icorner);

    // Note: icorner is endpoint 0 of diagonal icorner
    //       but endpoint 1 of the diagonal in the opposite direction.
    const CTYPE2 * end1 = cube.DiagonalEnd0Coord(icorner);
    const CTYPE2 * end0 = cube.DiagonalEnd1Coord(icorner);

    bool flag =
      iso_intersects_line_segment
      (end0, end1, coord, gradient, scalar, isovalue);

    return(flag);
  }

  // **************************************************
  // TEMPLATE ROUTINES TO COMPUTE DISTANCE
  // **************************************************

  /// Compute distance from point to plane defined by gradient field.
  /// @param gfield_gradient[] Gradient over entire gradient field.
  /// @pre Magnitude of gfield_gradient[] must be greater than zero.
  /// @param gfield_point[] Point in the gradient field.
  /// @param gfield_point_scalar Scalar value at gfield_point.
  /// @param point[] Compute distance from point to plane.
  /// @param plane_scalar Scalar value of plane in gradient field.
  /// @param[out] distance Distance from point[] to plane.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename STYPE0, typename STYPE1, typename DIST_TYPE>
  void compute_distance_to_gfield_plane
  (const CTYPE0 gfield_gradient[DIM3], const CTYPE1 gfield_point[DIM3],
   const STYPE0 gfield_point_scalar, const CTYPE2 point[DIM3],
   const STYPE1 plane_scalar, DIST_TYPE & distance)
  {
    CTYPE0 gradient_magnitude =
      gfield_gradient[0]*gfield_gradient[0] +
      gfield_gradient[1]*gfield_gradient[1] +
      gfield_gradient[2]*gfield_gradient[2];

    gradient_magnitude = std::sqrt(gradient_magnitude);

    distance =
      gfield_gradient[0]*(point[0]-gfield_point[0]) +
      gfield_gradient[1]*(point[1]-gfield_point[1]) +
      gfield_gradient[2]*(point[2]-gfield_point[2]) +
      gfield_point_scalar - plane_scalar;
    distance = distance/gradient_magnitude;
  }

};

#endif

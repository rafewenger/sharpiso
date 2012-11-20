
#include <cmath>

#include "sharpiso_types.h"

#ifndef _SHARPISO_DIST_
#define _SHARPISO_DIST_

/// Definitions
namespace SHARPISO {

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
      gfield_gradient[2]*(point[1]-gfield_point[2]) +
      gfield_point_scalar - plane_scalar;
    distance = distance/gradient_magnitude;
  }

}

#endif


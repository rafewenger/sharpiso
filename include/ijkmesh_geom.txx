/// \file ijkmesh_geom.txx
/// ijk templates for handling polyhedral mesh geometry
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

#ifndef _IJKMESH_GEOM_
#define _IJKMESH_GEOM_

#include "ijk.txx"
#include "ijkcoord.txx"
#include <vector>

namespace IJK {


  // **************************************************
  // TRIANGULATE POLYGONS
  // **************************************************

  /// Compute the cosine of an angle given by 3 points.
  /// Compute cosine of angle between (coord1-coord0) and (coord2-coord0).
  /// Cosine is represented as (cos_numerator/cos_denominator).
  /// @param cos_numerator Numerator of the cosine.
  /// @param cos_denominator Denominator of the cosine.
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename MTYPE>
  void compute_cos_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, CTYPE3 & cos_angle)
  {
    IJK::ARRAY<CTYPE0> u1(dimension);
    IJK::ARRAY<CTYPE1> u2(dimension);

    IJK::subtract_coord(dimension, coord1, coord0, u1.Ptr());
    IJK::subtract_coord(dimension, coord2, coord0, u2.Ptr());
    IJK::normalize_vector(dimension, u1.Ptr(), max_small_magnitude, u1.Ptr());
    IJK::normalize_vector(dimension, u2.Ptr(), max_small_magnitude, u2.Ptr());

    IJK::compute_inner_product(dimension, u1.Ptr(), u2.Ptr(), cos_angle);
  }

  /// Triangulate a single polygon.
  /// Split maximum polygon angle.
  /// Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// Add new triangles to vector tri_vert.
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_polygon_split_max_angle
  (const DTYPE dimension,
   const CTYPE * vert_coord,
   const NTYPE num_poly_vert,
   const VTYPE0 * poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    CTYPE min_cos;

    if (num_poly_vert < 3) { return; };

    NTYPE k_min_cos = 0;        // Index of the vertex with min_cos
    NTYPE iv0 = poly_vert[0];
    NTYPE iv1 = poly_vert[num_poly_vert-1];
    NTYPE iv2 = poly_vert[1];

    compute_cos_angle
      (dimension, vert_coord+iv0*dimension, vert_coord+iv1*dimension,
       vert_coord+iv2*dimension, max_small_magnitude, min_cos);

    for (NTYPE i = 1; i < num_poly_vert; i++) {
      iv0 = poly_vert[i];
      iv1 = poly_vert[i-1];
      iv2 = poly_vert[((i+1)%num_poly_vert)];
      CTYPE cos_v0;
      compute_cos_angle
        (dimension, vert_coord+iv0*dimension, vert_coord+iv1*dimension,
         vert_coord+iv2*dimension, max_small_magnitude, cos_v0);

      if (cos_v0 < min_cos) {
        k_min_cos = i;
        min_cos = cos_v0;
      }
    }

    triangulate_polygon(num_poly_vert, poly_vert, k_min_cos, tri_vert);
  }

  /// Triangulate a list of quadrilaterals.
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const CTYPE * vert_coord,
   const VTYPE0 * quad_vert,
   const NTYPE num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {

      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      triangulate_polygon_split_max_angle
        (dimension, vert_coord, NUM_VERT_PER_QUAD, quad_vert+k, 
         max_small_magnitude, tri_vert);
    }
  }

  /// Triangulate a list of quadrilaterals.
  /// C++ STL vector format for vert_coord and quad_vert.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vert_coord,
   const std::vector<VTYPE0> & quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_split_max_angle
      (dimension, IJK::vector2pointer(vert_coord),
       IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
       tri_vert);
  }

}

#endif

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
  // CONVERT QUADRILATERALS TO TRIANGLES
  // **************************************************

  /// Convert quadrilaterals to triangles.
  /// Split maximum quadrilateral angle.
  /// Quadrilateral vertices are listed in order:
  ///   Lower-left, Lower-right, Upper-Left, Upper-Right
  /// Add new triangles to vector tri_vert.
  /// @param max_small_magnitude Vectors with magnitude less than
  ///   max_small_magnitude are considered zero vectors.
  template <typename DTYPE, typename VTYPE0, typename CTYPE0, 
            typename MTYPE, typename VTYPE1>
  void convert_quad_to_tri_split_max_angle
  (const DTYPE dimension,
   const std::vector<CTYPE0> & vert_coord,
   const std::vector<VTYPE0> & quad_vert, 
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    VTYPE0 v[NUM_VERT_PER_QUAD];
    CTYPE0 w[NUM_VERT_PER_QUAD][dimension];
    MTYPE cos_w[NUM_VERT_PER_QUAD];
    SIZE_TYPE isplit;
    MTYPE cos_vsplit;
                                

    for (SIZE_TYPE i = 0; i < quad_vert.size(); i+=NUM_VERT_PER_QUAD) {

      // Get quadrilateral vertices in CYCLIC ORDER.
      v[0] = quad_vert[i];
      v[1] = quad_vert[i+1];
      v[2] = quad_vert[i+3]; // quad_vert[i+3] first for cyclic order
      v[3] = quad_vert[i+2];

      for (SIZE_TYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
        SIZE_TYPE j2 = (j+1)%NUM_VERT_PER_QUAD;

        VTYPE0 ivj = v[j];
        VTYPE0 ivj2 = v[j2];

        subtract_coord
          (dimension, &(vert_coord[ivj2*dimension]), 
           &vert_coord[ivj*dimension], w[j]);

        normalize_vector(dimension, w[j], max_small_magnitude, w[j]);
      }

      for (SIZE_TYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
        SIZE_TYPE j3 = (j+3)%NUM_VERT_PER_QUAD;

        compute_inner_product(dimension, w[j], w[j3], cos_w[j]);
      }

      isplit = 0;
      cos_vsplit = cos_w[0];
      for (SIZE_TYPE j = 1; j < NUM_VERT_PER_QUAD; j++) {

        if (cos_w[j] > cos_vsplit) {
          isplit = j;
          cos_vsplit = cos_w[j];
        }
      }

      // Add first triangle to tri_vert.
      tri_vert.push_back(v[isplit]);
      tri_vert.push_back(v[(isplit+1)%NUM_VERT_PER_QUAD]);
      tri_vert.push_back(v[(isplit+2)%NUM_VERT_PER_QUAD]);

      // Add second triangle to tri_vert.
      tri_vert.push_back(v[(isplit+2)%NUM_VERT_PER_QUAD]);
      tri_vert.push_back(v[(isplit+3)%NUM_VERT_PER_QUAD]);
      tri_vert.push_back(v[isplit]);
    }

  }

}

#endif

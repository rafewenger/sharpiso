/// \file sharpiso_intersect.cxx
/// Use sharp formulat to compute intersection of isosurface and grid edge.
/// Version v0.1.1

/*
 IJK: Isosurface Jeneration Code
 Copyright (C) 2012-2013 Rephael Wenger
 
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

#include "sharpiso_intersect.h"

#include "ijkcoord.txx"
#include "ijkinterpolate.txx"

namespace {

  using namespace SHARPISO;

  void compute_edge_intersection
  (const SCALAR_TYPE s, const GRADIENT_COORD_TYPE g,
   const SCALAR_TYPE isovalue, SCALAR_TYPE & t);

  bool select_t0
  (const SCALAR_TYPE s0, const SCALAR_TYPE s1,
   const GRADIENT_COORD_TYPE g0, const GRADIENT_COORD_TYPE g1,
   const COORD_TYPE t0, const COORD_TYPE t1);
}

// Compute intersection of edge and plane determined by gradient g, scalar s.
// Intersection point is v0+t*dir, 0 <= t <= 1.
// If plane does not intersect edge in a single point, then t < 0 or t > 1.
void SHARPISO::compute_edge_intersection
(const SCALAR_TYPE s, const GRADIENT_COORD_TYPE g,
 const SCALAR_TYPE isovalue, SCALAR_TYPE & t)
{
  SCALAR_TYPE sdiff = isovalue - s;

  if (sdiff > 0) {
    if (g > sdiff) { t = sdiff/g; }
    else if (g == sdiff) { t = 1; }
    else if (g <= 0) { t = -1; }
    else { t = 2; }
  }
  else if (sdiff < 0) {
    if (g < sdiff) { t = sdiff/g; }
    else if (g == sdiff) { t = 1; }
    else if (g >= 0) { t = -1; }
    else { t = 2; }
  }
  else { t = 0; }
}

// Compute intersection of isosurface and grid edge.
void SHARPISO::compute_isosurface_grid_edge_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
 COORD_TYPE p[DIM3])
{
  const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
  const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);
  const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
  const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
  COORD_TYPE t0, t1;
  COORD_TYPE coord0[DIM3], coord1[DIM3];

  compute_edge_intersection(s0, g0, isovalue, t0);
  compute_edge_intersection(s1, -g1, isovalue, t1);
  t1 = 1-t1;

  if (s0 == isovalue) { 
    scalar_grid.ComputeCoord(iv0, p);
    return; 
  }

  if (s1 == isovalue) { 
    scalar_grid.ComputeCoord(iv1, p);
    return;
  }

  scalar_grid.ComputeCoord(iv0, coord0);
  scalar_grid.ComputeCoord(iv1, coord1);

  if (0 <= t0 && t0 <= 1) {
    if (0 <= t1 && t1 <= 1) {
      if (select_t0(s0, s1, g0, g1, t0, t1)) {
        IJK::linear_interpolate_coord(DIM3, 1-t0, coord0, coord1, p);
      }
      else {
        IJK::linear_interpolate_coord(DIM3, 1-t1, coord0, coord1, p);
      }
    }
    else {
      // Use t0.
      IJK::linear_interpolate_coord(DIM3, 1-t0, coord0, coord1, p);
    }
  }
  else {
    if (0 <= t1 && t1 <= 1) {
      // Use t1.
      IJK::linear_interpolate_coord(DIM3, 1-t1, coord0, coord1, p);
    }
    else {
      // Use linear interpolation to compute intersection point.
      IJK::linear_interpolate_coord<float>
        (DIM3, s0, coord0, s1, coord1, isovalue, p);
    }
  }

}


// Compute intersection of isosurface and grid edge and 
//    normal at the intersection point.
void SHARPISO::compute_isosurface_grid_edge_intersection
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
 const GRADIENT_COORD_TYPE max_small_magnitude,
 COORD_TYPE p[DIM3],
 GRADIENT_COORD_TYPE normal[DIM3])
{
  const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
  const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);
  const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
  const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
  COORD_TYPE t0, t1;
  COORD_TYPE coord0[DIM3], coord1[DIM3];

  compute_edge_intersection(s0, g0, isovalue, t0);
  compute_edge_intersection(s1, -g1, isovalue, t1);
  t1 = 1-t1;

  if (s0 == isovalue) { 
    scalar_grid.ComputeCoord(iv0, p);
    IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv0), 
                          max_small_magnitude, normal);
    return; 
  }

  if (s1 == isovalue) { 
    scalar_grid.ComputeCoord(iv1, p);
    IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv1), 
                          max_small_magnitude, normal);
    return;
  }

  scalar_grid.ComputeCoord(iv0, coord0);
  scalar_grid.ComputeCoord(iv1, coord1);

  if (0 <= t0 && t0 <= 1) {
    if (0 <= t1 && t1 <= 1) {
      if (select_t0(s0, s1, g0, g1, t0, t1)) {
        IJK::linear_interpolate_coord(DIM3, 1-t0, coord0, coord1, p);
        IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv0), 
                              max_small_magnitude, normal);  
      }
      else {
        IJK::linear_interpolate_coord(DIM3, 1-t1, coord0, coord1, p);
        IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv1), 
                              max_small_magnitude, normal);
      }
    }
    else {
      // Use t0.
      IJK::linear_interpolate_coord(DIM3, 1-t0, coord0, coord1, p);
      IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv0), 
                            max_small_magnitude, normal);  
    }
  }
  else {
    if (0 <= t1 && t1 <= 1) {
      // Use t1.
      IJK::linear_interpolate_coord(DIM3, 1-t1, coord0, coord1, p);
      IJK::normalize_vector(DIM3, gradient_grid.VectorPtrConst(iv1), 
                            max_small_magnitude, normal);
    }
    else {
      // Use linear interpolation to compute intersection point
      //   and surface normal.
      IJK::linear_interpolate_coord<float>
        (DIM3, s0, coord0, s1, coord1, isovalue, p);
      IJK::linear_interpolate_coord<float>
        (DIM3, s0, gradient_grid.VectorPtrConst(iv0), 
         s1, gradient_grid.VectorPtrConst(iv1), isovalue, normal);
      IJK::normalize_vector(DIM3, normal, max_small_magnitude, normal);
    }
  }

}


namespace {

  using namespace SHARPISO;

  // Return true if t0 should be selected.
  // @pre 0 <= t0 <= 1 and 0 <= t1 <= 1.
  bool select_t0
  (const SCALAR_TYPE s0, const SCALAR_TYPE s1,
   const GRADIENT_COORD_TYPE g0, const GRADIENT_COORD_TYPE g1,
   const COORD_TYPE t0, const COORD_TYPE t1)
  {
    if (s0 < s1) {
      if (g0 > g1) {
        if (t0 < t1) { return(false); }
      }
      else {
        // g0 <= g1
        if (t0 > t1) { return(false); }
      }
    }
    else if (s0 > s1) {
      if (g0 > g1) {
        if (t0 > t1) { return(false); }
      }
      else {
        // g0 <= g1
        if (t0 < t1) { return(false); }
      }
    }

    return(true);
  }

}

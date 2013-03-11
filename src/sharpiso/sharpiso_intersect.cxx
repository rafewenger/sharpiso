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


void SHARPISO::intersect_isosurface_grid_edge_sharp3D
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue, 
 const VERTEX_INDEX iv0, const VERTEX_INDEX iv1, const int dir,
 VERTEX_INDEX & iv2, COORD_TYPE * coord)
{
  const SCALAR_TYPE s0 = scalar_grid.Scalar(iv0);
  const SCALAR_TYPE s1 = scalar_grid.Scalar(iv1);
  const GRADIENT_COORD_TYPE g0 = gradient_grid.Vector(iv0, dir);
  const GRADIENT_COORD_TYPE g1 = gradient_grid.Vector(iv1, dir);
  COORD_TYPE t0, t1;

  compute_edge_intersection(s0, g0, isovalue, t0);
  compute_edge_intersection(s1, -g1, isovalue, t1);
  t1 = 1-t1;

  iv2 = iv0;  // default

  if (s0 == isovalue) { 
    scalar_grid.ComputeCoord(iv0, coord);
    return; 
  }

  if (s1 == isovalue) { 
    iv2 = iv1; 
    scalar_grid.ComputeCoord(iv1, coord);
    return;
  }

  if (0 <= t1 && t1 <= 1) {
    if (t0 < 0 || t0 > 1) { iv2 = iv1; }
    else {
      if (!select_t0(s0, s1, g0, g1, t0, t1))
        { iv2 = iv1; }
    }
  }

  scalar_grid.ComputeCoord(iv0, coord);
  if (iv2 == iv0) {
    t0 = IJK::clamp01_coord(t0);
    coord[dir] += t0;
  }
  else {
    t1 = IJK::clamp01_coord(t1);
    coord[dir] += t1;
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

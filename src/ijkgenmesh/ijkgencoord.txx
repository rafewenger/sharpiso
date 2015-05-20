/// \file ijkgencoord.txx
/// routines for generating coordinates

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2014 Rephael Wenger

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

#ifndef _IJKGENCOORD_
#define _IJKGENCOORD_

#include <cmath>
#include <cstdlib>

#include "ijkcoord.txx"

/// Routines for generating coordinates.
namespace IJKGENCOORD {

  // **********************************************************************
  // Type conversion functions.  Round integers.
  // **********************************************************************

  /// Convert to given type.
  template<typename T>
  T convert2type(const double x)
  {
    return(T(x));
  }

  /// Convert to type int
  template<int>
  int convert2type(const double x)
  {
    double y = floor(x+0.5);
    return(int(y));
  }

  /// Convert to type char
  template<char>
  char convert2type(const double x)
  {
    double y = floor(x+0.5);
    return(char(y));
  }

  /// Convert to type char
  template<unsigned char>
  unsigned char convert2type(const double x)
  {
    double y = floor(x+0.5);
    return((unsigned char)(y));
  }

  // **********************************************************************
  // Vector operators.
  // **********************************************************************

  /// Compute 3D orthogonal basis.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void compute_orthogonal_basis3D
  (const VTYPE0 v0[], const VTYPE1 v1[],
   VTYPE2 xdir[], VTYPE2 ydir[], VTYPE2 zdir[])
  {
    const int DIM3 = 3;

    normalize_vector(DIM3, v0, xdir);
    IJK::compute_orthogonal_vector(DIM3, v1, xdir, ydir);
    normalize_vector(DIM3, ydir);
    IJK::compute_cross_product_3D(xdir, ydir, zdir);
  }

  /// Compute w=(v[0]*xdir[]+v[1]*ydir[]+v[2]*zdir[]).
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void transform_basis3D
  (const VTYPE0 v[], 
   const VTYPE1 xdir[], const VTYPE1 ydir[], const VTYPE1 zdir[],
   VTYPE2 w[])
  {
    const int DIM3 = 3;

    IJK::multiply_coord_3D(v[0], xdir, w);
    IJK::add_scaled_coord_3D(v[1], ydir, w, w);
    IJK::add_scaled_coord_3D(v[2], zdir, w, w);
  }

  /// Compute (L2) length of a vector.
  template <typename DTYPE, typename VTYPE, typename LTYPE>
  void compute_length
  (const DTYPE dimension, const VTYPE v[], LTYPE & length)
  {
    IJK::compute_magnitude(dimension, v, length);
  }

  /// Compute (L1) length of a vector.
  template <typename DTYPE, typename VTYPE, typename LTYPE>
  void compute_length_L1
  (const DTYPE dimension, const VTYPE v[], LTYPE & length)
  {
    length = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      length = length + std::abs(v[d]);
    }
  }

  template <typename DTYPE, typename VTYPE>
  void normalize_vector(const DTYPE dimension, VTYPE v[])
  {
    double length;

    compute_length(dimension, v, length);

    if (length == 0.0) {
      IJK::set_coord(dimension, 0, v);
      return;
    }

    for (DTYPE d = 0; d < dimension; d++) 
      { v[d] = v[d]/length; };
  }

  /// Normalize vector v[].
  /// @param[out] flag_zero_vector True if v[] is a zero vector.
  template <typename DTYPE, typename VTYPE>
  void normalize_vector(const DTYPE dimension, VTYPE v[],
                        bool & flag_zero_vector)
  {
    double length;

    flag_zero_vector = false;
    compute_length(dimension, v, length);

    if (length == 0.0) {
      flag_zero_vector = true;
      IJK::set_coord(dimension, 0, v);
      return;
    }

    for (DTYPE d = 0; d < dimension; d++) 
      { v[d] = v[d]/length; };
  }

  template <typename DTYPE, typename VTYPE0, typename VTYPE1>
  void normalize_vector
  (const DTYPE dimension, const VTYPE0 v0[], VTYPE1 v1[])
  {
    double length;

    compute_length(dimension, v0, length);

    if (length == 0.0) {
      IJK::set_coord(dimension, 0, v1);
      return;
    }

    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = v0[d]/length; };
  }

  /// Normalize vector v0[] and store in v1[].
  /// @param[out] flag_zero_vector True if v0[] is a zero vector.
  template <typename DTYPE, typename VTYPE0, typename VTYPE1>
  void normalize_vector
  (const DTYPE dimension, const VTYPE0 v0[], VTYPE1 v1[],
   bool & flag_zero_vector)
  {
    double length;

    flag_zero_vector = false;
    compute_length(dimension, v0, length);

    if (length == 0.0) {
      flag_zero_vector = true;
      IJK::set_coord(dimension, 0, v1);
      return;
    }

    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = v0[d]/length; };
  }

  /// Compute and normalize component of vector orthogonal to given direction.
  /// @pre dir is a unit vector.
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void compute_normalized_orthogonal_vector
  (const DTYPE dimension, const VTYPE0 v[], 
   const VTYPE1 dir[], VTYPE2 v_orth[])
  {
    IJK::compute_orthogonal_vector(dimension, v, dir, v_orth);
    normalize_vector(dimension, v_orth);
  }

  /// Compute and normalize component of vector orthogonal to given direction.
  /// @param flag_zero_vector True if component is the zero vector.
  /// @pre dir is a unit vector.
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void compute_normalized_orthogonal_vector
  (const DTYPE dimension, const VTYPE0 v[], 
   const VTYPE1 dir[], VTYPE2 v_orth[], bool & flag_zero_vector)
  {
    IJK::compute_orthogonal_vector(dimension, v, dir, v_orth);
    normalize_vector(dimension, v_orth, flag_zero_vector);
  }

  /// Set vector v to -v.  Set magnitude m to -m.
  template <typename DTYPE, typename CTYPE, typename MTYPE>
  void flip_vector
  (const DTYPE dimension, CTYPE v[], MTYPE & m)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { v[d] = -v[d]; }
    m = -m;
  }

  /// Compute angle (radians) between two vectors.
  /// @param dimension Dimension (vector length).
  /// @param v0[] Vector.
  /// @pre v0[] is a unit vector.
  /// @param v1[] Vector.
  /// @pre v1[] is a unit vector.
  /// @param[out] angle Angle in radians between v0[] and v1[].
  /// @remark If v0[] or v1[] are zero vectors, then the output is undefined.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename ANGLE_TYPE>
  void compute_angle
  (const DTYPE dimension, const CTYPE0 v0[], const CTYPE1 v1[],
   ANGLE_TYPE & angle)
  {
    double cos_angle;

    IJK::compute_inner_product(dimension, v0, v1, cos_angle);

    // Snap to interval [-1,1] to correct possible numerical error.
    if (cos_angle > 1) { cos_angle = 1; };
    if (cos_angle < -1) { cos_angle = -1; };

    angle = std::acos(cos_angle);
  }

  /// Compute rotation of v0 to v1 and apply to w1.
  /// @param[out] w1[] Result of rotation applied to w0.
  /// @pre Array w1[] is pre-allocated to length at least dimension.
  template<typename DTYPE, typename MTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, typename CTYPE3>
  void rotate_vector
  (const DTYPE dimension, const MTYPE max_small_magnitude, 
   const CTYPE0 v0[], const CTYPE1 v1[],
   const CTYPE2 w0[], CTYPE3 w1[])
  {
    IJK::ARRAY<double> normalized_v0(dimension), normalized_v1(dimension);
    IJK::ARRAY<double> v0_orth(dimension), v1_orth(dimension);
    IJK::ARRAY<double> u0(dimension), u1(dimension);
    double a0, a1;

    IJK::normalize_vector
      (dimension, v0, max_small_magnitude, normalized_v0.Ptr());
    IJK::normalize_vector
      (dimension, v1, max_small_magnitude, normalized_v1.Ptr());

    IJK::compute_inner_product(dimension, w0, normalized_v0.PtrConst(), a0);
    IJK::multiply_coord(dimension, a0, normalized_v1.PtrConst(), w1);

    if (dimension <= 1) { return; };

    // Compute vectors orthogonal to v0 and v1 in plane (v0,v1).
    IJK::compute_orthogonal_vector
      (dimension, normalized_v1.PtrConst(), normalized_v0.PtrConst(), 
       v0_orth.Ptr());
    IJK::normalize_vector
      (dimension, v0_orth.PtrConst(), max_small_magnitude, v0_orth.Ptr());
    IJK::compute_orthogonal_vector
      (dimension, normalized_v0.PtrConst(), normalized_v1.PtrConst(), 
       v1_orth.Ptr());
    IJK::normalize_vector
      (dimension, v1_orth.PtrConst(), max_small_magnitude, v1_orth.Ptr());

    // Multiply v1_orth by -1 for orientation of (v0_orth,v1_orth)
    //   to match orientation of (v0,v1).
    IJK::multiply_coord(dimension, -1, v1_orth.PtrConst(), v1_orth.Ptr());

    IJK::compute_orthogonal_vector
      (dimension, w0, normalized_v0.PtrConst(), u0.Ptr());
    IJK::compute_inner_product
      (dimension, u0.PtrConst(), v0_orth.PtrConst(), a1);
    IJK::add_scaled_coord(dimension, a1, v1_orth.PtrConst(), w1, w1);

    // Compute component of w0 orthogonal to v0 and v0_orth.
    IJK::compute_orthogonal_vector
      (dimension, u0.PtrConst(), v0_orth.PtrConst(), u1.Ptr());

    IJK::add_coord(dimension, u1.PtrConst(), w1, w1);
  }

  /// Compute rotation of e_d = (0,0,...,1) to v1 and apply to w1.
  /// @param[out] w1[] Result of rotation applied to w0.
  /// @pre Array w1[] is pre-allocated to length at least dimension.
  template<typename DTYPE, typename MTYPE,
           typename CTYPE1, typename CTYPE2, typename CTYPE3>
  void rotate_vector_ed
  (const DTYPE dimension, const MTYPE max_small_magnitude, 
   const CTYPE1 v1[], const CTYPE2 w0[], CTYPE3 w1[])
  {
    IJK::ARRAY<double> ed(dimension);

    if (dimension < 1) { return; }

    IJK::set_coord(dimension, 0, ed.Ptr());
    ed[dimension-1] = 1;

    rotate_vector(dimension, max_small_magnitude, ed.PtrConst(), v1, w0, w1);
  }

  // **********************************************************************
  // Compute distances
  // **********************************************************************

  /// Compute L1 distance between points.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE>
  void compute_L1_distance
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance)
  {
    distance = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      double diff = double(point1[d]-point0[d]);
      if (diff < 0) { diff = -diff; }
      distance += diff;
    }
  }

  /// Compute Linf distance between points.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance)
  {
    distance = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      double diff = double(point1[d]-point0[d]);
      if (diff < 0) { diff = -diff; }
      if (diff > distance)
        { distance = diff; };
    }
  }

  /// Compute Linf distance between points.
  /// @param[out] icoord Coordinate achieving distance.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE, typename ITYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance, ITYPE & icoord)
  {
    distance = 0;
    icoord = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      double diff = double(point1[d]-point0[d]);
      if (diff < 0) { diff = -diff; }
      if (diff > distance) {
        distance = diff; 
        icoord = d;
      }
    }
  }

  /// Compute L2 distance between points.
  /// @param[out] v[] Vector from point0[] to point1[].
  /// @pre v[] is pre-allocated to length at least dimension.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE, typename CTYPE2>
  void compute_L2_distance
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance, CTYPE2 v[])
  {
    IJK::subtract_coord(dimension, point1, point0, v);
    IJK::compute_magnitude(dimension, v, distance);
  }

  /// Compute L2 distance between points.
  /// @param[out] v[] Vector from point0[] to point1[].
  /// @pre v[] is pre-allocated to length at least dimension.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE>
  void compute_L2_distance
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance)
  {
    distance = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      DIST_TYPE diff = point1[d]-point0[d];
      distance += (diff*diff);
    }
    distance = std::sqrt(distance);
  }

  /// Compute distance from a point to a line.
  /// @param dimension Volume dimension.
  /// @param coord0[] Coordinates of a point on the line.
  /// @param dir0[] Line direction.
  /// @pre dir0[] is a unit vector
  /// @param coord1 Point coordinates.
  /// @param[out] distance Distance from point to the line.
  /// @param[out] vdiff[] Vector from line to point.
  /// @param[out] vproj[] Vector from coord0 to orthogonal projection 
  ///             of coord1 on line.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE, 
           typename CTYPE3, typename CTYPE4>
  void compute_dist2line_L2
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance, 
   CTYPE3 v_diff[], CTYPE4 v_proj[])
  {
    IJK::ARRAY<double> v(dimension);

    IJK::subtract_coord(dimension, coord1, coord0, v.Ptr());
    IJK::project_vector(dimension, v.PtrConst(), dir0, v_proj);
    IJK::subtract_coord(dimension, v.PtrConst(), v_proj, v_diff);
    compute_length(dimension, v_diff, distance);
  }

  /// Compute distance from a point to a line.
  /// @param dimension Volume dimension.
  /// @param coord0[] Coordinates of a point on the line.
  /// @param dir0[] Line direction.
  /// @pre dir0[] is a unit vector
  /// @param coord1 Point coordinates.
  /// @param[out] distance Distance from point to the line.
  /// @param[out] vdiff[] Vector from line to point.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE, typename CTYPE3>
  void compute_dist2line_L2
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance, CTYPE3 v_diff[])
  {
    IJK::ARRAY<double> v_proj(dimension);

    compute_dist2line_L2(dimension, coord0, dir0, coord1, distance, 
                         v_diff, v_proj.Ptr());
  }

  /// Compute distance from a point to a line.
  /// @param dimension Volume dimension.
  /// @param coord0[] Coordinates of a point on the line.
  /// @param dir0[] Line direction.
  /// @param coord1 Point coordinates.
  /// @param[out] distance Distance from point to the line.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE>
  void compute_dist2line_L2
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance)
  {
    IJK::ARRAY<double> vdiff(dimension);

    compute_dist2line_L2
      (dimension, coord0, dir0, coord1, distance, vdiff.Ptr());
  }

  /// Compute L1 distance from projection of point on plane.
  /// @pre dir0[] is a unit vector.
  /// @pre dir1[] is a unit vector orthogonal to dir0[].
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIR_TYPE, typename DIST_TYPE,
           typename UTYPE, typename CTYPE2>
  void compute_planar_dist_L1
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const DIR_TYPE dir0[], const DIR_TYPE dir1[], 
   const CTYPE1 coord1[], DIST_TYPE & distance,
   UTYPE & u0, UTYPE & u1, CTYPE2 vdiff[])
  {
    IJK::subtract_coord(dimension, coord1, coord0, vdiff);
    IJK::compute_inner_product(dimension, vdiff, dir0, u0);
    IJK::compute_inner_product(dimension, vdiff, dir1, u1);
    distance = std::abs(u0) + std::abs(u1);
  }

  /// Compute L1 distance from projection of point on plane.
  /// @pre dir0[] is a unit vector.
  /// @pre dir1[] is a unit vector orthogonal to dir0[].
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIR_TYPE, typename DIST_TYPE>
  void compute_planar_dist_L1
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const DIR_TYPE dir0[], const DIR_TYPE dir1[], 
   const CTYPE1 coord1[], DIST_TYPE & distance)
  {
    IJK::ARRAY<double> vdiff(dimension);
    double u0, u1;

    compute_planar_dist_L1
      (dimension, coord0, dir0, dir1, coord1, distance, u0, u1, vdiff.Ptr());
  }

  /// Compute L_inf distance from projection of point on plane.
  /// @pre dir0[] is a unit vector
  /// @param xdir[] X-direction in plane.
  /// @param ydir[] Y-direction in plane.
  /// @pre xdir[] and ydir[] are unit vectors.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIR_TYPE, typename DIST_TYPE>
  void compute_planar_dist_Linf
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const DIR_TYPE xdir[], const DIR_TYPE ydir[], 
   const CTYPE1 coord1[], DIST_TYPE & distance)
  {
    IJK::ARRAY<double> v(dimension);
    double ux, uy;

    IJK::subtract_coord(dimension, coord1, coord0, v.Ptr());
    IJK::compute_inner_product(dimension, v.PtrConst(), xdir, ux);
    IJK::compute_inner_product(dimension, v.PtrConst(), ydir, uy);
    ux = std::abs(ux);
    uy = std::abs(uy);
    distance = std::max(ux, uy);
  }

  /// Compute distance to a cylinder.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename RTYPE, typename CTYPE2, typename DIST_TYPE>
  void compute_dist2cylinder_L2
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 dir0[],
   const RTYPE radius, const CTYPE2 coord1[], DIST_TYPE & distance)
  {
    compute_dist2line_L2(dimension, coord0, dir0, coord1, distance);

    distance = distance - radius;
    if (distance < 0) { distance = -distance; }
  }

  /// Compute distance to a cylinder.
  /// @param[out] distance Distance from point to the cylinder.
  /// @param[out] vdir[] Unit vector pointing from cylinder to point.
  ///    If coord1[] equals coord0[], then vdir[] is zero.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename RTYPE, typename CTYPE2, typename DIST_TYPE,
           typename CTYPE3>
  void compute_dist2cylinder_L2
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 axis_dir[],
   const RTYPE radius, const CTYPE2 coord1[], DIST_TYPE & distance,
   CTYPE3 vdir[])
  {
    compute_dist2line_L2
      (dimension, coord0, axis_dir, coord1, distance, vdir);
    if (distance > 0) {
      double ratio = (distance-radius)/distance;
      IJK::multiply_coord(dimension, ratio, vdir, vdir);
    }

    distance = distance - radius;
    if (distance < 0) { distance = -distance; }
  }

  /// Compute unsigned distance to a plane.
  /// @param[out] vdiff[] Vector from plane to point
  /// @pre orth_dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE, typename CTYPE3>
  void compute_dist2plane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 orth_dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance, CTYPE3 vdiff[])
  {
    IJK::ARRAY<double> v(dimension);

    IJK::subtract_coord(dimension, coord1, coord0, v.Ptr());
    IJK::project_vector(dimension, v.PtrConst(), orth_dir0, vdiff);

    compute_length(dimension, vdiff, distance);
  }

  /// Compute unsigned distance to a plane.
  /// @pre orth_dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE>
  void compute_dist2plane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 orth_dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance)
  {
    IJK::ARRAY<double> vdiff(dimension);

    compute_dist2plane
      (dimension, coord0, orth_dir0, coord1, distance, vdiff.Ptr());
  }

  /// Compute signed distance to a plane.
  /// @param[out] vdiff[] Vector from plane to point
  /// @pre orth_dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE, typename CTYPE3>
  void compute_signed_dist2plane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 orth_dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance, CTYPE3 vdiff[])
  {

    compute_dist2plane(dimension, coord0, orth_dir0, coord1, distance, vdiff);

    double p;
    IJK::compute_inner_product(dimension, orth_dir0, vdiff, p);
    if (p < 0) { distance = -distance; };
  }

  /// Compute signed distance to a plane.
  /// @param[out] vdiff[] Vector from plane to point
  /// @pre orth_dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CTYPE2, typename DIST_TYPE>
  void compute_signed_dist2plane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 orth_dir0[], 
   const CTYPE2 coord1[], DIST_TYPE & distance)
  {
    IJK::ARRAY<double> vdiff(dimension);

    compute_signed_dist2plane
      (dimension, coord0, orth_dir0, coord1, distance, vdiff.Ptr());
  }

  /// Compute distance to annulus.
  /// @param center[] = Coordinates of annulus center.
  /// @param axis_dir[] = Axis direction.
  /// @pre axis_dir[] is a unit vector
  /// @param coord[] = Point coordinates.
  /// @param half_diff_height_width Half of difference 
  ///    of annulus height and width.
  ///    If negative, then width is greater than height.
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename DIFF_TYPE, typename RADIUS_TYPE, typename DIST_TYPE>
  void compute_annulus_dist
  (const DTYPE dimension, 
   const CTYPE0 center[], const CTYPE1 axis_dir[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, const DIFF_TYPE half_diff_height_width,
   DIST_TYPE & distance)
  {
    DIST_TYPE dist0, dist1;

    compute_dist2cylinder_L2(dimension, center, axis_dir, radius, coord, dist0);
    compute_dist2plane(dimension, center, axis_dir, coord, dist1);

    distance  = std::max(dist0, dist1-half_diff_height_width);
  }

  /// Compute distance to flange.
  /// @param half_diff_height_width Half of difference 
  ///    of annulus height and width.
  ///    If negative, then width is greater than height.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename RADIUS_TYPE, typename DIFF_TYPE, typename DIST_TYPE>
  void compute_flange_dist
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 dir0[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, const DIFF_TYPE half_diff_height_width,
   DIST_TYPE & distance)
  {
    DIST_TYPE dist0, dist1, min_dist, max_dist;

    compute_dist2cylinder_L2(dimension, coord0, dir0, radius, coord, dist0);
    compute_dist2plane(dimension, coord0, dir0, coord, dist1);

    dist1 = dist1 - half_diff_height_width;

    min_dist = std::min(dist0, dist1);
    max_dist = std::max(dist0, dist1);
    max_dist = max_dist/2;
    distance = std::max(max_dist, min_dist);
  }

  /// Compute distance to circle.
  /// @param center[] = Coordinates of circle.
  /// @param axis_dir[] = Axis direction.
  /// @pre axis_dir[] is a unit vector
  /// @param coord[] = Point coordinates.
  /// @param[out] vdiff[] Vector from circle to point.
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, 
           typename CTYPE2, typename CTYPE3,
           typename RADIUS_TYPE, typename DIST_TYPE>
  void compute_dist2circle
  (const DTYPE dimension, 
   const CTYPE0 center[], const CTYPE1 axis_dir[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, DIST_TYPE & distance, CTYPE3 vdiff[])
  {
    DIST_TYPE dist0, dist1;
    IJK::ARRAY<double> vdiff0(dimension), vdiff1(dimension);

    compute_dist2cylinder_L2
      (dimension, center, axis_dir, radius, coord, dist0, vdiff0.Ptr());
    compute_dist2plane(dimension, center, axis_dir, coord, dist1, vdiff1.Ptr());

    IJK::add_coord(dimension, vdiff0.PtrConst(), vdiff1.PtrConst(), vdiff);
    distance = std::sqrt(dist0*dist0+dist1*dist1);
  }

  /// Compute distance to circle.
  /// @param center[] = Coordinates of circle.
  /// @param axis_dir[] = Axis direction.
  /// @pre axis_dir[] is a unit vector
  /// @param coord[] = Point coordinates.
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename RADIUS_TYPE, typename DIST_TYPE>
  void compute_dist2circle
  (const DTYPE dimension, 
   const CTYPE0 center[], const CTYPE1 axis_dir[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, DIST_TYPE & distance)
  {
    DIST_TYPE dist0, dist1;

    compute_dist2cylinder_L2
      (dimension, center, axis_dir, radius, coord, dist0);
    compute_dist2plane(dimension, center, axis_dir, coord, dist1);

    distance = std::sqrt(dist0*dist0+dist1*dist1);
  }

  /// Compute signed distance to cone.
  /// @param apex0[] Coordinates of apex of cone.
  /// @param axis_dir[] Axis direction.
  /// @pre axis_dir[] is a unit vector.
  /// @param alpha Cone angle alpha in radians.
  /// @pre alpha > 0.
  /// @param coord[] Point coordinates.
  /// @param[out] v_diff[] Vector from axis to point.
  /// @param[out] v_proj[] Projection of vector from apex0 to point
  ///                      onto axis.
  template<typename DTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename CTYPE3, typename CTYPE4,
           typename LTYPE3, typename LTYPE4,
           typename ANGLE_TYPE, typename DIST_TYPE>
  void compute_dist2cone
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance,
   CTYPE3 v_diff[], LTYPE3 & v_diff_length,
   CTYPE4 v_proj[], LTYPE4 & v_proj_length)
  {
    DIST_TYPE x1, x2;    // v_diff_length = x1 + x2.
    double y;

    compute_dist2line_L2(dimension, apex, axis_dir, coord, v_diff_length,
                         v_diff, v_proj);
    compute_length(dimension, v_proj, v_proj_length);
    x1 = v_proj_length*std::tan(alpha);
    IJK::compute_inner_product(dimension, axis_dir, v_proj, y);
    if (y < 0) {
      x2 = v_diff_length - x1;
    }
    else {
      x2 = v_diff_length + x1;
    }
    distance = x2 * std::cos(alpha);
  }

  /// Compute signed distance to cone.
  /// @param apex0[] Coordinates of apex of cone.
  /// @param axis_dir[] Axis direction.
  /// @pre axis_dir[] is a unit vector.
  /// @param alpha Cone angle alpha in radians.
  /// @pre alpha > 0.
  /// @param coord[] Point coordinates.
  template<typename DTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename ANGLE_TYPE, typename DIST_TYPE>
  void compute_dist2cone
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance)
  {
    IJK::ARRAY<double> v_diff(dimension);
    IJK::ARRAY<double> v_proj(dimension);
    double v_diff_length, v_proj_length;

    compute_dist2cone(dimension, apex, axis_dir, coord, alpha, distance,
                      v_diff.Ptr(), v_diff_length, v_proj.Ptr(), v_proj_length);
  }

  /// Compute signed distance to smooth tipped cone.
  /// @param apex0[] Coordinates of apex of cone.
  /// @param axis_dir[] Axis direction.
  /// @pre axis_dir[] is a unit vector.
  /// @param alpha Cone angle alpha in radians.
  /// @pre alpha > 0.
  /// @param coord[] Point coordinates.
  /// @param[out] flag_dist2apex If true, distance = |coord-apex|.
  /// @param[out] v_diff[] Vector from axis to point.
  /// @param[out] v_proj[] Projection of vector from apex0 to point
  ///                      onto axis.
  template<typename DTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename CTYPE3, typename CTYPE4,
           typename LTYPE3, typename LTYPE4,
           typename ANGLE_TYPE, typename DIST_TYPE>
  void compute_dist2cone_smooth_tip
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance,
   bool & flag_dist2apex,
   CTYPE3 v_diff[], LTYPE3 & v_diff_length,
   CTYPE4 v_proj[], LTYPE4 & v_proj_length)
  {
    DIST_TYPE x1, x2;    // v_diff_length = x1 + x2.
    double y;

    flag_dist2apex = false;

    compute_dist2line_L2(dimension, apex, axis_dir, coord, v_diff_length,
                         v_diff, v_proj);
    compute_length(dimension, v_proj, v_proj_length);
    x1 = v_proj_length*std::tan(alpha);
    IJK::compute_inner_product(dimension, axis_dir, v_proj, y);

    if (y < 0) {
      x2 = v_diff_length - x1;
      distance = x2 * std::cos(alpha);
    }
    else {
      x2 = v_diff_length + x1;

      if (v_proj_length*v_proj_length <= v_diff_length*x1) {
        distance = x2 * std::cos(alpha);
      }
      else {
        // Compute distance to smooth tip of cone.
        IJK::compute_distance(dimension, coord, apex, distance);
        flag_dist2apex = true;
      }
    }
  }

  /// Compute signed distance to smooth tipped cone.
  /// @param apex0[] Coordinates of apex of cone.
  /// @param axis_dir[] Axis direction.
  /// @pre axis_dir[] is a unit vector.
  /// @param alpha Cone angle alpha in radians.
  /// @pre alpha > 0.
  /// @param coord[] Point coordinates.
  /// @param[out] flag_dist2apex If true, distance = |coord-apex|.
  template<typename DTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename ANGLE_TYPE, typename DIST_TYPE>
  void compute_dist2cone_smooth_tip
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance,
   bool & dist2apex)
  {
    IJK::ARRAY<double> v_diff(dimension);
    IJK::ARRAY<double> v_proj(dimension);
    double v_diff_length, v_proj_length;

    compute_dist2cone_smooth_tip
      (dimension, apex, axis_dir, coord, alpha, distance, dist2apex,
       v_diff.Ptr(), v_diff_length, v_proj.Ptr(), v_proj_length);
  }
  /// Compute signed distance to smooth tipped cone.
  /// @param apex0[] Coordinates of apex of cone.
  /// @param axis_dir[] Axis direction.
  /// @pre axis_dir[] is a unit vector.
  /// @param alpha Cone angle alpha in radians.
  /// @pre alpha > 0.
  /// @param coord[] Point coordinates.
  template<typename DTYPE,
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename ANGLE_TYPE, typename DIST_TYPE>
  void compute_dist2cone_smooth_tip
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance)
  {
    bool dist2apex;

    compute_dist2cone_smooth_tip
      (dimension, apex, axis_dir, coord, alpha, distance, dist2apex);
  }

}

#endif

/// \file ijkscalarfield.txx
/// scalar field generation and manipulation routines

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2014 Rephael Wenger

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

#ifndef _IJKSCALARFIELD_
#define _IJKSCALARFIELD_
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>

#include "ijkscalar_grid.txx"
#include "ijkcoord.txx"

/// Routines for generating scalar fields.
namespace IJKSCALARFIELD {

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

  /// Return the dot product of two vectors.
  template <typename DTYPE, typename VTYPE0, typename VTYPE1>
  VTYPE0 dot_product
  (const DTYPE dimension, const VTYPE0 v0[], const VTYPE1 v1[])
  {
    VTYPE0 result(0);
    for (DTYPE d = 0; d < dimension; d++) 
      { result = result + v0[d] * v1[d]; };

    return result;
  }

  /// Compute the determinant of a 2x2 matrix
  template <typename T>
  T determinant2D(const T a0, const T a1, const T b0, const T b1)
  {
    T result = a0*b1 - a1*b0;
    return result;
  }

  /// Compute the cross product of two 3D vectors.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void cross_product3D(const VTYPE0 v0[], const VTYPE1 v1[],
                       VTYPE2 v2[])
  {
    v2[0] = determinant2D(v0[1], v0[2], v1[1], v1[2]);
    v2[1] = -determinant2D(v0[0], v0[2], v1[0], v1[2]);
    v2[2] = determinant2D(v0[0], v0[1], v1[0], v1[1]);
  }

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
    cross_product3D(xdir, ydir, zdir);
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
    LTYPE L_squared = dot_product(dimension, v, v);
    length = sqrt(L_squared);
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
    u0 = dot_product(dimension, vdiff, dir0);
    u1 = dot_product(dimension, vdiff, dir1);
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
    ux = dot_product(dimension, v.PtrConst(), xdir);
    uy = dot_product(dimension, v.PtrConst(), ydir);
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

  // **********************************************************************
  // Transform and combine scalar fields.
  // **********************************************************************

  /// Compute min of two scalar fields.
  template<typename GRID_TYPE>
  void min_scalar(const GRID_TYPE & gridA, const GRID_TYPE & gridB,
                  GRID_TYPE & gridC)
  {
    IJK::PROCEDURE_ERROR error("min_scalar");

    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    if (!gridA.Check(gridB, error)) { throw error; };

    gridC.SetSize(gridA);

    for (VTYPE iv = 0; iv < gridA.NumVertices(); iv++) {
      STYPE s = min(gridA.Scalar(iv), gridB.Scalar(iv));
      gridC.Set(iv, s);
    }
  }

  /// Compute max of two scalar fields.
  template<typename GRID_TYPE>
  void max_scalar(const GRID_TYPE & gridA, const GRID_TYPE & gridB,
                  GRID_TYPE & gridC)
  {
    IJK::PROCEDURE_ERROR error("max_scalar");

    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    if (!gridA.Check(gridB, error)) { throw error; };

    gridC.SetSize(gridA);

    for (VTYPE iv = 0; iv < gridA.NumVertices(); iv++) {
      STYPE s = max(gridA.Scalar(iv), gridB.Scalar(iv));
      gridC.Set(iv, s);
    }
  }

  /// Invert scalar field.  Map x to (s0 - x).
  template<typename GRID_TYPE>
  void invert_scalar(const typename GRID_TYPE::SCALAR_TYPE s0, 
                     GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      STYPE s = s0 - grid.Scalar(iv);
      grid.Set(iv, s);
    }
  }

  /// Invert scalar field.  Change scalar field to f(p) + s0.
  template<typename GRID_TYPE>
  void add_to_scalar(const typename GRID_TYPE::SCALAR_TYPE s0, 
                     GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      STYPE s = s0 + grid.Scalar(iv);
      grid.Set(iv, s);
    }
  }

  /// Multiple scalar field by a scalar.
  template<typename GRID_TYPE, typename ATYPE>
  void multiply_scalar(const ATYPE a,
                       GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      double x = a*grid.Scalar(iv);
      STYPE s = convert2type<STYPE>(x);
      grid.Set(iv, s);
    }
  }


  /// Set values at boundary vertices to s.
  template<typename GRID_TYPE, typename STYPE>
  void set_boundary_values(const STYPE s, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    IJK::GRID_BOUNDARY_VERTEX_LIST<VTYPE> boundary_vlist(grid);

    for (VTYPE i = 0; i < boundary_vlist.NumVertices(); i++) {
      VTYPE iv = boundary_vlist.VertexIndex(i);
      grid.Set(iv, s);
    }
  }

  // **********************************************************************
  // Generate scalar fields.
  // **********************************************************************

  /// Generate a scalar field representing the L2 distance to a point
  template<typename GRID_TYPE, typename COORD_TYPE>
  void gen_dist2point_L2(const COORD_TYPE coord0[], GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      double distance;
      compute_L2_distance(dimension, coord0, coord.PtrConst(), distance);
      STYPE s = convert2type<STYPE>(distance);

      grid.Set(iv, s);
    }
  }
  
  /// Generate a scalar field representing the L1 distance to a point
  template<typename GRID_TYPE, typename COORD_TYPE>
  void gen_dist2point_L1(const COORD_TYPE coord0[], GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    double distance;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_L1_distance(dimension, coord0, coord.PtrConst(), distance);
      STYPE s = convert2type<STYPE>(distance);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field representing the L_infinity distance to a point
  template<typename GRID_TYPE, typename COORD_TYPE>
  void gen_dist2point_Linf(const COORD_TYPE coord0[], GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    double distance;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_Linf_distance
        (grid.Dimension(), coord0, coord.PtrConst(), distance);
      STYPE s = convert2type<STYPE>(distance);
      
      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field representing distance to a line.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE>
  void gen_dist2line_L2
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double distance;

    normalize_vector(grid.Dimension(), dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2line_L2
        (grid.Dimension(), coord0, normalized_dir0.PtrConst(), 
         coord.PtrConst(), distance);
      STYPE s = convert2type<STYPE>(distance);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces are closed cylinders
  /// @param diff_length_diameter Difference of cylinder length and diameter.
  ///    If negative, then diameter is greater than length.
  template<typename GRID_TYPE, typename COORD_TYPE, 
           typename DIR_TYPE, typename DIFF_TYPE>
  void gen_closed_cylinder
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const DIFF_TYPE diff_length_diameter, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_length_diameter/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2line_L2(dimension, coord0, normalized_dir0.PtrConst(), 
                           coord.PtrConst(), dist0);
      compute_dist2plane(dimension, coord0, normalized_dir0.PtrConst(), 
                         coord.PtrConst(), dist1);
      x = std::max(dist0, dist1-half_diff);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces are closed square cylinders.
  /// @param xdir = Direction of x-axis.
  /// @param ydir = Direction of y-axis.  Square lies in xy-plane.
  /// @param zdir = Direction of cylinder axis.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIFF_TYPE>
  void gen_closed_square_cylinder
  (const COORD_TYPE coord0[], const COORD_TYPE xdir[], 
   const COORD_TYPE ydir[], const COORD_TYPE zdir[],
   const DIFF_TYPE diff_height_width, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_xdir(dimension);
    IJK::ARRAY<double> normalized_zdir(dimension);
    IJK::ARRAY<double> ydir_orth(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, xdir, normalized_xdir.Ptr());
    compute_normalized_orthogonal_vector
      (dimension, ydir, normalized_xdir.PtrConst(), ydir_orth.Ptr());
    normalize_vector(dimension, zdir, normalized_zdir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_planar_dist_Linf
        (dimension, coord0, normalized_xdir.PtrConst(), ydir_orth.PtrConst(), 
         coord.PtrConst(), dist0);
      compute_dist2plane
        (dimension, coord0, normalized_zdir.PtrConst(), 
         coord.PtrConst(), dist1);
      x = std::max(dist0, dist1-half_diff);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces are cibes
  /// @param xdir = Direction of x-axis.
  /// @param ydir = Direction of y-axis.  Square lies in xy-plane.
  /// @param zdir = Direction of cylinder axis.
  template<typename GRID_TYPE, typename COORD_TYPE>
  void gen_cube
  (const COORD_TYPE coord0[], const COORD_TYPE xdir[], 
   const COORD_TYPE ydir[], const COORD_TYPE zdir[],
   GRID_TYPE & grid)
  {
    const COORD_TYPE diff_height_width = 0;

    gen_closed_square_cylinder
      (coord0, xdir, ydir, zdir, diff_height_width, grid);
  }

  /// Generate a scalar field whose isosurfaces are closed square cylinders
  ///   rotated 45 degrees around dir1.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIFF_TYPE>
  void gen_closed_square_cylinder_rot45
  (const COORD_TYPE coord0[], const COORD_TYPE dir0[], 
   const COORD_TYPE dir1[], const COORD_TYPE dir2[],
   const DIFF_TYPE diff_height_width,
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    IJK::ARRAY<double> normalized_dir2(dimension);
    IJK::ARRAY<double> dir1_orth(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());
    compute_normalized_orthogonal_vector
      (dimension, dir1, normalized_dir0.PtrConst(), dir1_orth.Ptr());
    normalize_vector(dimension, dir2, normalized_dir2.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_planar_dist_L1
        (dimension, coord0, normalized_dir0.PtrConst(), 
         dir1_orth.PtrConst(), coord.PtrConst(), dist0);
      dist0 = dist0/std::sqrt(2.0);
      compute_dist2plane
        (dimension, coord0, normalized_dir2.PtrConst(), 
         coord.PtrConst(), dist1);
      x = std::max(dist0, dist1-half_diff);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces are closed cones.
  /// @param dimension Volume dimension.
  /// @param apex0 Apex of cone with isovalue 0.
  /// @param axis_dir Cone axis direction.
  /// @param angle Cone angle.
  /// @param height0 Height of cone with isovalue 0.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename ANGLE_TYPE, typename DIST_TYPE>
   void gen_closed_cone
  (const COORD_TYPE apex0[], const DIR_TYPE axis_dir[],
   const ANGLE_TYPE angle, const DIST_TYPE height0, 
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double angle_radians = (angle*M_PI)/180.0;
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_axis_dir(dimension);
    double dist0, dist1, x;

    if (angle <= 0) {
      gen_closed_cylinder(apex0, axis_dir, height0, grid);
      return;
    }

    normalize_vector(dimension, axis_dir, normalized_axis_dir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2cone
        (dimension, apex0, normalized_axis_dir.PtrConst(), 
         coord.PtrConst(), angle_radians, dist0);
      compute_signed_dist2plane
        (dimension, apex0, normalized_axis_dir.PtrConst(), 
         coord.PtrConst(), dist1);
      dist1 = -dist1;
      x = std::max(dist0, dist1-height0);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces are frustra (truncated cones).
  /// @param dimension Volume dimension.
  /// @param apex0 Apex of cone with isovalue 0.
  /// @param axis_dir Cone axis direction.
  /// @param angle Cone angle.
  /// @param dist2near0 Distance from cone apex to near plane with isovalue 0.
  /// @param dist2far0 Distance from cone apex to far plane with isovalue 0.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename ANGLE_TYPE, typename DIST0_TYPE, typename DIST1_TYPE>
   void gen_frustrum
  (const COORD_TYPE apex0[], const DIR_TYPE axis_dir[],
   const ANGLE_TYPE angle, 
   const DIST0_TYPE dist2near0, const DIST1_TYPE dist2far0,
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double angle_radians = (angle*M_PI)/180.0;
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_axis_dir(dimension);
    double dist0, dist1, dist2, x;

    if (angle <= 0) {
      grid.SetAll(0);
      return;
    }

    normalize_vector(dimension, axis_dir, normalized_axis_dir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2cone
        (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
         angle_radians, dist0);
      compute_signed_dist2plane
        (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
         dist1);
      dist2 = -dist1;
      x = std::max(dist1+dist2near0, dist2-dist2far0);
      x = std::max(dist0, x);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces bound thickened annuli.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename RADIUS_TYPE, typename DIFF_TYPE>
  void gen_annulus
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, const DIFF_TYPE diff_height_width, 
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_annulus_dist
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         radius, half_diff, x);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces form a flange.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename RADIUS_TYPE, typename DIFF_TYPE>
  void gen_flange
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, const DIFF_TYPE diff_height_width, 
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_flange_dist
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         radius, half_diff, x);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces bound a torus.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename RADIUS_TYPE>
  void gen_torus
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2circle
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         radius, x);
      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field of maximum distance to planes.
  /// @param coord0[] All planes pass through coord0.
  /// @param normal[] normal[i*dim+j] is j'th coordinate of normal to plane i.
  ///        where dim is the grid dimension.
  template<typename GRID_TYPE, typename CTYPE0, typename CTYPE1,
           typename ITYPE>
  void gen_max_dist2planes
  (const CTYPE0 coord0[], const CTYPE1 normal[], const ITYPE num_planes,
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE0> coord(dimension);
    IJK::ARRAY<CTYPE1> unit_normal(dimension*num_planes);
    double x, distance;

    if (num_planes < 1) {
      grid.SetAll(0);
      return;
    }

    for (ITYPE i = 0; i < num_planes; i++) {
      normalize_vector
        (dimension, normal+i*dimension, unit_normal.Ptr()+i*dimension);
    }

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_signed_dist2plane
        (dimension, coord0, unit_normal.PtrConst(), coord.PtrConst(), x);

      for (ITYPE i = 1; i < num_planes; i++) {
        compute_signed_dist2plane
          (dimension, coord0, unit_normal.PtrConst()+i*dimension, 
           coord.PtrConst(), distance);

        if (distance > x) { x = distance; }
      }

      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field of minimum distance to planes.
  /// @param coord0[] All planes pass through coord0.
  /// @param normal[] normal[i*dim+j] is j'th coordinate of normal to plane i.
  ///        where dim is the grid dimension.
  template<typename GRID_TYPE, typename CTYPE0, typename CTYPE1,
           typename ITYPE>
  void gen_min_dist2planes
  (const CTYPE0 coord0[], const CTYPE1 normal[], const ITYPE num_planes,
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE0> coord(dimension);
    IJK::ARRAY<CTYPE1> unit_normal(dimension*num_planes);
    double x, distance;

    if (num_planes < 1) {
      grid.SetAll(0);
      return;
    }

    for (ITYPE i = 0; i < num_planes; i++) {
      normalize_vector
        (dimension, normal+i*dimension, unit_normal.Ptr()+i*dimension);
    }

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_signed_dist2plane
        (dimension, coord0, unit_normal.PtrConst(), coord.PtrConst(), x);

      for (ITYPE i = 1; i < num_planes; i++) {
        compute_signed_dist2plane
          (dimension, coord0, unit_normal.PtrConst()+i*dimension, 
           coord.PtrConst(), distance);

        if (distance < x) { x = distance; }
      }

      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field with constant gradient
  template<typename GRID_TYPE, typename COORD_TYPE, typename GCOORD_TYPE>
  void gen_constant_gradient
  (const COORD_TYPE coord0[], const GCOORD_TYPE grad0[], 
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> v(dimension);
    double x;

    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      IJK::subtract_coord(dimension, coord.PtrConst(), coord0, v.Ptr());
      IJK::compute_inner_product(dimension, v.PtrConst(), grad0, x);
      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field with constant gradient with unit magnitude
  /// Scalar value at coord0[] is 0.
  template<typename SCALAR_GRID_TYPE,
           typename COORD_TYPE, typename DIR_COORD_TYPE>
  void gen_constant_unit_gradient
  (const COORD_TYPE coord0[], const DIR_COORD_TYPE dir0[],
   SCALAR_GRID_TYPE & grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<double> normalized_dir0(dimension);

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());
    gen_constant_gradient(coord0, normalized_dir0.PtrConst(), grid);
  }

  /// Generate a random scalar field.
  /// Scalar values are uniformly distributed integers from 0 to maxval.
  /// @pre Random seed has already been set.
  template<typename GRID_TYPE, typename VAL_TYPE>
  void gen_random_int
  (const VAL_TYPE maxval, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      //DEBUG
	  //long int x = random();
	  long int x = rand();
      x = x%(maxval+1);
      STYPE s = convert2type<STYPE>(x);
      grid.Set(iv, s);
    }
  }

  /// Generate a random scalar field.
  /// @param seed Random seed.
  template<typename GRID_TYPE, typename VAL_TYPE,
           typename SEED_TYPE>
  void gen_random_int
  (const VAL_TYPE maxval, const SEED_TYPE seed, GRID_TYPE & grid)
  {
	//DEBUG
    srand(seed);
    gen_random_int(maxval, grid);
  }

  /// Generate scalar values for cube representing an isotable entry.
  /// @param cube_index Cube index.  
  ///    Set scalar values of vertices of cube cube_index.
  /// @param table_index Isosurface lookup table index.
  /// @param neg_value Scalar value at negative cube vertex.
  /// @param pos_value Scalar value at positive cube vertex.
  /// @pre cube_index is the index of the lower-left vertex of a grid cube.
  template <typename GRID_TYPE, typename VTYPE, typename ISOTABLE_INDEX_TYPE,
            typename T0, typename T1>
  void gen_isotable_entry
  (const VTYPE cube_index, const ISOTABLE_INDEX_TYPE isotable_index,
   const T0 neg_value, const T1 pos_value, GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = grid.Dimension();
    const NTYPE num_vertices = grid.NumCubeVertices();
    IJK::PROCEDURE_ERROR error("gen_isotable_entry");

    ISOTABLE_INDEX_TYPE mask = 1;
    for (NTYPE i = 0; i < num_vertices; i++) {
      VTYPE iv = grid.CubeVertex(cube_index, i);
      ISOTABLE_INDEX_TYPE x = (mask & isotable_index);

      if (x == 0) 
        { grid.Set(iv, neg_value); }
      else
        { grid.Set(iv, pos_value); }

      mask = (mask << 1);
    }
  }

};

#endif

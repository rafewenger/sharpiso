/// \file ijkgradientfield.txx
/// Gradient field generation and manipulation routines
/// Generate scalar and gradient at each grid vertex.

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

#ifndef _IJKGRADIENTFIELD_
#define _IJKGRADIENTFIELD_

#include <cmath>

#include "ijk.txx"
#include "ijkcoord.txx"
#include "ijkscalar_grid.txx"
#include "ijkscalarfield.txx"

namespace IJKSCALARFIELD {

  // **********************************************************************
  // Check routine.
  // **********************************************************************

  /// Check that gradient grid matches scalar grid.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename ERROR_TYPE>
  bool check_gradient_grid
  (const SCALAR_GRID_TYPE & scalar_grid, 
   const GRADIENT_GRID_TYPE & gradient_grid, ERROR_TYPE & error)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (!scalar_grid.Check
        (gradient_grid, "scalar grid", "gradient grid", error))
      { throw error; }

    if (scalar_grid.Dimension() < 0) {
      error.AddMessage("Illegal negative scalar grid dimension: ",
                       scalar_grid.Dimension(), ".");
      return(false);
    }

    if (gradient_grid.VectorLength() != gradient_grid.Dimension()) {
      error.AddMessage
        ("Incorrect gradient length ", gradient_grid.VectorLength(), ".");
      error.AddMessage("  Gradient length should equal volume dimension ",
                       gradient_grid.Dimension(), ".");
      return(false);
    }

    return(true);
  }

  // **********************************************************************
  // Compute min/max scalar
  // **********************************************************************

  /// Compute min scalar.
  template<typename DTYPE, 
           typename STYPE0, typename STYPE1, typename STYPE2,
           typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void min_scalar_tie_zero
  (const DTYPE dimension, const STYPE0 s0, const CTYPE0 v0[],
   const STYPE1 s1, const CTYPE1 v1[], STYPE2 & s2, CTYPE2 v2[])
  {
    if (s0 < s1) {
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
    else if (s0 > s1) {
      s2 = s1;
      IJK::copy_coord(dimension, v1, v2);
    }
    else {
      // s0 == s1
      s2 = s0;
      IJK::set_coord(dimension, 0, v2);
    }
  }

  /// Compute min scalar.
  template<typename DTYPE, 
           typename STYPE0, typename STYPE1, typename STYPE2,
           typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void min_scalar_tie_use_v0
  (const DTYPE dimension, const STYPE0 s0, const CTYPE0 v0[],
   const STYPE1 s1, const CTYPE1 v1[], STYPE2 & s2, CTYPE2 v2[])
  {
    if (s0 < s1) {
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
    else if (s0 > s1) {
      s2 = s1;
      IJK::copy_coord(dimension, v1, v2);
    }
    else {
      // s0 == s1
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
  }

  /// Compute max scalar.
  template<typename DTYPE, 
           typename STYPE0, typename STYPE1, typename STYPE2,
           typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void max_scalar_tie_zero
  (const DTYPE dimension, const STYPE0 s0, const CTYPE0 v0[],
   const STYPE1 s1, const CTYPE1 v1[], STYPE2 & s2, CTYPE2 v2[])
  {
    if (s0 > s1) {
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
    else if (s0 < s1) {
      s2 = s1;
      IJK::copy_coord(dimension, v1, v2);
    }
    else {
      // s0 == s1
      s2 = s0;
      IJK::set_coord(dimension, 0, v2);
    }
  }

  /// Compute max scalar.
  template<typename DTYPE, 
           typename STYPE0, typename STYPE1, typename STYPE2,
           typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void max_scalar_tie_use_v0
  (const DTYPE dimension, const STYPE0 s0, const CTYPE0 v0[],
   const STYPE1 s1, const CTYPE1 v1[], STYPE2 & s2, CTYPE2 v2[])
  {
    if (s0 > s1) {
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
    else if (s0 < s1) {
      s2 = s1;
      IJK::copy_coord(dimension, v1, v2);
    }
    else {
      // s0 == s1
      s2 = s0;
      IJK::copy_coord(dimension, v0, v2);
    }
  }

  // **********************************************************************
  // Compute distances and gradients
  // **********************************************************************

  /// Compute L1 distance and gradient from point0 to point1.
  /// @param[out] flag_discontinuity True if there is a discontinuity
  ///     in the gradient at the given point. Array gradient[] is arbitrarily
  ///     set to one of the gradients at the discontinuity.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE, typename GCOORD_TYPE>
  void compute_L1_gradient
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance, GCOORD_TYPE * gradient, bool & flag_discontinuity)
  {
    DTYPE icoord;
    compute_L1_distance(dimension, point0, point1, distance);

    IJK::set_coord(dimension, 0, gradient);
    flag_discontinuity = false;

    for (DTYPE d = 0; d < dimension; d++) {
      if (point1[d] > point0[d]) {
        gradient[d] = 1;
      }
      else if (point1[d] < point0[d]) {
        gradient[d] = -1;
      }
      else {
        gradient[d] = 1;  // arbitrarily set to 1
        flag_discontinuity = true;
      }
    }

  }

  /// Compute Linf distance and gradient from point0 to point1.
  /// @param[out] flag_discontinuity True if there is a discontinuity
  ///     in the gradient at the given point. Array gradient[] is arbitrarily
  ///     set to one of the gradients at the discontinuity.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIST_TYPE, typename GCOORD_TYPE>
  void compute_Linf_gradient
  (const DTYPE dimension, const CTYPE0 point0[], const CTYPE1 point1[], 
   DIST_TYPE & distance, GCOORD_TYPE * gradient, bool & flag_discontinuity)
  {
    DTYPE icoord;
    compute_Linf_distance(dimension, point0, point1, distance, icoord);

    IJK::set_coord(dimension, 0, gradient);
    flag_discontinuity = false;

    if (dimension < 1) { return; };

    // set gradient
    double diff = double(point1[icoord]-point0[icoord]);
    if (diff >= 0) { gradient[icoord] = 1; }
    else { gradient[icoord] = -1; }

    // check for discontuity
    for (DTYPE d = 0; d < dimension; d++) {
      double diff = double(point1[d]-point0[d]);
      if (diff < 0) { diff = -diff; };
      if (diff == distance && d != icoord) 
        { flag_discontinuity = true; }
    }

  }

  /// Compute L_inf gradient from projection of point on plane.
  /// @param xdir[] X-direction in plane.
  /// @param ydir[] Y-direction in plane.
  /// @pre xdir[] and ydir[] are unit vectors.
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename DIR_TYPE, typename DIST_TYPE, 
           typename GCOORD_TYPE>
  void compute_planar_Linf_gradient
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const DIR_TYPE xdir[], const DIR_TYPE ydir[], 
   const CTYPE1 coord1[], DIST_TYPE & distance,
   GCOORD_TYPE gradient[])
  {
    IJK::ARRAY<double> v(dimension);
    IJK::ARRAY<double> xdir2(dimension), ydir2(dimension);
    double ux, uy;

    IJK::subtract_coord(dimension, coord1, coord0, v.Ptr());
    IJK::copy_coord(dimension, xdir, xdir2.Ptr());
    IJK::copy_coord(dimension, ydir, ydir2.Ptr());
    ux = dot_product(dimension, v.PtrConst(), xdir2.Ptr());
    uy = dot_product(dimension, v.PtrConst(), ydir2.Ptr());

    if (ux < 0)
      { flip_vector(dimension, xdir2.Ptr(), ux); }
    if (uy < 0)
      { flip_vector(dimension, ydir2.Ptr(), uy); }

    max_scalar_tie_zero
      (dimension, ux, xdir2.PtrConst(), uy, ydir2.PtrConst(), 
       distance, gradient);
  }

  /// Compute distance and gradient for cylinder scalar field.
  /// @param diff_length_diameter Difference of cylinder length and diameter.
  ///    If negative, then diameter is greater than length.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename DIFF_TYPE, typename DIST_TYPE, typename GRADIENT_TYPE>
  void compute_cylinder_gradient
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 dir0[], const CTYPE2 coord[], 
   const DIFF_TYPE diff_length_diameter, 
   DIST_TYPE & distance, GRADIENT_TYPE * gradient)
  {
    const double half_diff = (diff_length_diameter/2.0);
    IJK::ARRAY<double> v0(dimension);
    IJK::ARRAY<double> v1(dimension);
    double dist0, dist1;

    compute_dist2line_L2
      (dimension, coord0, dir0, coord, dist0, v0.Ptr());
    compute_dist2plane
      (dimension, coord0, dir0, coord, dist1, v1.Ptr());

    // normalize gradient vectors
    normalize_vector(dimension, v0.Ptr(), v0.Ptr());
    normalize_vector(dimension, v1.Ptr(), v1.Ptr());

    max_scalar_tie_zero
      (dimension, dist0, v0.PtrConst(), dist1-half_diff, v1.PtrConst(), 
       distance, gradient);
  }

  /// Compute distance and gradient for square cylinder scalar field.
  /// @param xdir = X-direction in plane.
  /// @param ydir = Y-direction in plane.
  /// @param zdir = Direction of cylinder axis.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height
  /// @pre xdir[], ydir[] and zdir[] are unit vectors.
  /// @pre ydir[] is orthogonal to xdir[].
  template<typename DTYPE, typename CTYPE0, typename CTYPE1,
           typename CXTYPE, typename CYTYPE, typename CZTYPE,
           typename DIFF_TYPE, typename DIST_TYPE, typename GCOORD_TYPE>
  void compute_square_cylinder_gradient
  (const DTYPE dimension, const CTYPE0 coord0[],
   const CXTYPE xdir[], const CYTYPE ydir[], const CZTYPE zdir[],
   const CTYPE1 coord[], const DIFF_TYPE diff_height_width, 
   DIST_TYPE & distance, GCOORD_TYPE gradient[])
  {
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<double> v0(dimension);
    IJK::ARRAY<double> v1(dimension);
    double dist0, dist1;

    compute_planar_Linf_gradient
      (dimension, coord0, xdir, ydir, coord, dist0, v0.Ptr());
    compute_dist2plane
      (dimension, coord0, zdir, coord, dist1, v1.Ptr());

    // normalize gradient vector v1
    normalize_vector(dimension, v1.Ptr(), v1.Ptr());

    max_scalar_tie_zero
      (dimension, dist0, v0.PtrConst(), dist1-half_diff, v1.PtrConst(), 
       distance, gradient);
  }

  /// Compute distance and gradient for annulus scalar field.
  /// @param half_diff_height_width Half of difference 
  ///    of annulus height and width.
  ///    If negative, then width is greater than height.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename RADIUS_TYPE, typename DIFF_TYPE, 
           typename DIST_TYPE, typename GRADIENT_TYPE>
  void compute_annulus_gradient
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 dir0[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, const DIFF_TYPE half_diff_height_width,
   DIST_TYPE & distance, GRADIENT_TYPE * gradient)
  {
    IJK::ARRAY<double> v0(dimension);
    IJK::ARRAY<double> v1(dimension);
    double dist0, dist1;

    compute_dist2cylinder_L2
      (dimension, coord0, dir0, radius, coord, dist0, v0.Ptr());
    compute_dist2plane
      (dimension, coord0, dir0, coord, dist1, v1.Ptr());

    // normalize gradient vector v0
    normalize_vector(dimension, v0.Ptr(), v0.Ptr());

    // normalize gradient vector v1
    normalize_vector(dimension, v1.Ptr(), v1.Ptr());

    max_scalar_tie_zero
      (dimension, dist0, v0.PtrConst(), dist1-half_diff_height_width, 
       v1.PtrConst(), distance, gradient);
  }

  /// Compute distance and gradient for flange scalar field.
  /// @param half_diff_height_width Half of difference 
  ///    of annulus height and width.
  ///    If negative, then width is greater than height.
  /// @pre dir0[] is a unit vector
  template<typename DTYPE, 
           typename CTYPE0, typename CTYPE1, typename CTYPE2, 
           typename RADIUS_TYPE, typename DIFF_TYPE, 
           typename DIST_TYPE, typename GRADIENT_TYPE>
  void compute_flange_gradient
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 dir0[], const CTYPE2 coord[], 
   const RADIUS_TYPE radius, const DIFF_TYPE half_diff_height_width,
   DIST_TYPE & distance, GRADIENT_TYPE * gradient)
  {
    IJK::ARRAY<double> v0(dimension);
    IJK::ARRAY<double> v1(dimension);
    IJK::ARRAY<double> min_grad(dimension), max_grad(dimension);
    double dist0, dist1, min_dist, max_dist;

    compute_dist2cylinder_L2
      (dimension, coord0, dir0, radius, coord, dist0, v0.Ptr());
    compute_dist2plane
      (dimension, coord0, dir0, coord, dist1, v1.Ptr());

    dist1 = dist1 - half_diff_height_width;

    // normalize gradient vector v0
    normalize_vector(dimension, v0.Ptr(), v0.Ptr());

    // normalize gradient vector v1
    normalize_vector(dimension, v1.Ptr(), v1.Ptr());

    min_scalar_tie_zero
      (dimension, dist0, v0.PtrConst(), dist1, v1.PtrConst(), 
       min_dist, min_grad.Ptr());
    max_scalar_tie_zero
      (dimension, dist0, v0.PtrConst(), dist1, v1.PtrConst(), 
       max_dist, max_grad.Ptr());

    // rescale max scalar so that max scalar level sets 
    //   intersects min scalar level sets
    const float rescale_max = 0.5;
    max_dist = max_dist * rescale_max;
    IJK::multiply_coord(dimension, rescale_max, max_grad.PtrConst(), 
                        max_grad.Ptr());

    max_scalar_tie_zero
      (dimension, min_dist, min_grad.PtrConst(), max_dist, max_grad.PtrConst(),
       distance, gradient);
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
           typename ANGLE_TYPE, typename DIST_TYPE, typename GRADIENT_TYPE>
  void compute_dist2cone_gradient
  (const DTYPE dimension, const CTYPE0 apex[], const CTYPE1 axis_dir[], 
   const CTYPE2 coord[], const ANGLE_TYPE alpha, DIST_TYPE & distance,
   GRADIENT_TYPE gradient[])
  {
    IJK::ARRAY<double> v_diff(dimension);
    IJK::ARRAY<double> v_proj(dimension);
    IJK::ARRAY<double> u(dimension);
    double v_diff_length, v_proj_length;
    
    compute_dist2cone(dimension, apex, axis_dir, coord, alpha, distance,
                      v_diff.Ptr(), v_diff_length, v_proj.Ptr(), v_proj_length);

    if (alpha <= 0) {
      IJK::copy_coord(dimension, v_diff.PtrConst(), gradient);
    }
    else {
      double x = v_diff_length*tan(alpha);
      IJK::multiply_coord(dimension, x, axis_dir, u.Ptr());
      IJK::add_coord(dimension, v_diff.PtrConst(), u.PtrConst(), gradient);
    }

    normalize_vector(dimension, gradient, gradient);
  }


  // **********************************************************************
  // Generate gradient fields.
  // **********************************************************************

  /// Generate a scalar field representing the L2 distance to a point.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_L2
  (const COORD_TYPE coord0[], SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    double distance;
    IJK::PROCEDURE_ERROR error("gen_gradient_dist2point_L2");

    if (!check_gradient_grid(scalar_grid, gradient_grid, error)) 
      { throw error; }

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_L2_distance
        (dimension, coord0, coord.PtrConst(), distance, gradient.Ptr());
      normalize_vector(dimension, gradient.Ptr());
      STYPE s = convert2type<STYPE>(distance);

      scalar_grid.Set(iv, s);
      gradient_grid.Set(iv, gradient.PtrConst());
    }
  }

  /// Generate a scalar field representing the L_infinity distance to a point.
  ///   Set gradient to zero at gradient discontinuity.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_Linf_zero
  (const COORD_TYPE coord0[], SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    double distance;
    IJK::PROCEDURE_ERROR error("gen_gradient_dist2point_Linf");

    if (!check_gradient_grid(scalar_grid, gradient_grid, error)) 
      { throw error; }

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      bool flag_gradient_discontinuity;
      compute_Linf_gradient
        (dimension, coord0, coord.PtrConst(), distance, gradient.Ptr(), 
         flag_gradient_discontinuity);
      STYPE s = convert2type<STYPE>(distance);

      scalar_grid.Set(iv, s);

      if (flag_gradient_discontinuity) {
        // set gradient to zero
        IJK::set_coord(dimension, 0, gradient.Ptr());
      }

      gradient_grid.Set(iv, gradient.PtrConst());
    }
  }

  /// Generate a scalar field representing the L_infinity distance to a point.
  ///   Select gradient from one incident surface at gradient discontinuity.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_Linf_select
  (const COORD_TYPE coord0[], SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    double distance;
    IJK::PROCEDURE_ERROR error("gen_gradient_dist2point_Linf");

    if (!check_gradient_grid(scalar_grid, gradient_grid, error)) 
      { throw error; }

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      bool flag_gradient_discontinuity;
      compute_Linf_gradient
        (dimension, coord0, coord.PtrConst(), distance, gradient.Ptr(), 
         flag_gradient_discontinuity);
      STYPE s = convert2type<STYPE>(distance);

      scalar_grid.Set(iv, s);

      gradient_grid.Set(iv, gradient.Ptr());
    }
  }

  /// Generate a scalar field representing the L_infinity distance to a point.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_Linf
  (const COORD_TYPE coord0[], const bool tie_zero,
   SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    if (tie_zero) {
      gen_gradient_dist2point_Linf_zero(coord0, scalar_grid, gradient_grid); 
    }
    else {
      gen_gradient_dist2point_Linf_select(coord0, scalar_grid, gradient_grid); 
    }
  }

  /// Generate a scalar field representing the L1 distance to a point.
  ///   Set gradient to zero at gradient discontinuity.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_L1_zero
  (const COORD_TYPE coord0[], SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    double distance;
    IJK::PROCEDURE_ERROR error("gen_gradient_dist2point_L1");

    if (!check_gradient_grid(scalar_grid, gradient_grid, error)) 
      { throw error; }

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      bool flag_gradient_discontinuity;
      compute_L1_gradient
        (dimension, coord0, coord.PtrConst(), distance, gradient.Ptr(),
         flag_gradient_discontinuity);

      STYPE s = convert2type<STYPE>(distance);
      scalar_grid.Set(iv, s);

      if (flag_gradient_discontinuity) {
        // set gradient to zero
        IJK::set_coord(dimension, 0, gradient.Ptr());
      }

      gradient_grid.Set(iv, gradient.PtrConst());
    }
  }

  /// Generate a scalar field representing the L1 distance to a point.
  ///   Select gradient from one incident surface at gradient discontinuity.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_dist2point_L1_select
  (const COORD_TYPE coord0[], SCALAR_GRID_TYPE & scalar_grid,
   GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    double distance;
    IJK::PROCEDURE_ERROR error("gen_gradient_dist2point_L1");

    if (!check_gradient_grid(scalar_grid, gradient_grid, error)) 
      { throw error; }

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      bool flag_gradient_discontinuity;
      compute_L1_gradient
        (dimension, coord0, coord.PtrConst(), distance, gradient.Ptr(),
         flag_gradient_discontinuity);

      STYPE s = convert2type<STYPE>(distance);
      scalar_grid.Set(iv, s);

      gradient_grid.Set(iv, gradient.PtrConst());
    }
  }

  /// Generate a scalar field representing the distance to a line.
  /// Isosurfaces are open cylinders.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_COORD_TYPE>
  void gen_gradient_dist2line_L2
  (const COORD_TYPE coord0[], const DIR_COORD_TYPE dir0[],
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VITYPE;
    typedef typename GRADIENT_GRID_TYPE::VECTOR_COORD_TYPE VCTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    IJK::ARRAY<double> v(dimension);
    IJK::ARRAY<VCTYPE> grad(dimension);
    double distance;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VITYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2line_L2
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         distance, v.Ptr());
      normalize_vector(dimension, v.PtrConst(), grad.Ptr());

      STYPE s = convert2type<STYPE>(distance);
    
      scalar_grid.Set(iv, s);
      gradient_grid.Set(iv, grad.PtrConst());
    }
  }

  /// Generate a scalar field whose isosurfaces are closed cylinders.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE, typename DIFF_TYPE>
  void gen_gradient_closed_cylinder
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const DIFF_TYPE diff_length_diameter, 
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_cylinder_gradient
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         diff_length_diameter, x, gradient.VectorPtr(iv));

      STYPE s = convert2type<STYPE>(x);
    
      scalar_grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces bound thickened annuli.
  /// @param half_diff_height_width Half of difference 
  ///    of annulus height and width.
  ///    If negative, then width is greater than height.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE,
           typename RADIUS_TYPE, typename DIFF_TYPE>
  void gen_gradient_annulus
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, const DIFF_TYPE diff_height_width, 
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double x;

    normalize_vector(scalar_grid.Dimension(), dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_annulus_gradient
        (scalar_grid.Dimension(), coord0, normalized_dir0.PtrConst(), 
         coord.PtrConst(), radius, half_diff, x, gradient.VectorPtr(iv));

      STYPE s = convert2type<STYPE>(x);
    
      scalar_grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces form a flange.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE,
           typename RADIUS_TYPE, typename DIFF_TYPE>
  void gen_gradient_flange
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, const DIFF_TYPE diff_height_width,
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double x;

    normalize_vector(scalar_grid.Dimension(), dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_flange_gradient
        (scalar_grid.Dimension(), coord0, normalized_dir0.PtrConst(), 
         coord.PtrConst(), radius, half_diff, x, gradient.VectorPtr(iv));

      STYPE s = convert2type<STYPE>(x);
    
      scalar_grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces are closed square cylinders.
  /// @param xdir = X-direction in plane.
  /// @param ydir = Y-direction in plane.
  /// @param zdir = Direction of cylinder axis.
  /// @param diff_height_width Difference of cylinder height and width.
  ///    If negative, then width is greater than height
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIFF_TYPE>
  void gen_gradient_closed_square_cylinder
  (const COORD_TYPE coord0[], const COORD_TYPE xdir[], 
   const COORD_TYPE ydir[], const COORD_TYPE zdir[], 
   const DIFF_TYPE diff_height_width,
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_xdir(dimension);
    IJK::ARRAY<double> normalized_zdir(dimension);
    IJK::ARRAY<double> ydir_orth(dimension);
    double x;

    normalize_vector(dimension, xdir, normalized_xdir.Ptr());
    compute_normalized_orthogonal_vector
      (dimension, ydir, normalized_xdir.PtrConst(), ydir_orth.Ptr());
    normalize_vector(dimension, zdir, normalized_zdir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_square_cylinder_gradient
        (dimension, coord0, normalized_xdir.PtrConst(), 
         ydir_orth.PtrConst(), normalized_zdir.PtrConst(),
         coord.PtrConst(), diff_height_width, x, gradient.VectorPtr(iv));

      STYPE s = convert2type<STYPE>(x);
    
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces are cubes.
  /// @param xdir = X-direction in plane.
  /// @param ydir = Y-direction in plane.
  /// @param zdir = Direction of cylinder axis.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE>
  void gen_gradient_cube
  (const COORD_TYPE coord0[], const COORD_TYPE xdir[], 
   const COORD_TYPE ydir[], const COORD_TYPE zdir[],
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    const COORD_TYPE diff_height_width = 0;

    gen_gradient_closed_square_cylinder
      (coord0, xdir, ydir, zdir, diff_height_width, grid, gradient);
  }

  /// Generate a scalar field whose isosurfaces are closed square cylinders
  ///   rotated 45 degrees around dir1.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIFF_TYPE>
  void gen_gradient_closed_square_cylinder_rot45
  (const COORD_TYPE coord0[], const COORD_TYPE dir0[], 
   const COORD_TYPE dir1[], const COORD_TYPE dir2[],
   const DIFF_TYPE diff_height_width,
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double half_diff = (diff_height_width/2.0);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    IJK::ARRAY<double> normalized_dir2(dimension);
    IJK::ARRAY<double> dir1_orth(dimension);
    double dist0, dist1, distance;
    double u0, u1;
    IJK::ARRAY<double> vdiff(dimension), v1(dimension);
    IJK::ARRAY<double> plane_grad_array(4*dimension);
    double * plane_grad[4];
    int igrad;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());
    compute_normalized_orthogonal_vector
      (dimension, dir1, normalized_dir0.PtrConst(), dir1_orth.Ptr());
    normalize_vector(dimension, dir2, normalized_dir2.Ptr());

    // Set the four plane gradients.
    for (int i = 0; i < 4; i++) 
      { plane_grad[i] = plane_grad_array.Ptr() + i*dimension; }
    IJK::add_coord(dimension, normalized_dir0.PtrConst(), 
                   dir1_orth.PtrConst(), plane_grad[0]);
    IJK::subtract_coord(dimension, dir1_orth.PtrConst(), 
                        normalized_dir0.PtrConst(),  plane_grad[1]);
    IJK::multiply_coord(dimension, -1, plane_grad[0], plane_grad[2]);
    IJK::multiply_coord(dimension, -1, plane_grad[1], plane_grad[3]);
    for (int i = 0; i < 4; i++) {
      normalize_vector(dimension, plane_grad[i]);
    }

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_planar_dist_L1
        (dimension, coord0, normalized_dir0.PtrConst(), dir1_orth.PtrConst(), 
         coord.PtrConst(), dist0, u0, u1, vdiff.Ptr());
      dist0 = dist0/std::sqrt(2.0);
      compute_dist2plane
        (dimension, coord0, normalized_dir2.PtrConst(), coord.PtrConst(), 
         dist1, v1.Ptr());

      // Select plane gradient.
      if (u0 >= 0) {
        if (u1 >= 0) { igrad = 0; }
        else { igrad = 3; }
      }
      else {
        if (u1 >= 0) { igrad = 1; }
        else { igrad = 2; }
      }

      normalize_vector(dimension, v1.Ptr());

      max_scalar_tie_zero
        (dimension, dist0, plane_grad[igrad], dist1-half_diff, v1.PtrConst(),
         distance, gradient.VectorPtr(iv));

      STYPE s = convert2type<STYPE>(distance);
      grid.Set(iv, s);
    }

  }

  /// Generate a scalar field whose isosurfaces are closed cones.
  /// @param dimension Volume dimension.
  /// @param apex0 Apex of cone with isovalue 0.
  /// @param axis_dir Cone axis direction.
  /// @param angle Cone angle.
  /// @param height0 Height of cone with isovalue 0.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE,
           typename ANGLE_TYPE, typename DIST_TYPE>
  void gen_gradient_closed_cone
  (const COORD_TYPE apex0[], const DIR_TYPE axis_dir[],
   const ANGLE_TYPE angle, const DIST_TYPE height0, 
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double angle_radians = (angle*M_PI)/180.0;
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_axis_dir(dimension);
    IJK::ARRAY<double> v0(dimension), v1(dimension);
    double dist0, dist1, x;

    if (angle <= 0) {
      gen_gradient_closed_cylinder
        (apex0, axis_dir, height0, grid, gradient);
      return;
    }

    normalize_vector(dimension, axis_dir, normalized_axis_dir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2cone_gradient
        (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
         angle_radians, dist0, v0.Ptr());
      compute_signed_dist2plane
        (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
         dist1, v1.Ptr());
      dist1 = -dist1;
      normalize_vector(dimension, v1.Ptr());

      max_scalar_tie_zero
        (dimension, dist0, v0.PtrConst(), dist1-height0, v1.PtrConst(),
         x, gradient.VectorPtr(iv));
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
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE, typename ANGLE_TYPE,
           typename DIST0_TYPE, typename DIST1_TYPE>
   void gen_gradient_frustrum
  (const COORD_TYPE apex0[], const DIR_TYPE axis_dir[],
   const ANGLE_TYPE angle, 
   const DIST0_TYPE dist2near0, const DIST1_TYPE dist2far0,
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double angle_radians = (angle*M_PI)/180.0;
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_axis_dir(dimension);
    IJK::ARRAY<double> v0(dimension), v1(dimension);
    IJK::ARRAY<double> reversed_axis_dir(dimension);
    double dist0, dist1, dist2, x;

    if (angle <= 0) {
      grid.SetAll(0);
      return;
    }

    normalize_vector(dimension, axis_dir, normalized_axis_dir.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2cone_gradient
        (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
         angle_radians, dist0, v0.Ptr());
      compute_signed_dist2plane
        (dimension, apex0, normalized_axis_dir.PtrConst(), 
         coord.PtrConst(), dist1);
      dist2 = -dist1;
 
     IJK::multiply_coord(dimension, -1, normalized_axis_dir.PtrConst(), 
                         reversed_axis_dir.Ptr());
     max_scalar_tie_zero
       (dimension, dist1+dist2near0, normalized_axis_dir.PtrConst(), 
        dist2-dist2far0, reversed_axis_dir.PtrConst(), x, v1.Ptr());

      max_scalar_tie_zero
        (dimension, x, v1.PtrConst(), dist0, v0.PtrConst(), x, 
         gradient.VectorPtr(iv));
      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
    }
  }

  /// Generate a scalar field whose isosurfaces bound a torus.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_TYPE, typename RADIUS_TYPE>
  void gen_gradient_torus
  (const COORD_TYPE coord0[], const DIR_TYPE dir0[],
   const RADIUS_TYPE radius, 
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<GTYPE> gradient(dimension);
    IJK::ARRAY<double> normalized_dir0(dimension);
    double dist0, dist1, x;

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());

    for (VTYPE iv = 0; iv < scalar_grid.NumVertices(); iv++) {
      scalar_grid.ComputeScaledCoord(iv, coord.Ptr());

      compute_dist2circle
        (dimension, coord0, normalized_dir0.PtrConst(), coord.PtrConst(), 
         radius, x, gradient.Ptr());
      normalize_vector(dimension, gradient.Ptr());

      STYPE s = convert2type<STYPE>(x);
    
      scalar_grid.Set(iv, s);
      gradient_grid.Set(iv, gradient.PtrConst());
    }

  }

  /// Generate a scalar field of maximum distance to planes.
  /// @param coord0[] All planes pass through coord0.
  /// @param normal[] normal[i*dim+j] is j'th coordinate of normal to plane i.
  ///        where dim is the grid dimension.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename CTYPE0, typename CTYPE1, typename ITYPE>
  void gen_gradient_max_dist2planes
  (const CTYPE0 coord0[], const CTYPE1 normal[], const ITYPE num_planes,
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE0> coord(dimension);
    IJK::ARRAY<double> vertex_gradient(dimension);
    IJK::ARRAY<CTYPE1> unit_normal(dimension*num_planes);
    double x, distance;

    if (num_planes < 1) {
      grid.SetAll(0);
      gradient.SetAllCoord(0);
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
      IJK::copy_coord
        (dimension, unit_normal.PtrConst(), vertex_gradient.Ptr());

      for (ITYPE i = 1; i < num_planes; i++) {
        compute_signed_dist2plane
          (dimension, coord0, unit_normal.PtrConst()+i*dimension, 
           coord.PtrConst(), distance);

        max_scalar_tie_zero
          (dimension, x, vertex_gradient.Ptr(), distance, 
           unit_normal.PtrConst()+i*dimension, x, vertex_gradient.Ptr());
      }

      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
      gradient.Set(iv, vertex_gradient.PtrConst());
    }
  }

  /// Generate a scalar field of minimum distance to planes.
  /// @param coord0[] All planes pass through coord0.
  /// @param normal[] normal[i*dim+j] is j'th coordinate of normal to plane i.
  ///        where dim is the grid dimension.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename CTYPE0, typename CTYPE1, typename ITYPE>
  void gen_gradient_min_dist2planes
  (const CTYPE0 coord0[], const CTYPE1 normal[], const ITYPE num_planes,
   SCALAR_GRID_TYPE & grid, GRADIENT_GRID_TYPE & gradient)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE0> coord(dimension);
    IJK::ARRAY<double> vertex_gradient(dimension);
	//DEBUG
    //CTYPE1 unit_normal[dimension*num_planes];
	 IJK::ARRAY<CTYPE1> unit_normal(dimension*num_planes);
    double x, distance;

    if (num_planes < 1) {
      grid.SetAll(0);
      gradient.SetAllCoord(0);
      return;
    }

    for (ITYPE i = 0; i < num_planes; i++) {
      normalize_vector
		//DEBUG
        //(dimension, normal+i*dimension, unit_normal+i*dimension);
		 (dimension, normal+i*dimension, unit_normal.Ptr()+i*dimension);
    }

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());
	  //DEBUG
      //compute_signed_dist2plane
       // (dimension, coord0, unit_normal, coord.PtrConst(), x);
		
		compute_signed_dist2plane
        (dimension, coord0, unit_normal.Ptr(), coord.PtrConst(), x);
      //DEBUG
	  //IJK::copy_coord(dimension, unit_normal, vertex_gradient.Ptr());
	  IJK::copy_coord(dimension, unit_normal.Ptr(), vertex_gradient.Ptr());

      for (ITYPE i = 1; i < num_planes; i++) {
        compute_signed_dist2plane
          
		  //DEBUG
		  //(dimension, coord0, unit_normal+i*dimension, coord.PtrConst(), 
          //distance);
		  (dimension, coord0, unit_normal.PtrConst()+i*dimension,
		   coord.PtrConst(), distance);

		//DEBUG
        //min_scalar_tie_zero
        //  (dimension, x, vertex_gradient.Ptr(), distance, 
        //   unit_normal+i*dimension, x, vertex_gradient.Ptr());
		min_scalar_tie_zero
          (dimension, x, vertex_gradient.Ptr(), distance, 
           unit_normal.PtrConst()+i*dimension, x, vertex_gradient.Ptr());
      }

      STYPE s = convert2type<STYPE>(x);

      grid.Set(iv, s);
      gradient.Set(iv, vertex_gradient.PtrConst());
    }
  }

  /// Generate a scalar field with constant unit magnitude gradient
  /// Scalar value at coord0[] is 0.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE,
           typename COORD_TYPE, typename DIR_COORD_TYPE>
  void gen_constant_unit_gradient
  (const COORD_TYPE coord0[], const DIR_COORD_TYPE dir0[],
   SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE 
      VCOORD_TYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<VCOORD_TYPE> normalized_dir0(dimension);

    gen_constant_unit_gradient(coord0, dir0, scalar_grid);

    normalize_vector(dimension, dir0, normalized_dir0.Ptr());
    gradient_grid.SetAll(normalized_dir0.PtrConst());
  }

  // **********************************************************************
  // Transform and combine gradient fields.
  // **********************************************************************

  /// Select gradient from gradientA or gradientB
  ///   depending upon whether gridC matches gridA or gridB
  ///   Set gradient to zero if both match.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  void select_gradient_tie_zero
  (const SCALAR_GRID_TYPE & gridA, const GRADIENT_GRID_TYPE & gradientA,
   const SCALAR_GRID_TYPE & gridB, const GRADIENT_GRID_TYPE & gradientB,
   const SCALAR_GRID_TYPE & gridC, GRADIENT_GRID_TYPE & gradientC)
  {
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = gridA.Dimension();
    VTYPE numv = gridA.NumVertices();
    IJK::PROCEDURE_ERROR error("select_gradient");

    if (!gridB.Check(gridA, error)) { throw error; };
    if (!gridC.Check(gridA, error)) { throw error; };
    if (!gradientB.Check(gradientA, error)) { throw error; };
    if (!check_gradient_grid(gridA, gradientA, error)) 
      { throw error; }

    /// *** DOES THIS MAKE SENSE? ***
    gradientC.SetSize(gradientC);

    for (VTYPE iv = 0; iv < numv; iv++) {

      if (gridA.Scalar(iv) == gridC.Scalar(iv)) {

        if (gridB.Scalar(iv) == gridC.Scalar(iv)) {
          // Both gridA and gridB match.  Set gradient to zero.
          IJK::set_coord
            (dimension, 0, gradientC.VectorPtr(iv));
        }
        else {
          gradientC.Set(iv, gradientA.VectorPtrConst(iv));
        }

      }
      else if (gridB.Scalar(iv) == gridC.Scalar(iv)) {
        gradientC.Set(iv, gradientB.VectorPtrConst(iv));
      }
      else {
        IJK::set_coord
          (dimension, 0, gradientC.VectorPtr(iv));
      }
    }

  }

  /// Select gradient from gradientA or gradientB
  ///   depending upon whether gridC matches gridA or gridB
  ///   Set gradient to gradientA if both match.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  void select_gradient_tie_useA
  (const SCALAR_GRID_TYPE & gridA, const GRADIENT_GRID_TYPE & gradientA,
   const SCALAR_GRID_TYPE & gridB, const GRADIENT_GRID_TYPE & gradientB,
   const SCALAR_GRID_TYPE & gridC, GRADIENT_GRID_TYPE & gradientC)
  {
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRADIENT_GRID_TYPE::SCALAR_TYPE GTYPE;

    DTYPE dimension = gridA.Dimension();
    VTYPE numv = gridA.NumVertices();
    IJK::PROCEDURE_ERROR error("select_gradient");

    if (!gridB.Check(gridA, error)) { throw error; };
    if (!gridC.Check(gridA, error)) { throw error; };
    if (!gradientB.Check(gradientA, error)) { throw error; };
    if (!check_gradient_grid(gridA, gradientA, error)) 
      { throw error; }

    /// *** DOES THIS MAKE SENSE? ***
    gradientC.SetSize(gradientC);

    for (VTYPE iv = 0; iv < numv; iv++) {

      if (gridA.Scalar(iv) == gridC.Scalar(iv)) {
        gradientC.Set(iv, gradientA.VectorPtrConst(iv));
      }
      else if (gridB.Scalar(iv) == gridC.Scalar(iv)) {
        gradientC.Set(iv, gradientB.VectorPtrConst(iv));
      }
      else {
        IJK::set_coord
          (dimension, 0, gradientC.VectorPtr(iv));
      }

    }

  }

  /// Compute min of two scalar fields.
  /// Select gradient based on min.
  /// @param tie_zero If true, set gradient to zero wherever scalars are equal.
  ///            If false, select gradient from A wherever scalars are equal.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  void min_scalar_select_gradient
  (const SCALAR_GRID_TYPE & gridA, const GRADIENT_GRID_TYPE & gradientA,
   const SCALAR_GRID_TYPE & gridB, const GRADIENT_GRID_TYPE & gradientB,
   const bool tie_zero,
   SCALAR_GRID_TYPE & gridC, GRADIENT_GRID_TYPE & gradientC)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRADIENT_GRID_TYPE::VECTOR_COORD_TYPE GCTYPE;

    const DTYPE dimension = gridA.Dimension();
    const VTYPE numv = gridA.NumVertices();
    IJK::PROCEDURE_ERROR error("min_scalar_select_gradient");

    if (!gridB.Check(gridA, error)) { throw error; };
    if (!gradientB.Check(gradientA, error)) { throw error; };
    if (!check_gradient_grid(gridA, gradientA, error)) 
      { throw error; }

    gridC.SetSize(gridA);
    gradientC.SetSize(gradientA);

    if (tie_zero) {
      STYPE s;

      for (VTYPE iv = 0; iv < numv; iv++) {

        min_scalar_tie_zero
          (dimension, gridA.Scalar(iv), gradientA.VectorPtrConst(iv),
           gridB.Scalar(iv), gradientB.VectorPtrConst(iv),
           s, gradientC.VectorPtr(iv));
        gridC.Set(iv, s);
      }
    }
    else {
      STYPE s;

      for (VTYPE iv = 0; iv < numv; iv++) {

        min_scalar_tie_use_v0
          (dimension, gridA.Scalar(iv), gradientA.VectorPtrConst(iv),
           gridB.Scalar(iv), gradientB.VectorPtrConst(iv),
           s, gradientC.VectorPtr(iv));
        gridC.Set(iv, s);
      }
    }

  }

  /// Compute max of two scalar fields.
  /// Select gradient based on max.
  /// @param tie_zero If true, set gradient to zero wherever scalars are equal.
  ///            If false, select gradient from A wherever scalars are equal.
  template<typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  void max_scalar_select_gradient
  (const SCALAR_GRID_TYPE & gridA, const GRADIENT_GRID_TYPE & gradientA,
   const SCALAR_GRID_TYPE & gridB, const GRADIENT_GRID_TYPE & gradientB,
   const bool tie_zero,
   SCALAR_GRID_TYPE & gridC, GRADIENT_GRID_TYPE & gradientC)
  {
    typedef typename SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SCALAR_GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRADIENT_GRID_TYPE::VECTOR_COORD_TYPE GCTYPE;

    const DTYPE dimension = gridA.Dimension();
    const VTYPE numv = gridA.NumVertices();
    IJK::PROCEDURE_ERROR error("max_scalar_select_gradient");

    if (!gridB.Check(gridA, error)) { throw error; };
    if (!gradientB.Check(gradientA, error)) { throw error; };
    if (!check_gradient_grid(gridA, gradientA, error)) 
      { throw error; }

    gridC.SetSize(gridA);
    gradientC.SetSize(gradientA);

    if (tie_zero) {
      STYPE s;

      for (VTYPE iv = 0; iv < numv; iv++) {

        max_scalar_tie_zero
          (dimension, gridA.Scalar(iv), gradientA.VectorPtrConst(iv),
           gridB.Scalar(iv), gradientB.VectorPtrConst(iv),
           s, gradientC.VectorPtr(iv));
        gridC.Set(iv, s);
      }
    }
    else {
      STYPE s;

      for (VTYPE iv = 0; iv < numv; iv++) {

        max_scalar_tie_use_v0
          (dimension, gridA.Scalar(iv), gradientA.VectorPtrConst(iv),
           gridB.Scalar(iv), gradientB.VectorPtrConst(iv),
           s, gradientC.VectorPtr(iv));
        gridC.Set(iv, s);
      }
    }

  }

};

#endif

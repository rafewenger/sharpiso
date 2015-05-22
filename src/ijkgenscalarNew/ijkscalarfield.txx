/// \file ijkscalarfield.txx
/// scalar field generation and manipulation routines

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2015 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

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

#include <cmath>
#include <cstdlib>

#include "ijkscalar_grid.txx"
#include "ijkcoord.txx"
#include "ijkgencoord.txx"

/// Routines for generating scalar fields.
namespace IJKSCALARFIELD {

  using namespace IJKGENCOORD;

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

  /// Generate a scalar field whose isosurfaces are closed cones
  ///   with smooth tips.
  /// @param dimension Volume dimension.
  /// @param apex0 Apex of cone with isovalue 0.
  /// @param axis_dir Cone axis direction.
  /// @param angle Cone angle.
  /// @param height0 Height of cone with isovalue 0.
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename ANGLE_TYPE, typename DIST_TYPE>
   void gen_closed_cone_smooth_tip
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

      compute_dist2cone_smooth_tip
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

  /// Generate a scalar field whose isosurfaces are cannon shaped.
  /// Cannon shape is the union of a frustrum and a ball.
  /// @param dimension Volume dimension.
  /// @param apex0 Apex of cone with isovalue 0.
  /// @param axis_dir Cone axis direction.
  /// @param angle Cone angle.
  /// @param dist2near0 Distance from cone apex to near plane with isovalue 0.
  /// @param dist2center Distance from cone apex to ball center.
  /// @pre (dist2center-dist2near) > dist2center * sin(angle).
  template<typename GRID_TYPE, typename COORD_TYPE, typename DIR_TYPE,
           typename ANGLE_TYPE, typename DIST0_TYPE, typename DIST1_TYPE>
   void gen_cannon
  (const COORD_TYPE apex0[], const DIR_TYPE axis_dir[],
   const ANGLE_TYPE angle, 
   const DIST0_TYPE dist2near0, const DIST1_TYPE dist2center,
   GRID_TYPE & grid)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();
    const double angle_radians = (angle*M_PI)/180.0;
    const double sin_angle = std::sin(angle_radians);
    IJK::ARRAY<COORD_TYPE> coord(dimension);
    IJK::ARRAY<double> normalized_axis_dir(dimension);
    IJK::ARRAY<double> u(dimension);
    double u_length;
    double dist0, dist1, x, cos_axis_u;

    if (angle <= 0) {
      grid.SetAll(0);
      return;
    }

    normalize_vector(dimension, axis_dir, normalized_axis_dir.Ptr());

    IJK::ARRAY<double> ball_center(dimension);
    IJK::add_scaled_coord
      (dimension, -dist2center, normalized_axis_dir.PtrConst(), apex0,
       ball_center.Ptr());

    for (VTYPE iv = 0; iv < grid.NumVertices(); iv++) {
      grid.ComputeScaledCoord(iv, coord.Ptr());

      IJK::subtract_coord
        (dimension, coord.PtrConst(), ball_center.PtrConst(), u.Ptr());
      IJK::compute_magnitude(dimension, u.PtrConst(), u_length);
      normalize_vector(dimension, u.Ptr(), u.Ptr());
      IJK::compute_inner_product
        (dimension, u.PtrConst(), normalized_axis_dir.PtrConst(), cos_axis_u);

      if (cos_axis_u >= sin_angle) {
        compute_dist2cone
          (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
           angle_radians, dist0);
        compute_signed_dist2plane
          (dimension, apex0, normalized_axis_dir.PtrConst(), coord.PtrConst(), 
           dist1);
        x = std::max(dist0, dist1+dist2near0);
      }
      else {
        x = u_length - dist2center*sin_angle;
      }

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

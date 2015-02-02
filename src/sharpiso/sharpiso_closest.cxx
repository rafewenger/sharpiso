/// \file sharpiso_closest.cxx
/// Compute closest point to a given point.
/// Version v0.1.1


#include <iostream>
#include <cmath>
#include <limits>

#include "ijkcoord.txx"
#include "ijkgrid.txx"
#include "ijkinterpolate.txx"
#include "sharpiso_linear_alg.txx"

#include "sharpiso_closest.h"

using namespace std;
using namespace SHARPISO;

///  Compute the closest point to point p on a given line.
void SHARPISO::compute_closest_point_on_line
(const COORD_TYPE point[DIM3],
 const COORD_TYPE line_origin[DIM3],
 const COORD_TYPE line_direction[DIM3],
 const GRADIENT_COORD_TYPE zero_tolerance_squared,
 COORD_TYPE closest_point[DIM3])
{
  COORD_TYPE a[DIM3];
  COORD_TYPE c[DIM3];
  COORD_TYPE line_direction_normalized[DIM3];
  float dotpdt;
  bool is_zero = false;

  // Normalize the line direction
  normalize_vector_3D
    (line_direction, line_direction_normalized, 
     is_zero, zero_tolerance_squared);
  
  // the vector perpendicular is given by
  // c = a - (a dot b^)* b^
  //set a and b is given  by the ray direction
  IJK::subtract_coord(DIM3, point, line_origin, a);
  
  // dot product
  IJK::compute_inner_product(DIM3, a, line_direction_normalized, dotpdt);
  
  // compute c and invert it at same time
  for (int d=0; d<DIM3; d++)
    { c[d] = -(a[d]-dotpdt*line_direction_normalized[d]); }
  
  IJK::add_coord(DIM3, point, c, closest_point);
}


void update_t_minus
(const int X,const int Y, const COORD_TYPE ray_direction_normalized[DIM3],
 const COORD_TYPE point[DIM3],const COORD_TYPE line_origin[DIM3], 
 COORD_TYPE &t, bool &istTrue)
{
  float RxMinusRy=ray_direction_normalized[X]-ray_direction_normalized[Y];
  if (RxMinusRy!=0) {
    t = ((line_origin[Y]-point[Y]) - (line_origin[X]-point[X]))/ RxMinusRy;
    istTrue = true;
  }
  else
    istTrue=false;
}

void update_t_plus
(const int X,const int Y, const COORD_TYPE ray_direction_normalized[DIM3],
 const COORD_TYPE point[DIM3], const COORD_TYPE line_origin[DIM3], 
 COORD_TYPE &t, bool &istTrue)
{
  float RxPlusRy=(ray_direction_normalized[X]+ ray_direction_normalized[Y]);
  if (RxPlusRy !=0) {
    t =((point[X]-line_origin[X])-(line_origin[Y]-point[Y]))/ RxPlusRy;
    istTrue=true;
  }
  else
    istTrue=false;
  
}

void compute_linf_point
(const COORD_TYPE coord[DIM3],const COORD_TYPE ray_direction_normalized[DIM3],
 const float t_i, COORD_TYPE point[DIM3])
{
  for (int j=0; j<DIM3; j++)
    point[j]=coord[j]+t_i*ray_direction_normalized[j];
  
}

void compute_linf_dist
(const COORD_TYPE pointA[DIM3],
 const COORD_TYPE pointB[DIM3], float & new_dist)
{
  float max =abs(pointA[0]-pointB[0]);
  for (int d=0; d<DIM3; d++) {
    if (abs(pointA[d]-pointB[d]) > max) {
      max=abs(pointA[d]-pointB[d]) ;
    }
  }
  new_dist=max;
}

///  Compute the closest point under the Linf metric to point p on a given line.
void SHARPISO::compute_closest_point_on_line_linf
(const COORD_TYPE point[DIM3],
 const COORD_TYPE line_origin[DIM3],
 const COORD_TYPE line_direction[DIM3],
 const GRADIENT_COORD_TYPE zero_tolerance_squared,
 COORD_TYPE closest_point[DIM3])
{
  COORD_TYPE line_direction_normalized[DIM3];
  bool is_zero = false;

  // initialize closest_point to line_origin
  IJK::copy_coord(DIM3, line_origin, closest_point);

  // Normalize the line direction
  normalize_vector_3D
    (line_direction, line_direction_normalized, 
     is_zero, zero_tolerance_squared);

  int X=0,Y=1,Z=2;
  // store the values of t (2 for each)
  float t[6]={0.0};
  bool  istTrue[6];
  
  update_t_minus(X, Y, line_direction_normalized, point, 
                 line_origin, t[0], istTrue[0]);
  update_t_plus(X, Y, line_direction_normalized, point, 
                line_origin, t[1], istTrue[1]);
  
  update_t_minus(X, Z, line_direction_normalized, point, 
                 line_origin, t[2], istTrue[2]);
  update_t_plus(X, Z, line_direction_normalized, point, 
                line_origin, t[3], istTrue[3]);
  
  update_t_minus( Y,Z, line_direction_normalized, point, 
                  line_origin, t[4], istTrue[4]);
  update_t_plus( Y,Z, line_direction_normalized, point, 
                 line_origin, t[5], istTrue[5]);

  COORD_TYPE point_on_line[DIM3];
  COORD_TYPE dist(0);
  bool is_dist_set = false;
  COORD_TYPE new_dist;
  for (int i=0; i<6; i++) {
    if (istTrue[i]) {
      compute_linf_point
        (line_origin, line_direction_normalized, t[i], point_on_line);
      compute_linf_dist(point, point_on_line, new_dist); 
      if (!is_dist_set || (new_dist < dist)) {
        dist = new_dist;
        IJK::copy_coord_3D(point_on_line, closest_point);
        is_dist_set = true;
      }
    }
  }

}

///  Compute the closest point under the Linf metric to point p on a given line.
void SHARPISO::compute_closest_point_on_line_unscaled_linf
(const COORD_TYPE point[DIM3],
 const COORD_TYPE line_origin[DIM3],
 const COORD_TYPE line_direction[DIM3],
 const GRADIENT_COORD_TYPE zero_tolerance_squared,
 const COORD_TYPE spacing[DIM3],
 COORD_TYPE closest_point[DIM3])
{
  COORD_TYPE point2[DIM3], line_origin2[DIM3], line_direction2[DIM3];

  // Unscale point and line.
  for (int d = 0; d < DIM3; d++) {
    point2[d] = point[d]/spacing[d];
    line_origin2[d] = line_origin[d]/spacing[d];
    line_direction2[d] = line_direction[d]/spacing[d];
  }

  compute_closest_point_on_line_linf
    (point2, line_origin2, line_direction2, zero_tolerance_squared, 
     closest_point);

  // Rescale closest point.
  IJK::scale_coord(DIM3, spacing, closest_point, closest_point);
}

// compute the linf distance between a point and a ray
void SHARPISO::compute_closest_point_to_cube_center_linf
(const GRID_COORD_TYPE cube_coord[],
 const COORD_TYPE coord[],
 const COORD_TYPE ray_direction[],
 COORD_TYPE closest_point[DIM3])
{
  const GRADIENT_COORD_TYPE ZERO_TOLERANCE_SQUARED = 0.0001;

  // find the cube center
  COORD_TYPE center_offset[DIM3] = {0.5, 0.5, 0.5};
  COORD_TYPE cube_center[DIM3];
  IJK::add_coord(DIM3, center_offset, cube_coord, cube_center);

  compute_closest_point_on_line_linf
    (cube_center, coord, ray_direction, ZERO_TOLERANCE_SQUARED,
     closest_point);
}


// Compute the point on a plane closest (L2) to point p.
// @pre orth_dir[] is a unit vector.
void SHARPISO::compute_closest_point_on_plane
(const COORD_TYPE pcoord[DIM3],
 const COORD_TYPE qcoord[DIM3],
 const COORD_TYPE orth_dir[DIM3],
 COORD_TYPE closest_point[DIM3])
{
  COORD_TYPE v_diff[DIM3], v_proj[DIM3];

  IJK::subtract_coord_3D(qcoord, pcoord, v_diff);
  IJK::project_vector(DIM3, v_diff, orth_dir, v_proj);
  IJK::add_coord_3D(pcoord, v_proj, closest_point);
}

// **************************************************
// INTERSECT LINE AND SQUARE CYLINDER
// **************************************************

namespace {

  // Compute the intersection of a line and two parallel planes.
  // Returns t[0], t[1], where the intersection is
  //   (line_p0 + line_dir*t[0], line_p0 + line_dir*t[1])
  //   and t[0] <= t[1].
  void intersect_line_two_planes_t
  (const COORD_TYPE plane_p0[DIM3], 
   const NUM_TYPE plane_orth_dir,
   const COORD_TYPE distance_between_planes,
   const COORD_TYPE line_p0[DIM3],
   const COORD_TYPE line_dir[DIM3],
   const COORD_TYPE max_small_grad,
   COORD_TYPE t[2],
   bool & flag_intersect)
  {
    const NUM_TYPE dir = plane_orth_dir;

    if (-max_small_grad <= line_dir[plane_orth_dir] &&
        line_dir[plane_orth_dir] <= max_small_grad) {
      flag_intersect = false;
      return;
    }
    
    flag_intersect = true;
    COORD_TYPE point_diff = plane_p0[dir] - line_p0[dir];
    t[0] = point_diff/line_dir[dir];
    t[1] = (point_diff + distance_between_planes)/line_dir[dir];

    if (t[1] < t[0]) { std::swap(t[0], t[1]); };
  }

}

/// Intersect a line and a square cylinder
void SHARPISO::intersect_line_square_cylinder
(const COORD_TYPE square_p0[DIM3],
 const NUM_TYPE cylinder_axis_dir,
 const COORD_TYPE square_width,
 const COORD_TYPE half_cylinder_length,
 const COORD_TYPE line_p0[DIM3],
 const COORD_TYPE line_dir[DIM3],
 const COORD_TYPE max_small_grad,
 COORD_TYPE end[2][DIM3],
 bool & flag_intersect)
{
  COORD_TYPE t[2], new_t[2];
  COORD_TYPE cylinder_p0[DIM3];

  IJK::copy_coord_3D(square_p0, cylinder_p0);
  cylinder_p0[cylinder_axis_dir] -= half_cylinder_length;

  flag_intersect = false;
  for (NUM_TYPE d = 0; d < DIM3; d++) {

    bool flag;
    COORD_TYPE distance_between_planes;

    if (d == cylinder_axis_dir) 
      { distance_between_planes = 2*half_cylinder_length; }
    else
      { distance_between_planes = square_width; }

    intersect_line_two_planes_t
      (cylinder_p0, d, distance_between_planes, line_p0, line_dir,
       max_small_grad, new_t, flag);

    if (flag) {
      if (!flag_intersect) {
        t[0] = new_t[0];
        t[1] = new_t[1];
        flag_intersect = true;
      }
      else {
        if (t[0] < new_t[0]) { t[0] = new_t[0]; }
        if (t[1] > new_t[1]) { t[1] = new_t[1]; }
      }
    }
    else {
      if (line_p0[d] < cylinder_p0[d] ||
          line_p0[d] > cylinder_p0[d] + distance_between_planes) {
        // Line is not between planes.
        flag_intersect = false;
        return;
      }
    }
  }

  if (!flag_intersect) { return; }

  if (t[0] > t[1]) {
    flag_intersect = false;
    return;
  }


  for (NUM_TYPE i = 0; i < 2; i++) {
    for (NUM_TYPE d = 0; d < DIM3; d++) {
      end[i][d] = line_p0[d] + t[i]*line_dir[d];
    }
  }

}


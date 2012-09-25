// FindIntersect.h
// Find intersection of a ray and a cube.


#ifndef _SHARPISO_FINDINTERSECT_
#define _SHARPISO_FINDINTERSECT_

#include<iostream>
#include <Eigen/Dense>
#include "sharpiso_types.h"
#include "sharpiso_eigen.h"

namespace SHARPISO {
  
  const float THRESHOLD_CLAMP = 0.001; // clamp threshold
  const float EPSILON = 0.0001;
  
  /// Compute the midpoint of the intersection of the ray and the cube.
  /// Return true if the ray intersects the cube
  bool calculate_point_intersect
  (const COORD_TYPE * point, const COORD_TYPE *dir, COORD_TYPE *intersect);
  
  /// Separate version which translate back and forth
  ///   to find the intersection with the unit cube.
  bool calculate_point_intersect
  (const GRID_COORD_TYPE cube_coord[], const COORD_TYPE *p,
   const COORD_TYPE *dir, COORD_TYPE *intersect);
  
  // Calculate the intersection with an enlarged cube.
  bool calculate_point_intersect_complex
  (const GRID_COORD_TYPE cube_coord[], const COORD_TYPE *original_pt,
   const COORD_TYPE *dir, const float th,  COORD_TYPE *intersect);
  
  
  // Check if the ray intersects the cube (No intersection calculated)
  bool calculate_point_intersect_complex
  (const GRID_COORD_TYPE cube_coord[],
   const COORD_TYPE *original_pt,
   const COORD_TYPE *dir,
   const float th
   );
  
  ///  Compute the closest point to point p on a given line.
  void compute_closest_point_on_line
  (const COORD_TYPE point[DIM3],
   const COORD_TYPE line_origin[DIM3],
   const COORD_TYPE line_direction[DIM3],
   const GRADIENT_COORD_TYPE zero_tolerance_squared,
   COORD_TYPE closest_point[DIM3]);

  // Calculate the closest point to a cube center
  void compute_closest_point_to_cube_center
  (const GRID_COORD_TYPE cube_coord[DIM3],
   const COORD_TYPE coord[DIM3],
   const COORD_TYPE ray_direction[DIM3],
   const GRADIENT_COORD_TYPE zero_tolerance_squared,
   COORD_TYPE closest_point[DIM3]);
  
  // compute the linf dis
  void compute_closest_point_to_cube_center_linf
  (const GRID_COORD_TYPE cube_coord[],
   const COORD_TYPE coord[],
   const COORD_TYPE ray_direction[],
   COORD_TYPE closest_point[DIM3]);

  /// Intersect a line and a square cylinder
  /// @param[out] end[][] end[i][d] is coordinate d of endpoint i.
  /// @param[flag_intersect] True if line intersects cylinder.
  void intersect_line_square_cylinder
  (const COORD_TYPE square_p0[DIM3],
   const NUM_TYPE cylinder_axis_dir,
   const COORD_TYPE square_width,
   const COORD_TYPE cylinder_length,
   const COORD_TYPE line_p0[DIM3],
   const COORD_TYPE line_dir[DIM3],
   const COORD_TYPE max_small_grad,
   COORD_TYPE end[2][DIM3],
   bool & flag_intersect);

};

#endif

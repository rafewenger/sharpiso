/// \file sharpiso_types.h
/// Definitions for sharp isosurface processing.
/// Version v0.1.2

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011-2015 Rephael Wenger

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

#ifndef _SHARPISO_TYPES_
#define _SHARPISO_TYPES_

#include <utility>

/// Definitions for sharp isosurface processing.
namespace SHARPISO {

  // **************************************************
  // TYPES
  // **************************************************

  typedef float SCALAR_TYPE;
  typedef float COORD_TYPE;
  typedef float GRADIENT_COORD_TYPE;
  typedef int GRID_COORD_TYPE;
  typedef GRADIENT_COORD_TYPE ANGLE_TYPE;
  typedef double WEIGHT_TYPE;

  typedef int AXIS_SIZE_TYPE;
  typedef int GRADIENT_LENGTH_TYPE;
  typedef unsigned int BOUNDARY_BITS_TYPE;

  /// Coordinate type which is guaranteed to be signed.
  /// (Note: COORD_TYPE may be signed or unsigned.)
  typedef float SIGNED_COORD_TYPE;

  typedef int NUM_TYPE;
  typedef int VERTEX_INDEX;

  /// Index and grid coord difference types.  Must be signed.
  typedef int INDEX_DIFF_TYPE;
  typedef int GRID_COORD_DIFF_TYPE;

  typedef std::pair<VERTEX_INDEX, VERTEX_INDEX> VERTEX_PAIR;

  typedef enum { LOC_NONE, CENTROID, CUBE_CENTER, LOC_SVD } LOC_TYPE ;

  typedef enum { GRAD_3x3x3, GRAD_5x5x5, GRAD_7x7x7, GRAD_9x9x9,
                 GRAD_C, GRAD_N, GRAD_CS, GRAD_NS, GRAD_XS, GRAD_IE, GRAD_IES,
                 GRAD_IE_DIR, GRAD_CD, GRAD_CD_DUP, GRAD_NIE, GRAD_NIES,
                 GRAD_BIES,
                 GRAD_BIES_OLD, // Old, deprecated version.
                 GRAD_EDGEI_INTERPOLATE, GRAD_EDGEI_GRAD,
                 UNKNOWN_GRAD_SELECTION_METHOD }
    GRAD_SELECTION_METHOD;

  /// Isosurface vertex types
  typedef enum 
  { CORNER_ISOVERT, EDGE_ISOVERT, SHARP_ISOVERT, SMOOTH_ISOVERT, ALL_ISOVERT }
  ISOVERT_TYPE;

  /// Cube adjacency types
  typedef enum
  { FACET_ADJACENT, EDGE_ADJACENT, VERTEX_ADJACENT, NO_ADJACENT_CUBE }
  CUBE_ADJACENCY_TYPE;


  // **************************************************
  // CONSTANTS
  // **************************************************

  const NUM_TYPE DIM3 = 3;
  const NUM_TYPE NUM_CUBE_VERTICES3D = 8;
  const NUM_TYPE NUM_CUBE_EDGES3D = NUM_CUBE_VERTICES3D*DIM3/2;
  const NUM_TYPE NUM_CUBE_FACETS3D = 2*DIM3;
  const NUM_TYPE NUM_CUBE_FACET_VERTICES3D = NUM_CUBE_VERTICES3D/2;
  const NUM_TYPE NUM_CUBE_FACET_EDGES3D =
    (NUM_CUBE_EDGES3D-NUM_CUBE_FACET_VERTICES3D)/2;
  const NUM_TYPE NUM_CUBE_DIAGONALS3D = NUM_CUBE_FACET_VERTICES3D;
  const NUM_TYPE NUM_CUBE_NEIGHBORS3D = 
    (NUM_CUBE_VERTICES3D+NUM_CUBE_EDGES3D+NUM_CUBE_FACETS3D);

  // *** DEPRECATED ***
  /// Number of quadrilateral vertices.
  const NUM_TYPE NUM_QUAD_VERTICES = 4;

  /// Number of quadrilateral vertices.
  const NUM_TYPE NUM_VERT_PER_TRI = 3;

  /// Number of quadrilateral vertices.
  const NUM_TYPE NUM_VERT_PER_QUAD = 4;

  /// Number of vertices in two cubes sharing a facet.
  const NUM_TYPE NUM_TWO_CUBE_VERTICES3D = DIM3 * NUM_CUBE_VERTICES3D;


  // *** NOTE:  SHOULD BE REPLACED BY PARAMETER ***
  // clamp very small values to the cube.
  const SCALAR_TYPE clamp_threshold = 0.0001; 

};

#endif

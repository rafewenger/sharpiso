/// \file sharpiso_types.h
/// Definitions for sharp isosurface processing.
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011,2012 Rephael Wenger

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

  // Coordinate type which is guaranteed to be signed.
  // Note: COORD_TYPE may be signed or unsigned.
  typedef float SIGNED_COORD_TYPE;

  typedef int NUM_TYPE;
  typedef int VERTEX_INDEX;

  typedef enum { LOC_NONE, CENTROID, CUBE_CENTER, LOC_SVD } LOC_TYPE ;

  typedef enum { GRAD_C, GRAD_N, GRAD_CS, GRAD_NS, GRAD_IE, GRAD_IES,
                 GRAD_IE_DIR, GRAD_CD, GRAD_CD_DUP, GRAD_NIE, GRAD_NIES,
                 GRAD_EDGEI_INTERPOLATE, GRAD_EDGEI_GRAD,
                 UNKNOWN_GRAD_SELECTION_METHOD }
    GRAD_SELECTION_METHOD;


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

  /// Number of vertices in two cubes sharing a facet.
  const NUM_TYPE NUM_TWO_CUBE_VERTICES3D = DIM3 * NUM_CUBE_VERTICES3D;


  // *** NOTE:  SHOULD BE REPLACED BY PARAMETER ***
  // clamp very small values to the cube.
  const SCALAR_TYPE clamp_threshold = 0.0001; 

};

#endif

/// \file isodual3D_types.h
/// Type definitions for isodual3D.

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

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

#ifndef _ISODUAL3D_TYPES_
#define _ISODUAL3D_TYPES_

#include <string>

#include "ijk.txx"
#include "ijkisopoly.txx"
#include "ijkscalar_grid.txx"
#include "ijkmerge.txx"
#include "sharpiso_types.h"

using namespace SHARPISO;
namespace ISODUAL3D {


// **************************************************
// SCALAR TYPES
// **************************************************

  typedef float GRADIENT_TYPE;   ///< Gradient coordinate type.
  typedef int LENGTH_TYPE;       ///< Vector length type.
  typedef int ISO_VERTEX_INDEX;  ///< Isosurface vertex index type.
  typedef int MERGE_INDEX;       ///< Merge index type.

  /// Edge index type.
  /// Vertex and edge indices must have the same type.
  typedef VERTEX_INDEX EDGE_INDEX;

  /// Facet vertex index type.
  typedef unsigned char FACET_VERTEX_INDEX;

// **************************************************
// ARRAY TYPES
// **************************************************

  typedef std::vector<COORD_TYPE> COORD_ARRAY;   ///< Grid coordinate array.
  typedef std::vector<VERTEX_INDEX>              /// Vertex index array.
    VERTEX_INDEX_ARRAY;
  typedef std::vector<SCALAR_TYPE> SCALAR_ARRAY; ///< Scalar array.

// **************************************************
// ENUMERATED TYPES
// **************************************************

  /// Interpolation type.
  /// - LINEAR_INTERPOLATION: Determine the location of an isosurface vertex
  ///   using linear interpolation on the endpoints of the cube, pyramid
  ///   or simplex containing the isosurface vertex.
  /// - MULTILINEAR_INTERPOLATION: Determine the location of an
  ///   isosurface vertex using multilinear interpolation on the cube vertices.
  typedef enum { LINEAR_INTERPOLATION, MULTILINEAR_INTERPOLATION }
  INTERPOLATION_TYPE;

  /// Isosurface vertex position method.
  /// CUBE_CENTER: Position isosurface vertices at cube centers.
  /// CENTROID_EDGE_ISO: Position isosurface vertices at the centroid
  ///                    of the edge isosurface intersections.
  /// GRADIENT_POSITIONING: Position using gradients.
  /// EDGEI_INTERPOLATE: Interpolate isosurface-edge intersection.
  /// EDGEI_GRADIENT:  Compute isosurface-edge intersection using gradients.
  /// EDGEI_INPUT_DATA: Isosurface-edge intersections and normals are
  ///                   provided as input data.
  typedef enum { CUBECENTER, CENTROID_EDGE_ISO,
                 GRADIENT_POSITIONING, EDGEI_INTERPOLATE, EDGEI_GRADIENT,
                 EDGEI_INPUT_DATA } 
  VERTEX_POSITION_METHOD;

  /// Quadrilateral triangulation method.
  typedef enum { UNDEFINED_TRI, UNIFORM_TRI, SPLIT_MAX_ANGLE }
  QUAD_TRI_METHOD;

  /// Ambiguity status.
  typedef enum { AMBIGUITY_NOT_SET, NOT_AMBIGUOUS, SEPARATE_POS, 
                 SEPARATE_NEG, UNDECIDED_AMBIGUITY, 
                 CONFLICTING_SEPARATION } 
  AMBIGUITY_STATUS;

  typedef unsigned char AMBIGUITY_TYPE;

// **************************************************
// CLASSES
// **************************************************

  typedef IJK::BOX<VERTEX_INDEX> GRID_BOX;  ///< Grid box type.

  typedef IJK::CUBE_FACE_INFO<int,int,int> ISODUAL3D_CUBE_FACE_INFO;

}

#endif

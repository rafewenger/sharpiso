/// \file ijkgenscalar.h
/// generate a scalar field
/// Version v0.1.3

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2013 Rephael Wenger

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

/*!
  \mainpage IJKGENSCALAR: Generate scalar grid.

  IJKGENSCALAR is a program for generating a regular scalar grid
  representing various scalar fields.  It can also be used to generate
  gradient vectors at the grid vertices.
*/

#ifndef _IJKGENSCALAR_
#define _IJKGENSCALAR_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"

#include "ijkgenscalar.txx"


namespace IJKGENSCALAR {

  // ********************************************************
  // TYPES
  // ********************************************************

  typedef size_t AXIS_SIZE_TYPE;
  typedef float SCALAR_TYPE;
  typedef float GRADIENT_COORD_TYPE;
  typedef float COORD_TYPE;
  typedef float RADIUS_TYPE;
  typedef float DIFF_TYPE;
  typedef float ANGLE_TYPE;
  typedef int NUM_TYPE;
  typedef unsigned long ISOTABLE_INDEX_TYPE;
  typedef unsigned int SEED_TYPE;
  typedef unsigned int MIN_MAX_TYPE;
  typedef IJK::GRID_NEIGHBORS<int, AXIS_SIZE_TYPE, int, int, int> BASE_GRID;
  typedef IJK::GRID_SPACING<COORD_TYPE, BASE_GRID> GRID;
  typedef IJK::SCALAR_GRID<GRID, SCALAR_TYPE> SCALAR_GRID;
  typedef IJK::VECTOR_GRID<GRID, NUM_TYPE, GRADIENT_COORD_TYPE> GRADIENT_GRID;

  typedef OBJECT_PROPERTIES_T
  <int, COORD_TYPE, GRADIENT_COORD_TYPE, RADIUS_TYPE, DIFF_TYPE, 
   ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE>
  OBJECT_PROPERTIES;
  typedef FIELD_PARAM_T
  <int, AXIS_SIZE_TYPE, COORD_TYPE, GRADIENT_COORD_TYPE, RADIUS_TYPE, 
   DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, SEED_TYPE, MIN_MAX_TYPE, 
   ISOTABLE_INDEX_TYPE, NUM_TYPE> 
  FIELD_PARAM;

  typedef FIELD_INFO_T<OBJECT_PROPERTIES, SCALAR_GRID, GRADIENT_GRID> FIELD_INFO;
};

#endif

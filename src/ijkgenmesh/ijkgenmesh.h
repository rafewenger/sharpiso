/// \file ijkgenmesh.h
/// generate a mesh
/// Version v0.1.4

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

/*!
  \mainpage IJKGENMESH: Generate a mesh.

  IJKGENMESH is a program for generating a mesh.
*/

#ifndef _IJKGENMESH_
#define _IJKGENMESH_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkmesh.txx"

#include "ijkgenmesh.txx"


namespace IJKGENMESH {

  // ********************************************************
  // TYPES
  // ********************************************************

  typedef size_t AXIS_SIZE_TYPE;
  typedef float COORD_TYPE;
  typedef float RADIUS_TYPE;
  typedef float DIFF_TYPE;
  typedef float ANGLE_TYPE;
  typedef int NUM_TYPE;
  typedef unsigned int SEED_TYPE;

  typedef IJK::POLYMESH<int, int, int> POLYMESH_TYPE;

  typedef MESH_OBJECT_PROPERTIES_T
  <int, COORD_TYPE, COORD_TYPE, RADIUS_TYPE, DIFF_TYPE, 
   ANGLE_TYPE, COORD_TYPE, NUM_TYPE, COORD_TYPE>
  OBJECT_PROPERTIES;
  typedef IJKGENGEOM::GEOM_PARAM_T<OBJECT_PROPERTIES, AXIS_SIZE_TYPE, SEED_TYPE>
    GEOM_PARAM;
  typedef MESH_PARAM_T<GEOM_PARAM> MESH_PARAM;

  typedef MESH_INFO_T<int, COORD_TYPE, POLYMESH_TYPE, OBJECT_PROPERTIES>
  MESH_INFO;

  // ********************************************************
  // CONSTANTS
  // ********************************************************

  const int DIM3(3);
  const int NUM_VERT_PER_TRIANGLE(3);
  const int NUM_VERT_PER_QUAD(4);


  // ********************************************************
  // ROUTINES
  // ********************************************************

  /// Generate a cube.
  void gen_cube(const int dimension, const OBJECT_PROPERTIES & prop,
                std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);

  /// Generate an annulus.
  /// @param prop Annulus properties.  
  ///   If prop.flag_flange is true, then annulus has a flange.
  void gen_annulus
    (const int dimension, const OBJECT_PROPERTIES & prop,
     std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);

  /// Generate an annulus with no flange.
  void gen_annulus_no_flange
  (const int dimension, const OBJECT_PROPERTIES & prop,
   std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);

  /// Generate a flanged annulus.
  void gen_flanged_annulus
  (const int dimension, const OBJECT_PROPERTIES & prop,
   std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);

};

#endif

/// \file ijkgenmesh.txx
/// generate a mesh

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

#ifndef _IJKGENMESH_TXX_
#define _IJKGENMESH_TXX_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkmesh.txx"

#include "ijkgengeom.txx"


namespace IJKGENMESH {

  // ********************************************************
  // TEMPLATE CLASS MESH_INFO_T
  // ********************************************************

  /// Information about a mesh
  /// Field name, scalar field generator, gradient field generator, etc.
  template <typename DIM_TYPE, typename COORD_TYPE, typename MESH_TYPE,
            typename OBJECT_PROPERTIES_TYPE>
  class MESH_INFO_T:public IJKGENGEOM::GEOM_INFO_T<OBJECT_PROPERTIES_TYPE> {

  public:
    typedef void (*FUNCTION_PTR)
      (const DIM_TYPE dimension, const OBJECT_PROPERTIES_TYPE & prop, 
       std::vector<COORD_TYPE> & coord, MESH_TYPE & mesh);

  public:
    FUNCTION_PTR function_ptr;

    /// If true, object has distance to "center/centerline".
    /// Often equivalent to isovalue in corresponding scalar field.
    bool flag_distance;  

  public:
    MESH_INFO_T() {};

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr)
    {
      this->name = field_name;
      function_ptr = f_ptr;
      this->SetAllFlags(false);
    }
  };

  // ********************************************************
  // TEMPLATE CLASS MESHOBJECT_PROPERTIES_T
  // ********************************************************

  /// Object properties.
  template <typename DIM_TYPE, typename COORD_TYPE, typename DIR_TYPE,
            typename RADIUS_TYPE, typename DIFF_TYPE, typename ANGLE_TYPE,
            typename SCALAR_TYPE, typename NUM_TYPE, typename DIST_TYPE>
  class MESH_OBJECT_PROPERTIES_T:
    public IJKGENGEOM::OBJECT_PROPERTIES_T
    <DIM_TYPE,COORD_TYPE,DIR_TYPE,RADIUS_TYPE,
     DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> {
    
  protected:
    void Init();

  public:

    /// Distance to "center/centerline".
    /// Often equivalent to isovalue in corresponding scalar field.
    IJKGENGEOM::SET_VALUE<DIST_TYPE> distance;

  public:
    MESH_OBJECT_PROPERTIES_T(){ Init(); };
    

    DIST_TYPE Distance() const
    { return(distance.Value()); }

    /// Copy.
    void Copy(const MESH_OBJECT_PROPERTIES_T
              <DIM_TYPE, COORD_TYPE, DIR_TYPE, RADIUS_TYPE,
              DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE, DIST_TYPE> 
              & right)
    {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
    };

  };

  // ********************************************************
  // TEMPLATE CLASS MESH_PARAM_T
  // ********************************************************

  /// Input parameters determining field properties.
  template <typename GEOM_PARAM_TYPE>
  class MESH_PARAM_T:public GEOM_PARAM_TYPE {

  protected:

    void Init();

  public:
    typedef typename GEOM_PARAM_TYPE::NUMBER_TYPE NUMBER_TYPE;

  public:
    bool flag_triangulate;

  public:
    MESH_PARAM_T(){ Init(); };

    void SetMeshIndex(const NUMBER_TYPE field_index)
    { this->SetGeomInfoIndex(field_index); }

    NUMBER_TYPE MeshIndex() const
    { return(this->GeomInfoIndex()); }

  };

  // ********************************************************
  // MACROS FOR TEMPLATE MEMBER FUNCTIONS
  // ********************************************************

#define _MESH_INFO_T_HEADER_                                         \
  template <typename OBJECT_PROPERTIES_TYPE>                         \
  void MESH_INFO_T<OBJECT_PROPERTIES_TYPE>::


#define _MESH_OBJECT_PROPERTIES_T_HEADER_                            \
  template <typename DIM_TYPE,                                       \
            typename COORD_TYPE, typename DIR_TYPE,                  \
            typename RADIUS_TYPE, typename DIFF_TYPE,                \
            typename ANGLE_TYPE, typename SCALAR_TYPE,               \
            typename NUM_TYPE, typename DIST_TYPE>                   \
  void MESH_OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE, DIR_TYPE,      \
                           RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,       \
                                SCALAR_TYPE, NUM_TYPE, DIST_TYPE>::

#define _MESH_PARAM_T_HEADER_                                        \
  template <typename GEOM_PARAM_TYPE>                                \
  void MESH_PARAM_T<GEOM_PARAM_TYPE>::


  // ********************************************************
  // TEMPLATE CLASS MESH_OBJECT_PROPERTIES_T MEMBER FUNCTIONS
  // ********************************************************

  _MESH_OBJECT_PROPERTIES_T_HEADER_
  Init()
  {
    this->distance = 0;
  }

  // ********************************************************
  // TEMPLATE CLASS MESH_PARAM_T MEMBER FUNCTIONS
  // ********************************************************

  _MESH_PARAM_T_HEADER_
  Init()
  {
    flag_triangulate = true;
  }

};

#endif

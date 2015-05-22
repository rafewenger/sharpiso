/// \file ijkgenscalar.txx
/// generate a scalar field

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

#ifndef _IJKGENSCALAR_TXX_
#define _IJKGENSCALAR_TXX_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.txx"
#include "ijkgengeom.txx"

namespace IJKGENSCALAR {

  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO_T
  // ********************************************************

  /// Information about a field.
  /// Field name, scalar field generator, gradient field generator, etc.
  template <typename OBJECT_PROPERTIES_TYPE,
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  class FIELD_INFO_T:public IJKGENGEOM::GEOM_INFO_T<OBJECT_PROPERTIES_TYPE> {

  protected:
    void Init();

  public:
    typedef void (*FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop, SCALAR_GRID_TYPE & grid);
    typedef void (*GRADIENT_FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop,
       SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid);

  public:
    FUNCTION_PTR function_ptr;
    GRADIENT_FUNCTION_PTR gradient_function_ptr;
    bool flag_maxval;           ///< If true, field has a maximum scalar value.

  public:
    FIELD_INFO_T() { Init(); };

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr)
    {
      this->name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = NULL;
      SetAllFlags(false);
    }

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr,
             GRADIENT_FUNCTION_PTR grad_f_ptr)
    {
      this->name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = grad_f_ptr;
      SetAllFlags(false);
    }

    void SetAllFlags(const bool flag);
    void SetRandomFlags();

    /// Return true if gradient_function_ptr is not NULL.
    bool IsGradientImplemented() const
    {
      if (gradient_function_ptr == NULL) 
        { return(false); }
      else
        { return(true); }
    }

  };

  // ********************************************************
  // TEMPLATE CLASS FIELD_OBJECT_PROPERTIES_T
  // ********************************************************

  /// Object properties.
  template <typename DIM_TYPE, typename COORD_TYPE, typename DIR_TYPE,
            typename RADIUS_TYPE, typename DIFF_TYPE, typename ANGLE_TYPE,
            typename SCALAR_TYPE, typename NUM_TYPE>
  class FIELD_OBJECT_PROPERTIES_T:
    public IJKGENGEOM::OBJECT_PROPERTIES_T
    <DIM_TYPE,COORD_TYPE,DIR_TYPE,RADIUS_TYPE,
     DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> {
    
  protected:
    void Init();

  public:
    FIELD_OBJECT_PROPERTIES_T(){ Init(); };
    bool gradient_discontinuity_zero;

    /// Copy.
    void Copy(const FIELD_OBJECT_PROPERTIES_T
              <DIM_TYPE, COORD_TYPE, DIR_TYPE, RADIUS_TYPE,
              DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> 
              & right)
    {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
    };

  };

  // ********************************************************
  // TEMPLATE CLASS FIELD_PARAM_T
  // ********************************************************

  /// Input parameters determining field properties.
  template <typename GEOM_PARAM_TYPE, 
            typename MIN_MAX_TYPE, typename ISOTABLE_INDEX_TYPE>
  class FIELD_PARAM_T:public GEOM_PARAM_TYPE {

  public:
    typedef typename GEOM_PARAM_TYPE::NUMBER_TYPE NUMBER_TYPE;

  protected:
    void Init();

  public:
    FIELD_PARAM_T(){ Init(); };

    IJKGENGEOM::SET_VALUE<MIN_MAX_TYPE> maxval;      ///< Max scalar value.

    /// Isosurface lookup table index.
    IJKGENGEOM::SET_VALUE<ISOTABLE_INDEX_TYPE> isotable_index;
    
    /// If true, set all boundary vertices to zero.
    bool flag_set_boundary2zero;

    MIN_MAX_TYPE MaxVal() const
    { return(maxval.Value()); }
    ISOTABLE_INDEX_TYPE IsotableIndex() const
    { return(isotable_index.Value()); }

    void SetFieldIndex(const NUMBER_TYPE field_index)
    { this->SetGeomInfoIndex(field_index); }

    NUMBER_TYPE FieldIndex() const
    { return(this->GeomInfoIndex()); }

  };

  // ********************************************************
  // MACROS FOR TEMPLATE MEMBER FUNCTIONS
  // ********************************************************

#define _FIELD_INFO_T_HEADER_                                        \
  template <typename OBJECT_PROPERTIES_TYPE,                         \
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>  \
  void FIELD_INFO_T<OBJECT_PROPERTIES_TYPE,                          \
                    SCALAR_GRID_TYPE,GRADIENT_GRID_TYPE>::


#define _FIELD_OBJECT_PROPERTIES_T_HEADER_                           \
  template <typename DIM_TYPE,                                       \
            typename COORD_TYPE, typename DIR_TYPE,                  \
            typename RADIUS_TYPE, typename DIFF_TYPE,                \
            typename ANGLE_TYPE, typename SCALAR_TYPE,               \
            typename NUM_TYPE>                                       \
  void FIELD_OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE, DIR_TYPE,     \
                           RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,       \
                           SCALAR_TYPE, NUM_TYPE>::

#define _FIELD_PARAM_T_HEADER_                                       \
  template <typename GEOM_PARAM_TYPE,                                \
            typename MIN_MAX_TYPE, typename ISOTABLE_INDEX_TYPE>     \
  void FIELD_PARAM_T<GEOM_PARAM_TYPE,                                \
                     MIN_MAX_TYPE, ISOTABLE_INDEX_TYPE>::



  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_INFO_T_HEADER_
    Init()
    {
      SetAllFlags(false);
    }


  _FIELD_INFO_T_HEADER_
  SetRandomFlags()
  {
    this->flag_random_seed = true;
    flag_maxval = true;
    this->flag_allow_multi_centers = false;
    this->flag_allow_tilt = false;
  }

  _FIELD_INFO_T_HEADER_
  SetAllFlags(const bool flag)
  {
    this->IJKGENGEOM::GEOM_INFO_T<OBJECT_PROPERTIES_TYPE>::SetAllFlags(flag);
    this->flag_maxval = false;
  }

  // ********************************************************
  // TEMPLATE CLASS OBJECT_PROPERTIES_T MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_OBJECT_PROPERTIES_T_HEADER_
  Init()
  {
    this->gradient_discontinuity_zero = true;
  }

  // ********************************************************
  // TEMPLATE CLASS FIELD_PARAM_T MEMBER FUNCTIONS
  // ********************************************************

  template <typename T> 
  void convert_to_string(const T x, std::string & s)
  {
    std::ostringstream ostr("");

    ostr << x;
    s = ostr.str();
  }
    
	
  template <typename ITYPE> 
  void convert_to_string(const ITYPE list_begin, const ITYPE list_end,
                         std::string & s)
  {
      std::ostringstream ostr("");

      for (ITYPE iter = list_begin; iter != list_end; iter++) {
        ostr << *iter;
        if (iter+1 != list_end) { ostr << " "; }
      }
      s = ostr.str();
  }

  _FIELD_PARAM_T_HEADER_
  Init()
  {
    this->gradient_discontinuity_zero = true;
  }
};

#endif

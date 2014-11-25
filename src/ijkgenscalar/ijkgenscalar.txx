/// \file ijkgenscalar2.h
/// generate a scalar field

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

#ifndef _IJKGENSCALAR_TXX_
#define _IJKGENSCALAR_TXX_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.txx"

namespace IJKGENSCALAR {

  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO_T
  // ********************************************************

  /// Information about a field.
  /// Field name, scalar field generator, gradient field generator, etc.
  template <typename OBJECT_PROPERTIES_TYPE,
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  class FIELD_INFO_T {

  protected:
    void Init();

  public:
    typedef void (*FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop, SCALAR_GRID_TYPE & grid);
    typedef void (*GRADIENT_FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop,
       SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid);

  public:
    std::string name;
    std::vector<std::string> description;
    FUNCTION_PTR function_ptr;
    GRADIENT_FUNCTION_PTR gradient_function_ptr;
    bool flag_center;           ///< If true, field has a center.
    bool flag_direction;        ///< If true, field has a direction.
    bool flag_radius;           ///< If true, field has a radius.
    bool flag_angle;            ///< If true, field has an angle.
    bool flag_length_difference;  ///< If true, field has length difference.
    bool flag_cylinder;         ///< If true, field has cylinder properties.
    bool flag_clip;             ///< If true, field has clipping planes.
    bool flag_translate;        ///< If true, field has a translation vector.
    bool flag_random_seed;      ///< If true, field has a random seed.
    bool flag_maxval;           ///< If true, field has a maximum scalar value.

    /// If true, allow combining multiple copies of the field.
    bool flag_allow_multi_centers;

    /// If true, isosurface may have tilt parameters.
    bool flag_allow_tilt;

    /// If true, isosurfaces may have a flange.
    bool flag_allow_flange;

    /// If true, isosurfaces may be a wedge
    bool flag_allow_wedge;

  public:
    FIELD_INFO_T() { Init(); };

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr)
    {
      name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = NULL;
      SetAllFlags(false);
    }

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr,
             GRADIENT_FUNCTION_PTR grad_f_ptr)
    {
      name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = grad_f_ptr;
      SetAllFlags(false);
    }

    void AddToDescription(const char * s)
    { description.push_back(s); }

    void SetDescription(const char * s)
    {
      description.clear();
      AddToDescription(s);
    }

    void SetFlagCenter(const bool flag_center);
    void SetCubeFlags();
    void SetSphereFlags();
    void SetCylinderFlags();
    void SetCrossingCylinderFlags();
    void SetAnnulusFlags();
    void SetTorusFlags();
    void SetConeFlags();
    void SetFrustrumFlags();
    void SetOctahedronFlags();
    void SetRandomFlags();
    void SetAllFlags(const bool flag);

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
  // TEMPLATE CLASS OBJECT_PROPERTIES_T
  // ********************************************************

  // forward declaration of template class SET_VALUE
  template <typename T> class SET_VALUE;

  /// Object properties.
  template <typename DIM_TYPE, typename COORD_TYPE, typename DIR_TYPE,
            typename RADIUS_TYPE, typename DIFF_TYPE, typename ANGLE_TYPE,
            typename SCALAR_TYPE, typename NUM_TYPE>
  class OBJECT_PROPERTIES_T {

  protected:
    void Init();

  public:
    OBJECT_PROPERTIES_T(){ Init(); };

    SET_VALUE<DIM_TYPE> dimension;
    bool flag_tilt;
    bool flag_flange;
    bool flag_wedge;
    bool gradient_discontinuity_zero;

    std::vector<COORD_TYPE> center;
    std::vector<DIR_TYPE> direction;
    std::vector<DIR_TYPE> side_direction;
    std::vector<DIR_TYPE> rotation_direction;
    std::vector<DIR_TYPE> normal;
    std::vector<RADIUS_TYPE> radius;
    std::vector<DIFF_TYPE> length_difference;
    std::vector<COORD_TYPE> flange_width;
    std::vector<COORD_TYPE> flange_height;
    std::vector<DIFF_TYPE> dist2near0;
    std::vector<DIFF_TYPE> dist2far0;
    std::vector<ANGLE_TYPE> angle;
    std::vector<ANGLE_TYPE> wedge_angle;
    std::vector<SCALAR_TYPE> wedge_isovalue;

    // Get functions
    DIM_TYPE Dimension() const
    { return(dimension.Value()); }
    COORD_TYPE * CenterPtr(const NUM_TYPE i)
      { return(&(center[0]) + i*Dimension()); }
    const COORD_TYPE * CenterPtrConst(const NUM_TYPE i) const
      { return(&(center[0]) + i*Dimension()); }
    DIR_TYPE * DirectionPtr(const NUM_TYPE i)
      { return(&(direction[0]) + i*Dimension()); }
    const DIR_TYPE * DirectionPtrConst(const NUM_TYPE i) const
      { return(&(direction[0]) + i*Dimension()); }
    const DIR_TYPE * SideDirectionPtrConst(const NUM_TYPE i) const
      { return(&(side_direction[0]) + i*Dimension()); }
    const DIR_TYPE * RotationDirectionPtrConst() const
      { return(&(rotation_direction[0])); }
    const DIR_TYPE * NormalPtrConst() const
      { return(&(normal[0])); }

    NUM_TYPE NumCenters() const
    { return(center.size()/Dimension()); }
    NUM_TYPE NumDirections() const
    { return(direction.size()/Dimension()); }
    NUM_TYPE NumSideDirections() const
    { return(side_direction.size()/Dimension()); }
    NUM_TYPE NumRadii() const
    { return(radius.size()); }
    NUM_TYPE NumLengthDifferences() const
    { return(length_difference.size()); }
    NUM_TYPE NumNormals() const
    { return(normal.size()/Dimension()); }
    NUM_TYPE NumRotationDirections() const
    { return(rotation_direction.size()/Dimension()); }
    NUM_TYPE NumAngles() const
    { return(angle.size()); }
    NUM_TYPE NumFlangeWidth() const
    { return(flange_width.size()); }
    NUM_TYPE NumFlangeHeight() const
    { return(flange_height.size()); }
    NUM_TYPE NumWedgeAngles() const
    { return(wedge_angle.size()); }
    NUM_TYPE NumWedgeIsovalues() const
    { return(wedge_isovalue.size()); }

    /// Copy.
    void Copy(const OBJECT_PROPERTIES_T
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
  template <typename DIM_TYPE, typename AXIS_SIZE_TYPE,
            typename COORD_TYPE, typename DIR_TYPE,
            typename RADIUS_TYPE, typename DIFF_TYPE, typename ANGLE_TYPE,
            typename SCALAR_TYPE,
            typename SEED_TYPE,
            typename MIN_MAX_TYPE, 
            typename ISOTABLE_INDEX_TYPE,
            typename NUM_TYPE>
  class FIELD_PARAM_T:
    public OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE, DIR_TYPE, RADIUS_TYPE,
                               DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> {

  protected:
    void Init();

  public:
    FIELD_PARAM_T(){ Init(); };

    SET_VALUE<NUM_TYPE> field_index;     ///< Index into field info[].
    SET_VALUE<NUM_TYPE> num_objects;     ///< Number of cubes, tori, etc.
    SET_VALUE<SEED_TYPE> random_seed;    ///< Random seed.
    SET_VALUE<SEED_TYPE> randompos_seed; ///< Random positions seed.
    SET_VALUE<SEED_TYPE> randomdir_seed; ///< Random directions seed.
    SET_VALUE<MIN_MAX_TYPE> maxval;      ///< Max scalar value.

    /// Isosurface lookup table index.
    SET_VALUE<ISOTABLE_INDEX_TYPE> isotable_index;

    std::vector<AXIS_SIZE_TYPE> axis_size;
    std::vector<COORD_TYPE> spacing;
    std::vector<DIR_TYPE> translate;

    bool flag_multi_centers;
    bool flag_multi_normals;
    bool flag_stack;
    
    /// If true, set all boundary vertices to zero.
    bool flag_set_boundary2zero;

    /// If true, use grid center as center coordinates.
    bool flag_grid_center;

    /// If true, generate random center coordinates.
    bool flag_random_centers;

    /// If true, generate random directions.
    bool flag_random_directions;

    // Set functions
    void SetDimension(const DIM_TYPE dimension)
    { this->dimension.Set(dimension); }
    void SetFieldIndex(const NUM_TYPE field_index)
    { this->field_index.Set(field_index); }
    void SetAxisSize(const AXIS_SIZE_TYPE size);
    
    // Get functions
    DIM_TYPE Dimension() const
    { return(OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE,
             DIR_TYPE, RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,
             SCALAR_TYPE, NUM_TYPE>::Dimension()); }
    NUM_TYPE FieldIndex() const
    { return(field_index.Value()); }
    NUM_TYPE NumObjects() const
    { return(num_objects.Value()); }
    NUM_TYPE NumTranslate() const
    { return(this->translate.size()/Dimension()); }
    const SEED_TYPE RandomSeed() const
    { return(random_seed.Value()); }
    MIN_MAX_TYPE MaxVal() const
    { return(maxval.Value()); }
    ISOTABLE_INDEX_TYPE IsotableIndex() const
    { return(isotable_index.Value()); }
    const DIR_TYPE * TranslatePtrConst(const NUM_TYPE i) const
      { return(&(translate[0]) + i*Dimension()); }
    const COORD_TYPE * SpacingPtrConst() const
      { return(&(spacing[0])); }

    void GetCenterString(const NUM_TYPE i, std::string & s) const;
    void GetDirectionString(const NUM_TYPE i, std::string & s) const;
    void GetSideDirectionString(const NUM_TYPE i, std::string & s) const;
    void GetRotationDirectionString(const NUM_TYPE i, std::string & s) const;
    void GetNormalString(const NUM_TYPE i, std::string & s) const;
    void GetTranslateString(const NUM_TYPE i, std::string & s) const;
    void GetRadiusString(const NUM_TYPE i, std::string & s) const;
    void GetAngleString(const NUM_TYPE i, std::string & s) const;
    void GetLengthDifferenceString(const NUM_TYPE i, std::string & s) const;
    void GetDist2Near0String(const NUM_TYPE i, std::string & s) const;
    void GetDist2Far0String(const NUM_TYPE i, std::string & s) const;
    void GetFlangeWidthString(const NUM_TYPE i, std::string & s) const;
    void GetFlangeHeightString(const NUM_TYPE i, std::string & s) const;
    void GetWedgeAngleString(const NUM_TYPE i, std::string & s) const;
    void GetWedgeIsovalueString(const NUM_TYPE i, std::string & s) const;
  };

  // ********************************************************
  // TEMPLATE CLASS SET_VALUE
  // ********************************************************

  template <typename T>
  class SET_VALUE {

  protected:
    T v;
    bool is_set;

    void Init() 
    { is_set = false; }

    void Init(const T value)
    {
      v = value;
      is_set = false;
    }

  public:
    SET_VALUE(){ Init(); };
    SET_VALUE(const T value) { Init(value); };

    void Set(const T value)
    {
      v = value;
      is_set = true;
    }

    bool IsSet() const { return(is_set); };
    T Value() const { return(v); };
  };

  // ********************************************************
  // MACROS FOR TEMPLATE MEMBER FUNCTIONS
  // ********************************************************

#define _FIELD_INFO_T_HEADER_                                        \
  template <typename OBJECT_PROPERTIES_TYPE,                         \
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>  \
  void FIELD_INFO_T<OBJECT_PROPERTIES_TYPE,                          \
                    SCALAR_GRID_TYPE,GRADIENT_GRID_TYPE>::


#define _OBJECT_PROPERTIES_T_HEADER_                                 \
  template <typename DIM_TYPE,                                       \
            typename COORD_TYPE, typename DIR_TYPE,                  \
            typename RADIUS_TYPE, typename DIFF_TYPE,                \
            typename ANGLE_TYPE, typename SCALAR_TYPE,               \
            typename NUM_TYPE>                                       \
  void OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE, DIR_TYPE,           \
                           RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,       \
                           SCALAR_TYPE, NUM_TYPE>::

#define _FIELD_PARAM_T_HEADER_                                       \
  template <typename DIM_TYPE, typename AXIS_SIZE_TYPE,              \
            typename COORD_TYPE, typename DIR_TYPE,                  \
            typename RADIUS_TYPE, typename DIFF_TYPE,                \
            typename ANGLE_TYPE, typename SCALAR_TYPE,               \
            typename SEED_TYPE,                                      \
            typename MIN_MAX_TYPE,                                   \
            typename ISOTABLE_INDEX_TYPE,                            \
            typename NUM_TYPE>                                       \
  void FIELD_PARAM_T<DIM_TYPE, AXIS_SIZE_TYPE, COORD_TYPE, DIR_TYPE, \
                     RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,             \
                     SCALAR_TYPE, SEED_TYPE,                         \
                     MIN_MAX_TYPE, ISOTABLE_INDEX_TYPE,              \
                     NUM_TYPE>::

#define _GET_COORDINATE_STRING_(_function_name,_obj)                 \
  _FIELD_PARAM_T_HEADER_ _function_name                              \
  (const NUM_TYPE i, std::string & s) const                          \
  {                                                                  \
    convert_to_string(this->_obj.begin()+i*Dimension(),              \
                      this->_obj.begin()+(i+1)*Dimension(), s);      \
  }

#define _GET_SCALAR_STRING_(_function_name,_obj)                     \
  _FIELD_PARAM_T_HEADER_ _function_name                              \
  (const NUM_TYPE i, std::string & s) const                          \
  {                                                                  \
    convert_to_string(this->_obj.begin()+i,                          \
                      this->_obj.begin()+i+1, s);                    \
  }

  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_INFO_T_HEADER_
    Init()
    {
      SetAllFlags(false);
    }


  _FIELD_INFO_T_HEADER_
    SetFlagCenter(const bool flag_center)
    { 
      this->flag_center = flag_center;
    }

  _FIELD_INFO_T_HEADER_
    SetCubeFlags()
  {
    flag_center = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = true;
    flag_allow_flange = false;
    flag_allow_wedge = false;
  }

  _FIELD_INFO_T_HEADER_
    SetSphereFlags()
  {
    flag_center = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = false;
    flag_allow_wedge = false;
  }

  _FIELD_INFO_T_HEADER_
  SetCylinderFlags()
  { 
    flag_center = true;
    flag_direction = true;
    flag_length_difference = true;
    flag_cylinder = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = true;
    flag_allow_flange = true;
    flag_allow_wedge = true;
  }

  _FIELD_INFO_T_HEADER_
  SetCrossingCylinderFlags()
  {
    SetCylinderFlags();
    flag_allow_multi_centers = false;
    flag_allow_flange = false;
    flag_allow_wedge = false;
  }

  _FIELD_INFO_T_HEADER_
  SetAnnulusFlags()
  { 
    SetCylinderFlags();
    flag_radius = true;
  }

  _FIELD_INFO_T_HEADER_
  SetTorusFlags()
  { 
    flag_center = true;
    flag_direction = true;
    flag_radius = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = true;
    flag_allow_wedge = true;
  }

  _FIELD_INFO_T_HEADER_
  SetConeFlags()
  { 
    flag_center = true;
    flag_direction = true;
    flag_angle = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = true;
    flag_allow_flange = false;
    flag_allow_wedge = true;
  }

  _FIELD_INFO_T_HEADER_
  SetFrustrumFlags()
  { 
    SetConeFlags();
    flag_allow_flange = true;
    flag_clip = true;
  }

  _FIELD_INFO_T_HEADER_
  SetOctahedronFlags()
  {
    flag_center = true;
    flag_allow_multi_centers = true;
    flag_allow_tilt = false;
    flag_allow_wedge = false;
  }

  _FIELD_INFO_T_HEADER_
  SetRandomFlags()
  {
    flag_random_seed = true;
    flag_maxval = true;
    flag_allow_multi_centers = false;
    flag_allow_tilt = false;
  }

  _FIELD_INFO_T_HEADER_
  SetAllFlags(const bool flag)
  {
    flag_center = flag;
    flag_direction = flag;
    flag_radius = flag;
    flag_angle = flag;
    flag_length_difference = flag;
    flag_cylinder = flag;
    flag_clip = flag;
    flag_translate = flag;
    flag_random_seed = false;
    flag_maxval = false;
    flag_allow_multi_centers = flag;
    flag_allow_tilt = flag;
    flag_allow_flange = flag;
    flag_allow_wedge = flag;
  }

  // ********************************************************
  // TEMPLATE CLASS OBJECT_PROPERTIES_T MEMBER FUNCTIONS
  // ********************************************************

  _OBJECT_PROPERTIES_T_HEADER_
  Init()
  {
    flag_tilt = true;
    flag_flange = false;
    flag_wedge = false;
    gradient_discontinuity_zero = true;
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
    this->flag_tilt = true;
    this->flag_flange = false;
    this->flag_wedge = false;
    flag_multi_centers = false;
    flag_multi_normals = false;
    flag_set_boundary2zero = false;
    flag_grid_center = false;
    flag_random_centers = false;
    flag_random_directions = false;
    this->gradient_discontinuity_zero = true;
  }

  _FIELD_PARAM_T_HEADER_
  SetAxisSize(const AXIS_SIZE_TYPE size)
  {
    IJK::PROCEDURE_ERROR error("FIELD_PARAM_T.SetAxisSIze");

    if (!this->dimension.IsSet()) {
      error.AddMessage
        ("Programming error. Dimension must be set before axis size.");
      throw error;
    }

    for (DIM_TYPE d = 0; d < Dimension(); d++) 
      { axis_size.push_back(size); }
  }

  _GET_COORDINATE_STRING_(GetCenterString, center);
  _GET_COORDINATE_STRING_(GetDirectionString, direction);
  _GET_COORDINATE_STRING_(GetSideDirectionString, side_direction);
  _GET_COORDINATE_STRING_(GetRotationDirectionString, rotation_direction);
  _GET_COORDINATE_STRING_(GetNormalString, normal);
  _GET_COORDINATE_STRING_(GetTranslateString, translate);

  _GET_SCALAR_STRING_(GetLengthDifferenceString, length_difference);
  _GET_SCALAR_STRING_(GetRadiusString, radius);
  _GET_SCALAR_STRING_(GetAngleString, angle);
  _GET_SCALAR_STRING_(GetWedgeAngleString, wedge_angle);
  _GET_SCALAR_STRING_(GetWedgeIsovalueString, wedge_isovalue);
  _GET_SCALAR_STRING_(GetFlangeWidthString, flange_width);
  _GET_SCALAR_STRING_(GetFlangeHeightString, flange_height);
  _GET_SCALAR_STRING_(GetDist2Near0String, dist2near0);
  _GET_SCALAR_STRING_(GetDist2Far0String, dist2far0);
};

#endif

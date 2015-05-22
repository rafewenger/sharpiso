/// \file ijkgenscalarIO.cxx
/// Read/write/prompt for ijkgenscalar
/// Version v0.1.3

/*
  IJK: Isosurface Jeneration Code
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

#include <cmath>
#include <iomanip>
#include <iostream>

#include "ijkgenscalar.h"
#include "ijkgenscalarIO.h"

#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"
#include "ijkstring.txx"

// Types
typedef IJK::NRRD_DATA<int, int> NRRD_HEADER;

using namespace IJKGENSCALAR;
using namespace std;

// **************************************************
// WRITE NRRD FILE
// **************************************************

void make_string
(const string & s_in, const int i, string & s_out)
{
  string s;

  IJK::val2string(i,s);
  s_out = s_in + s;
}

void add_nrrd_key_values
(const FIELD_PARAM & field_param, NRRD_HEADER & nrrd_header,
 const std::vector<FIELD_INFO> & field_info)
{
  const int field_index = field_param.FieldIndex();
  string field_name = field_info[field_index].name;
  string s;

  if (field_info[field_index].flag_center) {

    if (field_param.NumCenters() > 0) {
      if (field_param.NumCenters() == 1 || !field_param.flag_multi_centers) {
        field_param.GetCenterString(0, s);
        nrrd_header.AddKeyValue("center", s);
      }
      else {
        for (int i = 0; (i < field_param.NumCenters()) && 
               (i < field_param.NumObjects()); i++) {
          string key;
          make_string("center", i, key);
          field_param.GetCenterString(i, s);
          nrrd_header.AddKeyValue(key, s);
        }
      }
    }
  }

  if (field_name == "cube" || 
      field_name == "squareX" || field_name == "square45X") {
    if (field_param.flag_tilt) {

      if (field_param.NumDirections() > 0) {
        if (field_param.NumDirections() == 1 || 
            !field_param.flag_multi_centers) {
          field_param.GetDirectionString(0, s);
          nrrd_header.AddKeyValue("firstAxisDirection", s);
        }
        else {
          for (int i = 0; (i < field_param.NumDirections()) && 
                 (i < field_param.NumDirections()); i++) {
            string key;
            make_string("firstAxisDirection", i, key);
            field_param.GetDirectionString(i, s);
            nrrd_header.AddKeyValue(key, s);
          }
        }
      }

      if (field_param.NumSideDirections() == 1) {
        field_param.GetSideDirectionString(0, s);
        nrrd_header.AddKeyValue("secondAxisDirection", s);
      }
    }
  }
  else if (field_info[field_index].flag_direction) {

    if (field_param.NumDirections() > 0) {
      if (field_param.NumDirections() == 1 || !field_param.flag_multi_centers) {
        field_param.GetDirectionString(0, s);
        nrrd_header.AddKeyValue("axisDirection", s);
      }
      else {
        for (int i = 0; (i < field_param.NumDirections()) && 
               (i < field_param.NumDirections()); i++) {
          string key;
          make_string("axisDirection", i, key);
          field_param.GetDirectionString(i, s);
          nrrd_header.AddKeyValue(key, s);
        }
      }
    }
  }

  if (field_name == "constant_unit_gradient") {

    if (field_param.NumDirections() == 1) {
      string coord_string;
      field_param.GetDirectionString(0, coord_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "gradientDirection", coord_string.c_str());
    }
  }

  if (field_info[field_index].flag_radius) {

    if (field_param.NumRadii() == 1) {
      field_param.GetRadiusString(0, s);
      nrrd_header.AddKeyValue("radius", s);
    }
  }

  if (field_info[field_index].flag_length_difference) {
    if (field_param.NumLengthDifferences() == 1) {
      field_param.GetLengthDifferenceString(0, s);
      nrrd_header.AddKeyValue("lengthDifference", s);
    }
  }

  if (field_name == "square") {

    if (field_param.NumSideDirections() == 1) {
      string coord_string;
      field_param.GetSideDirectionString(0, coord_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "sideDirection", coord_string.c_str());
    }

  }

  if (field_name == "square45") {

    if (field_param.NumSideDirections() == 1) {
      string coord_string;
      field_param.GetSideDirectionString(0, coord_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "squareCornerDirection", coord_string.c_str());
    }
  }

  if (field_name == "edge") {

    if (field_param.NumAngles() == 1) {
      string angle_string;
      field_param.GetAngleString(0, angle_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "angle", angle_string.c_str());
    }

    if (field_param.NumNormals() == 1) {
      string normal_string;
      field_param.GetNormalString(0, normal_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "normal0", normal_string.c_str());
    }

    if (field_param.NumRotationDirections() == 1) {
      string direction_string;
      field_param.GetNormalString(0, direction_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "rotationDirection", direction_string.c_str());
    }

  }

  if (field_info[field_index].flag_angle) {
    if (field_param.NumAngles() == 1) {
      string angle_string;
      field_param.GetAngleString(0, angle_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "angle", angle_string.c_str());
    }
  }

  if (field_name == "frustrum") {
    if (field_param.dist2near0.size() == 1) {
      field_param.GetDist2Near0String(0, s);
      nrrd_header.AddKeyValue("dist2near0", s);
    }

    if (field_param.dist2far0.size() == 1) {
      field_param.GetDist2Far0String(0, s);
      nrrd_header.AddKeyValue("dist2far0", s);
    }
  }

  if (field_param.flag_flange) {
    if (field_param.NumFlangeWidth() == 1) {
      field_param.GetFlangeWidthString(0, s);
      nrrd_header.AddKeyValue("flangeWidth", s);
    }

    if (field_param.NumFlangeHeight() == 1) {
      field_param.GetFlangeHeightString(0, s);
      nrrd_header.AddKeyValue("flangeHeight", s);
    }
  }

  if (field_param.flag_wedge) {

    if (field_param.NumSideDirections() == 1) {
      field_param.GetSideDirectionString(0, s);
      nrrd_header.AddKeyValue("wedgeNormal", s);
    }

    if (field_param.NumWedgeAngles() == 1) {
      field_param.GetWedgeAngleString(0, s);
      nrrd_header.AddKeyValue("wedgeAngle", s);
    }

    if (field_param.NumWedgeIsovalues() == 1) {
      field_param.GetWedgeIsovalueString(0, s);
      nrrd_header.AddKeyValue("wedgeIsovalue", s);
    }
  }

  if (field_param.flag_stack) {

    if (field_param.NumTranslate() == 1) {
      string coord_string;
      field_param.GetTranslateString(0, coord_string);
      nrrdKeyValueAdd
        (nrrd_header.DataPtr(), "translate", coord_string.c_str());
    }

    string num_obj_str;
    make_string("", field_param.NumObjects(), num_obj_str);
    nrrdKeyValueAdd
      (nrrd_header.DataPtr(), "numStackObj", num_obj_str.c_str());

  }

  if (field_info[field_index].flag_random_seed) {
    string seed_string;
    IJKGENSCALAR::convert_to_string(field_param.RandomSeed(), seed_string);
    nrrdKeyValueAdd
      (nrrd_header.DataPtr(), "randomSeed", seed_string.c_str());
  }

  if (field_info[field_index].flag_maxval) {
    string maxval_string;
    IJKGENSCALAR::convert_to_string(field_param.MaxVal(), maxval_string);
    nrrdKeyValueAdd
      (nrrd_header.DataPtr(), "maxValue", maxval_string.c_str());
  }

  if (field_name == "isotable_entry") {
    IJKGENSCALAR::convert_to_string(field_param.IsotableIndex(), s);
    nrrd_header.AddKeyValue("isotableIndex", s);
  }

}


void IJKGENSCALAR::write_scalar_grid
(const char * output_filename, const SCALAR_GRID & grid, 
 const FIELD_PARAM & field_param,  const std::vector<FIELD_INFO> & field_info,
 const bool flag_gzip)
{
  const int dimension = grid.Dimension();
  const int field_index = field_param.FieldIndex();
  NRRD_HEADER nrrd_header;
  IJK::ARRAY<double> grid_spacing(dimension, 1);

  // Store grid spacing in array of double.
  for (int d = 0; d < dimension; d++) 
    { grid_spacing[d] = grid.Spacing(d); }

  nrrd_header.SetSize(grid.Dimension(), grid.AxisSize());
  nrrdAxisInfoSet_nva(nrrd_header.DataPtr(), nrrdAxisInfoSpacing, 
                      grid_spacing.PtrConst());
  nrrdKeyValueAdd
    (nrrd_header.DataPtr(), "ijkDatasetType",
     field_info[field_index].name.c_str());

  add_nrrd_key_values(field_param, nrrd_header, field_info);

  if (flag_gzip) {
    write_scalar_grid_nrrd_gzip(output_filename, grid, nrrd_header);
  }
  else {
    write_scalar_grid_nrrd(output_filename, grid, nrrd_header);
  }
}

void IJKGENSCALAR::write_gradient_grid
(const string & output_filename, const GRADIENT_GRID & grid,
 const FIELD_PARAM & field_param, const std::vector<FIELD_INFO> & field_info,
 const bool flag_gzip)
{
  const int dimension = grid.Dimension();
  const int field_index = field_param.FieldIndex();
  const char * DISCONTINUITY_KEY = "handleGradientDiscontinuity";
  NRRD_HEADER nrrd_header;
  IJK::ARRAY<int> nrrd_axis_size(dimension+1);
  IJK::ARRAY<double> grid_spacing(dimension+1, 1);
  
  nrrd_axis_size[0] = grid.VectorLength();
  for (int d = 1; d <= dimension; d++) {
    nrrd_axis_size[d] = grid.AxisSize(d-1); 
    grid_spacing[d] = grid.Spacing(d-1);
  }

  nrrd_header.SetSize(dimension+1, nrrd_axis_size.PtrConst());
  nrrdAxisInfoSet_nva(nrrd_header.DataPtr(), nrrdAxisInfoSpacing, 
                      grid_spacing.PtrConst());

  std::string dataset_type = field_info[field_index].name;
  dataset_type += "Gradients";

  nrrdKeyValueAdd
    (nrrd_header.DataPtr(), "ijkDatasetType", dataset_type.c_str());

  add_nrrd_key_values(field_param, nrrd_header, field_info);

  if (field_param.gradient_discontinuity_zero) {
    nrrdKeyValueAdd
      (nrrd_header.DataPtr(), DISCONTINUITY_KEY, "setToZero");
  }
  else {
    nrrdKeyValueAdd
      (nrrd_header.DataPtr(), DISCONTINUITY_KEY, "selectOne");
  }

  if (flag_gzip) {
    write_vector_grid_nrrd_gzip(output_filename, grid, nrrd_header);
  }
  else {
    write_vector_grid_nrrd(output_filename, grid, nrrd_header);
  }
}

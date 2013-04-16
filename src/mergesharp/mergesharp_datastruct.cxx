/// \file mergesharp_datastruct.cxx
/// mergesharp data structures

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

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "mergesharp_types.h"
#include "mergesharp_datastruct.h"

using namespace IJK;
using namespace MERGESHARP;


// **************************************************
// DUAL CONTOURING ISOSURFACE CLASS
// **************************************************

void DUAL_ISOSURFACE::Clear()
{
  vertex_coord.clear();
  quad_vert.clear();
  tri_vert.clear();
}

// **************************************************
// CLASS MERGESHARP_PARAM
// **************************************************

// Initialize MERGESHARP_PARAM
void MERGESHARP_PARAM::Init()
{
  interpolation_type = LINEAR_INTERPOLATION;
  vertex_position_method = CENTROID_EDGE_ISO;
  flag_merge_sharp = true;
  quad_tri_method = UNDEFINED_TRI;
  flag_convert_quad_to_tri = false;
  use_only_cube_gradients = false;
  use_selected_gradients = true;
  grad_selection_cube_offset = 0;
  max_small_eigenvalue = 0.1;
  allow_multiple_iso_vertices = true;
  flag_resolve_ambiguous_facets = true;
  flag_separate_neg = true;
  flag_split_non_manifold = false;
  flag_delete_isolated_vertices = true;
  flag_store_isovert_info = false;
  flag_grad2hermite = false;
  flag_grad2hermiteI = false;
}

/// Set type of interpolation
void MERGESHARP_PARAM::SetInterpolationType
(const INTERPOLATION_TYPE interpolation_type)
{
  this->interpolation_type = interpolation_type;
}

/// Set isosurface vertex position method.
void MERGESHARP_PARAM::SetVertexPositionMethod
(const VERTEX_POSITION_METHOD vertex_position_method)
{
  this->vertex_position_method = vertex_position_method;
}

/// Set flag for using selected gradients
void MERGESHARP_PARAM::SetUseSelectedGradients
(const bool flag)
{
  this->use_selected_gradients = flag;
}

/// Set flag for using only cube gradients.
void MERGESHARP_PARAM::SetUseOnlyCubeGradients
(const bool flag)
{
  this->use_only_cube_gradients = flag;
}

/// Set flag for using gradients at endpoints of intersected edges.
void MERGESHARP_PARAM::SetUseIntersectedEdgeEndpointGradients
(const bool flag)
{
  this->use_intersected_edge_endpoint_gradients = flag;
}

/// Set
void MERGESHARP_PARAM::Set
(const MERGESHARP_PARAM & param)
{
  *this = param;
}

/// Return true if gradient data required.
bool MERGESHARP_PARAM::GradientsRequired() const
{
  if (VertexPositionMethod() == GRADIENT_POSITIONING ||
      VertexPositionMethod() == EDGEI_INTERPOLATE ||
      VertexPositionMethod() == EDGEI_GRADIENT ||
      flag_grad2hermite || flag_grad2hermiteI)
    { return(true); }
  else
    { return(false); }
}

/// Return true if edge-isosurface normal data required.
bool MERGESHARP_PARAM::NormalsRequired() const
{
  if (VertexPositionMethod() == EDGEI_INPUT_DATA) 
    { return(true); }
  else
    { return(false); }
}

// **************************************************
// CLASS MERGESHARP_DATA
// **************************************************

// Initialize MERGESHARP_DATA
void MERGESHARP_DATA::Init()
{
  is_scalar_grid_set = false;
  is_gradient_grid_set = false;
  are_edgeI_set = false;
}

void MERGESHARP_DATA::FreeAll()
{
  is_scalar_grid_set = false;
  is_gradient_grid_set = false;
}


// Copy scalar grid
void MERGESHARP_DATA::CopyScalarGrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid2)
{
  scalar_grid.Copy(scalar_grid2);
  scalar_grid.SetSpacing(scalar_grid2.SpacingPtrConst());
  is_scalar_grid_set = true;
}

// Copy gradient grid
void MERGESHARP_DATA::CopyGradientGrid
(const GRADIENT_GRID_BASE & gradient_grid2)
{
  gradient_grid.Copy(gradient_grid2);
  gradient_grid.SetSpacing(gradient_grid2.SpacingPtrConst());
  is_gradient_grid_set = true;
}

// Subsample scalar grid
void MERGESHARP_DATA::SubsampleScalarGrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid2, const int subsample_resolution)
{
  scalar_grid.Subsample(scalar_grid2, subsample_resolution);
  scalar_grid.SetSpacing(subsample_resolution, scalar_grid2.SpacingPtrConst());
  is_scalar_grid_set = true;
}

// Supersample scalar grid
void MERGESHARP_DATA::SupersampleScalarGrid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid2, 
 const int supersample_resolution)
{
  scalar_grid.Supersample(scalar_grid2, supersample_resolution);
  scalar_grid.SetSpacing(float(1.0/supersample_resolution),
                         scalar_grid2.SpacingPtrConst());
 is_scalar_grid_set = true;
}

/// Subsample gradient grid
/// Rescale gradients by multiplying them by the subsample_resolution.
void MERGESHARP_DATA::SubsampleGradientGrid
(const GRADIENT_GRID_BASE & gradient_grid2, const int subsample_resolution)
{
  gradient_grid.Subsample(gradient_grid2, subsample_resolution);
  gradient_grid.ScalarMultiply(subsample_resolution);
  gradient_grid.SetSpacing(subsample_resolution,
                           gradient_grid2.SpacingPtrConst());
  is_gradient_grid_set = true;
}

// Copy, subsample or supersample scalar grid.
void MERGESHARP_DATA::SetScalarGrid
(const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid,
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution)
{
  PROCEDURE_ERROR error("MERGESHARP_DATA::SetScalarGrid");

  if (flag_subsample && flag_supersample) {
    error.AddMessage("Scalar grid cannot both be subsampled and supersampled.");
    throw error;
  }

  if (flag_subsample) {
    // subsample grid
    SubsampleScalarGrid(full_scalar_grid, subsample_resolution);
  }
  else if (flag_supersample) {
    // supersample grid
    SupersampleScalarGrid(full_scalar_grid, supersample_resolution);
  }
  else {
    CopyScalarGrid(full_scalar_grid);
  };
}

// Copy, subsample or supersample scalar and gradient grids.
void MERGESHARP_DATA::SetGrids
(const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid,
 const GRADIENT_GRID_BASE & full_gradient_grid,
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution)
{
  PROCEDURE_ERROR error("MERGESHARP_DATA::SetGrids");

  if (flag_subsample && flag_supersample) {
    error.AddMessage("Scalar grid cannot both be subsampled and supersampled.");
    throw error;
  }

  if (flag_subsample) {
    SubsampleScalarGrid(full_scalar_grid, subsample_resolution);
    SubsampleGradientGrid(full_gradient_grid, subsample_resolution);
  }
  else if (flag_supersample) {
    error.AddMessage("Supersampling of gradient grid is not yet implemented.");
    throw error;
  }
  else {
    CopyScalarGrid(full_scalar_grid);
    CopyGradientGrid(full_gradient_grid);
  };

}

/// Set edge-isosurface intersections and normals.
void MERGESHARP_DATA::SetEdgeI
(const std::vector<COORD_TYPE> & edgeI_coord,
 const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord)
{
  IJK::PROCEDURE_ERROR error("MERGESHARP_DATA::SetEdgeI");

  if (edgeI_coord.size() != edgeI_normal_coord.size()) {
    error.AddMessage
      ("Programming error.  Number of normal coordinates does not");
    error.AddMessage
      ("  match number of edge-isosurface intersection coordinates.");
    error.AddMessage("  Number of edge-isosurface intersection coordinates = ",
                     edgeI_coord.size(), ".");
    error.AddMessage("  Number of normal coordinates = ",
                     edgeI_normal_coord.size(), ".");
    throw error;
  }

  this->edgeI_coord.resize(edgeI_coord.size());
  this->edgeI_normal_coord.resize(edgeI_normal_coord.size());
  std::copy(edgeI_coord.begin(), edgeI_coord.end(), this->edgeI_coord.begin());
  std::copy(edgeI_normal_coord.begin(), edgeI_normal_coord.end(), 
            this->edgeI_normal_coord.begin());

  are_edgeI_set = true;
}

/// Check data structure
bool MERGESHARP_DATA::Check(IJK::ERROR & error) const
{
  IJK::ERROR error2;

  if (!IsScalarGridSet()) {
    error.AddMessage("Scalar grid is not set.");
    return(false);
  }

  return(true);
}


// **************************************************
// MERGESHARP TIME
// **************************************************

MERGESHARP::MERGESHARP_TIME::MERGESHARP_TIME()
// constructor
{
  Clear();
}

void MERGESHARP::MERGESHARP_TIME::Clear()
{
  preprocessing = 0.0;
  extract = 0.0;
  merge_identical = 0.0;
  position = 0.0;
  merge_sharp = 0.0;
  total = 0.0;
}

void MERGESHARP::MERGESHARP_TIME::Add(const MERGESHARP_TIME & isodual_time)
{
  preprocessing += isodual_time.preprocessing;
  extract += isodual_time.extract;
  merge_identical += isodual_time.merge_identical;
  position += isodual_time.position;
  merge_sharp += isodual_time.merge_sharp;
  total += isodual_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

MERGESHARP::GRID_INFO::GRID_INFO()
{
  Clear();
}

void MERGESHARP::GRID_INFO::Clear()
{
  num_cubes = 0;
}

void MERGESHARP::SCALAR_INFO::Init(const int dimension)
{
  this->dimension = 0;
  SetDimension(dimension);

  Clear();
}

void MERGESHARP::SCALAR_INFO::FreeAll()
{
  dimension = 0;
}

void MERGESHARP::SCALAR_INFO::Clear()
{
  num_non_empty_cubes = 0;
  num_bipolar_edges = 0;
}

void MERGESHARP::SCALAR_INFO::SetDimension(const int dimension)
{
  FreeAll();

  this->dimension = dimension;

  Clear();
}

void MERGESHARP::SCALAR_INFO::Copy(const SCALAR_INFO & info)
{
  Init(info.Dimension());
  num_non_empty_cubes = info.num_non_empty_cubes;
  num_bipolar_edges = info.num_bipolar_edges;
}

/// Copy assignment.
const SCALAR_INFO &  MERGESHARP::SCALAR_INFO::operator =
(const SCALAR_INFO & right)
{
  if (&right != this) {
    FreeAll();
    Copy(right);
  }

  return *this;
}

MERGESHARP::SCALAR_INFO::~SCALAR_INFO()
{
  dimension = 0;
  Clear();
}

MERGESHARP::SHARPISO_INFO::SHARPISO_INFO()
{
  Clear();
}

MERGESHARP::MERGESHARP_INFO::MERGESHARP_INFO()
{
  Clear();
}

void MERGESHARP::SHARPISO_INFO::Clear()
{
  ISOVERT_INFO::Clear();

  /* OBSOLETE
  num_conflicts = 0;
  num_merged_iso_vertices = 0;
  num_sharp_corners = 0;
  num_sharp_edges = 0;
  num_smooth_vertices = 0;
  num_Linf_iso_vertex_locations = 0;
  */

  num_edge_collapses = 0;
  num_repositioned_vertices = 0;

  num_cube_not_ambiguous = 0;
  num_cube_separate_pos = 0;
  num_cube_separate_neg = 0;
  num_cube_unresolved_ambiguity = 0;

  num_cube_single_isov = 0;
  num_cube_multi_isov = 0;
  num_non_disk_isopatches = 0;
  num_non_manifold_split = 0;

  vertex_info.clear();
}

// Increment num_sharp_corners or num_sharp_edges or num_smooth_vertices
//   depending upon the number of large eigenvalues.
void MERGESHARP::SHARPISO_INFO::
IncrementIsoVertexNum(const int num_large_eigenvalues)
{
  if (num_large_eigenvalues >= 3) { num_sharp_corners++; }
  else if (num_large_eigenvalues == 2) { num_sharp_edges++; }
  else { num_smooth_vertices++; }
}


MERGESHARP::MERGESHARP_INFO::MERGESHARP_INFO
(const int dimension):scalar(dimension)
{
  Clear();
}

void MERGESHARP::MERGESHARP_INFO::Clear()
{
  grid.Clear();
  scalar.Clear();
  time.Clear();
  sharpiso.Clear();
}

// **************************************************
// DUAL_ISOVERT_INFO
// **************************************************

void MERGESHARP::DUAL_ISOVERT_INFO::Set(const MERGESHARP::DUAL_ISOVERT & isov)
{
  cube_index = isov.cube_index;
  patch_index = isov.patch_index;
  table_index = isov.table_index;
}

// **************************************************
// SET ISOVERT INFORMATION
// **************************************************

void MERGESHARP::set_isovert_info
(const std::vector<MERGESHARP::DUAL_ISOVERT> & iso_vlist,
 const std::vector<GRID_CUBE> & gcube_list,
 std::vector<DUAL_ISOVERT_INFO> & isovert_info)
{
  IJK::PROCEDURE_ERROR error("set_isovert_info");
  CUBE_HASH_TABLE cube_hash;

  isovert_info.resize(iso_vlist.size());
  for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {
    VERTEX_INDEX cube_index = gcube_list[i].cube_index;
    cube_hash.insert(CUBE_HASH_TABLE::value_type(cube_index, i));
  }

  for (NUM_TYPE i = 0; i < iso_vlist.size(); i++) {
    isovert_info[i].Set(iso_vlist[i]);

    VERTEX_INDEX cube_index = isovert_info[i].cube_index;

    CUBE_HASH_TABLE::iterator cube_iter = cube_hash.find(cube_index);

    if (cube_iter == cube_hash.end()) {
      error.AddMessage("Programming error. Illegal isosurface vertex cube.");
      error.AddMessage("  Isosurface vertex: ", i,
                       ".  Cube index: ", cube_index, ".");
      throw error;
    }

    VERTEX_INDEX gcube_index = cube_iter->second;

    isovert_info[i].num_eigenvalues = gcube_list[gcube_index].num_eigenvalues;
    isovert_info[i].flag_centroid_location =
      gcube_list[gcube_index].flag_centroid_location;
  }

}

// Delete vertices i where flag_keep[i] = false.
void MERGESHARP::delete_vertices
(const std::vector<MERGESHARP::DUAL_ISOVERT> & iso_vlist,
 const std::vector<NUM_TYPE> & new_isovert_index,
 const std::vector<bool> & flag_keep,
 std::vector<MERGESHARP::DUAL_ISOVERT> & new_iso_vlist)
{
  IJK::PROCEDURE_ERROR error("delete_vertices");
  new_iso_vlist.clear();

  NUM_TYPE numv = 0;
  for (NUM_TYPE i = 0; i < flag_keep.size(); i++) {
    if (flag_keep[i]) { numv++; }
  }

  IJK::ARRAY<bool> is_set(numv, false);
  new_iso_vlist.resize(numv);

  for (NUM_TYPE i = 0; i < flag_keep.size(); i++) {
    if (flag_keep[i]) {
      NUM_TYPE j = new_isovert_index[i];
      is_set[j] = true;
      new_iso_vlist[j] = iso_vlist[i];
    }
  }

  for (NUM_TYPE j = 0; j < new_iso_vlist.size(); j++) {
    if (!is_set[j]) {
      error.AddMessage
        ("Programming error.  Reference to new vertex ", j,
         " missing from new_isovert_index.");
      throw error;
    }
  }

}


// **************************************************
// MERGE DATA
// **************************************************

void MERGESHARP::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  this->num_obj_per_grid_vertex =
    dimension*num_obj_per_edge + num_obj_per_vertex;
  compute_num_grid_vertices(dimension, axis_size, num_vertices);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj =
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX,MERGE_INDEX>::Init(num_obj);
}

bool MERGESHARP::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxNumInt() <
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}


/// \file isodual3D_position.cxx
/// Position dual isosurface vertices

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

#include "ijkcoord.txx"
#include "ijkcube.txx"
#include "ijkgrid_macros.h"
#include "ijkinterpolate.txx"
#include "ijkisopoly.txx"

#include "ijkdualtable.h"
#include "ijkdualtable_ambig.h"

#include "isodual3D_position.h"

#include "sharpiso_types.h"
#include "sharpiso_feature.h"

#include <iostream>
#include <iomanip>


using namespace ISODUAL3D;
using namespace SHARPISO;
using namespace std;

// **************************************************
// Position at cube centers
// **************************************************

/// Position dual isosurface vertices in cube centers
void ISODUAL3D::position_dual_isovertices_cube_center
(const SHARPISO_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
{
  const int dimension = grid.Dimension();
  GRID_COORD_TYPE grid_coord[dimension];
  COORD_TYPE unit_cube_center[dimension];
  IJK::PROCEDURE_ERROR error("position_dual_isovertices_center");

  IJK::set_coord(dimension, 0.5, unit_cube_center);

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    grid.ComputeCoord(iv, grid_coord);
    IJK::add_coord(dimension, grid_coord, unit_cube_center, coord+i*dimension);
  }
}

/// Position dual isosurface vertices in cube centers
void ISODUAL3D::position_dual_isovertices_cube_center
(const SHARPISO_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_cube_center(grid, vlist, &(coord.front()));
}

// **************************************************
// Position at centroid of isosurface-edge intersections
// **************************************************

/// Position dual isosurface vertices in centroid
///   of isosurface-edge intersections
void ISODUAL3D::position_dual_isovertices_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    compute_isosurface_grid_edge_centroid
      (scalar_grid, isovalue, iv, coord+i*dimension);
  }
}


/// Position dual isosurface vertices in centroid
/// of isosurface-edge intersections
void ISODUAL3D::position_dual_isovertices_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)

{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_centroid
    (scalar_grid, isovalue, vlist, &(coord.front()));
}

/// Position dual isosurface vertices in centroid
///   of isosurface-edge intersections.
/// Allow multiple vertices in a grid cube.
void ISODUAL3D::position_dual_isovertices_centroid_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  ISODUAL3D_CUBE_FACE_INFO  cube(dimension);
  IJKDUALTABLE::TABLE_INDEX it;

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    it = iso_vlist[i].table_index;

    if (isodual_table.NumIsoVertices(it) == 1) {

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isovalue, icube, coord+i*dimension);
    }
    else {
      FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isodual_table, isovalue, icube, ipatch,
         it, cube, coord+i*DIM3);
    }
  }
}

/// Position dual isosurface vertices in centroid
///   of isosurface-edge intersections.
/// Allowing multiple vertices in a grid cube.
/// Version using std::vector for array coord[].
void ISODUAL3D::position_dual_isovertices_centroid_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_centroid_multi
    (scalar_grid, isodual_table, isovalue, iso_vlist, &(coord.front()));
}

// **************************************************
// Position using gradients and svd
// **************************************************

/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord,
 SHARPISO_INFO & sharp_info)
{

  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isodual_param.grad_selection_cube_offset;

  SVD_INFO svd_info;
  svd_info.location = LOC_NONE;

  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};

  VERTEX_INDEX num_large_eigenvalues;

  OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    svd_compute_sharp_vertex_for_cube
      (scalar_grid, gradient_grid, iv, isovalue,
      isodual_param, cube_111, sharp_coord+i*DIM3,
      eigenvalues, num_large_eigenvalues, svd_info);

    if (svd_info.flag_conflict) 
      { sharp_info.num_conflicts++; }
    else if (svd_info.flag_Linf_iso_vertex_location)
      { sharp_info.num_Linf_iso_vertex_locations++; }

    sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
  }

}


/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_using_gradients
    (scalar_grid, gradient_grid, isovalue, isodual_param,
     vlist, &(coord.front()), sharp_info);
}


/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 COORD_TYPE * sharp_coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isodual_param.grad_selection_cube_offset;
  ISODUAL3D_CUBE_FACE_INFO  cube(dimension);
  IJKDUALTABLE::TABLE_INDEX it;

  SVD_INFO svd_info;
  svd_info.location = LOC_NONE;

  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};

  VERTEX_INDEX num_large_eigenvalues;

  OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    VERTEX_INDEX it = iso_vlist[i].table_index;

    if (isodual_table.NumIsoVertices(it) == 1) {

      svd_compute_sharp_vertex_for_cube
        (scalar_grid, gradient_grid, icube, isovalue,
         isodual_param, cube_111, sharp_coord+i*DIM3,
         eigenvalues, num_large_eigenvalues, svd_info);

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      else if (svd_info.flag_Linf_iso_vertex_location)
        { sharp_info.num_Linf_iso_vertex_locations++; }

      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
    else {
      FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isodual_table, isovalue, icube, ipatch,
         it, cube, sharp_coord+i*DIM3);
    }
  }

}

/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<COORD_TYPE> & sharp_coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();

  sharp_coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_using_gradients_multi
    (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param,
     iso_vlist, &(sharp_coord.front()), sharp_info);
}

/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
 const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
 const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
 COORD_TYPE * sharp_coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  const SIGNED_COORD_TYPE grad_selection_cube_offset =
    isodual_param.grad_selection_cube_offset;
  ISODUAL3D_CUBE_FACE_INFO  cube(dimension);

  SVD_INFO svd_info;
  svd_info.location = LOC_NONE;

  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};

  VERTEX_INDEX num_large_eigenvalues;

  OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);

  for (VERTEX_INDEX i = 0; i < iso_vlist_cube.size(); i++) {
    VERTEX_INDEX icube = iso_vlist_cube[i];

    IJKDUALTABLE::TABLE_INDEX it;
    IJK::compute_isotable_index
      (scalar_grid.ScalarPtrConst(), isovalue, icube,
       scalar_grid.CubeVertexIncrement(), num_cube_vertices, it);

    if (iso_vlist_cube_ambig[i] == SEPARATE_POS) 
      // Complement table entry.
      { it = isodual_table.NumTableEntries()-1 - it; }

    if (isodual_table.NumIsoVertices(it) == 1) {

      svd_compute_sharp_vertex_for_cube
        (scalar_grid, gradient_grid, icube, isovalue,
         isodual_param, cube_111, sharp_coord+i*DIM3,
         eigenvalues, num_large_eigenvalues, svd_info);

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      else if (svd_info.flag_Linf_iso_vertex_location)
        { sharp_info.num_Linf_iso_vertex_locations++; }

      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
    else {
      FACET_VERTEX_INDEX ipatch= iso_vlist_patch[i];

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isodual_table, isovalue, icube, ipatch,
         it, cube, sharp_coord+i*DIM3);
    }
  }

}

/// Position dual isosurface vertices using gradients
/// Return list of vertex conflicts.
/// Version using std::vector for array coord[].
void ISODUAL3D::position_dual_isovertices_using_gradients
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
 const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
 const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
 std::vector<COORD_TYPE> & sharp_coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();

  sharp_coord.resize(iso_vlist_cube.size()*dimension);
  position_dual_isovertices_using_gradients
    (scalar_grid, gradient_grid, isodual_table, isovalue,
     isodual_param, iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig,
     &(sharp_coord.front()), sharp_info);
}


// **************************************************
// Position using isovert information
// **************************************************

// Position dual isosurface vertices using isovert information.
// Allows multiple vertices in a grid cube.
void ISODUAL3D::position_dual_isovertices_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 COORD_TYPE * isov_coord)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  ISODUAL3D_CUBE_FACE_INFO cube(dimension);
  IJKDUALTABLE::TABLE_INDEX it;

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {

    VERTEX_INDEX cube_index = iso_vlist[i].cube_index;
    VERTEX_INDEX gcube_index = isovert.sharp_ind_grid.Scalar(cube_index);

    if (isovert.gcube_list[gcube_index].flag != SELECTED_GCUBE) {

      it = iso_vlist[i].table_index;


      if (isodual_table.NumIsoVertices(it) == 1) {
        if (isovert.gcube_list[gcube_index].flag == SMOOTH_GCUBE) {

          IJK::copy_coord_3D(isovert.gcube_list[gcube_index].isovert_coord,
                             isov_coord+i*DIM3);
        }
        else {
          compute_isosurface_grid_edge_centroid
            (scalar_grid, isovalue, cube_index, isov_coord+i*DIM3);
        }
      }
      else {
        FACET_VERTEX_INDEX ipatch= iso_vlist[i].patch_index;

        compute_isosurface_grid_edge_centroid
          (scalar_grid, isodual_table, isovalue, cube_index, ipatch,
           it, cube, isov_coord+i*DIM3);
      }
    }
    else {
      IJK::copy_coord_3D(isovert.gcube_list[gcube_index].isovert_coord,
                         isov_coord+i*DIM3);
    }
  }

}

// Position dual isosurface vertices using isovert information.
// Allows multiple vertices in a grid cube.
void ISODUAL3D::position_dual_isovertices_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<COORD_TYPE> & isov_coord)
{
  const int dimension = scalar_grid.Dimension();

  isov_coord.resize(iso_vlist.size()*dimension);
  position_dual_isovertices_multi
    (scalar_grid, isodual_table, isovalue, isovert,
     iso_vlist, &(isov_coord.front()));
}

// ********************************************************
// Position vertices using SVD on grid edge-isosurface intersections.
// ********************************************************

/// Position vertices using SVD on grid edge-isosurface intersections.
void ISODUAL3D::position_dual_isovertices_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const VERTEX_POSITION_METHOD position_method,
 COORD_TYPE * sharp_coord,
 SHARPISO_INFO & sharp_info)
{
  SVD_INFO svd_info;
  VERTEX_INDEX num_large_eigenvalues;
  EIGENVALUE_TYPE eigenvalues[DIM3];

  if (position_method == EDGEI_GRADIENT) {
    for (VERTEX_INDEX i = 0; i<vlist.size(); i++) {
      VERTEX_INDEX iv = vlist[i];

      svd_compute_sharp_vertex_edgeI_sharp_gradient
        (scalar_grid, gradient_grid, iv, isovalue,
         isodual_param, sharp_coord+i*DIM3,
         eigenvalues, num_large_eigenvalues, svd_info);

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
  }
  else {

    for (VERTEX_INDEX i = 0; i<vlist.size(); i++) {
      VERTEX_INDEX iv = vlist[i];

      svd_compute_sharp_vertex_edgeI_interpolate_gradients
        (scalar_grid, gradient_grid, iv, isovalue,
         isodual_param, sharp_coord+i*DIM3,
         eigenvalues, num_large_eigenvalues, svd_info);

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
  }

}

// Position vertices using SVD on grid edge-isosurface intersections.
// Version using std::vector for array coord[].
void ISODUAL3D::position_dual_isovertices_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const VERTEX_POSITION_METHOD position_method,
 std::vector<COORD_TYPE> & coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();
  coord.resize(vlist.size()*dimension);

  position_dual_isovertices_edgeI
    (scalar_grid, gradient_grid, isovalue, isodual_param, vlist,
     position_method, &(coord.front()), sharp_info);
}


/// Position vertices using SVD on grid edge-isosurface intersections.
/// Use isodual_table to determine connections.
void ISODUAL3D::position_dual_isovertices_edgeI_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const VERTEX_POSITION_METHOD position_method,
 COORD_TYPE * coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  ISODUAL3D_CUBE_FACE_INFO  cube(dimension);
  IJKDUALTABLE::TABLE_INDEX it;

  for (VERTEX_INDEX i = 0; i < iso_vlist.size(); i++) {
    VERTEX_INDEX icube = iso_vlist[i].cube_index;
    it = iso_vlist[i].table_index;

    if (isodual_table.NumIsoVertices(it) == 1) {

      SVD_INFO svd_info;
      svd_info.location = LOC_NONE;

      EIGENVALUE_TYPE eigenvalues[DIM3];
      VERTEX_INDEX num_large_eigenvalues;

      if (position_method == EDGEI_GRADIENT) {
        svd_compute_sharp_vertex_edgeI_sharp_gradient
          (scalar_grid, gradient_grid, icube, isovalue, isodual_param,
           coord+i*dimension, eigenvalues,  num_large_eigenvalues, svd_info);
      }
      else {
        svd_compute_sharp_vertex_edgeI_interpolate_gradients
          (scalar_grid, gradient_grid, icube, isovalue, isodual_param,
           coord+i*dimension, eigenvalues,  num_large_eigenvalues, svd_info);
      }

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
    else {
      FACET_VERTEX_INDEX ipatch = iso_vlist[i].patch_index;

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isodual_table, isovalue, icube, ipatch,
         it, cube, coord+i*DIM3);
    }

  }

}

/// Position vertices using SVD on grid edge-isosurface intersections.
/// Use isodual_table to determine connections.
/// Version using std::vector for array coord[].
void ISODUAL3D::position_dual_isovertices_edgeI_multi
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<DUAL_ISOVERT> & iso_vlist,
 const VERTEX_POSITION_METHOD position_method,
 std::vector<COORD_TYPE> & coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist.size()*dimension);

  position_dual_isovertices_edgeI_multi
    (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param, 
     iso_vlist, position_method, &(coord.front()), sharp_info);
}

/// Position vertices using SVD on grid edge-isosurface intersections.
/// Use isodual_table to determine connections.
/// Use iso_vlist_cube_ambig[] to resolve ambiguities.
void ISODUAL3D::position_dual_isovertices_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
 const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
 const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
 const VERTEX_POSITION_METHOD position_method,
 COORD_TYPE * coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();
  const int num_cube_vertices = scalar_grid.NumCubeVertices();
  ISODUAL3D_CUBE_FACE_INFO  cube(dimension);

  for (VERTEX_INDEX i = 0; i < iso_vlist_cube.size(); i++) {
    VERTEX_INDEX icube = iso_vlist_cube[i];

    IJKDUALTABLE::TABLE_INDEX it;
    IJK::compute_isotable_index
      (scalar_grid.ScalarPtrConst(), isovalue, icube,
       scalar_grid.CubeVertexIncrement(), num_cube_vertices, it);

    if (iso_vlist_cube_ambig[i] == SEPARATE_POS) 
      // Complement table entry.
      { it = isodual_table.NumTableEntries()-1 - it; }

    if (isodual_table.NumIsoVertices(it) == 1 ||
        (iso_vlist_cube_ambig[i] != SEPARATE_POS &&
         iso_vlist_cube_ambig[i] != SEPARATE_NEG)) {

      SVD_INFO svd_info;
      svd_info.location = LOC_NONE;

      EIGENVALUE_TYPE eigenvalues[DIM3];
      VERTEX_INDEX num_large_eigenvalues;

      if (position_method == EDGEI_GRADIENT) {
        svd_compute_sharp_vertex_edgeI_sharp_gradient
          (scalar_grid, gradient_grid, icube, isovalue, isodual_param,
           coord+i*dimension, eigenvalues, num_large_eigenvalues, svd_info);
      }
      else {
        svd_compute_sharp_vertex_edgeI_interpolate_gradients
          (scalar_grid, gradient_grid, icube, isovalue, isodual_param,
           coord+i*dimension, eigenvalues, num_large_eigenvalues, svd_info);
      }

      if (svd_info.flag_conflict) 
        { sharp_info.num_conflicts++; }
      sharp_info.IncrementIsoVertexNum(num_large_eigenvalues);
    }
    else {
      FACET_VERTEX_INDEX ipatch= iso_vlist_patch[i];

      compute_isosurface_grid_edge_centroid
        (scalar_grid, isodual_table, isovalue, icube, ipatch,
         it, cube, coord+i*DIM3);
    }
  }

}

/// Position vertices using SVD on grid edge-isosurface intersections.
/// Use isodual_table to determine connections.
/// Use iso_vlist_cube_ambig[] to resolve ambiguities.
/// Version using std::vector for array coord[].
void ISODUAL3D::position_dual_isovertices_edgeI
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
 const std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
 const std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
 const VERTEX_POSITION_METHOD position_method,
 std::vector<COORD_TYPE> & coord,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(iso_vlist_cube.size()*dimension);

  position_dual_isovertices_edgeI
    (scalar_grid, gradient_grid, isodual_table, isovalue, isodual_param, 
     iso_vlist_cube, iso_vlist_patch, iso_vlist_cube_ambig, position_method,
     &(coord.front()), sharp_info);
}


// **************************************************
// Compute routines
// **************************************************

/// Compute centroid of intersections of isosurface and grid edges
void ISODUAL3D::compute_isosurface_grid_edge_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX iv,
 COORD_TYPE * coord)
{
  const int dimension = scalar_grid.Dimension();
  GRID_COORD_TYPE grid_coord[dimension];
  COORD_TYPE vcoord[dimension];
  COORD_TYPE coord0[dimension];
  COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];

  IJK::PROCEDURE_ERROR error("compute_isosurface_grid_edge_centroid");

  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
    for (int k = 0; k < scalar_grid.NumFacetVertices(); k++) {
      VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
      VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0, edge_dir);

      SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
      bool is_end0_positive = true;
      if (s0 < isovalue)
      { is_end0_positive = false; };

      SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);
      bool is_end1_positive = true;
      if (s1 < isovalue)
      { is_end1_positive = false; };

      if (is_end0_positive != is_end1_positive) {

        scalar_grid.ComputeCoord(iend0, coord0);
        scalar_grid.ComputeCoord(iend1, coord1);

        IJK::linear_interpolate_coord
          (dimension, s0, coord0, s1, coord1, isovalue, coord2);

        IJK::add_coord(dimension, vcoord, coord2, vcoord);

        num_intersected_edges++;
      }
    }

    if (num_intersected_edges > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCoord(iv, vcoord);
      for (int d = 0; d < dimension; d++)
      { vcoord[iv] += 0.5; };
    }

    IJK::copy_coord(dimension, vcoord, coord);
}

/// Compute centroid of intersections of isosurface and cube edges.
/// Use only grid cube edges associated with isosurface patch ipatch.
void ISODUAL3D::compute_isosurface_grid_edge_centroid
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX icube,
 const FACET_VERTEX_INDEX ipatch, const IJKDUALTABLE::TABLE_INDEX it,
 const ISODUAL3D_CUBE_FACE_INFO & cube,
 COORD_TYPE * coord)
{
  static COORD_TYPE vcoord[DIM3];
  static COORD_TYPE coord0[DIM3];
  static COORD_TYPE coord1[DIM3];
  static COORD_TYPE coord2[DIM3];

  int num_intersected_edges = 0;
  IJK::set_coord_3D(0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    if (isodual_table.IsBipolar(it, ie)) {
      if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
        int k0 = cube.EdgeEndpoint(ie, 0);
        int k1 = cube.EdgeEndpoint(ie, 1);
        VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
        VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);
        SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
        SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

        scalar_grid.ComputeCoord(iend0, coord0);
        scalar_grid.ComputeCoord(iend1, coord1);

        IJK::linear_interpolate_coord
          (DIM3, s0, coord0, s1, coord1, isovalue, coord2);

        IJK::add_coord_3D(vcoord, coord2, vcoord);

        num_intersected_edges++;
      }
    }
  }

  if (num_intersected_edges > 0) {
    IJK::multiply_coord_3D
      (1.0/num_intersected_edges, vcoord, coord);
  }
  else {
    scalar_grid.ComputeCubeCenterCoord(icube, coord);
  }

}

// **************************************************
// Split routine
// **************************************************

// Split dual isosurface vertices.
// @param isodual_table Dual isosurface lookup table.
void ISODUAL3D::split_dual_isovert_ambig
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & cube_list,
 const std::vector<AMBIGUITY_TYPE> & cube_ambig,
 const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
 const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 std::vector<ISO_VERTEX_INDEX> & iso_vlist_cube,
 std::vector<FACET_VERTEX_INDEX> & iso_vlist_patch,
 std::vector<AMBIGUITY_TYPE> & iso_vlist_cube_ambig,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 VERTEX_INDEX & num_split)
{
  const int dimension = scalar_grid.Dimension();
  const NUM_TYPE num_cube_vertices = scalar_grid.NumCubeVertices();
  const NUM_TYPE num_facet_vertices = scalar_grid.NumFacetVertices();
  std::vector<IJKDUALTABLE::TABLE_INDEX> cube_table_index;
  std::vector<FACET_VERTEX_INDEX> num_isov;
  std::vector<ISO_VERTEX_INDEX> first_cube_isov;
  num_isov.resize(cube_list.size());
  first_cube_isov.resize(cube_list.size());
  IJK::CUBE_FACE_INFO<int,NUM_TYPE,NUM_TYPE> cube(dimension);

  num_split = 0;

  cube_table_index.resize(cube_list.size());
  ISO_VERTEX_INDEX total_num_isov = 0;
  for (ISO_VERTEX_INDEX i = 0; i < cube_list.size(); i++) {
    IJKDUALTABLE::TABLE_INDEX it;
    IJK::compute_isotable_index
      (scalar_grid.ScalarPtrConst(), isovalue, cube_list[i],
       scalar_grid.CubeVertexIncrement(), num_cube_vertices, it);

    if (cube_ambig[i] == SEPARATE_POS) 
      // Complement table entry.
      { it = isodual_table.NumTableEntries() - it - 1; }
    cube_table_index[i] = it;
    num_isov[i] = isodual_table.NumIsoVertices(it);
    total_num_isov += num_isov[i];

    if (num_isov[i] > 1) { num_split++; }
  }

  iso_vlist_cube.resize(total_num_isov);
  iso_vlist_patch.resize(total_num_isov);
  iso_vlist_cube_ambig.resize(total_num_isov);

  ISO_VERTEX_INDEX k = 0;
  for (ISO_VERTEX_INDEX i = 0; i < cube_list.size(); i++) {
    first_cube_isov[i] = k;
    for (FACET_VERTEX_INDEX j = 0; j < num_isov[i]; j++) {
      iso_vlist_cube[k+j] = cube_list[i];
      iso_vlist_patch[k+j] = j;
      iso_vlist_cube_ambig[k+j] = cube_ambig[i];
    }
    k = k+num_isov[i];
  }

  isoquad_vert.resize(isoquad_cube.size());

  for (ISO_VERTEX_INDEX i = 0; i < isoquad_cube.size(); i++) {
    ISO_VERTEX_INDEX k = isoquad_cube[i];
    IJKDUALTABLE::TABLE_INDEX it = cube_table_index[k];

    // Compute index of facet vertex opposite to facet_vertex[i]
    int facet_vertex_i = facet_vertex[i];
    int ifacet = cube.FacetIndex(facet_vertex_i);
    int j = facet_vertex_i - ifacet*num_facet_vertices;
    int opposite_vertex = (num_facet_vertices-1) - j;
    opposite_vertex += (ifacet*num_facet_vertices);

    isoquad_vert[i] = first_cube_isov[k] + 
      isodual_table.IncidentIsoVertex(it, opposite_vertex);

  }

}

// Split dual isosurface vertices.
// Don't split isosurface vertices in cubes merged with other cubes.
// @param isodual_table Dual isosurface lookup table.
void ISODUAL3D::split_dual_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
 const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 const ISODUAL_PARAM & isodual_param,
 std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = isodual_table.Dimension();
  const NUM_TYPE num_gcube = isovert.gcube_list.size();
  std::vector<ISO_VERTEX_INDEX> cube_list(num_gcube);
  std::vector<bool> no_split(num_gcube,true);

  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) {
    cube_list[i] = isovert.gcube_list[i].cube_index;
    if (isovert.gcube_list[i].flag == SMOOTH_GCUBE ||
        isovert.gcube_list[i].flag == UNAVAILABLE_GCUBE ||
        isovert.gcube_list[i].flag == NON_DISK_GCUBE)
      { no_split[i] = false; }
  }

  int num_split;
  if (isodual_param.flag_split_non_manifold) {
    IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO ambig_info(dimension);
    int num_non_manifold_split;

    IJK::split_dual_isovert_manifold
      (scalar_grid, isodual_table, ambig_info, isovalue, 
       cube_list, no_split, isoquad_cube, facet_vertex,
       iso_vlist, isoquad_vert, num_split, num_non_manifold_split);
    sharp_info.num_non_manifold_split = num_non_manifold_split;
  }
  else {
    IJK::split_dual_isovert
      (scalar_grid, isodual_table, isovalue, cube_list, no_split, 
       isoquad_cube, facet_vertex,
       iso_vlist, isoquad_vert, num_split);
  }
  sharp_info.num_cube_multi_isov = num_split;
  sharp_info.num_cube_single_isov = num_gcube - num_split;
}


// Split dual isosurface vertices.
// Split all isosurface vertices as determined by the lookup table
//   of the non-manifold edges (if flag_split_non_manifold is true.)
// @param isodual_table Dual isosurface lookup table.
void ISODUAL3D::full_split_dual_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
 const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 const ISODUAL_PARAM & isodual_param,
 std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 std::vector<CUBE_ISOVERT_DATA> & cube_isovert_data,
 SHARPISO_INFO & sharp_info)
{
  const int dimension = isodual_table.Dimension();
  const NUM_TYPE num_gcube = isovert.gcube_list.size();

  std::vector<ISO_VERTEX_INDEX> cube_list(num_gcube);

  for (NUM_TYPE i = 0; i < isovert.gcube_list.size(); i++) 
    { cube_list[i] = isovert.gcube_list[i].cube_index; }

  int num_split;
  if (isodual_param.flag_split_non_manifold) {
    IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG_INFO ambig_info(dimension);
    int num_non_manifold_split;

    IJK::split_dual_isovert_manifold
      (scalar_grid, isodual_table, ambig_info, isovalue, 
       cube_list, isoquad_cube, facet_vertex,
       iso_vlist, isoquad_vert, cube_isovert_data,
       num_split, num_non_manifold_split);
    sharp_info.num_non_manifold_split = num_non_manifold_split;
  }
  else {
    IJK::split_dual_isovert
      (scalar_grid, isodual_table, isovalue, cube_list,
       isoquad_cube, facet_vertex, iso_vlist, isoquad_vert, 
       cube_isovert_data, num_split);
  }
  sharp_info.num_cube_multi_isov = num_split;
  sharp_info.num_cube_single_isov = num_gcube - num_split;
}

// Split dual isosurface vertices.
// Split all isosurface vertices as determined by the lookup table
//   of the non-manifold edges (if flag_split_non_manifold is true.)
// @param isodual_table Dual isosurface lookup table.
void ISODUAL3D::full_split_dual_isovert
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const IJKDUALTABLE::ISODUAL_CUBE_TABLE & isodual_table,
 const SCALAR_TYPE isovalue,
 const ISOVERT & isovert,
 const std::vector<ISO_VERTEX_INDEX> & isoquad_cube,     
 const std::vector<FACET_VERTEX_INDEX> & facet_vertex,
 const ISODUAL_PARAM & isodual_param,
 std::vector<DUAL_ISOVERT> & iso_vlist,
 std::vector<VERTEX_INDEX> & isoquad_vert,
 SHARPISO_INFO & sharp_info)
{
  std::vector<CUBE_ISOVERT_DATA> cube_isovert_data;

  full_split_dual_isovert
    (scalar_grid, isodual_table,isovalue, isovert, isoquad_cube, facet_vertex,
     isodual_param, iso_vlist, isoquad_vert, cube_isovert_data, sharp_info);
}


// **************************************************
// Copy routine
// **************************************************

/// Copy isovert position from data structure isovert
void ISODUAL3D::copy_isovert_positions
(const std::vector<GRID_CUBE> & gcube_list,
 COORD_ARRAY & vertex_coord)
{
  VERTEX_INDEX num_isovert = gcube_list.size();

  vertex_coord.clear();
  vertex_coord.resize(DIM3*num_isovert);
  for (VERTEX_INDEX i = 0; i < num_isovert; i++) {
    std::copy(gcube_list[i].isovert_coord, 
              gcube_list[i].isovert_coord+DIM3,
              vertex_coord.begin() + i*DIM3);
  }
}

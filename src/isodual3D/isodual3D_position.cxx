/// \file isodual3D_position.cxx
/// Position dual isosurface vertices

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

#include "isodual3D_position.h"
#include "ijkinterpolate.txx"
#include "ijkcoord.txx"
#include "ijkgrid_macros.h"

#include "sharpiso_types.h"
#include "sharpiso_feature.h"

#include <iostream>
#include <iomanip>


using namespace ISODUAL3D;
using namespace SHARPISO;
using namespace std;

// **************************************************
// Position routines
// **************************************************

/// Position dual isosurface vertices in cube centers
void ISODUAL3D::position_dual_isovertices_cube_center
(const ISODUAL_GRID & grid,
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
(const ISODUAL_GRID & grid,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_cube_center(grid, vlist, &(coord.front()));
}

/// Position dual isosurface vertices in centroid
///   of isosurface-edge intersections
void ISODUAL3D::position_dual_isovertices_centroid
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
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
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)

{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_centroid
    (scalar_grid, isovalue, vlist, &(coord.front()));
}


/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{
  const int dimension = scalar_grid.Dimension();
  SVD_INFO svd_info;
  svd_info.location = LOC_NONE;

  EIGENVALUE_TYPE eigenvalues[DIM3];

  GRADIENT_COORD_TYPE max_small_mag(0.0);
  EIGENVALUE_TYPE max_small_eigenvalue(0.1);
  VERTEX_INDEX num_large_eigenvalues;


  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];

    svd_compute_sharp_vertex_in_cube
      (scalar_grid, gradient_grid, iv, isovalue,
      max_small_mag, max_small_eigenvalue,
      cube_offset2, sharp_coord+i*dimension,
      eigenvalues, num_large_eigenvalues, svd_info);
  }
}


/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
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
  }

}

/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_using_gradients
    (scalar_grid, gradient_grid, isovalue, cube_offset2, vlist, &(coord.front()));
}



/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();

  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_using_gradients
    (scalar_grid, gradient_grid, isovalue, isodual_param,
    vlist, &(coord.front()));
}


// **************************************************
// Position using interpolation along grid edge
// **************************************************

/// Position dual isosurface vertices using SVD and edge intersection simple
//  This function is called from isodual.cxx
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();
  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_using_edge_intersection_simple
    (scalar_grid, gradient_grid, isovalue, vlist, &(coord.front()));
};

/// Position dual isosurface vertices using SVD and edge intersection simple
//  same as above but takes the isodual_param paramter
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const ISODUAL_PARAM & isodual_param,
 std::vector<COORD_TYPE> & coord)
{
  const int dimension = scalar_grid.Dimension();
  coord.resize(vlist.size()*dimension);
  position_dual_isovertices_using_edge_intersection_simple
    (scalar_grid, gradient_grid, isovalue, vlist, isodual_param, 
     &(coord.front()));
};

/// Position dual isosurface vertices using SVD and edge intersection simple.
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{
  const int dimension = scalar_grid.Dimension();
  for (VERTEX_INDEX i = 0; i<vlist.size(); i++)
  {
    VERTEX_INDEX iv = vlist[i];

    //set up svd info
    SVD_INFO svd_info;
    svd_info.location = LOC_NONE;

    EIGENVALUE_TYPE eigenvalues[DIM3];

    GRADIENT_COORD_TYPE max_small_mag(0.0);
    EIGENVALUE_TYPE max_small_eigenvalue(0.1);
    VERTEX_INDEX cube_index(0);
    VERTEX_INDEX num_large_eigenvalues;

    svd_compute_sharp_vertex_in_cube_edge_based_simple
      (scalar_grid, gradient_grid, iv, isovalue,
      max_small_mag, max_small_eigenvalue, sharp_coord+i*dimension,
      eigenvalues, num_large_eigenvalues, svd_info);
  }
};


/// Position dual isosurface vertices using SVD and edge intersection simple.
// this one has the isovert param
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const ISODUAL_PARAM & isodual_param,
 COORD_TYPE * sharp_coord)
{
  const int dimension = scalar_grid.Dimension();
  for (VERTEX_INDEX i = 0; i<vlist.size(); i++)
  {
    VERTEX_INDEX iv = vlist[i];

    //set up svd info
    SVD_INFO svd_info;
    svd_info.location = LOC_NONE;

    EIGENVALUE_TYPE eigenvalues[DIM3];

    GRADIENT_COORD_TYPE max_small_mag(0.0);
    EIGENVALUE_TYPE max_small_eigenvalue(0.1);
    VERTEX_INDEX cube_index(0);
    VERTEX_INDEX num_large_eigenvalues;

    svd_compute_sharp_vertex_in_cube_edge_based_simple
      (scalar_grid, gradient_grid, iv, isovalue,
       max_small_mag, max_small_eigenvalue,sharp_coord+i*dimension,
       eigenvalues, num_large_eigenvalues, isodual_param,svd_info);

  }
}

// **************************************************
// Position using edge endpoint gradients
// **************************************************

/// Position dual isosurface vertices using SVD and edge intersection complex
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord){
   const int dimension = scalar_grid.Dimension();

   coord.resize(vlist.size()*dimension);
   position_dual_isovertices_using_edge_intersection_complex
     (scalar_grid, gradient_grid, isovalue, vlist, &(coord.front()));
};


/// Position dual isosurface vertices using SVD and edge intersection complex
//  also has the isodual_param 
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const ISODUAL_PARAM & isodual_param,
 std::vector<COORD_TYPE> & coord){
   const int dimension = scalar_grid.Dimension();

   coord.resize(vlist.size()*dimension);

   position_dual_isovertices_using_edge_intersection_complex
     (scalar_grid, gradient_grid, isovalue, vlist, isodual_param, 
      &(coord.front()));

};



/// Position dual isosurface vertices using SVD and edge intersection complex.
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{

  const int dimension = scalar_grid.Dimension();

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];
    //set up svd info
    SVD_INFO svd_info;
    svd_info.location = LOC_NONE;

    EIGENVALUE_TYPE eigenvalues[DIM3];

    GRADIENT_COORD_TYPE max_small_mag(0.0);
    EIGENVALUE_TYPE max_small_eigenvalue(0.1);
    VERTEX_INDEX cube_index(0);
    VERTEX_INDEX num_large_eigenvalues;

    svd_compute_sharp_vertex_in_cube_edge_based_cmplx
      (scalar_grid, gradient_grid, iv, isovalue, 
       max_small_mag, max_small_eigenvalue, sharp_coord+i*dimension, 
       eigenvalues, num_large_eigenvalues, svd_info);
  }
};

/// Position dual isosurface vertices using SVD and edge intersection complex.
/// this has isodual param 
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 const ISODUAL_PARAM & isodual_param,
 COORD_TYPE * sharp_coord)
{
  const int dimension = scalar_grid.Dimension();

  for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
    VERTEX_INDEX iv = vlist[i];
    //set up svd info
    SVD_INFO svd_info;
    svd_info.location = LOC_NONE;

    EIGENVALUE_TYPE eigenvalues[DIM3];

    GRADIENT_COORD_TYPE max_small_mag(0.0);
    EIGENVALUE_TYPE max_small_eigenvalue(0.1);
    VERTEX_INDEX cube_index(0);
    VERTEX_INDEX num_large_eigenvalues;

    svd_compute_sharp_vertex_in_cube_edge_based_cmplx
      (scalar_grid, gradient_grid, iv, isovalue,
       max_small_mag, max_small_eigenvalue, 
       sharp_coord+i*dimension, eigenvalues,
       num_large_eigenvalues, isodual_param, svd_info);

  }
}

// **************************************************
// Compute routines
// **************************************************

/// Compute centroid of intersections of isosurface and grid edges
void ISODUAL3D::compute_isosurface_grid_edge_centroid
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
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

// **************************************************
// Reposition close isosurface vertices
// **************************************************

namespace {

  typedef IJK::SCALAR_GRID_BASE<SHARPISO_GRID, VERTEX_INDEX>
  INDEX_GRID_BASE;
  typedef IJK::SCALAR_GRID<SHARPISO_GRID, VERTEX_INDEX>
  INDEX_GRID;
  typedef IJK::BOOL_GRID_BASE<SHARPISO_GRID>
  SHARPISO_BOOL_GRID_BASE;
  typedef IJK::BOOL_GRID<SHARPISO_GRID>
  SHARPISO_BOOL_GRID;

  void set_vertex_locations
  (const std::vector<ISO_VERTEX_INDEX> & vlist,
   INDEX_GRID_BASE & vloc, SHARPISO_BOOL_GRID_BASE & isovert_flag)
  {
    isovert_flag.SetAll(false);
    for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
      VERTEX_INDEX iv = vlist[i];
      vloc.Set(iv, i);
      isovert_flag.Set(iv, true);
    }
  }

  /// Returns true if \a s is greater than min scalar value of facet vertices
  ///   and \a s is less then or equal to max scalar value of facet vertices.
  /// @param scalar_grid Scalar grid.
  /// @pre GRID_TYPE must have member function FacetVertex(iv0,d)
  /// @pre scalar_grid.Dimension() > 0 so scalar_grid.NumCubeVertices() > 0.
  /// @param s Scalar value
  template <typename GRID_TYPE, typename ITYPE, typename DTYPE, 
            typename STYPE>
  bool is_gt_facet_min_le_facet_max
  (const GRID_TYPE & scalar_grid, const ITYPE iv0, const DTYPE orth_dir,
   const STYPE s)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    if (scalar_grid.Scalar(iv0) < s) {
      for (NTYPE k = 1; k < scalar_grid.NumFacetVertices(); k++) {
        VTYPE iv1 = scalar_grid.FacetVertex(iv0, orth_dir, k);
        if (scalar_grid.Scalar(iv1) >= s)
          { return(true); }
      }
    }
    else {
      for (NTYPE k = 1; k < scalar_grid.NumCubeVertices(); k++) {
        VTYPE iv1 = scalar_grid.FacetVertex(iv0, orth_dir, k);
        if (scalar_grid.Scalar(iv1) < s)
          { return(true); }
      }
    }

    return(false);
  }

  void get_cubes_around_interior_edge
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const VERTEX_INDEX iv0, const int dir,
   VERTEX_INDEX & jv0, VERTEX_INDEX & jv1,
   VERTEX_INDEX & jv2, VERTEX_INDEX & jv3)
  {
    int dir1 = (dir+1)%DIM3;
    int dir2 = (dir+2)%DIM3;
    jv2 = iv0;
    jv3 = scalar_grid.PrevVertex(iv0, dir2);
    jv1 = scalar_grid.PrevVertex(iv0, dir1);
    jv0 = scalar_grid.PrevVertex(jv1, dir2);
  }

  void reposition_close_isovertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const COORD_TYPE sep_dist_squared,
   const ISODUAL_PARAM & isodual_param,
   const VERTEX_INDEX iv0, const VERTEX_INDEX jloc0,
   const VERTEX_INDEX iv1, const VERTEX_INDEX jloc1,
   COORD_TYPE * isovert_coord)
  {
    COORD_TYPE distance_squared;
    IJK::compute_distance_squared
      (DIM3, isovert_coord+jloc0*DIM3, isovert_coord+jloc1*DIM3, 
       distance_squared);

    if (distance_squared < sep_dist_squared) {
      ISODUAL3D::compute_isosurface_grid_edge_centroid
        (scalar_grid, isovalue, iv0, isovert_coord+jloc0*DIM3);
    }
  }

}

/// Reposition to separate isosurface vertices
void ISODUAL3D::reposition_dual_isovertices
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const ISODUAL_PARAM & isodual_param,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * isovert_coord)
{
  INDEX_GRID vloc;
  SHARPISO_BOOL_GRID isovert_flag;
  GRID_COORD_TYPE v0coord[DIM3];

  COORD_TYPE separation_distance = isodual_param.separation_distance;
  COORD_TYPE sep_dist_squared = separation_distance*separation_distance;
  
  vloc.SetSize(scalar_grid);
  isovert_flag.SetSize(scalar_grid);

  set_vertex_locations(vlist, vloc, isovert_flag);

  for (VERTEX_INDEX i0 = 0; i0 < vlist.size(); i0++) {
    VERTEX_INDEX iv0 = vlist[i0];

    scalar_grid.ComputeCoord(iv0, v0coord);

    for (int d = 0; d < DIM3; d++) {

      if (v0coord[d]+1 < scalar_grid.AxisSize(d)) {
        VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, d);

        if (isovert_flag.Scalar(iv1)) {
          if (is_gt_facet_min_le_facet_max
              (scalar_grid, iv1, d, isovalue)) {
            VERTEX_INDEX i1 = vloc.Scalar(iv1);

            reposition_close_isovertices
              (scalar_grid, gradient_grid, isovalue, sep_dist_squared, 
               isodual_param, iv0, i0, iv1, i1, isovert_coord);
          }
        }
      }
    }
  }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iv0, dir, scalar_grid, VERTEX_INDEX) {
    VERTEX_INDEX iv1 = scalar_grid.NextVertex(iv0, dir);

    if (IJK::is_gt_min_le_max(scalar_grid, iv0, iv1, isovalue)) {
      VERTEX_INDEX jv0, jv1, jv2, jv3;
      get_cubes_around_interior_edge
        (scalar_grid, iv0, dir, jv0, jv1, jv2, jv3);
      VERTEX_INDEX j0 = vloc.Scalar(jv0);
      VERTEX_INDEX j1 = vloc.Scalar(jv1);
      VERTEX_INDEX j2 = vloc.Scalar(jv2);
      VERTEX_INDEX j3 = vloc.Scalar(jv3);

      reposition_close_isovertices
        (scalar_grid, gradient_grid, isovalue, sep_dist_squared, 
         isodual_param, jv0, j0, jv2, j2, isovert_coord);
      reposition_close_isovertices
        (scalar_grid, gradient_grid, isovalue, sep_dist_squared, 
         isodual_param, jv1, j1, jv3, j3, isovert_coord);
    }
  }
  
}


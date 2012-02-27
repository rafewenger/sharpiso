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
    svd_info.ray_intersect_cube = false;
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

/* OBSOLETE
/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const bool use_selected_gradients,
 const bool use_only_cube_gradients,
 const SIGNED_COORD_TYPE grad_selection_cube_offset,
 const SCALAR_TYPE ray_intersection_cube_offset,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{
  const int dimension = scalar_grid.Dimension();
  SVD_INFO svd_info;
  svd_info.ray_intersect_cube = false;
  svd_info.location = LOC_NONE;

  EIGENVALUE_TYPE eigenvalues[DIM3]={0.0};
    
  GRADIENT_COORD_TYPE max_small_mag(0.0);
  EIGENVALUE_TYPE max_small_eigenvalue(0.1);
  VERTEX_INDEX num_large_eigenvalues;
    
  if (use_selected_gradients) {
        
    OFFSET_CUBE_111 cube_111(grad_selection_cube_offset);
        
    if (use_only_cube_gradients) {
            
      for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
        VERTEX_INDEX iv = vlist[i];
                
        svd_compute_sharp_vertex_in_cube_S
          (scalar_grid, gradient_grid, iv, isovalue,
           max_small_mag, max_small_eigenvalue, ray_intersection_cube_offset,
           sharp_coord+i*dimension, 
           eigenvalues, num_large_eigenvalues, svd_info, cube_111);
      }
            
    }
    else {
            
      for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
        VERTEX_INDEX iv = vlist[i];
                
        svd_compute_sharp_vertex_neighborhood_S
          (scalar_grid, gradient_grid, iv, isovalue,
           max_small_mag, max_small_eigenvalue, ray_intersection_cube_offset,
           sharp_coord+i*dimension, 
           eigenvalues, num_large_eigenvalues, svd_info, cube_111);
      }
            
    }
  }
  else {
    if (use_only_cube_gradients) {
            
      for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
        VERTEX_INDEX iv = vlist[i];
                
        svd_compute_sharp_vertex_in_cube
          (scalar_grid, gradient_grid, iv, isovalue,
           max_small_mag, max_small_eigenvalue, ray_intersection_cube_offset,
           sharp_coord+i*dimension, 
           eigenvalues, num_large_eigenvalues, svd_info);
      }
    }
    else {
            
      for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
        VERTEX_INDEX iv = vlist[i];
                
        svd_compute_sharp_vertex_neighborhood
          (scalar_grid, gradient_grid, iv, isovalue,
           max_small_mag, max_small_eigenvalue, ray_intersection_cube_offset,
           sharp_coord+i*dimension, 
           eigenvalues, num_large_eigenvalues, svd_info);
      }
    }
  }
    
}
*/


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
  svd_info.ray_intersect_cube = false;
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


/// Position dual isosurface vertices using SVD and edge intersection simple.
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{
    
    const int dimension = scalar_grid.Dimension();
    for (VERTEX_INDEX i = 0; i<vlist.size(); i++)
    {
        VERTEX_INDEX iv = vlist[i];
       
        //set up svd info
        SVD_INFO svd_info;
        svd_info.ray_intersect_cube = false;
        svd_info.location = LOC_NONE;
        
        EIGENVALUE_TYPE eigenvalues[DIM3];
        
        GRADIENT_COORD_TYPE max_small_mag(0.0);
        EIGENVALUE_TYPE max_small_eigenvalue(0.1);
        VERTEX_INDEX cube_index(0);
        VERTEX_INDEX num_large_eigenvalues;
        
        svd_compute_sharp_vertex_in_cube_edge_based_simple
        (scalar_grid, gradient_grid, iv, isovalue,
         max_small_mag, max_small_eigenvalue,cube_offset2, sharp_coord+i*dimension,
         eigenvalues, num_large_eigenvalues, svd_info);
        
    }
};

/// Position dual isosurface vertices using SVD and edge intersection complex.
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 COORD_TYPE * sharp_coord)
{
    const int dimension = scalar_grid.Dimension();
    
    for (VERTEX_INDEX i = 0; i < vlist.size(); i++) {
        VERTEX_INDEX iv = vlist[i];
        //set up svd info
        SVD_INFO svd_info;
        svd_info.ray_intersect_cube = false;
        svd_info.location = LOC_NONE;
        
        EIGENVALUE_TYPE eigenvalues[DIM3];
        
        GRADIENT_COORD_TYPE max_small_mag(0.0);
        EIGENVALUE_TYPE max_small_eigenvalue(0.1);
        VERTEX_INDEX cube_index(0);
        VERTEX_INDEX num_large_eigenvalues;
        
        svd_compute_sharp_vertex_in_cube_edge_based_cmplx
        (scalar_grid, gradient_grid, iv, isovalue,
         max_small_mag, max_small_eigenvalue, cube_offset2, sharp_coord+i*dimension, eigenvalues,
         num_large_eigenvalues, svd_info);
    }
};



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

/* OBSOLETE
/// Position dual isosurface vertices using gradients
void ISODUAL3D::position_dual_isovertices_using_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const bool use_selected_gradients,
 const bool use_only_cube_gradients,
 const SIGNED_COORD_TYPE cube_offset,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
    const int dimension = scalar_grid.Dimension();
    
    coord.resize(vlist.size()*dimension);
    position_dual_isovertices_using_gradients
    (scalar_grid, gradient_grid, isovalue, 
     use_selected_gradients, use_only_cube_gradients, cube_offset, cube_offset2,
     vlist, &(coord.front()));
}
*/

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

/// Position dual isosurface vertices using SVD and edge intersection simple
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_simple
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord)
{
    
    const int dimension = scalar_grid.Dimension();
    coord.resize(vlist.size()*dimension);
    
    position_dual_isovertices_using_edge_intersection_simple
    (scalar_grid, gradient_grid, isovalue, cube_offset2, vlist, &(coord.front()));
};

/// Position dual isosurface vertices using SVD and edge intersection complex
void ISODUAL3D::position_dual_isovertices_using_edge_intersection_complex
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE cube_offset2,
 const std::vector<ISO_VERTEX_INDEX> & vlist,
 std::vector<COORD_TYPE> & coord){
    const int dimension = scalar_grid.Dimension();
    
    coord.resize(vlist.size()*dimension);
    position_dual_isovertices_using_edge_intersection_complex
    (scalar_grid, gradient_grid, isovalue, cube_offset2, vlist, &(coord.front()));
};

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

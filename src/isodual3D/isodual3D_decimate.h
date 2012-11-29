/// \file isodual3D_decimate.h
/// Decimate dual isosurface.

/*
IJK: Isosurface Jeneration Kode
Copyright (C) 2012 Rephael Wenger

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


#ifndef _ISODUAL3D_DECIMATE_
#define _ISODUAL3D_DECIMATE_

#include <vector>

#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"


namespace ISODUAL3D {


  // **************************************************
  // Merge some isosurface vertices
  // **************************************************

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & quad_vert,
   std::vector<VERTEX_INDEX> & gcube_map, SHARPISO_INFO & sharpiso_info);

  /// Merge isosurface vertices in cubes adjacent to selected sharp cubes.
  void merge_sharp_iso_vertices
  (const ISODUAL_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const ISOVERT & isovert, 
   const SHARP_ISOVERT_PARAM & sharp_isovert_param,
   std::vector<VERTEX_INDEX> & quad_vert,
   SHARPISO_INFO & sharpiso_info);


  // **************************************************
  // Function class: IS_ISOPATCH_DISK
  // **************************************************

  /// Function class for determining if isopatch is a disk.
  class IS_ISOPATCH_DISK {

  public:
    static const AXIS_SIZE_TYPE num_vert_along_region_axis = 4;
    static const AXIS_SIZE_TYPE region_edge_length = 
      num_vert_along_region_axis-1;

  protected:
    AXIS_SIZE_TYPE region_axis_size[DIM3];

    SHARPISO_SCALAR_GRID regionScalar;
    SHARPISO_BOOL_GRID cubeFlag;

    /// Indicates vertices on region boundary.
    SHARPISO_BOOL_GRID regionBoundary;

    /// Indicates vertices on boundary of region formed by selected cubes.
    SHARPISO_BOOL_GRID selectedCubeBoundary;

    SHARPISO_GRID_NEIGHBORS neighbor_grid;

    /// List of boundary cubes in region.
    std::vector<VERTEX_INDEX> region_boundary_cube;

    /// Increments for region vertices.
    /// Vertex k in region corresponds to vertex iv+region_vertex_increment[k]
    ///   around vertex iv.
    INDEX_DIFF_TYPE * region_vertex_increment;

    bool * visited;

    void SetCubeFlag
    (const VERTEX_INDEX cube_index,
     const ISOVERT & isovert,
     const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map);
    void SetScalar
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
     const VERTEX_INDEX cube_index);
    void SetSelectedCubeBoundary();

    /// Set region_vertex_increment.
    void SetRegionVertexIncrement(const SHARPISO_GRID & grid);

    void GetBoundaryVertices
    (const SCALAR_TYPE isovalue, const bool flag_pos,
     std::vector<int> & vlist) const;
    void GetBoundaryPosVertices
    (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const;
    void GetBoundaryNegVertices
    (const SCALAR_TYPE isovalue, std::vector<int> & vlist) const;
    void GetBoundaryEdges
    (const SCALAR_TYPE isovalue, const bool flag_pos,
     std::vector<int> & elist) const;
    void GetBoundaryPosEdges
    (const SCALAR_TYPE isovalue, std::vector<int> & elist) const;
    void GetBoundaryNegEdges
    (const SCALAR_TYPE isovalue, std::vector<int> & elist) const;

  public:
    IS_ISOPATCH_DISK(const SHARPISO_GRID & grid);
    ~IS_ISOPATCH_DISK();

    /// Return true is isopatch is a disk.
    bool IsIsopatchDisk
    (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
     SCALAR_TYPE isovalue,
     const VERTEX_INDEX cube_index,
     const ISOVERT & isovert,
     const std::vector<SHARPISO::VERTEX_INDEX> & gcube_map);

    /// Reverse merges to isosurface vertex at cube_index.
    void UnmapAdjacent
    (const NUM_TYPE cube_index, const ISODUAL3D::ISOVERT & isovert, 
     std::vector<SHARPISO::VERTEX_INDEX> & gcube_map) const;
  };

};


#endif

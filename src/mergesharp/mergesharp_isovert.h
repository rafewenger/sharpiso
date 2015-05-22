/// \file mergesharp_isovert.h
/// Data structures for creating and processing sharp isosurface vertices.

/*
Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _MERGESHARP_ISOVERT_H_
#define _MERGESHARP_ISOVERT_H_

#include "ijkdualtable.h"

#include "sharpiso_grids.h"
#include "sharpiso_feature.h"
#include "mergesharp_types.h"

#include <vector>

namespace MERGESHARP {

// **************************************************
// GRID CUBE DATA
// **************************************************

typedef enum {
	AVAILABLE_GCUBE,       ///< Cube is available for selection.
	SELECTED_GCUBE,        ///< Cube contains a sharp vertex.
	COVERED_A_GCUBE,       ///< Cube is near a cube containing a sharp vertex.
	COVERED_B_GCUBE,	     ///< Cube is covered by extended mapping.
	COVERED_CORNER_GCUBE,  ///< Covers a cube with >3 eigen value
	COVERED_POINT,         ///< The sharp vertex is in a covered cube.
	UNAVAILABLE_GCUBE,     ///< Cube is within 3x3 of a 2 covering.
	NON_DISK_GCUBE,        ///< Merging cube creates non-disk patch.
	SMOOTH_GCUBE		       ///< Cube contains smooth isosurface patch.
} GRID_CUBE_FLAG;


class GRID_CUBE_DATA {

protected:

  void Init();

public:
  GRID_COORD_TYPE cube_coord[DIM3];  ///< Cube coordinates (unscaled).
  COORD_TYPE isovert_coord[DIM3];    ///< Location of the sharp isovertex.
  COORD_TYPE isovert_coordB[DIM3];   ///< Substitute location.
  unsigned char num_eigenvalues;     ///< Number of eigenvalues.
  GRID_CUBE_FLAG flag;               ///< Type for this cube.
  BOUNDARY_BITS_TYPE boundary_bits;  ///< Boundary bits for the cube
  VERTEX_INDEX cube_index;           ///< Index of cube in scalar grid.

  /// Type of cube cover.
  CUBE_ADJACENCY_TYPE cover_type;

  /// If num_eigenvalues == 2, then direction = direction of isosurface edge.
  /// If num_eigenvalues == 1, then 
  ///   direction = direction orthogonal to isosurface.
  COORD_TYPE direction[DIM3];         

  /// Linf-dist from isovert_coord[] to cube-center.
  COORD_TYPE linf_dist;

  /// If true, location is centroid of (grid edge)-isosurface intersections.
  bool flag_centroid_location;

  /// If true, some other non-empty cube contains the isovert coord.
  bool flag_conflict;

  /// If true, cube is near corner cube.
  bool flag_near_corner;

  /// Grid index of cube containing the isovert_coord.
  VERTEX_INDEX cube_containing_isovert;

  /// If true, isovert_coord[] determined by an adjacent cube.
  bool flag_coord_from_other_cube;

  /// If true, isovert_coord[] determined by a grid vertex.
  bool flag_coord_from_vertex;

  /// If true, isovert_coord[] determined by a grid edge.
  bool flag_coord_from_edge;

  /// If true, using replacement coordinate.
  bool flag_using_substitute_coord;

  /// If true, coordinates have been recomputed.
  bool flag_recomputed_coord;

  /// If true, coordinates have been recomputed
  ///   with min gradient cube offset.
  bool flag_recomputed_coord_min_offset;

  /// If true, location is recomputed from location of adjacent vertices.
  bool flag_recomputed_using_adjacent;

  /// If true, svd coord were farther than max_dist.
  bool flag_far;

  /// If true, cube is selected despite mismatch.
  bool flag_ignore_mismatch;

  /// Index of cube configuration is isosurface lookup table.
  IJKDUALTABLE::TABLE_INDEX table_index;

  /// Grid index of cube which covered this cube.
  VERTEX_INDEX covered_by;

  /// Grid index of cube which this cube maps to.
  /// Currently, only used for output information.
  VERTEX_INDEX maps_to_cube;

  /// Return true if cube is covered or selected.
  bool IsCoveredOrSelected() const;

  GRID_CUBE_DATA() { Init(); }
};

typedef std::vector<GRID_CUBE_DATA> GRID_CUBE_DATA_ARRAY;


// **************************************************
// GCUBE_COMPARE
// **************************************************

class GCUBE_COMPARE {

public:
  const std::vector<GRID_CUBE_DATA> * gcube_list;

  GCUBE_COMPARE(const std::vector<GRID_CUBE_DATA> & gcube_list)
  { this->gcube_list = &gcube_list; };

  bool operator () (int i,int j)
  {
    int num_eigen_i = gcube_list->at(i).num_eigenvalues;
    int num_eigen_j = gcube_list->at(j).num_eigenvalues;

    if (num_eigen_i == num_eigen_j) {

      /* OBSOLETE
      int flag_i = int(gcube_list->at(i).flag_coord_from_other_cube);
      int flag_j = int(gcube_list->at(j).flag_coord_from_other_cube);
      */

      return ((gcube_list->at(i).linf_dist) < (gcube_list->at(j).linf_dist)); 

      /* OBSOLETE
      if (flag_i == flag_j) {
        return ((gcube_list->at(i).linf_dist) < (gcube_list->at(j).linf_dist)); 
      }
      else {
        // Process first cubes which generated their own iso vertices.
        return ((flag_i < flag_j));
      }
      */
    }
    else {
      return ((num_eigen_i > num_eigen_j));
    }
  }
};


// **************************************************
// ISOSURFACE VERTEX DATA
// **************************************************

class ISOVERT {
public:

	/// gcube_list containing the active cubes and their vertices.
	std::vector<GRID_CUBE_DATA> gcube_list;

	static const int NO_INDEX = -1;       ///< Flag for no index.

	/// Grid containing the index to the gcube_list.
	/// If cube is not active, then it is defined as NO_INDEX.
	SHARPISO_INDEX_GRID index_grid;

  /// Grid neighbor information.
  SHARPISO_GRID_NEIGHBORS grid;

  /// Return true if cube is active.
	bool isActive(const int cube_index) const;

  /// Return true if cube flag equals flag.
	bool isFlag(const int cube_index, GRID_CUBE_FLAG flag) const; 

  /// Return cube index.
  VERTEX_INDEX CubeIndex(const int gcube_index) const
  { return(gcube_list[gcube_index].cube_index); }

  /// Return gcube index or NO_INDEX.
  INDEX_DIFF_TYPE GCubeIndex(const int cube_index) const
  { return(index_grid.Scalar(cube_index)); }

  /// Return gcube index or NO_INDEX.
  /// Set error message if cube not active (gcube_index = NO_INDEX).
  INDEX_DIFF_TYPE GCubeIndex(const int cube_index, IJK::ERROR & error) const;

  /// Return number of eigenvalues.
  NUM_TYPE NumEigenvalues(const int gcube_index) const
  { return(gcube_list[gcube_index].num_eigenvalues); }

  /// Return pointer to isovert_coord[].
  const COORD_TYPE * IsoVertCoord(const int gcube_index) const
  { return(gcube_list[gcube_index].isovert_coord); }
};


// **************************************************
// ISOVERT INFO
// **************************************************

class ISOVERT_INFO {

public:

	int num_sharp_corners;
	int num_sharp_edges;
	int num_smooth_vertices;
  int num_merged_iso_vertices;
  int num_conflicts;
  int num_Linf_iso_vertex_locations;

  void Clear();                          ///< Clear all values.

  /// Constructor.
  ISOVERT_INFO()
  { Clear(); };

  /// Set.
  void Set(const ISOVERT_INFO & info)
  { *this = info; }
};


// **************************************************
// ROUTINES
// **************************************************

/// Compute dual isosurface vertices.
void compute_dual_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   const VERTEX_POSITION_METHOD vertex_position_method,
   ISOVERT & isovert);

/// Compute dual isosurface vertices.
void compute_dual_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const std::vector<COORD_TYPE> & edgeI_coord,
   const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT & isovert);

/// Recompute isosurface vertex positions for cubes 
///   which are not selected or covered.
/// also takes isovert_info as parameter
void recompute_isovert_positions(
	const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID_BASE & gradient_grid,
	const SCALAR_TYPE isovalue,
	const SHARP_ISOVERT_PARAM & isovert_param,
	ISOVERT & isovertData);


/// Recompute isosurface vertex positions for cubes 
///   which are not selected or covered.
void recompute_isovert_positions(
    const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const SHARP_ISOVERT_PARAM & isovert_param,
    ISOVERT & isovertData);


/// Recompute isosurface vertex positions for cubes 
///   which are not selected or covered.
/// Version for hermite data.
void recompute_isovert_positions
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<COORD_TYPE> & edgeI_coord,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert);

/// Recompute isovert positions for cubes containing covered points.
void recompute_covered_point_positions
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 ISOVERT & isovert);

/// Recompute isovert positions for cubes containing covered points.
/// Use voxel for gradient cube offset, 
///   not isovert_param.grad_selection_cube_offset.
/// @param flag_min_offset If true, voxel uses minimum gradient cube offset.
void recompute_covered_point_positions
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const GRADIENT_GRID_BASE & gradient_grid,
 const SHARPISO_BOOL_GRID & covered_grid,
 const SCALAR_TYPE isovalue,
 const SHARP_ISOVERT_PARAM & isovert_param,
 const OFFSET_VOXEL & voxel,
 const bool flag_min_offset,
 ISOVERT & isovert);

/// Recompute using adjacent isosurface vertex locations.
void recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, ISOVERT & isovert);

/// Recompute isovert in cube cube_index
///   using adjacent isosurface vertex locations.
void recompute_using_adjacent
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, 
 const VERTEX_INDEX cube_index, ISOVERT & isovert,
 bool & flag_loc_recomputed);

/// Set cover type of all covered cubes.
void set_cover_type(ISOVERT & isovert);

/// Set grid containing locations of edges in edgeI_coord[].
void set_edge_index(const std::vector<COORD_TYPE> & edgeI_coord,
                    SHARPISO_EDGE_INDEX_GRID & edge_index);

/// Count number of vertices on sharp corners or sharp edges.
/// Count number of smooth vertices.
void count_vertices
(const ISOVERT & isovert, ISOVERT_INFO & isovert_info);

/// Select all grid cubes which are not smooth.
void select_non_smooth(ISOVERT & isovert);

/// Get list of grid cubes from isovert.
void get_cube_list
  (const ISOVERT & isovert, std::vector<VERTEX_INDEX> & cube_list);

/// Transform GRID_CUBE_FLAG into a string
void convert2string(const GRID_CUBE_FLAG & flag, std::string & s);

// Transform CUBE_ADJACENCY_TYPE into a string
void convert2string
(const CUBE_ADJACENCY_TYPE & adj_type, std::string & s);

/// Copy isovert position from from_gcube to to_gcube.
void copy_isovert_position
(const SHARPISO_GRID & grid, 
 const NUM_TYPE from_gcube, const NUM_TYPE to_gcube, ISOVERT & isovert);

/// Compute overlap of 3x3x3 regions
bool find_3x3x3_overlap
(const GRID_COORD_TYPE cubeA_coord[DIM3],
 const GRID_COORD_TYPE cubeB_coord[DIM3],
 GRID_COORD_TYPE rmin[DIM3], GRID_COORD_TYPE rmax[DIM3],
 int & overlap_dim);

/// Compute overlap of 3x3x3 regions.
/// Version with grid and cube indices input.
template <typename GRID_TYPE>
bool find_3x3x3_overlap
(const GRID_TYPE & grid,
 const VERTEX_INDEX cubeA_index, const VERTEX_INDEX cubeB_index,
 GRID_COORD_TYPE rmin[DIM3], GRID_COORD_TYPE rmax[DIM3],
 int & overlap_dim);

// **************************************************
// APPLY TEMPLATES
// **************************************************

/// Return (b & (1 << i))
template <typename BTYPE, typename NTYPE>
inline BTYPE and_bit(const BTYPE b, const NTYPE i)
{
  return(b & (BTYPE(1) << i));
}

/// Convert cube corner index to vertex bit mask
template <typename DTYPE, typename ITYPE, typename BTYPE>
inline void convert2vertex_bit_mask
(const DTYPE dimension, const ITYPE corner_index, BTYPE & bit_mask)
{
  bit_mask = 0;
  for (DTYPE d = 0; d < dimension; d++) {
    if (and_bit(corner_index, d) == 0)
      { bit_mask = (bit_mask | (BTYPE(1) << (2*d))); }
    else
      { bit_mask = (bit_mask | (BTYPE(1) << (2*d+1))); }
  }
}

/// Convert cube corner index to edge bit mask
template <typename DTYPE, typename DTYPE2, typename ITYPE, typename BTYPE>
inline void convert2edge_bit_mask
(const DTYPE dimension, const DTYPE2 edge_dir,
 const ITYPE corner_index, BTYPE & bit_mask)
{
  bit_mask = 0;
  for (DTYPE d = 0; d < dimension; d++) {
    if (d != edge_dir) {
      if (and_bit(corner_index, d) == 0)
        { bit_mask = (bit_mask | (BTYPE(1) << (2*d))); }
      else
        { bit_mask = (bit_mask | (BTYPE(1) << (2*d+1))); }
    }
  }
}

/// Compute cube which shares i'th corner (and no other corners)
///   with icubeA.
template <typename GTYPE, typename VTYPE, typename ITYPE>
inline VTYPE compute_vertex_adjacent_cube
(const GTYPE & grid, const VTYPE icubeA, const ITYPE i)
{
  typedef typename GTYPE::DIMENSION_TYPE DTYPE;

  VTYPE icubeB = icubeA;
  for (DTYPE d = 0; d < grid.Dimension(); d++) {

    if (and_bit(i, d) == 0)
      { icubeB = grid.PrevVertex(d, icubeB); }
    else
      { icubeB = grid.NextVertex(d, icubeB); }
  }

  return(icubeB);
}

/// Compute cube which shares i'th edge in direction edge_dir 
///   (and no other edges) with icubeA.
template <typename GTYPE, typename VTYPE, typename DTYPE2, typename ITYPE>
inline VTYPE compute_edge_adjacent_cube
(const GTYPE & grid, const VTYPE icubeA, const DTYPE2 edge_dir,
 const ITYPE i)
{
  typedef typename GTYPE::DIMENSION_TYPE DTYPE;

  VTYPE icubeB = icubeA;
  for (DTYPE d = 0; d < grid.Dimension(); d++) {

    if (d != edge_dir) {
      if (and_bit(i, d) == 0)
        { icubeB = grid.PrevVertex(d, icubeB); }
      else
        { icubeB = grid.NextVertex(d, icubeB); }
    }
  }

  return(icubeB);
}

/// Apply to each cube which is facet adjacent to icubeA. (4+8 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8>
void apply_to_cubes_facet_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }
  else {

    for (DTYPE d = 0; d < grid.Dimension(); d++) {

      if (and_bit(boundary_bits, 2*d) == 0) {
        const VTYPE icubeB = grid.PrevVertex(d, icubeA);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }

      if (and_bit(boundary_bits, 2*d+1) == 0) {
        const VTYPE icubeB = grid.NextVertex(d, icubeA);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }
    }
  }
}

/// Apply to each cube which is facet adjacent to icubeA. (4+9 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8, typename ATYPE9>
void apply_to_cubes_facet_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
      const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }
  else {

    for (DTYPE d = 0; d < grid.Dimension(); d++) {

      if (and_bit(boundary_bits, 2*d) == 0) {
        const VTYPE icubeB = grid.PrevVertex(d, icubeA);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
      }

      if (and_bit(boundary_bits, 2*d+1) == 0) {
        const VTYPE icubeB = grid.NextVertex(d, icubeA);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
      }
    }
  }
}

/// Apply to each cube which is edge adjacent to icubeA. (4+8 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8>
void apply_to_cubes_edge_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
      const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }
  else {

    for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

      BTYPE mask_d = (BTYPE(1) << edge_dir);
      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) == 1) { continue; }

        BTYPE mask;
        convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

        if ((j & mask) == 0) {
          VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
      }
    }
  }
}

/// Apply to each cube which is edge adjacent to icubeA. (4+9 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8, typename ATYPE9>
void apply_to_cubes_edge_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
      const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }
  else {

    for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

      BTYPE mask_d = (BTYPE(1) << edge_dir);
      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) == 1) { continue; }

        BTYPE mask;
        convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

        if ((j & mask) == 0) {
          VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, 
            arg8, arg9);
        }
      }
    }
  }
}

/// Apply to each cube which is vertex adjacent to icubeA. (4+8 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8>
void apply_to_cubes_vertex_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
      const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }
  else {

    for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

      BTYPE mask;
      convert2vertex_bit_mask(grid.Dimension(), j, mask);

      if ((j & mask) == 0) {
        VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }
    }
  }
}

/// Apply to each cube which is vertex adjacent to icubeA. (4+9 arguments.)
template <typename FTYPE, typename GRID_TYPE, 
          typename VTYPE, typename BTYPE,
          typename ATYPE1, typename ATYPE2, typename ATYPE3,
          typename ATYPE4, typename ATYPE5, typename ATYPE6, 
          typename ATYPE7, typename ATYPE8, typename ATYPE9>
void apply_to_cubes_vertex_adjacent_to
(FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
{
  typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

  if (boundary_bits == 0) {
    // Cube icubeA is an interior cube.

    for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
      const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
      f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }
  else {

    for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

      BTYPE mask;
      convert2vertex_bit_mask(grid.Dimension(), j, mask);

      if ((j & mask) == 0) {
        VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
      }
    }
  }
}

/// Apply to each cube which is facet adjacent to a cube in gcube_list.
/// (3 + 7 arguments.)
template <typename FTYPE, typename ATYPE1, typename ATYPE2,
          typename ATYPE3,typename ATYPE4, typename ATYPE5,
          typename ATYPE6, typename ATYPE7>
void apply_to_cubes_facet_adjacent_to_list
(FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7)
{
  for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

    NUM_TYPE to_gcube = gcube_list[i];
    VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[to_gcube].boundary_bits;

    apply_to_cubes_facet_adjacent_to
      (f, isovert.grid, to_cube, boundary_bits, isovert,
       arg1, arg2, arg3, arg4, arg5, arg6, arg7);
  }
}

/// Apply to each cube which is edge adjacent to a in gcube_list.
template <typename FTYPE, typename ATYPE1, typename ATYPE2,
          typename ATYPE3,typename ATYPE4, typename ATYPE5,
          typename ATYPE6, typename ATYPE7>
void apply_to_cubes_edge_adjacent_to_list
(FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7)
{
  for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

    NUM_TYPE to_gcube = gcube_list[i];
    VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[to_gcube].boundary_bits;

    apply_to_cubes_edge_adjacent_to
      (f, isovert.grid, to_cube, boundary_bits, isovert,
       arg1, arg2, arg3, arg4, arg5, arg6, arg7);
  }
}

/// Apply to each cube which is vertex adjacent to a in gcube_list.
template <typename FTYPE, typename ATYPE1, typename ATYPE2,
          typename ATYPE3,typename ATYPE4, typename ATYPE5,
          typename ATYPE6, typename ATYPE7>
void apply_to_cubes_vertex_adjacent_to_list
(FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
 ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
 ATYPE6 & arg6, ATYPE7 & arg7)
{
  for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

    NUM_TYPE to_gcube = gcube_list[i];
    VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

    BOUNDARY_BITS_TYPE boundary_bits = 
      isovert.gcube_list[to_gcube].boundary_bits;

    apply_to_cubes_vertex_adjacent_to
      (f, isovert.grid, to_cube, boundary_bits, isovert,
       arg1, arg2, arg3, arg4, arg5, arg6, arg7);
  }
}


// **************************************************
// BIN_GRID ROUTINES
// **************************************************

/// Initialize bin_grid.
/// @param bin_width = number of cubes along each axis.
void init_bin_grid
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 BIN_GRID<VERTEX_INDEX> & bin_grid);

/// Insert cube cube_index into the bin_grid.
void bin_grid_insert
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid);

/// Remove cube cube_index into the bin_grid.
void bin_grid_remove
(const SHARPISO_GRID & grid, const AXIS_SIZE_TYPE bin_width,
 const VERTEX_INDEX cube_index, BIN_GRID<int> & bin_grid);

/// Get the selected vertices around iv.
void get_selected
(const SHARPISO_GRID & grid,
 const VERTEX_INDEX iv,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 std::vector<VERTEX_INDEX> & selected_list);


// **************************************************
// SUBROUTINES
// **************************************************

/// Return true if this vertex creates a triangle with a large angle.
/// @param check_triangle_angle If true, check it triangle has large angles.
/// @param bin_grid Contains the already selected vertices.
/// @param[out] v1,v2 vertex indices which form a triangle with iv.
bool creates_triangle
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const bool check_triangle_angle,
 const ISOVERT & isovertData,
 const VERTEX_INDEX iv,
 const SCALAR_TYPE isovalue,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 VERTEX_INDEX & v1,
 VERTEX_INDEX & v2);

/// Return true if selecting this vertex creates a triangle with a large angle.
/// @param bin_grid Contains the already selected vertices.
/// @param[out] iv1, iv2 vertex indices which form a triangle with iv0
bool creates_triangle_new
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
 const ISOVERT & isovert,
 const VERTEX_INDEX iv0,
 const SCALAR_TYPE isovalue,
 const BIN_GRID<VERTEX_INDEX> & bin_grid,
 const AXIS_SIZE_TYPE bin_width,
 VERTEX_INDEX & iv1,
 VERTEX_INDEX & iv2);

/// Select and sort cubes with more than one eigenvalue.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to sharp cubes
///    sorted by number of large eigenvalues and by distance 
///    of sharp coord from cube center.


/// @param[out] gcube_index_list Sorted list of references to grid cubes
///      with more than one eigenvalue.
///    List is sorted by decreasing number of eigenvalues and then by
///      increasing linf_dist.
void get_corner_or_edge_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list);

/// Get selected cubes.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to selected cubes
///    sorted by number of large eigenvalues and then by distance 
///    of sharp coord from cube center.
void get_selected_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list);

/// Get selected corner cubes.
///   Store references to cubes sorted by number of large eigenvalues
///     and then by increasing distance from isovert_coord to cube center.
/// @param gcube_index_list List of references to selected corner cubes
///    sorted by distance of sharp coord from cube center.
void get_selected_corner_cubes
(const std::vector<GRID_CUBE_DATA> & gcube_list,
 std::vector<NUM_TYPE> & gcube_index_list);

/// Store boundary bits for each cube in gcube_list.
void store_boundary_bits
(const SHARPISO_GRID & grid, GRID_CUBE_DATA_ARRAY & gcube_list);

/// Store isosurface lookup table index in gcube_list.
void store_table_index
(const std::vector<IJKDUALTABLE::TABLE_INDEX> & table_index,
 GRID_CUBE_DATA_ARRAY & gcube_list);

/// Return true if line through sharp edge in cube0 passes near cube1.
/// Return false if cube0 has no sharp edge.
bool does_sharp_edge_point_to_cube
(const SHARPISO_GRID & grid, const ISOVERT & isovert,
 const VERTEX_INDEX cube0_index, const VERTEX_INDEX cube1_index,
 const COORD_TYPE min_distance);

}

#endif /* _MERGESHARP_ISOVERT_H_ */

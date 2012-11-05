/// \file isodual3D_isovert.h
/// Data structures for creating and processing sharp isosurface vertices.
/// Created on: Oct 29, 2012
/// Author: arindam


#ifndef _ISODUAL3D_ISOVERT_H_
#define _ISODUAL3D_ISOVERT_H_

#include "isodual3D_types.h"
#include "sharpiso_grids.h"
#include "sharpiso_feature.h"

#include <vector>

namespace ISODUAL3D {

// **************************************************
// GRID CUBES INFORMATION
// **************************************************

typedef enum{
	AVAILABLE_GCUBE,   // available
	SELECTED_GCUBE,    // cube contains a sharp vertex
	COVERED_GCUBE,     // cube is within 3x3 of a cube containing a sharp vertex
	UNAVAILABLE_GCUBE,  // cube is within 3x3 of a 2 covering
	SMOOTH_GCUBE		// smooth
} GRID_CUBE_FLAG;


class GRID_CUBE{
public:
	COORD_TYPE isovert_coord[DIM3]; ///< Location of the sharp isovertex.
  unsigned char num_eigen;        ///< Number of eigenvalues.
	GRID_CUBE_FLAG flag;            ///< Type for this cube.
	SCALAR_TYPE linf_dist;          /// Linf-dist from sharp point to cube-center.
	int boundary_bits;             /// boundary bits for the cube
	VERTEX_INDEX cube_index;          /// Index of cube in scalar grid.
};


// **************************************************
// ISOSURFACE VERTEX DATA
// **************************************************

class ISOVERT {
public:

  /// gcube_list containing the active cubes and their vertices.
	std::vector<GRID_CUBE> gcube_list;

  static const int NO_INDEX = -1;       ///< Flag for no index.

  /// Grid containing the index to the gcube_list.
  /// If cube is not active, then it is defined as NO_INDEX.
	SHARPISO_INDEX_GRID sharp_ind_grid;

	SCALAR_TYPE linf_dist_to_sharp_pt;    ///< L-inf distance variable.
};


// **************************************************
// ISOVERT INFO
// **************************************************

class ISOVERT_INFO {

 public:

  int num_sharp_corners;
  int num_sharp_edges;
  int num_smooth_vertices;
  int num_vertex_collapses;
};


// **************************************************
// ROUTINES
// **************************************************

  // Compute the ISOVERT object.
  void compute_dual_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT &isovertData);

  /// Count number of vertices on sharp corners or sharp edges.
  /// Count number of smooth vertices.
  void count_vertices
  (const ISOVERT & isovert, ISOVERT_INFO & isovert_info);

}

#endif /* _ISODUAL3D_ISOVERT_H_ */

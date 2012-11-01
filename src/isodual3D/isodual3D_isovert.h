/*
 * isodual3D_isovert.h
 *
 *  Created on: Oct 29, 2012
 *      Author: arindam
 */

#ifndef ISODUAL3D_ISOVERT_H_
#define ISODUAL3D_ISOVERT_H_

#include "ijkgrid.txx"
#include "ijkscalar_grid.txx"
#include "sharpiso_types.h"
#include "sharpiso_grids.h"
#include "isodual3D_types.h"
#include "isodual3D_datastruct.h"



#include <vector>
namespace SHARPISO{

typedef enum{
	AVAILABLE_GCUBE,   // available
	SELECTED_GCUBE,    // cube contains a sharp vertex
	COVERED_GCUBE,     // cube is within 3x3 of a cube containing a sharp vertex
	UNAVAILABLE_GCUBE,  // cube is within 3x3 of a 2 covering
	SMOOTH_GCUBE		// smooth
} GRID_CUBE_FLAG;


class GRID_CUBE{
public:
	COORD_TYPE isovert_coord[DIM3];//location of the sharp isovertex
	unsigned char num_eigen;
	GRID_CUBE_FLAG flag; // defines the type for this cube
	SCALAR_TYPE l2dist;  // this is the l2 dist between the sharp point and the cube-center

};



// grid contains the index to the gcube_list defined in ISOVERT
// If tube is not sharp then it has [-1]
typedef IJK::SCALAR_GRID<SHARPISO_GRID, NUM_TYPE>SHARPISO_INDEX_GRID;

class ISOVERT{
public:
	std::vector<GRID_CUBE> gcube_list; // vector containing the cubes which has sharp vertex
	SHARPISO_INDEX_GRID sharp_ind_grid;
	SCALAR_TYPE linf_dist_to_sharp_pt;/// linf distance variable
};




}

namespace ISODUAL3D{

// compute the ISOVERT object
void computeDualIsovert(
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID_BASE & gradient_grid,
		const SCALAR_TYPE isovalue,
		const ISODUAL_PARAM & isodual_param,
		ISOVERT &isovertData
);
}

#endif /* ISODUAL3D_ISOVERT_H_ */

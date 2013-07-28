
#ifndef _FINDEDGECOUNT_
#define _FINDEDGECOUNT_
#include "countdegree_types.h"

void  count_edge_degrees
	(const int dim, const COORD_TYPE * coord, const int numv,
	const VERTEX_INDEX * edge_vert, const int nume, vector <int> &vert_degree);

#endif

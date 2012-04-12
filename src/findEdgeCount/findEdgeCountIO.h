// print statements and other help routines for findEdgeCount 
#ifndef _FINDEDGECOUNT_IO_
#define _FINDEDGECOUNT_IO_
#include "findEdgeCountTypes.h"
// **************************************************
// Output routines
// **************************************************
void output_edges
(const int dim, const COORD_TYPE * coord, const int numv,
 const VERTEX_INDEX * edge_vert, const int nume);
void output_vert_degree
( const int numv, vector <int> &vert_degree);
void output_vert_degree_2_file 
(const int numv, vector <int> &vert_degree,const string output_fn);

#endif


// *************************************
// Compute the degrees of each edges
// **************************************
#include <iostream>
#include "findEdgeCount.h"
using namespace std;
void  count_edge_degrees
	(const int dim, const COORD_TYPE * coord, const int numv,
	const VERTEX_INDEX * edge_vert, const int nume, vector <int> &vert_degree)
{
	for (int i=0;i<nume;i++)
	{
		VERTEX_INDEX iv0 = edge_vert[2*i];
		VERTEX_INDEX iv1 = edge_vert[2*i+1];
		vert_degree[iv0]++;
		vert_degree[iv1]++;
	}
}

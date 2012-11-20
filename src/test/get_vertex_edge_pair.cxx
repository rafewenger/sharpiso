
#include "sharpiso_types.h"
#include "sharpiso_grids.h"

using namespace SHARPISO;

// Get neighboring vertex and facet edge.
void get_vertex_edge_pair
(const SHARPISO_GRID & grid, const VERTEX_INDEX iv0, const int facet_index,
 const int facet_vertex_index, const int edge_index,
 VERTEX_INDEX & iw0, VERTEX_INDEX & iw1, VERTEX_INDEX & iw2,
 int & edir01, int & eorient01, int & edir12, int & eorient12)
{
  edir01 = facet_index % DIM3;
  edir12 = (edir01 + 1 + edge_index)%DIM3;

  iw1 = grid.FacetVertex(iv0, edir01, facet_vertex_index);

  if (facet_index < 3) {
    eorient01 = -1;
    iw0 = grid.PrevVertex(iw1, edir01);
  }
  else {
    eorient01 = 1;
    iw1 = grid.NextVertex(iw1, edir01);
    iw0 = grid.NextVertex(iw1, edir01);
  }

  int mask = (1 << edge_index);

  if ((facet_vertex_index ^ mask) == 0) {
    iw2 = grid.NextVertex(iw1, edir12);
    eorient12 = 1;
  }
  else {
    iw2 = grid.PrevVertex(iw1, edir12);
    eorient12 = -1;
  }
}

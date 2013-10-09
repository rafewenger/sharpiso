
#include "sharpiso_types.h"
#include "sharpiso_get_gradients.h"

#include "ijkprint.txx"

using namespace std;
using namespace IJK;
using namespace SHARPISO;

int main()
{
  VOXEL voxel;
  COORD_TYPE min_coord[DIM3];
  COORD_TYPE max_coord[DIM3];
  COORD_TYPE spacing[DIM3];

  for (int d = 0; d < DIM3; d++) {
    min_coord[d] = d+1;
    max_coord[d] = d+11;
    spacing[d] = (d+1)/2.0;
  }

  voxel.SetVertexCoord(min_coord, max_coord);

  for (int i = 0; i < voxel.NumVertices(); i++) {
    print_coord3D(cout, voxel.VertexCoord(i));
    cout << endl;
  }

  voxel.SetOffsetVertexCoord(spacing, 0.5);

  for (int i = 0; i < voxel.NumVertices(); i++) {
    print_coord3D(cout, voxel.VertexCoord(i));
    cout << endl;
  }

  return 0;
}

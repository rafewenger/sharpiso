
#include "sharpiso_grids.h"
#include "sharpiso_types.h"

using namespace std;
using namespace IJK;
using namespace SHARPISO;

int main()
{
  int axis_size3[DIM3] = { 10, 10, 10 };
  SHARPISO_SCALAR_GRID scalar_grid;
  COORD_TYPE coord[DIM3] = { 5.5, 1.1, 1 };
  GRID_COORD_TYPE cube_coord[DIM3];
  int cube_index = 115;

  scalar_grid.SetSize(DIM3, axis_size3);

  scalar_grid.ComputeCoord(cube_index, cube_coord);
  
  if (scalar_grid.CubeContainsPoint(cube_index, coord)) {
    cerr << "Cube " << cube_index << " ";
    ijkgrid_output_coord(cerr, DIM3, cube_coord);
    cerr << " contains coord ";
    ijkgrid_output_coord(cerr, DIM3, coord);
    cerr << endl;
  }
  else {
    cerr << "Cube " << cube_index << " ";
    ijkgrid_output_coord(cerr, DIM3, cube_coord);
    cerr << " does not contain coord ";
    ijkgrid_output_coord(cerr, DIM3, coord);
    cerr << endl;
  }
  
}

// Test scalar grid functions

#include "ijkgrid_macros.h"

#include "sharpiso_grids.h"
#include "sharpiso_types.h"

using namespace SHARPISO;

// routines
void output_scalar_grid(const SHARPISO_SCALAR_GRID_BASE & sgrid);
void check_scalar_grid(const SHARPISO_SCALAR_GRID_BASE & sgrid);
void set_scalar_grid(SHARPISO_SCALAR_GRID_BASE & sgrid);
void output_grid_cubes_containing
(const SHARPISO_SCALAR_GRID_BASE & sgrid, const SCALAR_TYPE s);

int main()
{
  const int DIM2 = 2;
  const int DIM3 = 3;
  SHARPISO_SCALAR_GRID sgrid;
  int axis_size2[DIM2] = { 3, 4 };
  int axis_size3[DIM3] = { 3, 4, 3 };

  using namespace std;

  try {

    sgrid.SetSize(DIM2, axis_size2);
    set_scalar_grid(sgrid);
    output_scalar_grid(sgrid);
    output_grid_cubes_containing(sgrid, 3.5);

    sgrid.SetSize(DIM3, axis_size3);
    set_scalar_grid(sgrid);
    output_scalar_grid(sgrid);
    output_grid_cubes_containing(sgrid, 16.5);
  }
  catch (IJK::ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

  return 0;
}

// **************************************************
// Set routines
// **************************************************

void set_scalar_grid(SHARPISO_SCALAR_GRID_BASE & sgrid)
{
  for (int iv = 0; iv < sgrid.NumVertices(); iv++) {
    SCALAR_TYPE s = iv + 0.1;
    sgrid.Set(iv, s);
  }
}

// **************************************************
// I/O routines
// **************************************************

void output_scalar_grid
(const SHARPISO_SCALAR_GRID_BASE & sgrid)
{
  std::cout << "dimension "<< sgrid.Dimension() << std::endl;
  output_scalar_grid(std::cout, sgrid);
}

void output_grid_cubes_containing
(const SHARPISO_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE s)
{
  using namespace std;

  typedef SHARPISO_SCALAR_GRID_BASE::VERTEX_INDEX_TYPE VTYPE;
  typedef SHARPISO_SCALAR_GRID_BASE::NUMBER_TYPE NTYPE;

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VTYPE) {

    cout << "Cube: " << icube;
    cout << "  Scalars:";
    for (NTYPE k = 0; k < scalar_grid.NumCubeVertices(); k++) {
      VTYPE iv = scalar_grid.CubeVertex(icube, k);
      cout << " " << scalar_grid.Scalar(iv);
    }

    if (is_gt_cube_min_le_cube_max(scalar_grid, icube, s)) {
      cout << ". Contains " << s << "." << endl;
    }
    else {
      cout << ". Does not contain " << s << "." << endl;
    }
  }
  cout << endl;
  cout << endl;
}




/// Test sharpiso_scalar routines

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "sharpiso_grids.h"
#include "sharpiso_types.h"

#include "sharpiso_scalar.txx"
#include "ijkcube.txx"

using namespace std;
using namespace IJK;
using namespace SHARPISO;

// type definitions
typedef UNIT_CUBE<int, int, int> TEST_CUBE;

// global variables

// output routines
void out_scalar_field_at_cube_vertices
(const TEST_CUBE & cube,
 const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar);
void out_intersects_diagonals
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE grad[DIM3], 
 const SCALAR_TYPE scalar, const SCALAR_TYPE isovalue);
void out_intersects_cube
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue);


// check routines
void check_intersects_cube
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue);

// miscellaneous routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{

  GRADIENT_COORD_TYPE gradA[DIM3] = { 1, 1, 1 };
  GRADIENT_COORD_TYPE gradB[DIM3] = { 1, 2, 3 };
  GRADIENT_COORD_TYPE gradC[DIM3] = { -3, 2, -1 };
  COORD_TYPE coord0[DIM3] = { 0, 0, 0 };
  COORD_TYPE coord1[DIM3] = { 1, 1, 1 };


  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);


    TEST_CUBE cube;

    cube.SetDimension(DIM3);
    out_scalar_field_at_cube_vertices(cube, coord0, gradA, 0);
    out_intersects_diagonals(cube, coord0, gradA, 0, 0.5);
    out_intersects_cube(cube, coord0, gradA, 0, 0.5);
    out_intersects_diagonals(cube, coord0, gradA, 0, 4);
    out_intersects_cube(cube, coord0, gradA, 0, 4);

    out_scalar_field_at_cube_vertices(cube, coord0, gradB, 0);
    out_intersects_diagonals(cube, coord0, gradB, 0, 1);
    out_intersects_cube(cube, coord0, gradB, 0, 1);
    out_intersects_diagonals(cube, coord0, gradB, 0, 2);
    out_intersects_cube(cube, coord0, gradB, 0, 2);

    out_scalar_field_at_cube_vertices(cube, coord1, gradB, 10);
    out_intersects_diagonals(cube, coord1, gradB, 10, 11);
    out_intersects_cube(cube, coord1, gradB, 10, 11);
    out_intersects_diagonals(cube, coord1, gradB, 10, 8);
    out_intersects_cube(cube, coord1, gradB, 10, 8);

    check_intersects_cube(cube, coord0, gradA, 0, 0.5);
    check_intersects_cube(cube, coord0, gradA, 0, 4);
    check_intersects_cube(cube, coord0, gradB, 0, 1);
    check_intersects_cube(cube, coord0, gradB, 0, 2);
    check_intersects_cube(cube, coord1, gradB, 10, 11);
    check_intersects_cube(cube, coord1, gradB, 8, 8);

    cout << "Passed check_intersects_cube." << endl;
  }
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(30);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// **************************************************
// I/O routines
// **************************************************

template <typename CTYPE>
void output_coord
(ostream & output, const int dimension, CTYPE * coord)
{
  cout << "(";
  for (int d = 0; d < dimension; d++) {
    cout << coord[d];
    if (d+1 < dimension) { cout << ","; }
  }
  cout << ")";
}

void out_scalar_field_at_cube_vertices
(const TEST_CUBE & cube,
 const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar)
{
  typedef TEST_CUBE::NUMBER_TYPE NTYPE;
  typedef TEST_CUBE::COORD_TYPE CTYPE;

  CTYPE vertex_coord[DIM3];

  cout << "Scalar field F(p) determined by coord ";
  output_coord(cout, DIM3, coord);
  cout << " grad ";
  output_coord(cout, DIM3, gradient);
  cout << " scalar " << scalar;
  cout << endl;

  for (NTYPE iv = 0; iv < cube.NumVertices(); iv++) {

    SCALAR_TYPE s =
      compute_gradient_based_scalar
      (cube.VertexCoord(iv), coord, gradient, scalar);

    cout << "F[";
    output_coord(cout, DIM3, cube.VertexCoord(iv));
    cout << "] = " << s << endl;
  }
  cout << endl;
}

void out_intersects_diagonal
(const TEST_CUBE & cube, const int idiagonal, 
 const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue)
{
  typedef TEST_CUBE::COORD_TYPE CTYPE;
  
  const CTYPE * end0_coord = cube.DiagonalCoord(idiagonal, 0);
  const CTYPE * end1_coord = cube.DiagonalCoord(idiagonal, 1);
  
  cout << "Coord ";
  output_coord(cout, DIM3, coord);
  cout << " grad ";
  output_coord(cout, DIM3, gradient);
  cout << " scalar " << scalar;
  cout << " isoval " << isovalue;

  if (iso_intersects_line_segment
      (end0_coord, end1_coord, coord, gradient, scalar, isovalue)) {
    cout << " intersects ";
  }
  else {
    cout << " does not intersect ";
  }

  cout << "diagonal ";
  output_coord(cout, DIM3, end0_coord);
  cout << " ";
  output_coord(cout, DIM3, end1_coord);
  cout << endl;
}

void out_intersects_diagonals
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue)
{
  typedef TEST_CUBE::NUMBER_TYPE NTYPE;

  for (NTYPE i = 0; i < cube.NumVertices(); i++) {
    out_intersects_diagonal
      (cube, i, coord, gradient, scalar, isovalue);
  }
  cout << endl;
}

void out_intersects_cube
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue)
{
  bool flag = 
    iso_intersects_cube(cube, coord, gradient, scalar, isovalue);

  cout << "Coord ";
  output_coord(cout, DIM3, coord);
  cout << " grad ";
  output_coord(cout, DIM3, gradient);
  cout << " scalar " << scalar;
  cout << " isoval " << isovalue;

  if (flag) { cout << " intersects "; }
  else { cout << " does not intersect "; }

  cout << "cube ";
  output_coord(cout, DIM3, cube.VertexCoord(0));
  cout << " ";
  cout << " edge length " << cube.EdgeLength() << endl;
  cout << endl;

}

// **************************************************
// Check routines
// **************************************************

void check_intersects_cube
(const TEST_CUBE & cube, const COORD_TYPE coord[DIM3],
 const GRADIENT_COORD_TYPE gradient[DIM3], 
 const SCALAR_TYPE scalar,
 const SCALAR_TYPE isovalue)
{
  typedef TEST_CUBE::DIMENSION_TYPE DTYPE;
  typedef TEST_CUBE::NUMBER_TYPE NTYPE;
  typedef TEST_CUBE::COORD_TYPE CTYPE;

  IJK::PROCEDURE_ERROR error("check_intersects_cube");

  bool flag =
    iso_intersects_cube(cube, coord, gradient, scalar, isovalue);

  bool flag2(false);
  NTYPE idiagonal(0);
  for (NTYPE k = 0; k < cube.NumDiagonals(); k++) {

    const CTYPE * end0_coord = cube.DiagonalCoord(k, 0);
    const CTYPE * end1_coord = cube.DiagonalCoord(k, 1);

    if (iso_intersects_line_segment
        (end0_coord, end1_coord, coord, gradient, scalar, isovalue)) {
      flag2 = true; 
      idiagonal = k;
      break;
    }
  }

  if (flag != flag2) {
    error.AddMessage("Error in computing cube-isosurface intersection.");
    error.AddMessage("  Isovalue: ", isovalue);
    if (flag) {
      error.AddMessage
        ("  Isosurface intersects cube but not any cube diagonal.");
    }
    else {
      error.AddMessage
        ("  Isosurface does not intersect cube but intersects cube diagonal ",
         idiagonal, ".");
    };
    throw error;
  }

}

// **************************************************
// Miscellaneous routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_msg()
{
  cerr << "Usage: testscalar" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void parse_command_line(int argc, char **argv)
{

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = string(argv[iarg]);
    usage_error();

    iarg++;
  }

}

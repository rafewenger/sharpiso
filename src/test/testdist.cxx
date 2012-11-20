/// Test compute_distance_to_gfield_plane routines

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012 Rephael Wenger

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

#include <cstdlib>
#include <new>
#include <iostream>

#include "sharpiso_scalar.txx"

using namespace std;
using namespace SHARPISO;

// global variables

// output routines
void output_distance_to_gfield_plane
(const GRADIENT_COORD_TYPE gfield_gradient[DIM3], 
 const COORD_TYPE gfield_point[DIM3],
 const SCALAR_TYPE gfield_point_scalar,
 const COORD_TYPE point[DIM3], const SCALAR_TYPE plane_scalar,
 const COORD_TYPE distance);
void output_distance_vary_gfield_point_scalar
(const GRADIENT_COORD_TYPE gfield_gradient[DIM3], 
 const COORD_TYPE gfield_point[DIM3],
 const COORD_TYPE point[DIM3], const SCALAR_TYPE plane_scalar);

// miscellaneous routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{
  const GRADIENT_COORD_TYPE sqrt3 = std::sqrt(3);

  GRADIENT_COORD_TYPE gradA[DIM3] = { 1, 0, 0 };
  GRADIENT_COORD_TYPE gradB[DIM3] = { 0, 1, 0 };
  GRADIENT_COORD_TYPE gradC[DIM3] = { 0, 0, 1 };
  GRADIENT_COORD_TYPE gradD[DIM3] = { 1/sqrt3, 1/sqrt3, 1/sqrt3 };
  GRADIENT_COORD_TYPE gradE[DIM3] = { 2/sqrt3, 2/sqrt3, 2/sqrt3 };

  COORD_TYPE coord0[DIM3] = { 0, 0, 0 };
  COORD_TYPE coord1[DIM3] = { 1, 1, 1 };
  COORD_TYPE coord2[DIM3] = { 2/sqrt3, 2/sqrt3, 2/sqrt3 };


  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    output_distance_vary_gfield_point_scalar(gradA, coord0, coord1, 1);
    cout << endl;

    output_distance_vary_gfield_point_scalar(gradB, coord0, coord1, 1);
    cout << endl;

    output_distance_vary_gfield_point_scalar(gradC, coord0, coord1, 1);
    cout << endl;

    output_distance_vary_gfield_point_scalar(gradD, coord0, coord2, 1);
    cout << endl;

    output_distance_vary_gfield_point_scalar(gradE, coord0, coord2, 1);
    cout << endl;

    output_distance_vary_gfield_point_scalar(gradE, coord0, coord2, 2);
    cout << endl;
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// **************************************************
// I/O routines
// **************************************************

void output_distance_vary_gfield_point_scalar
(const GRADIENT_COORD_TYPE gfield_gradient[DIM3], 
 const COORD_TYPE gfield_point[DIM3],
 const COORD_TYPE point[DIM3], const SCALAR_TYPE plane_scalar)
{
  COORD_TYPE distance;
  for (SCALAR_TYPE s = 0; s <= 2; s += 0.5) {

    compute_distance_to_gfield_plane
      (gfield_gradient, gfield_point, s, point, plane_scalar, distance);

    output_distance_to_gfield_plane
      (gfield_gradient, gfield_point, s, point, plane_scalar, distance);
  };

}

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

void output_distance_to_gfield_plane
(const GRADIENT_COORD_TYPE gfield_gradient[DIM3], 
 const COORD_TYPE gfield_point[DIM3],
 const SCALAR_TYPE gfield_point_scalar,
 const COORD_TYPE point[DIM3], const SCALAR_TYPE plane_scalar,
 const COORD_TYPE distance)
{
  cout << "Distance from ";
  output_coord(cout, DIM3, point);
  cout << " to {x : (x-";
  output_coord(cout, DIM3, gfield_point);
  cout << ") dot ";
  output_coord(cout, DIM3, gfield_gradient);
  cout << ") + " << gfield_point_scalar << " = " << plane_scalar
       << "} is " << distance << endl;
}

// **************************************************
// Check routines
// **************************************************

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
  cerr << "Usage: testdist" << endl;
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

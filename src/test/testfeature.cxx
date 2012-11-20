/// Test sharpiso_feature routines

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

#include "sharpiso_feature.h"

using namespace std;
using namespace IJK;
using namespace SHARPISO;

// global variables

// output routines
void output_cube(const OFFSET_CUBE_111 & cube);

// miscellaneous routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{
  try {

    std::set_new_handler(memory_exhaustion);

    OFFSET_CUBE_111 cube(0.1);

    output_cube(cube);

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

void output_cube(const OFFSET_CUBE_111 & cube)
{
  for (int i = 0; i < cube.NumVertices(); i++) {
    output_coord(cout, SHARPISO::DIM3, cube.VertexCoord(i));
    cout << endl;
  }
  cout << endl;
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
  cerr << "Usage: testfeature" << endl;
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

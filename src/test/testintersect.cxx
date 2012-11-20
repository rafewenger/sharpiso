/// Test intersect routines

#include <cstdlib>
#include <iostream>

#include "ijk.txx"
#include "ijkprint.txx"

#include "sharpiso_findIntersect.h"

using namespace std;
using namespace IJK;
using namespace SHARPISO;

// global variables

// output routines
void output_three_facet_line_location
(const COORD_TYPE facet_p0[DIM3],
 const COORD_TYPE line_p0[DIM3], const COORD_TYPE line_dir[DIM3]);
void output_facet_line_location
(const COORD_TYPE facet_p0[DIM3], const NUM_TYPE facet_orth_dir,
 const COORD_TYPE line_p0[DIM3], const COORD_TYPE line_dir[DIM3]);

// miscellaneous routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{
  COORD_TYPE origin[DIM3] = { 0, 0, 0};

  COORD_TYPE pointX[DIM3] = { 1, 0, 0 };
  COORD_TYPE pointY[DIM3] = { 0, 1, 0 };
  COORD_TYPE pointZ[DIM3] = { 0, 0, 1 };

  COORD_TYPE line_pA[DIM3] = { 0.1, 0.2, 0.3 };

  COORD_TYPE line_dirX[DIM3] = { 1, 0, 0 };
  COORD_TYPE line_dirY[DIM3] = { 0, 1, 0 };
  COORD_TYPE line_dirZ[DIM3] = { 0, 0, 1 };
  COORD_TYPE line_dir110[DIM3] = { 1, 1, 0 };

  try {

    output_three_facet_line_location(origin, line_pA, line_dirX);
    output_three_facet_line_location(origin, line_pA, line_dirY);
    output_three_facet_line_location(origin, line_pA, line_dirZ);

    cout << endl;
    output_three_facet_line_location(origin, line_pA, line_dir110);


  }
  catch (IJK::ERROR error) {
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

// Output intersection of three facets and a line
void output_three_facet_line_location
(const COORD_TYPE facet_p0[DIM3],
 const COORD_TYPE line_p0[DIM3], const COORD_TYPE line_dir[DIM3])
{
  for (NUM_TYPE d = 0; d < DIM3; d++)
    { output_facet_line_location(facet_p0, d, line_p0, line_dir); };
}

void output_facet_line_location
(const COORD_TYPE facet_p0[DIM3], const NUM_TYPE facet_orth_dir,
 const COORD_TYPE line_p0[DIM3], const COORD_TYPE line_dir[DIM3])
{
  COORD_TYPE end[2][DIM3];
  COORD_TYPE max_small_grad = 0.0001;
  bool flag_intersect;

  intersect_line_square_cylinder
    (facet_p0, facet_orth_dir, 1, 50, line_p0, line_dir, max_small_grad,
     end, flag_intersect);

  cout << "  Line: ";
  print_coord3D(cout, line_p0, " + t");
  print_coord3D(cout, line_dir, " ");
  if (flag_intersect) {
    cout << " intersects ";
  }
  else {
    cout << " does not intersect ";
  }

  cout << "facet square: ";
  print_coord3D(cout, facet_p0, " ");
  cout << "orth dir: " << facet_orth_dir;

  if (flag_intersect) {
    cout << " at (";
    print_coord3D(cout, end[0], ",");
    print_coord3D(cout, end[1], ").");
    cout << endl;
  }
  else {
    cout << "." << endl;
  }
  
}


// Test bin


#include "sharpiso_grids.h"

typedef SHARPISO::BIN_GRID <int> BGRID;

std::ostream & operator << 
  (std::ostream & output, const SHARPISO::BIN<int> & bin) 
{
  output << "(";
  for (int i = 0; i < bin.list.size(); i++) {
    output << bin.list[i];
    if (i+1 < bin.list.size()) { output << ","; }
  }
  output << ")";
  return(output);
}


// routines
void set_bin_grid(BGRID & bgrid);
void output_bin_grid(const BGRID & bgrid);


using namespace IJK;
using namespace SHARPISO;
using namespace std;


int main()
{
  const int DIM1 = 1;
  const int DIM2 = 2;
  const int DIM3 = 3;
  int axis_size1[DIM1] = { 4 };
  int axis_size2[DIM2] = { 7, 5 };
  int axis_size3[DIM3] = { 7, 3, 2 };
  int axis_size3B[DIM3] = { 5, 6, 7 };
  int min_coord2[DIM2] = { 2, 1 };
  int max_coord2[DIM2] = { 5, 3 };

  BGRID bgrid;

  using namespace std;

  try {

    bgrid.SetSize(DIM2, axis_size2);
    set_bin_grid(bgrid);
    output_bin_grid(bgrid);
  }
  catch (ERROR error) {
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

void set_bin_grid(BGRID & bgrid)
{
  for (int iv = 0; iv < bgrid.NumVertices(); iv++) {
    bgrid.Insert(iv, 10*iv);
    bgrid.Insert(iv, 10*iv+1);
    bgrid.Insert(iv, 10*iv+2);
  }
}

// **************************************************
// I/O routines
// **************************************************

void output_bin_grid
(const BGRID & bgrid)
{
  std::cout << "dimension "<< bgrid.Dimension() << std::endl;
  output_object_grid(std::cout, bgrid);
  std::cout << std::endl;
}


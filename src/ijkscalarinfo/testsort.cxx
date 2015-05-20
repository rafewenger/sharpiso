
#include <iostream>

#include "ijk.txx"
#include "ijksortgrid.txx"

using namespace std;

int main()
{
  const int NUM_SCALAR = 20;
  float scalar[NUM_SCALAR];
  int index[NUM_SCALAR];

  for (int i = 0; i < NUM_SCALAR; i++) 
    { scalar[i] = (NUM_SCALAR - i)%5; }

  IJK::sort_grid_vertices(scalar, NUM_SCALAR, index);

  cout << "Scalar values: ";
  for (int i = 0; i < NUM_SCALAR; i++)
    { cout << scalar[i] << " "; }
  cout << endl << endl;

  cout << "Sorted values: ";
  for (int i = 0; i < NUM_SCALAR; i++)
    { cout << scalar[index[i]] << " "; }
  cout << endl << endl;
}


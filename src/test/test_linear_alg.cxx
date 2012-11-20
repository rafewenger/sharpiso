
#include <iostream>

#include "sharpiso_linear_alg.txx"

using namespace std;
using namespace SHARPISO;

int main()
{
  const int DIM3 = 3;
  float v[DIM3];
  float w[DIM3];
  float zero_tolerance;
  bool is_v_zero;

  cout << "Enter vector (3 values): ";
  cin >> v[0] >> v[1] >> v[2];

  cout << "Enter tolerance: ";
  cin >> zero_tolerance;

  normalize_vector_3D(v, w, is_v_zero, zero_tolerance*zero_tolerance);

  if (is_v_zero) {
    cout << "Zero vector: ";
  }
  for (int d = 0; d < DIM3; d++) 
    { cout << " " << w[d]; }
  cout << endl;
}

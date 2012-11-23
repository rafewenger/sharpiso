// Test time for sorting.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "ijktime.txx"

// Data structure
class DATA {
public:
  int index;
  float x;
};

class DATA_COMPARE {

public:
  std::vector<DATA> data;

	DATA_COMPARE(std::vector<DATA> & data )
  { this->data = data; }

	bool operator () (int i,int j)
	{
    if (data[i].index == data[j].index)
      { return (data[i].x < data[j].x); }
    else
      { return (data[i].index > data[j].index); }
	}
  

};

class DATA_COMPARE2 {

public:
  const DATA * data;

	DATA_COMPARE2(const DATA * data)
  { this->data = data; }

	bool operator () (int i,int j)
	{
    if (data[i].index == data[j].index)
      { return (data[i].x < data[j].x); }
    else
      { return (data[i].index > data[j].index); }
	}
  

};

// routines
void set_int_list(const int n, std::vector<int> & int_list);
void set_data_list(const int n, std::vector<DATA> & data_list);
void set_float_list(const int n, std::vector<float> & float_list);
void init_index_sorted(const int n, std::vector<int> & index_sorted);


using namespace IJK;
using namespace std;

int main(int argc, char ** argv)
{
  vector<int> int_list;
  vector<float> float_list;
  vector<DATA> data_list;
  vector<int> index_sorted;
  clock_t t0, t1;
  float seconds;
  int n;

  cout << "Enter number of integers: ";
  cin >> n;

  srandom(100);

  set_int_list(n, int_list);
  set_float_list(n, float_list);
  set_data_list(n, data_list);

  t0 = clock();

  sort(int_list.begin(), int_list.end());

  t1 = clock();
  clock2seconds(t1-t0, seconds);

  cout << "int sort time (sec): " << seconds << endl;


  t0 = clock();

  sort(float_list.begin(), float_list.end());

  t1 = clock();
  clock2seconds(t1-t0, seconds);

  cout << "float sort time (sec): " << seconds << endl;

  init_index_sorted(n, index_sorted);
  DATA_COMPARE data_compare(data_list);

  t0 = clock();

  sort(index_sorted.begin(), index_sorted.end(), data_compare);

  t1 = clock();
  clock2seconds(t1-t0, seconds);

  cout << "data sort time (sec): " << seconds << endl;


  init_index_sorted(n, index_sorted);
  DATA_COMPARE2 data_compare2(&(data_list[0]));

  t0 = clock();

  sort(index_sorted.begin(), index_sorted.end(), data_compare2);

  t1 = clock();
  clock2seconds(t1-t0, seconds);

  cout << "data sort time using DATA_COMPARE2 (sec): " << seconds << endl;

  return 0;
}

void set_int_list(const int n, vector<int> & int_list)
{
  int_list.resize(n);

  for (int i = 0; i < n; i++) {
    int k = random() % (10*n);
    int_list.push_back(k);
  }
}

void set_float_list(const int n, vector<float> & float_list)
{
  float_list.resize(n);

  for (int i = 0; i < n; i++) {
    int k = float(random()) / RAND_MAX;
    float_list.push_back(k);
  }
}

void set_data_list(const int n, vector<DATA> & data_list)
{
  DATA data;

  data_list.resize(n);
  for (int i = 0; i < n; i++) {
    data.index = (random() % 3);
    data.x = float(random()) / RAND_MAX;
    data_list.push_back(data);
  }
}

void init_index_sorted(const int n, vector<int> & index_sorted)
{
  index_sorted.resize(n);

  for (int i = 0; i < n; i++) 
    { index_sorted[i] = i; }
}

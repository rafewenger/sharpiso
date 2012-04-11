// print statements and other help routines for findEdgeCount 

// **************************************************
// Output routines
// **************************************************
#include <iostream>

#include"findEdgeCountIO.h"

void output_edges
(const int dim, const COORD_TYPE * coord, const int numv,
 const VERTEX_INDEX * edge_vert, const int nume)
{
  for (int i = 0; i < nume; i++) {
    VERTEX_INDEX iv0 = edge_vert[2*i];
    VERTEX_INDEX iv1 = edge_vert[2*i+1];

    cout << "[";
    print_list(cout, coord+iv0*dim, dim);
    cout << " ";
    print_list(cout, coord+iv1*dim, dim);
    cout << "]" << endl;
  }
}

void output_vert_degree
( const int numv, vector <int> &vert_degree)
{
  vector <int> deg_one;
  vector <int> deg_two;
  vector <int> deg_threeandmore; 
  int cnt1(0),cnt2(0),cnt3(0); // number of vertices with degree one or  3 or more  
  for (int i=0;i<numv;i++)
  {
    if (vert_degree[i] == 1){
      deg_one.push_back(i);
      cnt1 ++;
    }
    else if (vert_degree[i] == 2){
      deg_two.push_back(i);
      cnt2++;
    }
    else{
      deg_threeandmore.push_back(i);
      cnt3++;
    }
  }

  cout <<"Vertices with degree 1 or 3 or more ["<<cnt1 + cnt3<<"]" << endl;
  cout <<"Vertices with degree 1              ["<<cnt1<<"]" << endl;
  cout <<"Vertices with degree 2              ["<<cnt2<<"]" << endl;
  cout <<"Total number of vertices            ["<<numv<<"]" << endl;  
}



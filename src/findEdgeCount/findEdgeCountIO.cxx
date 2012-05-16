// print statements and other help routines for findEdgeCount 

// **************************************************
// Output routines
// **************************************************
#include <iostream>
#include <cstdlib>     // function exit() is in cstdlib
#include <fstream>     // class ofstream() is in fstream
#include <string>  
#include"findEdgeCountIO.h"

// global variables 
int cnt0(0),cnt1(0),cnt2(0),cnt3(0); // number of vertices with degree one or  3 or more
int cntMoreThan3(0); 

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
( const int numv, vector <int> &vert_degree, const string &fname)
{
  vector <int> deg_one;
  vector <int> deg_two;
  vector <int> deg_threeandmore;  
  for (int i=0;i<numv;i++)
  {
    if (vert_degree[i] == 0)
      cnt0++;
    else if (vert_degree[i] == 1){
      //deg_one.push_back(i);
      cnt1++;
    }
    else if (vert_degree[i] == 2){
      //deg_two.push_back(i);
      cnt2++;
    }
    else if (vert_degree[i] == 3){
      //deg_threeandmore.push_back(i);
      cnt3++;
    }
    else 
    {
      cntMoreThan3++;
    }
  }

  cout <<"file name : "<<fname <<endl;
  cout <<"Vertices with degree 0              ["<<cnt0<<"]" << endl;
  cout <<"Vertices with degree 1              ["<<cnt1<<"]" << endl;
  cout <<"Vertices with degree 2              ["<<cnt2<<"]" << endl;
  cout <<"Vertices with degree 3              ["<<cnt3<<"]" << endl;
  cout <<"Vertices with degree > 3            ["<<cntMoreThan3<<"]" << endl;
  cout <<"Vertices with degree 1 or 3 or more ["<<cnt1 + cnt3 + cntMoreThan3<<"]" << endl;
  cout <<"Total number of non 0 vertices      ["<<cnt1+cnt2+cnt3+cntMoreThan3<<"]" 
    << endl;  
  cout <<"Total number of vertices            ["<<numv<<"]" << endl;  
}

void output_vert_degree_2_file 
(const int numv, vector <int> &vert_degree)
{

 vector <int> deg_one;
  vector <int> deg_two;
  vector <int> deg_threeandmore;  
  for (int i=0;i<numv;i++)
  {
    if (vert_degree[i] == 0)
      cnt0++;
    else if (vert_degree[i] == 1){
      //deg_one.push_back(i);
      cnt1++;
    }
    else if (vert_degree[i] == 2){
      //deg_two.push_back(i);
      cnt2++;
    }
    else if (vert_degree[i] == 3){
      //deg_threeandmore.push_back(i);
      cnt3++;
    }
    else 
    {
      cntMoreThan3++;
    }
  }

  cout << cnt0 <<" " <<  cnt1 << " "
	   << cnt2 <<" " << cnt3 <<" "
	   << cntMoreThan3 <<" "
	   << cnt1 + cnt3 + cntMoreThan3<<" "
	   << cnt1+cnt2+cnt3+cntMoreThan3<<" "
	   << numv<<endl; 
}


void output_vert_degree_2_file 
(const int numv, vector <int> &vert_degree, const string &fname )
{

 vector <int> deg_one;
  vector <int> deg_two;
  vector <int> deg_threeandmore;  
  for (int i=0;i<numv;i++)
  {
    if (vert_degree[i] == 0)
      cnt0++;
    else if (vert_degree[i] == 1){
      //deg_one.push_back(i);
      cnt1++;
    }
    else if (vert_degree[i] == 2){
      //deg_two.push_back(i);
      cnt2++;
    }
    else if (vert_degree[i] == 3){
      //deg_threeandmore.push_back(i);
      cnt3++;
    }
    else 
    {
      cntMoreThan3++;
    }
  }
 /*
  cout <<fname<<" "<< cnt0 <<" " <<  cnt1 << " "
	   << cnt2 <<" " << cnt3 <<" "
	   << cntMoreThan3 <<" "
	   << cnt1 + cnt3 + cntMoreThan3<<" "
	   << cnt1+cnt2+cnt3+cntMoreThan3<<" "
	   << numv<<endl; 
	   */
	   cout <<fname<<" "<< cnt1 + cnt3 + cntMoreThan3;
}


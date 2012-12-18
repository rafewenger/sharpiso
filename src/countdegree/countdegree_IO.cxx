// print statements and other help routines for findEdgeCount 

// **************************************************
// Output routines
// **************************************************
#include <iostream>
#include <cstdlib>     // function exit() is in cstdlib
#include <fstream>     // class ofstream() is in fstream
#include <string>  
#include <iomanip>
#include <cmath>
#include"countdegree_IO.h"

// global variables 
int cnt0(0),cnt1(0),cnt2(0),cnt3(0); // number of vertices with degree one or  3 or more
int cntMoreThan3(0); 

extern int real_degree_3_verts;
extern int real_degree_1_verts;

void output_edges
(const int dim, 
		const COORD_TYPE * coord,
		const int numv,
		const VERTEX_INDEX * edge_vert,
		const int nume,
		vector<int> vert_degree)
{
	for (int i = 0; i < nume; i++) {
		VERTEX_INDEX iv0 = edge_vert[2*i];
		VERTEX_INDEX iv1 = edge_vert[2*i+1];

		cout << "[";
		print_list(cout, coord+iv0*dim, dim);
		cout << " ";
		print_list(cout, coord+iv1*dim, dim);
		cout << "]"<< endl;
	}
}


void print_edge_info(int numv, const COORD_TYPE * coord, vector <int> vert_degree)
{
	// header
	cout <<"\n"<<setw(4)<<"vert"<<setw(28)<<"position"<<setw(5)<<"dg1"<<setw(5)
										  <<"dg3"<<setw(5)<<"dg>3"<<endl;
	float crd[3]={0.0};

	for (int i=0;i<numv;i++)
	{
		crd[0]= coord[3*i];crd[1]= coord[3*i+1];crd[2]= coord[3*i+2];

		if (vert_degree[i]==1){
			print_point(coord, vert_degree,  i, 1);
		}
		else if (vert_degree[i]==3){
			print_point(coord, vert_degree,  i, 2);
		}
		else if (vert_degree[i]>3){
			print_point(coord, vert_degree,  i, 3);
		}
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
			cnt1++;
		}
		else if (vert_degree[i] == 2){
			cnt2++;
		}
		else if (vert_degree[i] == 3){
			cnt3++;
		}
		else
		{
			cntMoreThan3++;
		}
	}

	// actual vertices with degree 3 and degree 1
	cnt3 = abs(cnt3-real_degree_3_verts);
	cnt1 = abs(cnt1-real_degree_1_verts);

	cout <<"File name : "<< fname <<endl;
	cout << "Vertices with degree 0              [" << cnt0 <<"]" << endl;
	cout << "Vertices with degree 1              [" << cnt1 <<"]" << endl;
	cout << "Vertices with degree 2              [" << cnt2 <<"]" << endl;
	cout << "Vertices with degree 3              [" << cnt3 <<"]" << endl;
	cout << "Vertices with degree > 3            [" << cntMoreThan3 <<"]" << endl;
	cout << "Vertices with degree 1 or 3 or more [" << cnt1 + cnt3 + cntMoreThan3 << "]" << endl;
	cout << "Total number of non 0 vertices      [" << cnt1+cnt2+cnt3+cntMoreThan3 << "]"
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


void output_short_info
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
			cnt1++;
		}
		else if (vert_degree[i] == 2){
			cnt2++;
		}
		else if (vert_degree[i] == 3){
			cnt3++;
		}
		else
		{
			cntMoreThan3++;
		}
	}
	// actual vertices with degree 3 and degree 1
	cnt3 = abs(cnt3-real_degree_3_verts);
	cnt1 = abs(cnt1-real_degree_1_verts);
	cout <<fname<<" , "<< cnt1 + cnt3 + cntMoreThan3<<endl;
}
/// HELPER FUNCTIONS
void  print_point (const COORD_TYPE * coordList,
		vector <int> vert_degree,
		const int vertId,
		const int deg)
{
	cout <<setw(5)<<vertId;
	int vert = 3*vertId;
	cout << fixed;
	cout << setprecision (5) ;
	cout << " ["<<setw(6)<<coordList[vert]<<","<<setw(6)<<coordList[vert+1]<<","<<setw(6)<<coordList[vert+2]<<"]";
	int w = deg*4;
	cout <<setw(w)<< vert_degree[vertId]<<endl;
}


/// \file find edge.cxx
/// Example of find edge

/*
IJK: Isosurface Jeneration Code
Copyright (C) 2008 Arindam Bhattacharya

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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <time.h>


#include "ijk.txx"
#include "ijkIO.txx"

using namespace std;
using namespace IJK;

typedef float COORD_TYPE;

// global variables
int dimension = 3;
int mesh_dimension;
int num_vertices = 0;
int num_simplices = 0;
int num_edges = 0;
COORD_TYPE * vertex_coord = NULL;
int * simplex_vert = NULL;
vector<int> new_simplex_vert;
vector<int> edge_vert;
vector<int> vec;
vector<COORD_TYPE> new_vertex;
vector<int> L;
string input_filename;
string out_filename;
bool output_specified = false;
double input_angle;
// output routine 
void output_mesh_info();
bool edge_equals( int a, int b);
//---------------
class point3D
{
public:

	double x,y,z;
	point3D()
	{
		x=0;
		y=0;
		z=0;
	}
	point3D(double a , double b , double c)
	{
		x=a;
		y=b;
		z=c;
	}
	point3D(int t)
	{
		x=vertex_coord[(t*3+0)];
		y=vertex_coord[(t*3+1)];
		z=vertex_coord[(t*3+2)];
	}
	void set(point3D m)
	{
		x=m.x;
		y=m.y;
		z=m.z;
	}

};
int x=0;
//---------------

// misc routines
void memory_exhaustion();
void parse_command_line(int argc, char **argv);
void usage_error();

struct myclass {
	bool operator() (const int& i0, const int& i1) { 
		if( L[2*i0] < L[2*i1])
		{ return true; }
		else if ( L[2*i0] > L[2*i1])
		{ return false; }
		else {
			if(L[2*i0+1] < L[2*i1+1])
				return true;
			else
				return false;
		}
	}
} myobject;

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
	std::set_new_handler(memory_exhaustion);
	time_t start,end;

	parse_command_line(argc, argv);

	try {
		ifstream in(input_filename.c_str(), ios::in);
		if (!in.good()) {
			cerr << "Unable to open file " << input_filename << "." << endl;
			exit(30);
		};

		ijkinOFF(in, dimension, mesh_dimension,
			vertex_coord, num_vertices, simplex_vert, num_simplices);
		in.close();
		// checj mesh dimension 
		if (mesh_dimension != 2)
		{
		  cout <<"The input must be a triangle mesh."<<endl;
		  exit (10);
		}
		output_mesh_info();
		cerr <<"File output to "<< out_filename<<endl;
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
	}

	delete [] simplex_vert;;
	delete [] vertex_coord;

	return(0);
}

void addtovector_simplex(int vert[])
{
	new_simplex_vert.push_back(vert[0]);
	new_simplex_vert.push_back(vert[1]);
	new_simplex_vert.push_back(vert[0]);
}

void addtovetor_edge(int vert[])
{
	edge_vert.push_back(vert[0]);
	edge_vert.push_back(vert[1]);

}

//class defintiions

void otherVert(int t,int e, int edge[])
{
	point3D n1,n2;
	int ov1,ov2;
	for(int i=0;i<3;i++)
	{
		if(simplex_vert[3*t+i]!=edge[0] && simplex_vert[3*t+i]!=edge[1])
		{ov1=simplex_vert[3*t+i];
		}

		if(simplex_vert[3*e+i]!=edge[0] && simplex_vert[3*e+i]!=edge[1])
		{ov2=simplex_vert[3*e+i];
		}
	}

	point3D p3=point3D(ov1);point3D p4=point3D(ov2);
	point3D p1=point3D(edge[0]);point3D p2=point3D(edge[1]);

	point3D vec_a=point3D(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
	point3D vec_b=point3D(p4.x-p1.x,p4.y-p1.y,p4.z-p1.z);
	point3D vec_c=point3D(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);

	double norm1x=(vec_a.y*vec_b.z)-(vec_a.z*vec_b.y);
	double norm1y=(vec_a.z*vec_b.x)-(vec_a.x*vec_b.z);
	double norm1z=(vec_a.x*vec_b.y)-(vec_a.y*vec_b.x);
	double norm1=sqrt(norm1x*norm1x + norm1y*norm1y + norm1z*norm1z);

	double norm2x=(vec_c.y*vec_a.z) - (vec_c.z*vec_a.y);
	double norm2y=(vec_c.z*vec_a.x) - (vec_c.x*vec_a.z);
	double norm2z=(vec_c.x*vec_a.y) - (vec_c.y*vec_a.x);
	double norm2=sqrt(norm2x*norm2x + norm2y*norm2y + norm2z*norm2z);

	n1.x=norm1x; n1.y=norm1y; n1.z=norm1z;
	n2.x=norm2x; n2.y=norm2y; n2.z=norm2z;

	double angle= acos ( (  n1.x*n2.x +n1.y*n2.y+n1.z*n2.z) /(norm1*norm2)) * (180.0/M_PI);

	if(  angle > (180-input_angle) ){
		addtovector_simplex(edge);
		addtovetor_edge(edge);
	}
}


void findCommonEdges(vector<int> I,vector<int> E,int numl)
{ 
	int e;
	vector<int > c1(numl,0);
	vector<int > Etri(numl,0);
	int op1,op2; //other vertex other than the common edge.
	int edge[2];point3D n1,n2;

	for (int t=0;t<num_simplices;t++)
		for (int j=0;j<3;j++)
		{
			e=I[3*t+j];
			if(c1[e]==0)
				Etri[e]=t;
			else
			{

				edge[0]=E[2*e];edge[1]=E[2*e+1];
				otherVert(t,Etri[e],edge);
			}

			c1[e]++;
		}
}

// **************************************************
// OUTPUT_ROUTINE
// **************************************************

void output_mesh_info()
{
	const int num_simplex_vertices = mesh_dimension+1;

	for (int i=0;i<num_vertices*3;i++)
		new_vertex.push_back(vertex_coord[i]);

	//algorrithm to construct arrays E[] and I[] from T[]
	int numl=3*num_simplices;

	//construct a list of all triangle edges
	int m=0;
	for (int t=0;t<num_simplices;t++)
		for(int j=0;j<=2;j++)
			for(int k=0 ;k<=2;k++)
				if(j!=k)
					L.push_back(simplex_vert[3*t+ k]);




	for(int e=0;e<numl;e++)
	{
		if(L[2*e]>L[2*e+1])

		{
			L[2*e]=L[2*e]^L[2*e+1];
			L[2*e+1]=L[2*e]^L[2*e+1];
			L[2*e]=L[2*e]^L[2*e+1];
		}

	}

	vector<int> index_sorted;
	for (int i =0;i<numl;i++)
	{
		index_sorted.push_back(i);
	}

	sort(index_sorted.begin(),index_sorted.end(),myobject);


	vector <int> C(num_simplices, 0);
	vector <int> E;


	vector <int> I(0);
	I.resize(numl,0);

	int i=0,k=0,t,e,c;
	int num_edges=0;


	while (i < numl)
	{
		k=i+1;
		while (k<numl && edge_equals(index_sorted[i],index_sorted[k]) )
		{
			k++;
		}//endwhile
		e=index_sorted[i];
		E.push_back(L[2*e]);
		E.push_back(L[2*e+1]);
		//store reference to edge in incident array 
		for (int j=i;j<k;j++)
		{
			e=index_sorted[j];
			t=e/3;
			c=C[t];
			I[3*t+c]=num_edges;
			C[t]++;
		}

		num_edges++;
		i=k;
	}



	findCommonEdges(I,E,numl);




	if (!output_specified)
	{
		size_t found;
		found=input_filename.find_last_of(".");
		out_filename = input_filename.substr(0,found);
		out_filename += ".line";

	}

	ofstream output_file;

	output_file.open(out_filename.c_str(), ios::out);
	int numv = new_vertex.size()/3;
	int nume = edge_vert.size()/2;
	int color[4]={1, 0, 0, 1};

	//ijkoutLINE(output_file, dimension, new_vertex,numv, edge_vert, nume);

	//ijkoutColorLINE(output_file, dimension, &(new_vertex[0]), numv, &(edge_vert[0]), nume, color);
	// out = output stream
	// dim = dimension
	// coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
	// numv = number of vertices
	// edge_vert[dim*j+k] = k'th vertex index of edge j
	// nums = number of simplices
	// rgba[] = array of R,G,B,A values
	//

	ijkoutColorLINE (output_file, dimension, &(new_vertex[0]), numv, &(edge_vert[0]), nume, color);


	output_file.close();


};



//edge_equals(L,index_sorted[i],index_sorted[k])
bool edge_equals(int a, int b)
{
	if(    (L[2*a]==L[2*b])&& (L[2*a+1]==L[2*b+1]) )
		return true;
	else
		return false;
}


// function to calculate angle given two vertex indexes


// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

void parse_command_line(int argc, char **argv)
{
	if (argc!=3) {usage_error();}
	int iarg=1;
	while (iarg < argc && argv[iarg][0]=='-')
	{
		string s = argv[iarg];
		if (s=="-o")
			iarg++;
		output_specified = true;
		out_filename=argv[iarg];
		iarg++;
	}

	input_filename = argv[iarg+1];
	input_angle=atoi(argv[iarg]);
}

void usage_error()
{
	cerr << "Usage: findedge -o {outputfile name ending in .line} {angle} {input filename}" << endl;
	cerr << "help : angle and input_filename  are required input"<<endl;
	exit(10);
}

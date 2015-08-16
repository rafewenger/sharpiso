/// \file find sharp edges and corners from a triangle mesh 
/// Version: v0.1.0

/*
IJK: Isosurface Jeneration Code
Copyright (C) 2008-2015 Arindam Bhattacharya

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

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <time.h>
#include <stdlib.h>

#include "findsharp_eigen_info.h"


#include "ijk.txx"
#include "ijkIO.txx"

using namespace std;
using namespace IJK;

typedef float COORD_TYPE;
typedef float ANGLE_TYPE;

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
bool output_specified = false;
const std::string VERSION = "v0.1.0";

const ANGLE_TYPE DEFAULT_ANGLE(140);
ANGLE_TYPE input_angle(DEFAULT_ANGLE);


// output routine 
void output_mesh_info();
bool edge_equals( int a, int b);

//If manually set to true then in DEBUG MODE
bool debugMode = true;

//eigen information
EIGEN_INFO eigen_info;

// file names
string input_filename;
string out_base_fname;
// colors
int red[4]={1, 0, 0, 1};
int blue[4]={0, 0, 1, 0.5};
int yellow[4]={1, 0, 1, 1.0};
// store the edge info
enum EdgeType {SHARP, SMOOTH, DEGEN};
class edgeInfo{
public:
	vector<int> sharp;
	vector<int> smooth;
	vector<int> degen;
	string out_sharp;
	string out_smooth;
	string out_degen;
};
edgeInfo ei;

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
void help(), usage_error();

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

		if (debugMode)
		{
			cout <<" num_vertices " << num_vertices << " num_simplices " << num_simplices << endl;
		}

		// check mesh dimension
		if (mesh_dimension != 2)
		{
			cout <<"The input must be a triangle mesh."<<endl;
			exit (10);
		}

		output_mesh_info();
		cout <<"output files: " << endl;

		cout <<"\tsharp: " <<ei.out_sharp<<endl;
		cout <<"\tsmooth: " <<ei.out_smooth<<endl;
		cout <<"\tdegen: " <<ei.out_degen<<endl;

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

void addtovetor_edge(int vertEdge[], const EdgeType  e)
{

	if (e == SHARP){
		ei.sharp.push_back(vertEdge[0]);
		ei.sharp.push_back(vertEdge[1]);
	}
	else if (e == SMOOTH){
		ei.smooth.push_back(vertEdge[0]);
		ei.smooth.push_back(vertEdge[1]);
	}
	else if (e == DEGEN){
		ei.degen.push_back(vertEdge[0]);
		ei.degen.push_back(vertEdge[1]);
	}

}

//class defintiions

void otherVert(int t,int e, int edge[])
{
	point3D n1,n2;
	int ov1=0,ov2=0;


	for(int i=0;i<3;i++)
	{
		if(simplex_vert[3*t+i]!=edge[0] && simplex_vert[3*t+i]!=edge[1])
		{
			ov1=simplex_vert[3*t+i];
		}

		if(simplex_vert[3*e+i]!=edge[0] && simplex_vert[3*e+i]!=edge[1])
		{
			ov2=simplex_vert[3*e+i];
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
	double angle=0;

	if ((abs(norm1-0.0) < 0.0000001)||(abs(norm2-0.0) <  0.0000001))
	{
		addtovetor_edge(edge, DEGEN);
		eigen_info.num_degen_edges++;
		//cout <<"There are degenerate edges." << endl;
	}
	else if (eigen_info.flag_eigen_based)
	{
	
		if ((eigen_info.num_eigen[edge[0]] > 1 && eigen_info.num_eigen[edge[1]] > 1)&&
				(eigen_info.flag_centroid[edge[0]] == false && eigen_info.flag_centroid[edge[1]] == false))
		{
			angle = acos ( (  n1.x*n2.x +n1.y*n2.y+n1.z*n2.z) /(norm1*norm2)) * (180.0/M_PI);
			if(angle > (180-input_angle)){
				addtovector_simplex(edge);
				addtovetor_edge(edge, SHARP);
			}
			else
				addtovetor_edge(edge, SMOOTH);
		}
		else
			addtovetor_edge(edge, SMOOTH);
	}
	else{
		angle = acos ( (  n1.x*n2.x +n1.y*n2.y+n1.z*n2.z) /(norm1*norm2)) * (180.0/M_PI);
		if(angle > (180-input_angle)){
			addtovector_simplex(edge);
			addtovetor_edge(edge, SHARP);
		}
		else
		{
			addtovetor_edge(edge, SMOOTH);
		}
	}


}

void findCommonEdges(vector<int> I,vector<int> E,int numl)
{ 
	int e;
	vector<int > c1(numl,0);
	vector<int > Etri(numl,0);
	int op1,op2; //other vertex other than the common edge.
	int edge[2];
	point3D n1,n2;

	for (int t=0; t<num_simplices; t++)
	{
		for (int j=0;j<3;j++)
		{
			e=I[3*t+j];
			if(c1[e] == 0)
				Etri[e]=t;
			else
			{
				edge[0]=E[2*e];
				edge[1]=E[2*e+1];
				otherVert(t,Etri[e],edge);
			}
			c1[e]++;
		}
	}

}

// **************************************************
// OUTPUT_ROUTINE
// **************************************************

void output_mesh_info()
{
	if (debugMode){
		cout <<"Starting findSharp computations"<<endl;
		cout <<"size of eigen Info "<< eigen_info.num_eigen.size() <<endl;
	}
	const int num_simplex_vertices = mesh_dimension+1;

	for (int i=0; i<num_vertices*3; i++){
		new_vertex.push_back(vertex_coord[i]);
	}

	//algorithm to construct arrays E[] and I[] from T[]
	int numl=3*num_simplices;

	//construct a list of all triangle edges
	int m=0;
	for (int t=0;t<num_simplices;t++)
		for(int j=0;j<=2;j++)
			for(int k=0 ;k<=2;k++)
				if(j!=k)
					L.push_back(simplex_vert[3*t+ k]);



	for(int e=0; e<numl; e++)
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
		out_base_fname = input_filename.substr(0,found);
		ei.out_sharp = out_base_fname + ".line";
		ei.out_smooth = out_base_fname + ".smooth.line";
		ei.out_degen = out_base_fname + ".degen.line";

	}

	ofstream output_file;
	int numv = new_vertex.size()/3;
	int nume = ei.sharp.size()/2;

	output_file.open(ei.out_sharp.c_str(), ios::out);

	ijkoutColorLINE (output_file, dimension, &(new_vertex[0]), numv, &(ei.sharp[0]), nume, red);
	output_file.close();

	output_file.open(ei.out_smooth.c_str(), ios::out);
	nume = ei.smooth.size()/2;
	ijkoutColorLINE (output_file, dimension, &(new_vertex[0]), numv, &(ei.smooth[0]), nume, blue);
	output_file.close();

	output_file.open(ei.out_degen.c_str(), ios::out);
	nume = ei.degen.size()/2;
	if(nume > 0 )
	{ ijkoutColorLINE (output_file, dimension, &(new_vertex[0]), numv,
	&(ei.degen[0]), nume, yellow);}
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
	if (argc < 2) {usage_error();}

  if (argc == 2 && string(argv[1]) == "-version") {
    cout << "Version: " << VERSION << endl;
    exit(0);
  }

  input_angle = DEFAULT_ANGLE;

	int iarg = 1;
	while (iarg < argc && argv[iarg][0]=='-')
	{
		string s = argv[iarg];
    if (s == "-angle") {
      iarg++;
      if (iarg >= argc) { usage_error(); }
      input_angle = atoi(argv[iarg]);
    }
		else if (s == "-o")
		{
			iarg++;
			output_specified = true;
			out_base_fname=argv[iarg];
		}
		else if (s == "-debug")
		{
			debugMode = true;
		}
		else if (s == "-eigen_info")
		{
			iarg++;
      if (iarg >= argc) { usage_error(); }
			eigen_info.flag_eigen_based = true;
			eigen_info.eigen_info_filename = string(argv[iarg]);
			eigen_info.read_file();
		}
    else if (s == "-version")
      { cout << "Version: " << VERSION << endl; }
    else if (s == "-help")
      { help(); }
		else
		{
      cerr << "Error.  Illegal option: " << s << endl;
      cerr << endl;
			usage_error();
		}
		iarg++;
	}

	if (argc != iarg+1) {
		usage_error();
		exit(10);
	}

  input_filename = argv[iarg];

  cout <<" input angle " << input_angle << endl;
  cout <<" input file name "<< input_filename <<endl;
}

void usage_msg(std::ostream & out)
{
	out << "Usage: findsharp [OPTIONS] {input filename}" << endl;
}

void usage_error()
{
  usage_msg(cerr);
  cerr << "OPTIONS:" << endl;
  cerr << "  [-angle {A}] [-o {output_filename}] [-eigen_info {isov info file}]" << endl;
  cerr << "  [-help] [-version]" << endl;

	exit(10);
}

void help()
{
  usage_msg(cout);
  cout << endl;
  cout << "findsharp - Output sharp mesh edges." << endl;
  cout << "            Input is a geomview .off file." << endl;
  cout << "            Output geomview .line file of sharp edges." << endl;
  cout << "            Sharp edges have dihedral angle less than {angle}."
       << endl;
  cout << "            Input file can only contain triangles." << endl;
  cout << endl;

  cout << "OPTIONS:" << endl;
  cout << "  -angle {A}:  Dihedral angle." 
       << "  (Default: " << DEFAULT_ANGLE << ".)" << endl
       << "     Edges with dihedral angle at most {A} are sharp." << endl;
  cout << "  -o {output_filename}:  Output to file geomview .line file output_filename." << endl;
  cout << "  -eigen_info {isov info file}:" << endl;
  cout << "     Output only mesh edges incident on isosurface vertices"
       << endl
       << "     associated with 2 or 3 eigenvalues." << endl;
  cout << "     Read number of eigenvalues for each isosurface vertex" << endl
       << "     from {isov info file}.  File {isov info file} can be " << endl
       << "     produced using shrec with option -write_isov_info." << endl;
  cout << "  -version: Print version." << endl;
  cout << "  -help:    Print this help message." << endl;

  exit(0);
}


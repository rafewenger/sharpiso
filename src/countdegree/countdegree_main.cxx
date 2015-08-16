// Countdegree main

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

#include "countdegree_types.h"
#include "countdegree_IO.h"
#include "countdegree.h"

// global variables
char * input_filename(NULL);
string output_fn;
COORD_TYPE * vertex_coord(NULL);
VERTEX_INDEX * edge_endpoint(NULL);
int dimension(3);
int num_vertices(0);
int num_edges(0);

// decides which functions to call for output

bool  flag_op_to_file_short = false;
bool  flag_op_to_file_long = false;
bool flag_print_edges = false;

// miscellaneous routines
void usage_error();
void parse_command_line(int argc, char **argv);
void compute_output_fn ();

int real_degree_3_verts = 0;
int real_degree_1_verts = 0;

// main function
int main(int argc, char **argv)
{
	try {
		parse_command_line(argc, argv);
    
		ifstream in(input_filename, ios::in);
		compute_output_fn ();
		if (!in.good()) {
			cerr << "Unable to open file " << input_filename << "." << endl;
			exit(30);
		};
    
		ijkinLINE(in, dimension, vertex_coord, num_vertices,
              edge_endpoint, num_edges);
		in.close();
    
		// count the degree of each vertex
		vector <int> vert_degree(num_vertices,0);
    
		// count the degrees of the different edge points
		count_edge_degrees
		(dimension, vertex_coord, num_vertices, edge_endpoint, num_edges, vert_degree);
    
		if ( flag_op_to_file_short){
      //compute the output file name
			output_short_info(num_vertices, vert_degree, output_fn);
		}
    else if (flag_op_to_file_long)
      {
      output_vert_degree_2_file
      (num_vertices,vert_degree);
      }
		else{
			// output the edge information
			output_vert_degree( num_vertices, vert_degree, output_fn);
		}
    
		// print the edges degrees for analysis
		if (flag_print_edges){
			print_edge_info
			(num_vertices, vertex_coord,  vert_degree);
    }
    
	}
	catch (ERROR & error) {
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
// Miscellaneous routines
// **************************************************

void usage_msg()
{
	cerr << "Usage: countdegree <options> <*.>" << endl;
	cerr <<"\t-e prints the vertices which have deg 1 , 3 or >3"<<endl;
  cerr <<"\t-fshort: short ouput, the sum of deg1, deg3 or deg>3"<<endl;
  cerr <<"\t-flong: long output,degree 0, for degree 1, degree 3, degree > 3, degree 1 or 3 or more, total num non-o ver and total num vert."<<endl;
	//cerr <<" example : countdegree -deg3 32 cube.off "<<endl;
}

void usage_error()
{
	usage_msg();
	exit(10);
}
// parse command line
void parse_command_line(int argc, char **argv)
{
	if (argc == 1)  {usage_error();}
	int iarg=1;
	while (iarg<argc && argv[iarg][0]=='-')
    {
		string s = argv[iarg];
		if (s=="-deg3"){
			real_degree_3_verts=atoi(argv[++iarg]);
			iarg++;
		}
		else if (s=="-deg1"){
			real_degree_1_verts=atoi(argv[++iarg]);
			iarg++;
		}
		else if (s=="-fshort")
      {
			flag_op_to_file_short = true;
			iarg++;
      }
  	else if (s=="-flong")
      {
			flag_op_to_file_long = true;
			iarg++;
      }
		else if (s=="-e")
      {
			flag_print_edges=true;
			iarg++;
      }
    else if (s=="-help" || s=="=h")
      {
      usage_msg();
      iarg++;
      exit(0);
      }
		else
      {
			cout <<"There is no option called ["<<argv[iarg]<<"] here are the possible options."<<endl;
      usage_msg();
			iarg++;
			exit(0);
      }
    }
	input_filename = argv[iarg];
}




void compute_output_fn ()
{
	size_t found;
	string infile = input_filename;
	found=infile.find_last_of(".");
	output_fn = infile.substr(0,found);
}

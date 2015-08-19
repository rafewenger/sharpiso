// Countdegree main
// Version: 0.1.0

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
const char * VERSION("0.1.0");

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
			cerr << "Error.  Unable to open file " << input_filename << "." << endl;
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

void usage_msg(std::ostream & out)
{
	out << "Usage: countdegree [OPTIONS] <input filename>" << endl;
}

void usage_error()
{
	usage_msg(cerr);
	cerr << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  [-e | -fshort | -flong] [-deg1 <N>] [-deg3 <N>] [-help] [-version]" << endl;

	exit(10);
}

void help()
{
	usage_msg(cout);
	cout << endl;
	cout << "countdegree - Count vertex degrees in graph composed of line segments." << endl;
	cout << "              Input is a .line file (usually produced by findsharp.)"
		<< endl;
	cout  << "              Output is the number of vertices with each degree."
		<< endl;
	cout << "  -e:      Print vertices which have degrees other than zero or two."
		<< endl;
	cout << "  -fshort: Short output format." << endl
		<< "           Print the number of vertices with degree other than zero or two."
		<< endl;
	cout << "  -flong:  Long (condensed) output format." << endl
       << "           Print comma separated list of number of vertices with:" << endl
       << "             1) degree 0; 2) degree 1; 3) degree 2; 4) degree 3;"
       << endl
       << "             5) degree > 3; 6) degree not 0 or 2; 7) degree not 0;"
       << endl
       << "             8) total number of vertices." << endl;
  cout << "  -deg1 <N>:  Expected number of degree 1 vertices." << endl;
  cout << "    When set, countdegree reports difference between N and"
       << endl
       << "    number of degree 1 vertices." << endl;
  cout << "  -deg3 <N>:  Expected number of degree 3 vertices." << endl;
  cout << "    When set, countdegree reports difference between N and"
       << endl
       << "    number of degree 3 vertices." << endl;
	cout << "  -version: Print version." << endl;
	cout << "  -help:    Print this help message." << endl;
	exit(0);
}

float get_option_float
(const char * option, const char * value_string)
{
  float x;
  std::istringstream v_string;

  v_string.str(value_string);

  v_string >> x;

  if (!v_string.eof()) {
    cerr << "Error in argument for option: " << option << endl;
    cerr << "Non-numeric character in string: " << value_string << endl;

    exit(50);
  }

  return(x);
}

// parse command line
void parse_command_line(int argc, char **argv)
{
	if (argc == 1)  {usage_error();}

	if (argc == 2 && std::string(argv[1]) == "-version") {
		cout << "Version: " << VERSION << endl;
		exit(0);
	}

	int iarg=1;
	while (iarg<argc && argv[iarg][0]=='-') {
		string s = argv[iarg];
		if (s == "-deg3") {
			iarg++;
			if (iarg >= argc) { usage_error(); }
			real_degree_3_verts = 
        get_option_float(argv[iarg-1], argv[iarg]);
		}
		else if (s == "-deg1") {
			iarg++;
			if (iarg >= argc) { usage_error(); }
			real_degree_1_verts =
        get_option_float(argv[iarg-1], argv[iarg]);
		}
		else if (s=="-fshort") 
		{ flag_op_to_file_short = true; }
		else if (s=="-flong")
		{ flag_op_to_file_long = true; }
		else if (s=="-e") 
		{ flag_print_edges=true; }
		else if (s == "-version") 
		{ cout << "Version: " << VERSION << endl; }
		else if (s == "-help" || s=="=h") 
		{ help(); }
		else {
			cerr <<"Error.  Illegal option " << argv[iarg] << "." << endl;
			cerr << endl;
			usage_error();
			exit(10);
		}

		iarg++;
	}

	if (iarg >= argc) {
		cerr << "Error.  Missing input filename." << endl;
		cerr << endl;
		usage_error();
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

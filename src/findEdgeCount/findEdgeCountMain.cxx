// Test ijkIO.txx

#include <cstdlib>
#include <cmath>
#include <fstream>

#include "findEdgeCountTypes.h"
#include "findEdgeCountIO.h"
#include "findEdgeCount.h"

// global variables
char * input_filename(NULL);
COORD_TYPE * vertex_coord(NULL);
VERTEX_INDEX * edge_endpoint(NULL);
int dimension(3);
int num_vertices(0);
int num_edges(0);
// miscellaneous routines
void usage_error();
void parse_command_line(int argc, char **argv);

int main(int argc, char **argv)
{
	try {
		parse_command_line(argc, argv);

		ifstream in(input_filename, ios::in);
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
        count_edge_degrees(dimension, vertex_coord, num_vertices,
			edge_endpoint, num_edges, vert_degree);
        // output the edge information 
		output_vert_degree( num_vertices, vert_degree);

        // DEBUG 
		/*             
		output_edges(dimension, vertex_coord, num_vertices,
		edge_endpoint, num_edges);
		*/ 
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
	cerr << "Usage: findedgecount <.line file>" << endl;
}

void usage_error()
{
	usage_msg();
	exit(10);
}

void parse_command_line(int argc, char **argv)
{
	int iarg = 1;

	if (argc != 2) { usage_error(); }

	input_filename = argv[1];
}


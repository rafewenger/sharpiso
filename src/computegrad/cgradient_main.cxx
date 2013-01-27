/// \file cgradient_main.cxx
/// compute gradient vectors from scalar data
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

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


#include <iostream>
#include <cmath>


#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"

#include "cgradient.h"
#include "isodual3D_datastruct.h"
#include "cgradient_inputinfo.h"

using namespace IJK;
using namespace ISODUAL3D;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
bool report_time_flag = false;
bool flag_gzip = false;
bool flag_out_param = false;

using namespace std;

// local subroutines
void memory_exhaustion();
void usage_error();
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info);
void output_param (INPUT_INFO & io_info);
void print_unreliable_grad_info
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid, const INPUT_INFO & input_info);

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
	time_t start_time;
	time(&start_time);
	INPUT_INFO input_info;
	input_info.set_defaults();
	IJK::ERROR error;

	try {

		std::set_new_handler(memory_exhaustion);

		parse_command_line(argc, argv, input_info);
		output_param(input_info);

		ISODUAL_SCALAR_GRID full_scalar_grid;
		GRID_NRRD_IN<int,int> nrrd_in;
		NRRD_DATA<int,int> nrrd_header;

		nrrd_in.ReadScalarGrid
		(scalar_filename, full_scalar_grid, nrrd_header, error);
		if (nrrd_in.ReadFailed()) { throw error; }
		GRADIENT_GRID  vertex_gradient_grid;
		// compute central difference
		if (input_info.flag_cdiff){
			compute_gradient_central_difference
			(full_scalar_grid, vertex_gradient_grid, input_info);

		}
		else if (input_info.flag_reliable_grad){
			compute_reliable_gradients
			(full_scalar_grid, vertex_gradient_grid, input_info);
			// print info
			cout <<"Total number of vertices "<<full_scalar_grid.NumVertices() << endl;
			cout <<"Number of vertices with reliable grads " << input_info.out_info.num_reliable << endl;
			cout <<"Number of vertices with un-reliable grads " << input_info.out_info.num_unreliable << endl;
			cout <<"Number of boundary vertices "<<input_info.out_info.boundary_verts <<endl;
			//
			if (input_info.flag_print_grad_loc){
				print_unreliable_grad_info (full_scalar_grid, input_info);
			}
		}
		else if (input_info.flag_reliable_grad_far){
			compute_reliable_gradients_far
			(full_scalar_grid, vertex_gradient_grid, input_info);
			// print info
			cout <<"Total number of vertices "<<full_scalar_grid.NumVertices() << endl;
			cout <<"Number of vertices with reliable grads " << input_info.out_info.num_reliable << endl;
			cout <<"Number of vertices with un-reliable grads " << input_info.out_info.num_unreliable << endl;
			cout <<"Number of boundary vertices "<<input_info.out_info.boundary_verts <<endl;
		}
		else
		{
			cerr <<"No gradients were computed"<<endl;
		}

		if (flag_gzip) {
			write_vector_grid_nrrd_gzip(gradient_filename, vertex_gradient_grid);
		}
		else {
			write_vector_grid_nrrd(gradient_filename, vertex_gradient_grid);
			cerr << "file "<< gradient_filename << " created."<<endl;
		}

		if (report_time_flag) {

			time_t end_time;
			time(&end_time);
			double total_elapsed_time = difftime(end_time, start_time);
			cout << endl;
			cout << "Elapsed time = " << total_elapsed_time << " seconds."
					<< endl;
		};

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
	};

}


// **************************************************
// MISC ROUTINES
// **************************************************

void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

/// print out the parameter informations
void output_param (INPUT_INFO & io_info){
	/*
	 * bool flag_cdiff;	 		// compute central difference
	bool flag_reliable_grad;  // reliable grad
	bool print_info;          // print info of the vertex
	bool flag_print_grad_loc;      // prints the location of the unreliable grads
	bool flag_reliable_grad_far; // compare the cdiff gradient with immediate neighbors or
								 // neighbors at a certain distance
	int print_info_vertex;
	float param_angle;
	float min_gradient_mag;   // minimum gradient
	float min_cos_of_angle;
	 */
	if (flag_out_param){
		cout <<" OUTPARAM \n ";
		if (io_info.flag_cdiff)
		{
			cout <<" flag_cdiff, ";
		}
		if (io_info.flag_reliable_grad)
		{
			cout <<" reliable_grad, ";
		}
		if (io_info.flag_reliable_grad_far)
		{
			cout <<" reliable grad far,";
			cout <<" reliable_grad_dist : "<< io_info.reliable_grad_far_dist <<endl;
		}
		cout << "min cos angle    " << (acos(io_info.min_cos_of_angle)*180.0)/M_PI << endl;
		cout << "min num agree    " << io_info.min_num_agree << endl;
		cout << "min_gradient_mag " << io_info.min_gradient_mag << endl;
	}
}
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info)
{
	int iarg = 1;

	while (iarg < argc && argv[iarg][0] == '-') {
		if (string(argv[iarg]) == "-time")
		{ report_time_flag = true;   }
		else if (string(argv[iarg]) == "-gzip")
		{ flag_gzip = true; }
		else if (string(argv[iarg])== "-cdiff")
		{
			io_info.flag_cdiff = true; // compute central difference
		}
		else if (string(argv[iarg])== "-reliable_grad")
		{
			io_info.flag_reliable_grad = true; // compute reliable gradient
		}
		else if (string(argv[iarg])== "-out_param")
		{
			flag_out_param = true;
		}
		else if (string(argv[iarg])== "-min_gradient_mag")
		{
			iarg++;
			io_info.min_gradient_mag= (float)atof(argv[iarg]);
		}
		else if (string(argv[iarg])== "-angle")
		{
			iarg++;
			io_info.param_angle= (float)atof(argv[iarg]);
			io_info.min_cos_of_angle = cos((io_info.param_angle*M_PI/180.0));
		}
		else if (string(argv[iarg])== "-min_num_agree")
		{
			iarg++;
			io_info.min_num_agree=(int)atoi(argv[iarg]);
		}
		else if (string(argv[iarg])== "-print_info")
		{
			iarg++;
			io_info.print_info = true;
			io_info.print_info_vertex= (int)atof(argv[iarg]);
		}
		else if (string(argv[iarg])== "-print_grad_loc")
		{
			io_info.flag_print_grad_loc = true;
		}
		else if (string (argv[iarg])== "-reliable_grad_far")
		{
			io_info.flag_reliable_grad_far = true;
		}
		else if (string(argv[iarg])=="-reliable_grad_far_dist"){
			iarg++;
			io_info.flag_reliable_grad_far = true;
			io_info.reliable_grad_far_dist=(int)atof(argv[iarg]);
		}
		else
		{ usage_error(); }
		iarg++;
	}

	if (iarg+2 != argc) { usage_error(); };

	scalar_filename = argv[iarg];
	gradient_filename = argv[iarg+1];
}

void usage_msg()
{
	cerr << "Usage: computegrad [options] {scalar nrrd file} {gradient nrrd file}"<<endl;
	cerr << "Options \n \t\t-cdiff : compute the central difference {default}"<<endl;
	cerr << "\t\t-reliable_grad: reliable gradient."<<endl;
	cerr << "\t\t-out_param"<< endl;
	cerr << "\t\t-min_gradient_mag [float] : min gradient magnitude." <<endl;
	cerr << "\t\t-print_info [int] : print_info of the vertex." <<endl;
	cerr << "\t\t-print_grad_loc : prints the location of all the vertices with unreliable grads" <<endl;
	cerr << "\t\t-angle : angle" <<endl;
	cerr << "\t\t-min_num_agree: default set to 4" <<endl;
}

void usage_error()
{
	usage_msg();
	exit(100);
}

/// print information of the unreliable gradients
void print_unreliable_grad_info
(const ISODUAL3D::ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const INPUT_INFO & input_info){
	COORD_TYPE coord0[3];
	for (int i=0; i<input_info.out_info.un_reliable_grads_vert_info.size();i++){
		VERTEX_INDEX v = input_info.out_info.un_reliable_grads_vert_info[i];
		scalar_grid.ComputeCoord(v,coord0);
		cout <<" iv " << v;
		cout <<" ["<<coord0[0]<<" "<<coord0[1]<<" "<<coord0[2]<<"]"<<endl;
	}
}

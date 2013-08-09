///Compute reliable gradients directly from the scalar data



#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"
#include "religrad_computations.h"
#include "religrad_inputIO.h"

using namespace IJK;
using namespace std;
using namespace RELIGRADIENT;

// global variables
char * scalar_filename = NULL;
char * gradient_filename = NULL;
bool report_time_flag = false;
bool flag_gzip = false;
bool flag_out_param = false;


// local subroutines
void memory_exhaustion();
void usage_error();
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info);
void output_param (INPUT_INFO & io_info);


timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

using namespace SHARPISO;

// Main program

int main(int argc, char **argv)
{

	time_t start_time;
	time(&start_time);
	INPUT_INFO input_info;
	input_info.set_defaults();
	IJK::ERROR error;
	bool OptChosen = false;
	try {

		std::set_new_handler(memory_exhaustion);
		//parse input
		parse_command_line(argc, argv, input_info);
		//print the input params
		output_param(input_info);

		
		RELIGRADIENT_SCALAR_GRID full_scalar_grid;
		GRID_NRRD_IN<int,int> nrrd_in;
		NRRD_DATA<int,int> nrrd_header;
		SCALAR_TYPE zero_vector [3] = { 0.0, 0.0, 0.0};

		nrrd_in.ReadScalarGrid
		(scalar_filename, full_scalar_grid, nrrd_header, error);
		if (nrrd_in.ReadFailed()) { throw error; }
		GRADIENT_GRID  vertex_gradient_grid;

    std::vector<COORD_TYPE> grid_spacing;
    nrrd_header.GetSpacing(grid_spacing);

    for (int d = 0; d < full_scalar_grid.Dimension(); d++) {
      full_scalar_grid.SetSpacing(d, grid_spacing[d]); 
      vertex_gradient_grid.SetSpacing(d, grid_spacing[d]);
    }

    // *** DEBUG ***
    using namespace std;
    cout << "Grid Spacing: ";
    for (int d = 0; d < full_scalar_grid.Dimension(); d++) {
      cout << " " << grid_spacing[d]; 
    }
    cout << endl;

    IJK::BOOL_GRID<RELIGRADIENT_GRID> reliable_grid;
		reliable_grid.SetSize(full_scalar_grid);
		reliable_grid.SetAll(true);

		// compute central difference
		if (input_info.flag_cdiff){
		
		compute_gradient_central_difference
			(full_scalar_grid, vertex_gradient_grid, input_info);
			OptChosen = true;
		}
  
		//check neighboring gradients
		//angle based
		if (input_info.angle_based)
		{
			cerr << "perform angle based reliability test.\n"
					<< "\tparameters \n\t[angle_based_dist] " << input_info.angle_based_dist
					<< "\t[min_num_agree] "<<input_info.min_num_agree
					<< "\t[min_cos_of_angle] "<<(acos(input_info.min_cos_of_angle)*180.0)/M_PI<<endl;

			time_t begin,end;
			time (&begin);
    
			compute_reliable_gradients_angle
			(full_scalar_grid, vertex_gradient_grid, reliable_grid, input_info);
			OptChosen = true;
			time (&end);

			cout << "num unreliable [" << input_info.out_info.num_unreliable
					<< "]\t num reliable [" << input_info.out_info.num_reliable
					<<"]\t grad mag > min ["<<input_info.out_info.grad_mag_zero
					<< "]\t total [" << full_scalar_grid.NumVertices() << "]"
					<< endl;
			//cout <<"time " <<difftime(end,begin) <<endl;
		
		}


    
		if (input_info.flag_reliable_scalar_prediction)
		{
			input_info.out_info.num_unreliable = 0;
			time_t begin,end;
			time (&begin);

			compute_reliable_gradients_SBP
			(full_scalar_grid, vertex_gradient_grid, reliable_grid, input_info);
			OptChosen = true;
			time (&end);

			cout <<"\tparameters\n\t[scalar_prediction_dist] "<<input_info.scalar_prediction_dist
					<<"\t[scalar_prediction_err] "<<input_info.scalar_prediction_err
					<<endl;
			cout <<"Num unreliable" << input_info.out_info.num_unreliable << endl;

		}
		
		if (!OptChosen)
		{
			cerr <<"No gradients were computed"<<endl;
			exit(0);
		}
		// set the correct gradients
		// set gradients
		int num_unreliable = 0;
		for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++)
		{
			if (!reliable_grid.Scalar(iv)){
				vertex_gradient_grid.Set(iv, zero_vector);
				num_unreliable++;
			}
		}

		if (flag_gzip)
		{
			write_vector_grid_nrrd_gzip(gradient_filename, vertex_gradient_grid);
		}
		else {
			write_vector_grid_nrrd(gradient_filename, vertex_gradient_grid);
			cerr << "file "<< gradient_filename << " created."<<endl;
		}

		if (report_time_flag)
		{

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
	if (flag_out_param){
		cout <<" OUTPARAM \n ";
		if (io_info.flag_cdiff)
		{
			cerr <<"[flag_cdiff] ";
		}
		if (io_info.flag_reliable_grad)
		{
			cout <<"reliable_grad, ";
		}
		if (io_info.angle_based)
		{
			cout << "Reliable grad far"<<endl;
			cout << "reliable_grad_dist : "<< io_info.angle_based_dist <<endl;
			cout << "min cos angle    " << (acos(io_info.min_cos_of_angle)*180.0)/M_PI << endl;
			cout << "min num agree    " << io_info.min_num_agree << endl;
			cout << "min_gradient_mag " << io_info.min_gradient_mag << endl;
			cout <<"\n";
		}
		if (io_info.flag_reliable_scalar_prediction)
		{
			cout <<"Scalar based prediction"<<endl;
			cout <<"scalar error tolerance: "<<io_info.scalar_prediction_err<<endl;
			cout <<"scalar error distance to neighbors: "<<io_info.scalar_prediction_dist<<endl;
		}

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
		else if (string(argv[iarg])== "-draw_vert")
		{
			iarg++;
			io_info.draw = true;
			io_info.draw_vert= atoi(argv[iarg]);

		}
		else if (string(argv[iarg])== "-print_grad_loc")
		{
			io_info.flag_print_grad_loc = true;
		}

		else if (string (argv[iarg])== "-angle_based")
		{
			io_info.angle_based = true;
		}

		else if (string(argv[iarg])=="-angle_based_dist"){
			iarg++;
			io_info.angle_based = true;
			io_info.angle_based_dist=(int)atof(argv[iarg]);
		}
		else if (string(argv[iarg])=="-reliable_scalar_pred_dist"){
			iarg++;
			io_info.flag_reliable_scalar_prediction = true;
			io_info.scalar_prediction_dist = atoi(argv[iarg]);
		}
		else if (string(argv[iarg])=="-scalar_pred_err"){
			iarg++;
			io_info.flag_reliable_scalar_prediction = true;
			io_info.scalar_prediction_err = atof(argv[iarg]);
		}
		else
		{
			cout <<"Error in  "<< string(argv[iarg]) << endl;
			usage_error();
		}
		iarg++;
	}

	if (iarg+2 != argc) { usage_error(); };

	scalar_filename = argv[iarg];
	gradient_filename = argv[iarg+1];
}

void usage_msg()
{
	cerr << "Usage: religrad [options] {scalar nrrd file} {gradient nrrd file}"<<endl;
	cerr << "Options \n \t\t-cdiff : compute the central difference {default}"<<endl;
	cerr << "\t\t-out_param"<< endl;
	cerr << "\t\t-min_gradient_mag [float] : min gradient magnitude." <<endl;
	cerr << "\t\t-print_info [int] : print_info of the vertex." <<endl;
	cerr << "\t\t-print_grad_loc : prints the location of all the vertices with unreliable grads" <<endl;
	cerr << "\t\t-angle_based : angle based reliable grad computations" <<endl;
	cerr << "\t\t\t-angle : angle for above" <<endl;
	cerr << "\t\t\t-min_num_agree: default set to 4" <<endl;
	cerr << "\t\t\t-angle_based_dist: <how far to look, e.g. 2>" <<endl;
	cerr << "\t\t-reliable_scalar_pred_dist: <how far to look, e.g. 2>" <<endl;
	cerr << "\t\t\t-scalar_pred_err: <threhold for error in prediction default 0.15>" <<endl;

}

void usage_error()
{
	usage_msg();
	exit(100);
}


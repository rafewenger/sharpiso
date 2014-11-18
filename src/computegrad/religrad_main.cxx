///Compute reliable gradients directly from the scalar data

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <time.h>


#include "ijkNrrd.h"
#include "ijkgrid_nrrd.txx"
#include "ijkprint.txx"

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
void usage_error(), help_msg();
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info);
void output_param(INPUT_INFO & io_info);
/*
timespec diff(timespec start, timespec end) {
timespec temp;
if ((end.tv_nsec - start.tv_nsec) < 0) {
temp.tv_sec = end.tv_sec - start.tv_sec - 1;
temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
} else {
temp.tv_sec = end.tv_sec - start.tv_sec;
temp.tv_nsec = end.tv_nsec - start.tv_nsec;
}
return temp;
}
*/

using namespace SHARPISO;

/*
* Output info
*/
void out_stats(const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
			   INPUT_INFO io_info)
{
	OUTPUT_INFO out = io_info.out_info;
	cout <<"********************************\n"
		<<"\tOutput stats"
		<<"\n********************************\n";
	cout <<"Num unreliable "<< out.num_unreliable 
		<<" [" << (out.num_unreliable * 1.0 / scalar_grid.NumVertices())*100
		<<"%] of total vertices.\n"
		<<"Num reliable "<< out.num_reliable << endl;
	
}

// Main program

int main(int argc, char **argv) {
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
		GRID_NRRD_IN<int, int> nrrd_in;
		NRRD_DATA<int, int> nrrd_header;
		GRADIENT_COORD_TYPE zero_vector[3] = { 0.0, 0.0, 0.0 };

		nrrd_in.ReadScalarGrid(scalar_filename, full_scalar_grid, nrrd_header,
			error);
		if (nrrd_in.ReadFailed()) {
			throw error;
		}
		GRADIENT_GRID vertex_gradient_grid;
		//store the magnitudes
		GRADIENT_MAGNITUDE_GRID magnitude_grid;

		const int dimension = full_scalar_grid.Dimension();

		std::vector<COORD_TYPE> grid_spacing;
		nrrd_header.GetSpacing(grid_spacing);

		full_scalar_grid.SetSpacing(&(grid_spacing[0]));

		IJK::BOOL_GRID<RELIGRADIENT_GRID> reliable_grid;
		reliable_grid.SetSize(full_scalar_grid);
		reliable_grid.SetAll(true);

		// compute central difference
		if (input_info.flag_cdiff) {
			compute_gradient_central_difference_normalized(full_scalar_grid,
				vertex_gradient_grid, magnitude_grid, input_info);
			OptChosen = true;
		}

		//check neighboring gradients
		//angle based
		if (input_info.angle_based) {
			cerr << "\nPerform Angle Based Reliability Test.\n"
				<< "PARAMETERS: \nAngle_based_dist: "
				<< input_info.angle_based_dist << ",\nMin_Num_Agree: "
				<< input_info.min_num_agree << ",\nMin_Cos_of_Angle: "
				<< (acos(input_info.min_cos_of_angle) * 180.0) / M_PI
				<< endl;

			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);

			compute_reliable_gradients_angle
				(full_scalar_grid, vertex_gradient_grid, magnitude_grid, reliable_grid,
				input_info);
			OptChosen = true;
			time(&end);
			finish = clock();

			cout << "num unreliable(only from angle computations): " << input_info.out_info.num_unreliable
				<< ", num reliable: " << input_info.out_info.num_reliable
				<< ",\n(grad mag < min): "
				<< input_info.out_info.grad_mag_zero << endl;
			cout << "Time estimates angle based: " << difftime(end, begin)
				<< ","<< (double(finish - start) / CLOCKS_PER_SEC) << endl;

		}

		if (input_info.flag_reliable_scalar_prediction) {
			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);
			
			compute_reliable_gradients_SBP
				(full_scalar_grid, vertex_gradient_grid, magnitude_grid, 
				reliable_grid, input_info);

			OptChosen = true;
			time(&end);
			finish = clock();

			cout << "\nPARAMETERS for  Scalar based predictions:\nScalar_prediction_dist: "
				<< input_info.scalar_prediction_dist
				<< ", Scalar_prediction_err: "
				<< input_info.scalar_prediction_err << endl;
			cout << "Num unreliable(only from scalar_test computations): " << input_info.out_info.num_unreliable;
			cout << " Num reliable: " << input_info.out_info.num_reliable
				<< endl;
			cout << "\nTime estimates scalar based : " << difftime(end, begin)
				<< ","<< (double(finish - start) / CLOCKS_PER_SEC) << endl;
			cout<<"\nTotal number of vertices: "
				<< full_scalar_grid.NumVertices() << endl;

		}
		if(input_info.adv_angle_based){
			//reset num unreliables. 
			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);

			compute_reliable_gradients_advangle
			(full_scalar_grid, vertex_gradient_grid, magnitude_grid, 
			reliable_grid, input_info);
			out_stats(full_scalar_grid, input_info);
			
			OptChosen = true;
			time(&end);
			finish = clock();			
		}
		if(input_info.adv_angle_based_v2){
			//reset num unreliables. 
			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);

			/*compute_reliable_gradients_advangle
			(full_scalar_grid, vertex_gradient_grid, magnitude_grid, 
			reliable_grid, input_info);
			out_stats(full_scalar_grid, input_info);
			*/


			compute_reliable_gradients_advangle_version2
				(full_scalar_grid, vertex_gradient_grid, magnitude_grid, 
				reliable_grid, input_info);
			out_stats(full_scalar_grid, input_info);


			OptChosen = true;
			time(&end);
			finish = clock();			
		}

		if (!OptChosen) {
			cerr << "No gradients were computed" << endl;
			exit(0);
		}
		// set the correct gradients
		// set gradients
		int num_unreliable = 0;
		for (VERTEX_INDEX iv = 0; iv < full_scalar_grid.NumVertices(); iv++) 
		{
			if (!reliable_grid.Scalar(iv)) {
				vertex_gradient_grid.Set(iv, zero_vector);
				num_unreliable++;
			}
			else
			{
				GRADIENT_COORD_TYPE grad_temp[DIM3] = { 0.0, 0.0, 0.0 };
				std::copy(vertex_gradient_grid.VectorPtrConst(iv),
					vertex_gradient_grid.VectorPtrConst(iv) + DIM3, grad_temp);
				
				GRADIENT_COORD_TYPE mag = magnitude_grid.Scalar(iv);
				
				for (int l = 0; l < DIM3; l++)
					{
						grad_temp[l] = grad_temp[l]*mag;
				}
				//set the gradient
				vertex_gradient_grid.Set(iv,grad_temp);
			}
		}

		NRRD_DATA<int, int> gradient_nrrd_header;
		std::vector<AXIS_SIZE_TYPE> gradient_nrrd_axis_size(dimension + 1, 1);
		std::vector<double> gradient_nrrd_spacing(dimension + 1);
		gradient_nrrd_axis_size[0] = vertex_gradient_grid.VectorLength();
		gradient_nrrd_spacing[0] = 1;
		for (int d = 1; d <= dimension; d++) {
			gradient_nrrd_axis_size[d] = vertex_gradient_grid.AxisSize(d - 1);
			gradient_nrrd_spacing[d] = vertex_gradient_grid.Spacing(d - 1);
		}

		gradient_nrrd_header.SetSize(dimension + 1,
			&(gradient_nrrd_axis_size[0]));
		nrrdAxisInfoSet_nva(gradient_nrrd_header.DataPtr(), nrrdAxisInfoSpacing,
			&(gradient_nrrd_spacing[0]));

		if (flag_gzip) {
			write_vector_grid_nrrd_gzip(gradient_filename, vertex_gradient_grid,
				gradient_nrrd_header);
		} else {
			write_vector_grid_nrrd(gradient_filename, vertex_gradient_grid,
				gradient_nrrd_header);
		}
		cerr << "file " << gradient_filename << " created." << endl;

		if (report_time_flag) {

			time_t end_time;
			time(&end_time);
			double total_elapsed_time = difftime(end_time, start_time);
			cout << endl;
			cout << "Elapsed time = " << total_elapsed_time << " seconds."
				<< endl;
		};
	} catch (ERROR error) {
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		} else {
			error.Print(cerr);
		}
		cerr << "Exiting." << endl;
		exit(20);
	} catch (...) {
		cerr << "Unknown error." << endl;
		exit(50);
	};

}

// **************************************************
// MISC ROUTINES
// **************************************************

void memory_exhaustion() {
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

/// print out the parameter informations
void output_param(INPUT_INFO & io_info) {
	if (flag_out_param) {
		cout << "*********************************\n";
		cout <<"\tOut Parameters.\n*********************************\n";
		if (io_info.flag_cdiff) {
			cerr << "Central Difference for computing gradients.\n";
		}
		if (io_info.flag_reliable_grad) {
			cout << "reliable_grad, ";
		}
		if (io_info.angle_based) {
			cout << "Reliable grad far" << endl;
			cout << "reliable_grad_dist : " << io_info.angle_based_dist << endl;
			cout << "min cos angle    "
				<< (acos(io_info.min_cos_of_angle) * 180.0) / M_PI << endl;
			cout << "min num agree    " << io_info.min_num_agree << endl;
			cout << "min_gradient_mag " << io_info.min_gradient_mag << endl;
			cout << "\n";
		}
		if (io_info.flag_reliable_scalar_prediction) {
			cout << "Scalar based prediction" << endl;
			cout << "scalar error tolerance: " << io_info.scalar_prediction_err
				<< endl;
			cout << "scalar error distance to neighbors: "
				<< io_info.scalar_prediction_dist << endl;
		}
		if(io_info.adv_angle_based)
		{
			cout << "*********************************\n";
			cout <<"Advanced Angle based Reliability criteria."<<endl;
			cout << "*********************************\n";
			cout <<"Angle threshold for gradients to agree "<<io_info.param_angle << endl;
			cout <<"Neighbor angle "<<io_info.neighbor_angle_parameter<<endl;
			cout <<"\t{angle between grad at v and vector vv' where v' is an edge neighbor}"<<endl;
			
		}

	}
}
//Parse command line
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info) {
	int iarg = 1;

	while (iarg < argc && argv[iarg][0] == '-') {

		string s = string(argv[iarg]);

		if (s == "-time") {
			report_time_flag = true;
		} 
		else if (s == "-cdiff") {
			io_info.flag_cdiff = true; // compute central difference
		}
		else if (s == "-angle_test") {
			io_info.flag_cdiff = true;
			io_info.angle_based = true;
		}
		else if (s == "-scalar_test") {
			io_info.flag_cdiff = true;
			io_info.flag_reliable_scalar_prediction = true;
		}
		else if (s == "-out_param") {
			flag_out_param = true;
		} 
		else if (s == "-min_gradient_mag") {
			iarg++;
			io_info.min_gradient_mag = (float) atof(argv[iarg]);
		} 
		else if (s == "-angle") {
			iarg++;
			io_info.param_angle = (float) atof(argv[iarg]);
			io_info.min_cos_of_angle = cos(
				(io_info.param_angle * M_PI / 180.0));
		} 
		else if (s == "-min_num_agree") {
			iarg++;
			io_info.min_num_agree = (int) atoi(argv[iarg]);
		} 
		else if (s == "-print_info") {
			iarg++;
			io_info.print_info = true;
			io_info.print_info_vertex = (int) atof(argv[iarg]);
		} 
		else if (s == "-draw_vert") {
			iarg++;
			io_info.draw = true;
			io_info.draw_vert = atoi(argv[iarg]);
		} 
		else if (s == "-print_grad_loc") {
			io_info.flag_print_grad_loc = true;
		}
		else if (s == "-angle_based_dist") {
			iarg++;
			io_info.angle_based = true;
			io_info.angle_based_dist = (int) atof(argv[iarg]);
		}
		else if (s == "-reliable_scalar_pred_dist") {
			iarg++;
			io_info.flag_reliable_scalar_prediction = true;
			io_info.scalar_prediction_dist = atoi(argv[iarg]);
		} 
		else if (s == "-scalar_pred_err") 
		{
			iarg++;
			io_info.flag_reliable_scalar_prediction = true;
			io_info.scalar_prediction_err = atof(argv[iarg]);
		} 
		else if(s == "-advangle")
		{
			io_info.flag_cdiff = true;
			io_info.adv_angle_based = true;
		}
		else if(s == "-advanglev2")
		{
			io_info.flag_cdiff = true;
			io_info.adv_angle_based_v2 = true;
		}
		else if (s == "-neighbor_angle")
		{
			iarg++;
			io_info.flag_cdiff = true;
			io_info.adv_angle_based = true;
			io_info.neighbor_angle_parameter = atof(argv[iarg]);
		}
		else if (s == "-gzip") {
			flag_gzip = true;
		} 
		else if (s == "-help") {
			help_msg();
		}
		else {
			cout << "Error in  " << s << endl;
			usage_error();
		}
		iarg++;
	}

	if (iarg + 2 != argc) {
		usage_error();
	};

	scalar_filename = argv[iarg];
	gradient_filename = argv[iarg + 1];
}

void usage_msg() {
	cerr << "Usage: religrad [OPTIONS] {scalar nrrd file} {gradient nrrd file}"
		<< endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  [-cdiff] [-angle_test] [-scalar_test] [-advangle]" << endl;
	cerr << "  [min_gradient_mag {M}] [-angle {A}] [-min_num_agree {N}]" << endl;
	cerr << "  [-angle_based_dist {D}] [-reliable_scalar_pred_dist {D}]" << endl;
	cerr << "  [-neighbor_angle {A}]"<< endl;
	cerr << "  [-scalar_pred_err {E}]" << endl;
	cerr << "  [-gzip]" << endl;
	cerr << "  [-out_param] [-print_info {V}] [-print_grad_loc] [-help]" << endl;
}

void help_msg() {
	cerr << "Usage: religrad [OPTIONS] {scalar nrrd file} {gradient nrrd file}"
		<< endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -cdiff: Compute the central difference (default)." << endl;
	cerr << "  -angle_test: Apply angle test." << endl;
	cerr << "  -scalar_test: Apply scalar test." << endl;
	cerr << "  -min_gradient_mag {M}:  Set min gradient magnitude to {M} (float)." << endl;
	cerr << "  -angle {A}: Set angle to {A} (float)." << endl;
	cerr << "  -min_num_agree {N}: Set minimum agree number to {N}." << endl;
	cerr << "     {N} gradient directions must agree to pass the angle test."
		<< endl;
	cerr << "     (Default 4.)" << endl;
	cerr << "  -angle_based_dist {D}: Distance (integer) to neighboring vertices"
		<< endl
		<< "     in angle test.  (Default 1.)" << endl;
	cerr << "  -reliable_scalar_pred_dist: Distance (integer) to neighboring vertices" << endl
		<< "     in scalar test.  (Default 2.)" << endl;
	cerr <<"\t-neighbor_angle {A} Angle between grad at v and vector vv'\n\t"
		<<"Where v' is an edge neighbor. Default is 30." << endl;
	cerr << "  -scalar_pred_err {E}:  Error threshold for scalar test." 
		<< endl;
	cerr << "     Errors above the threshold fail the test. (Default 0.4.)" 
		<< endl;
	cerr << "  -gzip: Store gradients in compressed (gzip) format." << endl;
	cerr << "  -out_param:  Print parameters." << endl;
	cerr << "  -print_info {V} : Print information about vertex {IV}." << endl;
	cerr << "  -print_grad_loc : Print location of vertices with unreliable gradients." << endl;
	cerr << "  -help: Print this help message." << endl;

	exit(0);
}

void usage_error() {
	usage_msg();
	exit(100);
}


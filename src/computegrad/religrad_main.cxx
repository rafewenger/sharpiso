/// Compute reliable gradients directly from the scalar data
/// Version: 0.1.0

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
const char * VERSION = "0.1.0";

// local subroutines
void memory_exhaustion();
void usage_msg(const bool flag_list_all_options);
void usage_error(const bool flag_list_all_options);
void help(const bool flag_list_all_options);
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info);
void output_param(INPUT_INFO & io_info);


using namespace SHARPISO;

/*
* Output info
* If bool fio is true, then write single line output.
*/
void out_stats
	(const string method, const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO io_info, bool & fio)
{

	OUTPUT_INFO out = io_info.out_info;
		int num_unreliable=0, num_reliable=0;
	for (int i = 0; i < reliable_grid.NumVertices(); i++)
	{	
		if ( !reliable_grid.Scalar(i))
		{
			num_unreliable++;
		}
		else
		{
			num_reliable++;
		}
	}
	if (!fio){
	cout <<"********************************\n"
		<<"\tOutput stats ["<< method<<"]"
		<<"\n********************************\n";

	cout <<"TOTAL number of unreliable vertices: "<< num_unreliable << endl;
	cout <<"TOTAL number of reliable vertices: "<< num_reliable << endl;
	cout <<"TOTAL number of vertices: "<< num_unreliable + num_reliable<< endl;
	}
	else{
		cout << num_unreliable*100.0/reliable_grid.NumVertices()<<","
			<< num_reliable*100.0/reliable_grid.NumVertices()<< endl;}

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
		reliable_grid.SetAll(false);
		//set to false if any other option is chosen 
		bool only_cdiff = true;

		// compute central difference
		if (input_info.flag_cdiff) {
			
			compute_gradient_central_difference_normalized(full_scalar_grid,
				vertex_gradient_grid, magnitude_grid, input_info);
			OptChosen = true;
		}

		//check neighboring gradients
		//angle based
		if (input_info.angle_based) {
			only_cdiff = false;
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
			
			only_cdiff = false;

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
			only_cdiff = false;

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
			bool fileIO = true;
			out_stats("adv_angle_based", full_scalar_grid, reliable_grid, input_info, fileIO);

			OptChosen = true;
			time(&end);
			finish = clock();			
		}
		if(input_info.adv_angle_based_v2){
			
			only_cdiff = false;

			//reset num unreliables. 
			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);

			compute_reliable_gradients_advangle_version2
				(full_scalar_grid, vertex_gradient_grid, magnitude_grid, 
				reliable_grid, input_info);
			bool fileIO = true;
			out_stats("adv_angle_based_v2", full_scalar_grid, reliable_grid, input_info, fileIO);


			OptChosen = true;
			time(&end);
			finish = clock();			
		}

		if (input_info.algo12)
		{
			only_cdiff = false; 
			//create and fill boundary grid
			RELIGRADIENT::BOOL_GRID boundary_grid;
			boundary_grid.SetSize(full_scalar_grid);
			compute_boundary_grid(boundary_grid);

			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);	
			compute_reliable_gradients_curvature_based_algo12
				(full_scalar_grid, boundary_grid, vertex_gradient_grid, magnitude_grid, 
				reliable_grid, input_info);
			bool fileIO = true;
			out_stats("Algorithm 1 or 2", full_scalar_grid, reliable_grid, input_info, fileIO);

			OptChosen = true;
			time(&end);
			finish = clock();	

		}

		if(input_info.curv_based)
		{
			only_cdiff = false;

			//boundary grid
			//create and fill boundary grid
			RELIGRADIENT::BOOL_GRID boundary_grid;
			boundary_grid.SetSize(full_scalar_grid);
			compute_boundary_grid(boundary_grid);

			input_info.out_info.num_unreliable = 0;
			input_info.out_info.num_reliable = 0;
			time_t begin, end;
			clock_t start, finish;
			start = clock();
			time(&begin);	
			//curvature based computation version B 
			compute_reliable_gradients_curvature_basedB
				(full_scalar_grid, boundary_grid, vertex_gradient_grid, magnitude_grid, 
				reliable_grid, input_info);
			bool fileIO = true;
			out_stats("curvature based", full_scalar_grid, reliable_grid, input_info, fileIO);

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
			if (only_cdiff)
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
			else{
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

		if (io_info.flag_reliable_grad) {
			cout << "reliable_grad, ";
		}
		else if (io_info.angle_based) {
			cout << "Reliable grad far" << endl;
			cout << "reliable_grad_dist : " << io_info.angle_based_dist << endl;
			cout << "min cos angle    "
				<< (acos(io_info.min_cos_of_angle) * 180.0) / M_PI << endl;
			cout << "min num agree    " << io_info.min_num_agree << endl;
			cout << "min_gradient_mag " << io_info.min_gradient_mag << endl;
			cout << "\n";
		}
		else if (io_info.flag_reliable_scalar_prediction) {
			cout << "Scalar based prediction" << endl;
			cout << "scalar error tolerance: " << io_info.scalar_prediction_err
				<< endl;
			cout << "scalar error distance to neighbors: "
				<< io_info.scalar_prediction_dist << endl;
		}
		else if(io_info.adv_angle_based)
		{

			cout <<"Advanced Angle based Reliability criteria."<<endl;

			cout <<"Angle threshold for gradients to agree "<<io_info.param_angle << endl;
			cout <<"Neighbor angle "<<io_info.neighbor_angle_parameter<<endl;
			cout <<"\t{angle between grad at v and vector vv' where v' is an edge neighbor}"<<endl;

		}
		else if(io_info.curv_based)
		{

			cout <<"Curvature  based."<<endl;

			cout <<"Neighbor angle "<<io_info.neighbor_angle_parameter<<endl;
			cout <<"Angle threshold for gradients to agree "<<io_info.param_angle << endl;
		}
		if (io_info.flag_cdiff) {
			cerr << "Central Difference for computing gradients.\n";
		}

	}
}

void main_options_msg() {
	cerr << "OPTIONS:" << endl;
	cerr << "  [-cdiff] [-curvature_based] [-cdist {D}] [-extended_curv]"   
       << endl;
	cerr << "  [-min_gradient_mag {M}] [-angle {A}]" << endl;
	cerr << "  [-neighbor_angle {A}]" << endl;
	cerr << "  [-gzip] [-out_param] [-print_info {V}] [-print_grad_loc]"
       << endl;
  cerr << "  [-help] [-version] [-list_all_options]" << endl;
}

void testing_options_msg() {
	cerr << "TESTING OPTIONS:" << endl;
  cerr << "  [-angle_test] [-scalar_test] [-advangle]" << endl;
  cerr << "  [-min_num_agree {N}]" << endl;
	cerr << "  [-angle_based_dist {D}] [-reliable_scalar_pred_dist {D}]" << endl;
  cerr << "  [-scalar_pred_err {E}]" << endl;
  cerr << "  [-help_testing]" << endl;
}

void usage_msg(std::ostream & out) 
{
  out << "Usage: religrad [OPTIONS] {scalar nrrd file} {gradient nrrd file}"
       << endl;
}
  

void help_main_options()
{
  cout << "MAIN OPTIONS:" << endl;
	cout << "  -cdiff: Compute the central difference (default)." << endl;
	cout << "  -curvature_based: two parameters, alpha set by -angle and -neighbor_angle."<<endl;
	cout <<	"  -cdist {D} : distance associated with option -curvature_based."<< endl;
	cout <<	"     How check reliability at distance D (use 1 or 2) from vertex v" << endl;
	cout << "  -extended_curv: extended version of curvature based reliable gradients."<< endl
		 << "     Takes parameters -angle and -neighbor-angle"<<endl;
	cout << "  -min_gradient_mag {M}:  Set min gradient magnitude to {M} (float)." << endl;
	cout << "  -angle {A}: Set angle to {A} (float)." << endl;
	cout <<	"  -neighbor_angle {A} Angle between grad at v and vector vv'\n\t"
		<<	"Where v' is an edge neighbor. Default is 30." << endl;
	cout << "  -gzip: Store gradients in compressed (gzip) format." << endl;
	cout << "  -out_param:  Print parameters." << endl;
	cout << "  -print_info {V} : Print information about vertex {IV}." << endl;
	cout << "  -print_grad_loc : Print location of vertices with unreliable gradients." << endl;
  cout << "  -version: Print version." << endl;
  cout << "  -list_all_options: List all options." << endl;
	cout << "  -help:    Print this help message." << endl;
}

void help_testing_options()
{
  cout << "TESTING OPTIONS:" << endl;
	cout << "  -angle_test: Apply angle test." << endl;
	cout << "  -scalar_test: Apply scalar test." << endl;
	cout << "  -angle {A}: Set angle to {A} (float)." << endl;
	cout << "  -min_num_agree {N}: Set minimum agree number to {N}." << endl;
	cout << "     {N} gradient directions must agree to pass the angle test."<< endl;
	cout << "     (Default 4.)" << endl;
	cout << "  -angle_based_dist {D}: Distance (integer) to neighboring vertices"<< endl
		<<	"     in angle test.  (Default 1.)" << endl;
	cout << "  -reliable_scalar_pred_dist: Distance (integer) to neighboring vertices" << endl
		<<	"     in scalar test.  (Default 2.)" << endl;
	cout <<	"  -neighbor_angle {A} Angle between grad at v and vector vv'\n\t"
		<<	"Where v' is an edge neighbor. Default is 30." << endl;
	cout << "  -scalar_pred_err {E}:  Error threshold for scalar test." 
		<< endl;
	cout << "     Errors above the threshold fail the test. (Default 0.4.)" 
		<< endl;
	cout << "  -help_testing:    Print this help message." << endl;
}

void help_testing()
{
  usage_msg(cout);
  help_testing_options();

  exit(0);
}

void help(const bool flag_list_all_options = false) 
{
  usage_msg(cout);
  help_main_options();

  if (flag_list_all_options) {
    cout << endl;
    help_testing_options();
  }

	exit(0);
}

void usage_error(const bool flag_list_all_options = false) 
{
	usage_msg(cerr);
  main_options_msg();

  if (flag_list_all_options) {
    cerr << endl;
    testing_options_msg();
  }

	exit(100);
}

//Parse command line
void parse_command_line(int argc, char **argv, INPUT_INFO & io_info) {
	int iarg = 1;

  if (argc == 2 && std::string(argv[1]) == "-version") {
    cout << "Version: " << VERSION << endl;
    exit(0);
  }

  bool flag_list_all_options(false);
  bool flag_output_help(false);
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
			io_info.neighbor_angle_parameter = atof(argv[iarg]);
		}
		else if (s == "-curvature_based")
		{
			io_info.flag_cdiff = true;
			io_info.curv_based = true;
			//default parameters
			io_info.param_angle =  20; // alpha
			io_info.neighbor_angle_parameter = 20;
		}

		else if(s == "-cdist")
		{
			iarg++;
			io_info.cdist = atoi(argv[iarg]);
			if (io_info.curv_based != true)
			{
				cerr <<"USE the -curvature_based option."<< endl;
			}
		}
		else if (s == "-extended_curv")
		{
			io_info.flag_cdiff = true;
			io_info.curv_based = true;
			//default parameters
			io_info.param_angle = 20;
			io_info.neighbor_angle_parameter = 20;
			io_info.extended_curv_based = true;
		}
		else if (s == "-extend_max")
		{
			iarg++;
			io_info.extend_max = atoi(argv[iarg]);
		}
		else if(s == "-algo12")
		{
			io_info.algo12 = true; 
			iarg++;
			io_info.cdist = atoi(argv[iarg]);
			cerr <<"algorithm : " << io_info.cdist << ".\n";
		}
		else if (s == "-start_grads")
		{
			//read a file name
		}
		else if (s == "-gzip") {
			flag_gzip = true;
		} 
    else if (s == "-version") {
      cout << "Version: " << VERSION << endl;
    }
		else if (s == "-help") {
      flag_output_help = true;
		}
		else if (s == "-help_testing") {
      help_testing();
		}
    else if (s == "-list_all_options") {
      flag_list_all_options = true;
    }
		else {
			cerr << "Error in  " << s << endl;
      cerr << endl;
			usage_error();
		}
		iarg++;
	}

  if (flag_output_help) 
    { help(flag_list_all_options); }
  else if (flag_list_all_options) 
    { usage_error(true); }

	if (iarg + 2 != argc) {	usage_error(); };

	scalar_filename = argv[iarg];
	gradient_filename = argv[iarg + 1];
}

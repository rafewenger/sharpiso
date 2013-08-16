
#ifndef _RELIGRADIENT_INPUTINFO_H_
#define _RELIGRADIENT_INPUTINFO_H_
#include <cmath>
#include <vector>

/*
 * Things that you can print using the print_info function
 */
typedef enum output {
		CURR_VERTEX, PREV_VERTEX, NEXT_VERTEX, GRADIENT
	} out;

class OUTPUT_INFO{
public:
	unsigned int num_reliable; // total number of reliable vertices
	unsigned long int num_unreliable; // total number of un-reliable vertices
	unsigned  int boundary_verts;
	unsigned int grad_mag_zero;
	std::vector<unsigned int> un_reliable_grads_vert_info;
	void set_defaults(){
		num_reliable=0;
		num_unreliable=0;
		boundary_verts=0;
		grad_mag_zero=0;
	}
};
class INPUT_INFO {
public:
	bool flag_cdiff;	 		// compute central difference
	bool flag_reliable_grad;  // reliable grad
	bool print_info;          // print info of the vertex
	bool flag_print_grad_loc;      // prints the location of the unreliable grads

	//Compare the cdiff gradient with immediate neighbors or
	//neighbors at a certain distance
	bool angle_based;

	bool flag_reliable_scalar_prediction; // check how good the gradient predicts the scalar
	// of the neighborhood grid vertices
	int scalar_prediction_dist;
	float scalar_prediction_err;

	int angle_based_dist;
	int min_num_agree;
	int print_info_vertex;
	float param_angle;
	float min_gradient_mag;   // minimum gradient
	float min_cos_of_angle;
	bool draw;
	int draw_vert;
	int num_vertices_mag_grt_zero;


	OUTPUT_INFO	 out_info;    // generate some output_info
	//Set the default values
	void set_defaults(){
	  //central diff flags
		flag_cdiff = false;
		//angle based 
		angle_based = false;
		angle_based_dist = 2;
		//scalar_based
		flag_reliable_scalar_prediction = false;
		
		
		flag_reliable_grad = false;
		flag_print_grad_loc = false;
		print_info = false;
		min_gradient_mag = 0.001;
		param_angle=20;
		draw = false;
		
		scalar_prediction_dist = 2;
		scalar_prediction_err = .15;
		min_num_agree = 4;
		min_cos_of_angle = cos((param_angle*M_PI/180.0)); // 20 degrees=0.34906585, 30 degrees=0.523598776;
		out_info.set_defaults();
		num_vertices_mag_grt_zero=0;
	}

};


#endif 

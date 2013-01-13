/*
 * cgradient_inputinfo.h
 *
 *  Created on: Jan 5, 2013
 *      Author: arindam
 */

#ifndef CGRADIENT_INPUTINFO_H_
#define CGRADIENT_INPUTINFO_H_
#include <cmath>
#include <vector>

class OUTPUT_INFO{
public:
	unsigned int num_reliable; // total number of reliable vertices
	unsigned long int num_unreliable; // total number of un-reliable vertices
	unsigned  int boundary_verts;
	std::vector<unsigned int> un_reliable_grads_vert_info;
	void set_defaults(){
		num_reliable=0;
		num_unreliable=0;
		boundary_verts=0;
	}
};
class INPUT_INFO {
public:
	bool flag_cdiff;	 		// compute central difference
	bool flag_reliable_grad;  // reliable grad
	bool print_info;          // print info of the vertex
	bool flag_print_grad_loc;      // prints the location of the unreliable grads
	int print_info_vertex;
	float min_gradient_mag;   // minimum gradient
	float min_cos_of_angle;

	OUTPUT_INFO	 out_info;    // generate some output_info
	//Set the default values
	void set_defaults(){
		flag_cdiff = false;
		flag_reliable_grad = false;
		flag_print_grad_loc=false;
		print_info = false;
		min_gradient_mag = 0.001;
		float param=30;
		min_cos_of_angle = cos((param*M_PI/180.0));// 20 degrees=0.34906585, 30 degrees=0.523598776;
		out_info.set_defaults();
	}

};


#endif /* CGRADIENT_INPUTINFO_H_ */

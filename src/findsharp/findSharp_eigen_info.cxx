/*
 * findSharp_eigen_info.cxx
 *
 *  Created on: Mar 17, 2013
 *      Author: arindam
 */
#include "findSharp_eigen_info.h"
#include <fstream>
#include <iostream>
using namespace std;
void EIGEN_INFO::read_file()
{
	ifstream infile(eigen_info_fame.c_str());
	int n1,n2;
	int v_index , cube_index , lkup_table_index , patch_index ,
					 num_eigen_vals , flag_use_centroid;
		while(infile >> v_index >> cube_index >> lkup_table_index >> patch_index >>
				 num_eigen_vals >> flag_use_centroid)
		{
			//store the vertex index and the numb eigens in the vector
			vertex_index.push_back(v_index);
			num_eigen.push_back(num_eigen_vals);
			flag_centroid.push_back(flag_use_centroid);
		}
		cout <<"read complete"<< endl;
}




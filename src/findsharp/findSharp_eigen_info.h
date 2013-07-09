/*
 * findSharp_eigen_info.h
 *
 *  Created on: Mar 17, 2013
 *      Author: arindam
 */

#ifndef FINDSHARP_EIGEN_INFO_H_
#define FINDSHARP_EIGEN_INFO_H_
#include <string>
#include <vector>
using namespace std;

class EIGEN_INFO{
public:
	bool flag_eigen_based;
	string eigen_info_fame;
	void read_file();
	int num_degen_edges;
	vector <int> vertex_index;
	vector <int> num_eigen;
	vector <bool> flag_centroid;
	EIGEN_INFO()
	{
		flag_eigen_based = false;
		num_degen_edges = 0;

	}
};

#endif /* FINDSHARP_EIGEN_INFO_H_ */

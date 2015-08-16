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

class EIGEN_INFO{
public:
	bool flag_eigen_based;
  std::string eigen_info_filename;
	void read_file();
	int num_degen_edges;
  std::vector <int> vertex_index;
  std::vector <int> num_eigen;
  std::vector <bool> flag_centroid;
	EIGEN_INFO()
	{
		flag_eigen_based = false;
		num_degen_edges = 0;
	}
};

#endif /* FINDSHARP_EIGEN_INFO_H_ */

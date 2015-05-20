/// \file religrad_computations.cxx
///compute reliable gradients from scalar data
#include "religrad_computations.h"
#include "ijkscalar_grid.txx"
#include "ijkgrid.txx"
#include "ijkcoord.txx"
#include "sharpiso_scalar.txx"
#include "ijkgrid_macros.h"
#include "sharpiso_types.h"
#include <algorithm>
#include <vector>
using namespace RELIGRADIENT;
using namespace std;
using namespace SHARPISO;


/**
* Does orthogonal directions match
**/
bool does_ortho_match(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info);
bool does_ortho_matchA
	(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info);
bool does_ortho_matchB
	(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info);

//debug print 
void print_vertex(const int v,
				  const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid)
{
	COORD_TYPE c[DIM3];
	scalar_grid.ComputeCoord(v, c);
	cout <<"vertex "<< v <<" coord "<< c[0]<<","<<c[1]<<","<<c[2]<< endl;
}

/// Compute central difference per vertex
void compute_gradient_central_difference(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient,
	const GRADIENT_COORD_TYPE & min_gradient_mag) {
		const int dimension = scalar_grid.Dimension();
		for (int d = 0; d < dimension; d++) {
			VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
			VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
			COORD_TYPE dist = scalar_grid.Spacing(d) * 2.0;
			gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))
				/ dist;
		}
		// set min gradient to )
		SCALAR_TYPE mag = 0.0;
		IJK::compute_magnitude_3D(gradient, mag);
		if (mag < min_gradient_mag) {
			IJK::set_coord(DIM3, 0.0, gradient);
		}

}
//Compute central difference per vertex , normalized version
void compute_gradient_central_difference_normalized(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient,
	GRADIENT_COORD_TYPE & grad_mag,
	const GRADIENT_COORD_TYPE & min_gradient_mag) {
		const int dimension = scalar_grid.Dimension();
		for (int d = 0; d < dimension; d++) {
			VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
			VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
			COORD_TYPE dist = scalar_grid.Spacing(d) * 2.0;
			gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))
				/ dist;
		}
		// set min gradient to )

		IJK::compute_magnitude_3D(gradient, grad_mag);
		if (grad_mag < min_gradient_mag) {
			IJK::set_coord(DIM3, 0.0, gradient);
		} else {
			for (int i = 0; i < dimension; i++) {
				gradient[i] = gradient[i] / grad_mag;
			}
		}
}

/// Compute_boundary_gradients
void compute_boundary_gradient(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient) {
		const int dimension = scalar_grid.Dimension();
		//GRID_COORD_TYPE coord[dimension];
		GRID_COORD_TYPE * coord = new GRID_COORD_TYPE [dimension];

		scalar_grid.ComputeCoord(iv1, coord);

		for (int d = 0; d < dimension; d++) {
			if (coord[d] > 0) {
				VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
				if (coord[d] + 1 < scalar_grid.AxisSize(d)) {
					VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
					// use central difference
					COORD_TYPE dist = scalar_grid.Spacing(d) * 2.0;
					gradient[d] =
						(scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))
						/ (2.0 * dist);
				} else {
					COORD_TYPE dist = scalar_grid.Spacing(d);
					gradient[d] =
						(scalar_grid.Scalar(iv1) - scalar_grid.Scalar(iv0))
						/ dist;
				}
			} else if (coord[d] + 1 < scalar_grid.AxisSize(d)) {
				VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
				COORD_TYPE dist = scalar_grid.Spacing(d);
				gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv1))
					/ dist;
			} else {
				gradient[d] = 0;
			}
		}
		delete [] coord;
}

// compute the boundary gradient the normalized version/
void compute_boundary_gradient_normalized(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX iv1, GRADIENT_COORD_TYPE * gradient,
	GRADIENT_COORD_TYPE & grad_mag,
	const GRADIENT_COORD_TYPE & min_gradient_mag) {
		//call the generic boundary gradient
		compute_boundary_gradient(scalar_grid, iv1, gradient);

		const int dimension = scalar_grid.Dimension();
		IJK::compute_magnitude_3D(gradient, grad_mag);
		if (grad_mag < min_gradient_mag) {
			IJK::set_coord(DIM3, 0.0, gradient);
		} else {

			for (int i = 0; i < dimension; i++) {
				gradient[i] = gradient[i] / grad_mag;
			}
		}
}

/*
* Compute gradients using central difference
*/
void compute_gradient_central_difference(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	GRADIENT_GRID & gradient_grid, const INPUT_INFO & io_info) {
		const int dimension = scalar_grid.Dimension();
		gradient_grid.SetSize(scalar_grid, dimension);
		gradient_grid.SetSpacing(scalar_grid);
		BOOL_GRID boundary_grid;
		boundary_grid.SetSize(scalar_grid);
		compute_boundary_grid(boundary_grid);
		for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
			if (boundary_grid.Scalar(iv)) {
				compute_boundary_gradient(scalar_grid, iv,
					gradient_grid.VectorPtr(iv));
			} else {
				compute_gradient_central_difference(scalar_grid, iv,
					gradient_grid.VectorPtr(iv), io_info.min_gradient_mag);
			}
		}
}

//compute gradient using central difference
//also computed the normalized gradients and the gradient magnitudes
void compute_gradient_central_difference_normalized(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	GRADIENT_GRID & normalized_grad_grid,
	GRADIENT_MAGNITUDE_GRID & grad_magnitude_grid,
	const INPUT_INFO & io_info) {
		const int dimension = scalar_grid.Dimension();
		normalized_grad_grid.SetSize(scalar_grid, dimension);
		normalized_grad_grid.SetSpacing(scalar_grid);
		grad_magnitude_grid.SetSize(scalar_grid);

		BOOL_GRID boundary_grid;
		boundary_grid.SetSize(scalar_grid);
		compute_boundary_grid(boundary_grid);
		for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
			GRADIENT_COORD_TYPE grad_mag = 0.0;
			if (boundary_grid.Scalar(iv)) {
				compute_boundary_gradient_normalized(scalar_grid, iv,
					normalized_grad_grid.VectorPtr(iv), grad_mag,
					io_info.min_gradient_mag);
			} else {
				compute_gradient_central_difference_normalized(scalar_grid, iv,
					normalized_grad_grid.VectorPtr(iv), grad_mag,
					io_info.min_gradient_mag);
			}
			grad_magnitude_grid.Set(iv, grad_mag);
		}
}

/* 
*  Check if two gradients are within min_cos_angle.
*  @param vert2, vertex to compare with.
*  @param gradient1. gradient of vert1.
*  @param vert1 is redundant.
*  
*/

bool gradients_agree(
	const GRADIENT_GRID & gradient_grid,
	const GRADIENT_COORD_TYPE * gradient1,
	const VERTEX_INDEX vert2,
	const INPUT_INFO & io_info
	) 
{
	GRADIENT_COORD_TYPE mag1 = 0.0, mag2 = 0.0; // gradient magnitude of vert2
	GRADIENT_COORD_TYPE gradient2[DIM3] = { 0.0, 0.0, 0.0 };

	std::copy(gradient_grid.VectorPtrConst(vert2),
		gradient_grid.VectorPtrConst(vert2) + DIM3, &(gradient2[0]));
	IJK::compute_magnitude_3D(gradient2, mag2);
	if (mag2 > io_info.min_gradient_mag) {
		// compute unit vertex gradient
		//gradient2[0] = gradient2[0] / mag2;
		//gradient2[1] = gradient2[1] / mag2;
		//gradient2[2] = gradient2[2] / mag2;
		//IJK::normalize_vector(DIM3, gradient2,io_info.min_gradient_mag, gradient2);
		IJK::multiply_coord_3D(1.0/mag2, gradient2, gradient2);
		GRADIENT_COORD_TYPE inn_pdt = 0;
		IJK::compute_inner_product(DIM3, gradient1, gradient2, inn_pdt);

		if (inn_pdt > io_info.min_cos_of_angle) {
			return true;
		}
	}
	return false;
}
/* 
*  Check if two gradients are within min_cos_angle.
*  @param vert1, 
*  @param vert2, vertex to compare with. 
*/

bool gradients_agree(
	const GRADIENT_GRID & gradient_grid,
	const VERTEX_INDEX vert1,
	const VERTEX_INDEX vert2,
	const INPUT_INFO & io_info
	) 
{
	GRADIENT_COORD_TYPE mag1 = 0.0, mag2 = 0.0; // gradient magnitude of vert2

	GRADIENT_COORD_TYPE gradient2[DIM3] = { 0.0, 0.0, 0.0 };
	GRADIENT_COORD_TYPE gradient1[DIM3] = { 0.0, 0.0, 0.0 };

	std::copy(gradient_grid.VectorPtrConst(vert2),
		gradient_grid.VectorPtrConst(vert2) + DIM3, &(gradient2[0]));

	std::copy(gradient_grid.VectorPtrConst(vert1),
		gradient_grid.VectorPtrConst(vert1) + DIM3, &(gradient1[0]));

	IJK::compute_magnitude_3D(gradient2, mag2);
	IJK::compute_magnitude_3D(gradient1, mag1);
	if (mag1 > io_info.min_gradient_mag && mag2 > io_info.min_gradient_mag) {
		// compute unit vertex gradient
		IJK::multiply_coord_3D(1.0/mag2, gradient2, gradient2);
		IJK::multiply_coord_3D(1.0/mag1, gradient1, gradient1);

		GRADIENT_COORD_TYPE inn_pdt = 0;
		IJK::compute_inner_product(DIM3, gradient1, gradient2, inn_pdt);

		if (inn_pdt > io_info.min_cos_of_angle) {
			return true;
		}
	}
	return false;
}

/*
* Compute reliable gradients by comparing with neighboring gradients
* neighboring gradients are defined by "angle_based_dist"
*/
void compute_reliable_gradients_angle(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	GRADIENT_GRID & gradient_grid,
	GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	BOOL_GRID & reliable_grid,
	INPUT_INFO & io_info) 
{
	int numAgree = 0;
	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		numAgree = 0;
		GRADIENT_COORD_TYPE gradient_iv[DIM3] = { 0.0, 0.0, 0.0 };
		GRADIENT_COORD_TYPE gradient_iv_mag = grad_mag_grid.Scalar(iv);

		std::copy(gradient_grid.VectorPtrConst(iv),
			gradient_grid.VectorPtrConst(iv) + DIM3, &(gradient_iv[0]));

		if (gradient_iv_mag > io_info.min_gradient_mag) {
			COORD_TYPE coord_iv[DIM3];
			scalar_grid.ComputeCoord(iv, coord_iv);

			for (int d = 0; d < DIM3; d++) 
			{
				int k = min(io_info.angle_based_dist, int(coord_iv[d]));
				VERTEX_INDEX prev_vertex = 
					iv - k * scalar_grid.AxisIncrement(d);

				if (gradients_agree
					(gradient_grid, gradient_iv, prev_vertex, 
					io_info)){
						numAgree++;
				}

				k = min(io_info.angle_based_dist,
					(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));

				VERTEX_INDEX next_vertex = iv
					+ k * scalar_grid.AxisIncrement(d);

				if ( gradients_agree
					(gradient_grid, gradient_iv, next_vertex, 
					io_info))
				{
					numAgree++;
				}
			}

			if (numAgree < io_info.min_num_agree) 
			{
				reliable_grid.Set(iv, false);
				io_info.out_info.num_unreliable++;
			} 
			else 
			{
				io_info.out_info.num_reliable++;
			}
		} 
		else 
		{
			io_info.out_info.grad_mag_zero++;
		}
	}
}

// scalar based 
void compute_plane_point_dist(
	const GRADIENT_COORD_TYPE * normal,
	const COORD_TYPE * pt_on_plane,
	const COORD_TYPE * far_point,
	SCALAR_TYPE & dist) 
{
	COORD_TYPE d = 0; //  n.pt on plane
	COORD_TYPE temp = 0; // n.far_point
	IJK::compute_inner_product_3D(normal, pt_on_plane, d);
	IJK::compute_inner_product_3D(normal, far_point, temp);
	dist = temp - d;
}

// Scalar based prediction
void compute_reliable_gradients_SBP
	(const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info) 
{
	GRADIENT_COORD_TYPE grad_iv[DIM3];  // unscaled gradient vector
	GRADIENT_COORD_TYPE normalized_grad_iv[DIM3]; // normalized grad_iv
	COORD_TYPE coord_iv[DIM3]; // Coordinates of vertex iv.
	COORD_TYPE coord_nv[DIM3]; // Coordinates of vertex nv.
	GRADIENT_COORD_TYPE mag1;
	COORD_TYPE err_distance;


	bool debug = false;
	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {


		// only run the test if gradient at vertex iv is reliable
		if (reliable_grid.Scalar(iv)) {

			mag1 = grad_mag_grid.Scalar(iv);
			if (mag1 > io_info.min_gradient_mag) {

				// run test on unscaled, uniform grid (all edge lengths 1)
				scalar_grid.ComputeCoord(iv, coord_iv);

				// grad_iv[] is an unscaled gradient
				for (int d = 0; d < DIM3; d++) 
				{
					COORD_TYPE spacing_d = scalar_grid.Spacing(d);
					grad_iv[d] = mag1*spacing_d*gradient_grid.Vector(iv,d);
				}

				IJK::normalize_vector
					(DIM3, grad_iv, io_info.min_gradient_mag, normalized_grad_iv);
				if(debug)
				{
					cout <<"gradiv " << grad_iv[0]<<" "<< grad_iv[1]<<" "<<grad_iv[2]<<endl;
					cout <<"normalized " << normalized_grad_iv[0]<<" "<< normalized_grad_iv[1]
					<<" "<< normalized_grad_iv[2]<<endl;
				}
				int num_agree = 0;
				// set up a vector to keep track of the distances
				vector<SCALAR_TYPE> vec_scalar_dists;

				io_info.num_vertices_mag_grt_zero++;
				// find the normalized gradient
				// point on the plane
				// find neighbor vertices
				vector<VERTEX_INDEX> near_vertices;
				scalar_grid.GetVertexNeighbors
					(iv, io_info.scalar_prediction_dist, near_vertices);

				// for all the neighboring points find the distance of points to plane
				bool flag_correct = true;
				for (int nv_ind = 0; nv_ind < near_vertices.size(); nv_ind++) {

					VERTEX_INDEX nv = near_vertices[nv_ind];
					scalar_grid.ComputeCoord(nv, coord_nv);

					// compute distance to plane

					COORD_TYPE dist_to_plane = 0.0;
					// compute distance in unscaled grid (uniform, unit edge lengths)
					compute_plane_point_dist
						(normalized_grad_iv, coord_iv, coord_nv, dist_to_plane);

					// if dist_to_plane is within the threshold
					if (abs(dist_to_plane) <= 0.5) {

						// compute the distance between prediction and observed scalar value

						// compute distance in unscaled grid (uniform, unit edge lengths)
						compute_distance_to_gfield_plane
							(grad_iv, coord_iv, scalar_grid.Scalar(iv), coord_nv,
							scalar_grid.Scalar(nv), err_distance);

						if (err_distance > io_info.scalar_prediction_err) {
							flag_correct = false;
							break;
						}
					}
				}

				if (!flag_correct) {
					reliable_grid.Set(iv, false);
					io_info.out_info.num_unreliable++;
				}
				else {
					io_info.out_info.num_reliable++;
				}
			}
		}
	}
}

// Compute vector between two grid coords
// parameters
// scalar_grid, 
// cube_index1, 
// cube_index2,
//return vector between cube_index2-cube_index1
void compute_vector_between_grid_vertex(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const VERTEX_INDEX ci1,
	const VERTEX_INDEX ci2,
	COORD_TYPE *vec)
{
	COORD_TYPE coord1[DIM3], coord2[DIM3];
	scalar_grid.ComputeCoord(ci1, coord1);
	scalar_grid.ComputeCoord(ci2, coord2);
	for (int i = 0; i < DIM3; i++)
	{
		vec[i] = coord2[i]-coord1[i];
	}
}

/* 
* Update tangent_vertex_list.
* by checking of the angle between the vertex gradient and 
* coordinates.
* @param vertex.  The vertex being considered to be inserted to the list
* @param inn_pdt. Inner product between vector(original vertex and @param1)
*				  and gradient at original vertex.
* @param tangent_vertex_list. Updated tangent_vertex
* @Return  true, if inserted.
*/
bool update_tangent_vertex_list( 
	const VERTEX_INDEX vertex, 
	const float inn_pdt,
	const int degree_threshold,  // threshold for grid neighbors
	vector<VERTEX_INDEX> &tangent_vertex_list)
{
	float inn_pdt_degree = acos(inn_pdt)*57.2957795;
	if(inn_pdt >= 0.0)
	{
		if(inn_pdt_degree > degree_threshold)
		{
			tangent_vertex_list.push_back(vertex);
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		// inn pdt is negative
		if(inn_pdt_degree  < (180-degree_threshold))
		{
			tangent_vertex_list.push_back(vertex);
			return true;
		}
		else
		{
			return false;
		}
	}
}

/* 
* Advanced angle based reliable gradients computations
* @return grad_mag_grid. Gradient magnitude grid. 
* @return reliable_grid. Reliable gradients  grid.
*/
void compute_reliable_gradients_advangle(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> &reliable_grid,
	INPUT_INFO & io_info)
{
	bool debug = false;

	//float degree_param = 30*M_PI/180.0;
	float degree_param = io_info.neighbor_angle_parameter*M_PI/180.0;
	//DEBUG
	/*cout <<"degree param cos "<< cos(30*M_PI/180.0) <<endl;
	cout <<"num vertices "<<scalar_grid.NumVertices()<<endl;
	cout <<"axis size "<< scalar_grid.AxisSize(0) 
	<<", "<<scalar_grid.AxisSize(1)
	<<", "<<scalar_grid.AxisSize(2);

	cout <<" axis increment  "<<  scalar_grid.AxisIncrement(0)
	<<", " <<scalar_grid.AxisIncrement(1) <<"," <<scalar_grid.AxisIncrement(2);*/
	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) 
	{
		COORD_TYPE coord_iv[DIM3] = {0.0,0.0,0.0};
		COORD_TYPE coord1[DIM3] = {0.0,0.0,0.0};
		scalar_grid.ComputeCoord(iv, coord_iv);
		if(debug){
			cout <<"\n---------------------\nvertex "<< iv <<
				" coord ["<<coord_iv[0]<<","<<coord_iv[1]<<","<<coord_iv[2]<<"]"<<endl;
		}
		int numAgree = 0;

		vector<VERTEX_INDEX> tangent_vertex_list;

		GRADIENT_COORD_TYPE gradient_iv[DIM3] = { 0.0, 0.0, 0.0 };
		GRADIENT_COORD_TYPE gradient_iv_mag = grad_mag_grid.Scalar(iv);

		std::copy(gradient_grid.VectorPtrConst(iv),
			gradient_grid.VectorPtrConst(iv) + DIM3, &(gradient_iv[0]));

		if(debug){
			cout <<"Grad "<< gradient_iv[0]<<","<< gradient_iv[1]<<","<< gradient_iv[2]
			<<" Mag "<< gradient_iv_mag<< endl;
		}
		if (gradient_iv_mag > io_info.min_gradient_mag) 
		{
			for (int d = 0; d < DIM3; d++) 
			{
				//cout <<"D "<< d <<"\t";
				int k = min(io_info.angle_based_dist, int(coord_iv[d]));
				VERTEX_INDEX prev_vertex = 
					iv - k * scalar_grid.AxisIncrement(d);
				//scalar_grid.ComputeCoord(prev_vertex, &(coord1[0]));				

				if(debug){
					cout <<"previous vertex {from axis increment} "
						<<prev_vertex<<"  "<< coord1[0]<<","<<coord1[1]<<","<<coord1[2];
				}
				COORD_TYPE edge_vec[DIM3] = {0, 0, 0};
				compute_vector_between_grid_vertex
					(scalar_grid, iv, prev_vertex, &(edge_vec[0]));
				IJK::normalize_vector(DIM3, &(edge_vec[0]), 0.0, &(edge_vec[0]));

				//DEBUG
				//cout <<" ["<< edge_vec[0]<<","<< edge_vec[1]<<","<< edge_vec[2]<<"] ";

				float inn_pdt = 0.0;
				IJK::compute_inner_product(DIM3, gradient_iv, edge_vec, inn_pdt);

				//debug
				//cout <<" inn pdt " << inn_pdt <<" angle "<< acos(inn_pdt)*57.2957795;

				bool insert_flag = update_tangent_vertex_list
					(prev_vertex, inn_pdt, degree_param, tangent_vertex_list);

				//debug 
				/*if(insert_flag)
				{
				cout <<" inserted.  "<< endl;
				}
				else
				{
				cout <<" NOT inserted."<< endl;
				}*/
				k = min(io_info.angle_based_dist,
					(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));

				VERTEX_INDEX next_vertex = iv
					+ k * scalar_grid.AxisIncrement(d);
				//scalar_grid.ComputeCoord(next_vertex, &(coord1[0]));

				//cout <<"next vertex "<< next_vertex<<" "
				//	<< coord1[0]<<","<<coord1[1]<<","<<coord1[2];

				compute_vector_between_grid_vertex
					(scalar_grid, iv, next_vertex, &(edge_vec[0]));
				IJK::normalize_vector(DIM3, &(edge_vec[0]), 0.0, &(edge_vec[0]));

				//DEBUG
				//cout <<" ["<< edge_vec[0]<<","<< edge_vec[1]<<","<< edge_vec[2]<<"] ";
				IJK::compute_inner_product(DIM3, gradient_iv, edge_vec, inn_pdt);

				//debug
				//cout <<" inn pdt " << inn_pdt <<" angle "<< acos(inn_pdt)*57.2957795;
				insert_flag = update_tangent_vertex_list
					(next_vertex, inn_pdt, degree_param, tangent_vertex_list);

				//debug 
				/*if(insert_flag)
				{
				cout <<" inserted." << endl;
				}
				else
				{
				cout <<" NOT inserted." << endl;
				}*/
			}
			//cout <<"tangent neighbor sizes " << tangent_vertex_list.size() << endl;
			for (int v = 0; v < tangent_vertex_list.size(); v++)
			{
				if(gradients_agree(gradient_grid, iv, tangent_vertex_list[v], io_info))
				{
					numAgree++;
				}
			}

			if(numAgree != tangent_vertex_list.size())
			{
				//cout <<"not reliable "<< tangent_vertex_list.size()<<" vs "<< numAgree<<endl;

				reliable_grid.Set(iv, false);
				io_info.out_info.num_unreliable++;
			}
			else
			{
				//cout <<"reliable\n";
				io_info.out_info.num_reliable++;
			}

		}// if end 
	}
}



void compute_reliable_gradients_advangle_version2(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info) 
{
	GRADIENT_COORD_TYPE grad_iv[DIM3];  // unscaled gradient vector
	GRADIENT_COORD_TYPE normalized_grad_iv[DIM3]; // normalized grad_iv
	COORD_TYPE coord_iv[DIM3]; // Coordinates of vertex iv.
	COORD_TYPE coord_nv[DIM3]; // Coordinates of vertex nv.
	GRADIENT_COORD_TYPE mag1;
	COORD_TYPE err_distance;


	bool debug = false;
	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		int numAgree = 0;
		// only run the test if gradient at vertex iv is reliable
		if (reliable_grid.Scalar(iv)) {

			mag1 = grad_mag_grid.Scalar(iv);
			if (mag1 > io_info.min_gradient_mag) {

				// run test on unscaled, uniform grid (all edge lengths 1)
				scalar_grid.ComputeCoord(iv, coord_iv);

				// grad_iv[] is an unscaled gradient
				for (int d = 0; d < DIM3; d++) 
				{
					COORD_TYPE spacing_d = scalar_grid.Spacing(d);
					grad_iv[d] = mag1*spacing_d*gradient_grid.Vector(iv,d);
				}

				IJK::normalize_vector
					(DIM3, grad_iv, io_info.min_gradient_mag, normalized_grad_iv);
				if(debug)
				{
					cout <<"iv "<< iv << endl;
					cout <<"gradiv " << grad_iv[0]<<" "<< grad_iv[1]<<" "<<grad_iv[2]<<endl;
					cout <<"normalized " << normalized_grad_iv[0]<<" "<< normalized_grad_iv[1]
					<<" "<< normalized_grad_iv[2]<<endl;
				}
				int num_agree = 0;
				// set up a vector to keep track of the distances
				vector<SCALAR_TYPE> vec_scalar_dists;

				io_info.num_vertices_mag_grt_zero++;
				// find the normalized gradient
				// point on the plane
				// find neighbor vertices
				vector<VERTEX_INDEX> near_vertices;
				scalar_grid.GetVertexNeighbors
					(iv, io_info.scalar_prediction_dist, near_vertices);

				// for all the neighboring points find the distance of points to plane
				bool flag_correct = true;
				vector<VERTEX_INDEX> tangent_vertex_list;
				for (int nv_ind = 0; nv_ind < near_vertices.size(); nv_ind++) 
				{
					VERTEX_INDEX nv = near_vertices[nv_ind];
					scalar_grid.ComputeCoord(nv, coord_nv);

					// compute distance to plane

					COORD_TYPE dist_to_plane = 0.0;
					// compute distance in unscaled grid (uniform, unit edge lengths)
					compute_plane_point_dist
						(normalized_grad_iv, coord_iv, coord_nv, dist_to_plane);

					// if dist_to_plane is within the threshold
					if (abs(dist_to_plane) <= 0.3) 
					{
						tangent_vertex_list.push_back(nv);
					}
				}//neighbor loop ends


				for (int v = 0; v < tangent_vertex_list.size(); v++)
				{
					if(gradients_agree(gradient_grid, iv, tangent_vertex_list[v], io_info))
					{
						numAgree++;
					}
				}

				if(numAgree != tangent_vertex_list.size())
				{
					if(debug)
						cout <<"not reliable "<< tangent_vertex_list.size()<<" vs "<< numAgree<<endl;

					reliable_grid.Set(iv, false);
					io_info.out_info.num_unreliable++;
				}
				else
				{
					if(debug)
						cout <<"reliable\n"<< tangent_vertex_list.size()<<" vs "<< numAgree<<endl;
					io_info.out_info.num_reliable++;
				}
			}
		}
	}
}






//DEBUG

bool debug = false;
int test_iv = 145;

/*
* check if neighbor v of iv can be put into the
* tangent neighbor set or orthogonal set.
*/
void check_v_for_insert2_Nt_and_No_iv(
	const VERTEX_INDEX iv, 
	GRADIENT_COORD_TYPE *norm_grad_iv,
	const VERTEX_INDEX v, 
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	vector <VERTEX_INDEX> &tangent_neighbor_set,
	vector <VERTEX_INDEX> &ortho_neighbor_set, INPUT_INFO & io_info)
{
	COORD_TYPE edge_vec[DIM3] = {0, 0, 0};
	float degree_param = io_info.neighbor_angle_parameter;
	compute_vector_between_grid_vertex
		(scalar_grid, iv, v, &(edge_vec[0]));
	IJK::normalize_vector(DIM3, &(edge_vec[0]), 0.0, &(edge_vec[0]));
	float inn_pdt = 0.0;
	IJK::compute_inner_product(DIM3, norm_grad_iv, edge_vec, inn_pdt);
	bool insert_flag = update_tangent_vertex_list
		(v, inn_pdt, degree_param, tangent_neighbor_set);
	if (!insert_flag)
	{
		ortho_neighbor_set.push_back(v);
	}
}
/*
* Compute tangent neighors set
* @param iv , find tangent neighbors of iv
* Takes care of the boundary
*/
void compute_Nt_and_No_iv
	(const VERTEX_INDEX iv, GRADIENT_COORD_TYPE *norm_grad_iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	vector <VERTEX_INDEX> &tangent_neighbor_set,
	vector <VERTEX_INDEX> &ortho_neighbor_set, INPUT_INFO & io_info)
{
	//insert_flag = update_tangent_vertex_list
	//(next_vertex, inn_pdt, degree_param, tangent_vertex_list);
	COORD_TYPE coord_iv[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(iv, coord_iv);
	for (int d = 0; d < DIM3; d++) 
	{
		int k = min(1, int(coord_iv[d]));
		VERTEX_INDEX prev_vertex = 
			iv - k * scalar_grid.AxisIncrement(d);
		check_v_for_insert2_Nt_and_No_iv(iv, norm_grad_iv, prev_vertex, scalar_grid, 
			tangent_neighbor_set, ortho_neighbor_set, io_info);

		k = min(1,(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));
		VERTEX_INDEX next_vertex = 
			iv + k * scalar_grid.AxisIncrement(d);
		check_v_for_insert2_Nt_and_No_iv(iv, norm_grad_iv, next_vertex, scalar_grid, 
			tangent_neighbor_set, ortho_neighbor_set, io_info);
	}
}
/*
* Compute the neighbor set of a vertex iv
*/
void compute_neighbor_set(
	const VERTEX_INDEX iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	vector<VERTEX_INDEX> &neighbor_set)
{
	COORD_TYPE coord_iv[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(iv, coord_iv);

	for (int d = 0; d < DIM3; d++) 
	{
		int k = min(1, int(coord_iv[d]));
		VERTEX_INDEX prev_vertex = 
			iv - k * scalar_grid.AxisIncrement(d);
		neighbor_set.push_back(prev_vertex);
		k = min(1,(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));
		VERTEX_INDEX next_vertex = 
			iv + k * scalar_grid.AxisIncrement(d);
		neighbor_set.push_back(next_vertex);
	}
}


/**
* Compute phi between v and v1, there is a new version B 
* Version B Works for vertices at distance 2
**/
void compute_phi(const VERTEX_INDEX v, const VERTEX_INDEX v1,
				 const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
				 const GRADIENT_GRID & gradient_grid,
				 const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
				 float * phi_vector)
{
	COORD_TYPE gradv[DIM3]={0.0,0.0,0.0};
	COORD_TYPE gradv1[DIM3]={0.0,0.0,0.0};
	float min_grad_mag = 0.0;
	GRADIENT_COORD_TYPE magv = grad_mag_grid.Scalar(v);
	GRADIENT_COORD_TYPE magv1 = grad_mag_grid.Scalar(v1);

	// grad_iv[] is an unscaled gradient
	for (int d = 0; d < DIM3; d++) 
	{
		COORD_TYPE spacing_d = scalar_grid.Spacing(d);
		gradv[d] = magv*spacing_d*gradient_grid.Vector(v,d);
	}
	IJK::normalize_vector
		(DIM3, gradv, min_grad_mag, gradv);
	//v1
	for (int d = 0; d < DIM3; d++) 
	{
		COORD_TYPE spacing_d = scalar_grid.Spacing(d);
		gradv1[d] = magv1*spacing_d*gradient_grid.Vector(v1,d);
	}
	IJK::normalize_vector
		(DIM3, gradv1, min_grad_mag, gradv1);

	float innpdt = 0.0;
	IJK::compute_inner_product(DIM3, gradv, gradv1, innpdt);
	for (int d = 0; d < DIM3; d++)
	{
		phi_vector[d] =  2*(innpdt)*gradv1[d] - gradv[d];
		//phi_vector[d] = phi_vector[d] + 2*(innpdt)*gradv1[d]-gradv[d];
		//phi_vector[d]=phi_vector[d]+gradv[d]+2*(innpdt)*gradv1[d];
	}
	IJK::normalize_vector(DIM3, phi_vector, 0.0, phi_vector);
	//DEBUG
	if(debug)
		cout <<"phi "<< phi_vector[0]<<" "<< phi_vector[1]<<" "<<phi_vector[2]<<endl;
	//DEBUG END
}

/*
* lambda computations
* Returns an angle in degrees
* THIS IS NOT CORRECT, spacing is alredy accounted for. 
*/
float lambda_computations(const VERTEX_INDEX iv, 
						  const VERTEX_INDEX v1, const VERTEX_INDEX v2,
						  const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
						  const GRADIENT_GRID & gradient_grid,
						  const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid)
{
	float  phi_vector[DIM3] = {0.0,0.0,0.0};

	compute_phi (iv, v1, scalar_grid, gradient_grid, grad_mag_grid,
		&(phi_vector[0]));

	GRADIENT_COORD_TYPE magv2 = grad_mag_grid.Scalar(v2);
	COORD_TYPE gradv2[DIM3] = {0.0,0.0,0.0};
	for (int d = 0; d < DIM3; d++) 
	{
		COORD_TYPE spacing_d = scalar_grid.Spacing(d);
		gradv2[d] = magv2*spacing_d*gradient_grid.Vector(v2,d);
	}
	IJK::normalize_vector (DIM3, gradv2, 0.0, gradv2);

	//
	float inn_pdt = 0.0;
	IJK::compute_inner_product(DIM3, gradv2, phi_vector, inn_pdt);

	//DEBUG
	if(debug){
		cout <<"gradV2 "<< gradv2[0]<<","<<gradv2[1]<<","<<gradv2[2] <<" inn pdt "<< inn_pdt 
			<<" angle "<<acos (inn_pdt) * 180.0 / M_PI<<endl;}
	//DEBUG end
	return acos (inn_pdt) * 180.0 / M_PI;
}
/*
* Check if vectors iv,v1 and v1,v2 are collinear.
*/
bool are_collinear(const VERTEX_INDEX iv,
				   const VERTEX_INDEX v1,
				   const VERTEX_INDEX v2,
				   const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid)
{
	//compute vectors and dot pdt for collinear
	COORD_TYPE edge_vec1[DIM3]={0.0,0.0,0.0};
	COORD_TYPE edge_vec2[DIM3]={0.0,0.0,0.0};
	compute_vector_between_grid_vertex
		(scalar_grid, iv, v1, &(edge_vec1[0]));
	IJK::normalize_vector(DIM3, &(edge_vec1[0]), 0.0, &(edge_vec1[0]));
	compute_vector_between_grid_vertex
		(scalar_grid, v1, v2, &(edge_vec2[0]));
	IJK::normalize_vector(DIM3, &(edge_vec2[0]), 0.0, &(edge_vec2[0]));

	float inn_pdt = 0.0;
	IJK::compute_inner_product(DIM3, edge_vec1, edge_vec2, inn_pdt);

	//DEBUG
	if(debug){
		using namespace std;
		COORD_TYPE c[DIM3]={0.0,0.0,0.0};
		scalar_grid.ComputeCoord(v2,c);
		cout <<"neighbor of {"<<v1<<"}, is "<< v2 
			<<" coordv2 {"<<c[0]<<" "<<c[1]<<" "<<c[2]<<"} "
			<<" col "<<int(abs(inn_pdt))<<"\n";
	}
	//DEBUG end 

	//collinear
	if((int(abs(inn_pdt))-1)==0)
		return true;
	else
		return false;
}
/*
* DEPRECATED
* Tangent neighbor based computations old
*/
bool Nt_based_computations_for_iv_old(
	const VERTEX_INDEX iv,
	float alpha,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	vector <VERTEX_INDEX> &tangent_neighbor_set
	)
{
	for (int v = 0; v < tangent_neighbor_set.size(); v++)
	{
		VERTEX_INDEX v1 = tangent_neighbor_set[v];
		vector<VERTEX_INDEX> neighbor_set_v1;
		compute_neighbor_set(v1, scalar_grid, neighbor_set_v1);
		for (int j = 0; j < neighbor_set_v1.size(); j++)
		{
			VERTEX_INDEX v2 = neighbor_set_v1[j];

			bool flag_collinear = are_collinear(iv,v1,v2, scalar_grid);
			if(flag_collinear)
			{
				if(lambda_computations
					(iv,v1,v2, scalar_grid, gradient_grid, grad_mag_grid) > alpha)
				{
					if(debug)
						cout <<" **fail** ";
					return false;
				}
			}
		}
	}
	return true;
}
/*
* Orthogonal neighbor based computations  OLD
*/
bool No_based_computations_for_iv_old(
	const VERTEX_INDEX iv,
	float alpha,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	vector <VERTEX_INDEX> &ortho_neighbor_set)
{
	int m=0;
	for (int i = 0; i < ortho_neighbor_set.size(); i++)
	{
		VERTEX_INDEX v1 = ortho_neighbor_set[i];
		vector<VERTEX_INDEX> neighbor_set_v1;
		compute_neighbor_set(v1, scalar_grid, neighbor_set_v1);
		for (int j = 0; j < neighbor_set_v1.size(); j++)
		{
			VERTEX_INDEX v2 = neighbor_set_v1[j];
			bool flag_collinear = are_collinear(iv,v1,v2, scalar_grid);
			if(flag_collinear)
			{
				if(lambda_computations
					(iv,v1,v2, scalar_grid, gradient_grid, grad_mag_grid) > alpha)
				{
					if(debug)
						cout <<" **failNo** ";
					m++;
				}
			}
		}
	}
	if(m>1)
		return false;
	else
		return true;
}

/**
* Compute phi version B.
* @param distance from vertex iv,
*/
void compute_phi_k_B(
	const VERTEX_INDEX iv, 
	const VERTEX_INDEX v1,
	const int k,
	GRADIENT_COORD_TYPE *p,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const GRADIENT_MAGNITUDE_GRID & grad_mag_grid)
{
	GRADIENT_COORD_TYPE giv[DIM3], ngiv[DIM3], gv1[DIM3],
		ngv1[DIM3];
	const GRADIENT_COORD_TYPE magiv = grad_mag_grid.Scalar(iv);
	const GRADIENT_COORD_TYPE magv1 = grad_mag_grid.Scalar(v1);

	for (int d = 0; d < DIM3; d++)
	{
		giv[d]  = magiv*gradient_grid.Vector(iv,d);
		//ngiv[d] = gradient_grid.Vector(iv,d);

		gv1[d]  = magv1*gradient_grid.Vector(v1,d);
		//ngv1[d] = gradient_grid.Vector(v1,d);
	}
	IJK::normalize_vector(DIM3,giv,0.001,ngiv);
	IJK::normalize_vector(DIM3,gv1,0.001,ngv1);

	for (int t = 0; t < k; t++)
	{

		GRADIENT_COORD_TYPE tempg[DIM3]={0.0,0.0,0.0};
		float inn_pdt = 0.0;
		IJK::compute_inner_product(DIM3, ngiv, ngv1, inn_pdt);
		IJK::multiply_coord_3D(2.0*inn_pdt, ngv1, tempg);

		IJK::subtract_coord_3D(tempg, ngiv, p);

		IJK::copy_coord_3D(ngv1,ngiv);
		IJK::copy_coord_3D(p, ngv1);

		IJK::normalize_vector(DIM3,ngiv,0.001,ngiv);
		IJK::normalize_vector(DIM3,ngv1,0.001,ngv1);
	}
}

/**
* Lambda computations B
* @param k Using vertices at distance k from iv
* @return An angle in degrees between phi k and grad v2
**/
float lambda_computationsB(
	const VERTEX_INDEX iv, 
	const VERTEX_INDEX v1, 
	const VERTEX_INDEX v2,
	const int k,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const GRADIENT_MAGNITUDE_GRID & grad_mag_grid)
{
	GRADIENT_COORD_TYPE p[DIM3], gv2[DIM3];

	compute_phi_k_B(iv, v1, k, &(p[0]),
		scalar_grid, gradient_grid, grad_mag_grid);

	for (int d = 0; d < DIM3; d++)
	{
		float x = grad_mag_grid.Scalar(v2)* gradient_grid.Vector(v2,d);
		gv2[d] = x;
	}

	IJK::normalize_vector(DIM3,gv2, 0.001, gv2);
	IJK::normalize_vector(DIM3, p, 0.001, p);
	float inn_pdt = 0.0;
	IJK::compute_inner_product(DIM3, gv2, p, inn_pdt);
	float deg = acos (inn_pdt) * 180.0 / M_PI;
	return deg;
}

/*
* predict gradient at vertex [predict_at_v], using gradients at vertices, 
* using_vertex_v1, 
* and, using_vertex_v2
*/
float lambda_computationsC(
	const VERTEX_INDEX predict_at_v,
	const VERTEX_INDEX using_vertex_v1,
	const VERTEX_INDEX using_vertex_v2,
	const int predict_at_distance,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const GRADIENT_MAGNITUDE_GRID & grad_mag_grid
	)
{
	GRADIENT_COORD_TYPE p[DIM3], gradient_predict_at_v[DIM3];
	compute_phi_k_B(using_vertex_v1, using_vertex_v2,
		predict_at_distance, p, scalar_grid, gradient_grid, grad_mag_grid);
	for (int i = 0; i < DIM3; i++)
	{
		gradient_predict_at_v[i]= gradient_grid.Vector(predict_at_v,i);
	}
	IJK::normalize_vector(DIM3, gradient_predict_at_v, 0.0001, gradient_predict_at_v);
	IJK::normalize_vector(DIM3, p, 0.0001, p);
	float inn_pdt = 0.0;
	IJK::compute_inner_product(DIM3, gradient_predict_at_v, p, inn_pdt);
	float deg = acos (inn_pdt) * 180.0 / M_PI;
	return deg;
}

/*
* Compute v2 from v1
* Checks for boundary conditions
*/
void compute_v2(
	const VERTEX_INDEX iv,
	const VERTEX_INDEX v1,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	VERTEX_INDEX & v2
	)
{
	COORD_TYPE civ[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(iv, civ);

	COORD_TYPE cv1[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(v1, cv1);

	COORD_TYPE cdiff[DIM3] = {0.0,0.0,0.0};
	IJK::subtract_coord_3D(cv1,civ,cdiff);

	COORD_TYPE cv2[DIM3] = {0.0,0.0,0.0};
	IJK::add_coord_3D(cdiff,cv1,cv2);
	//

	for (int d = 0; d < DIM3; d++)
	{
		if(cv2[d] < 0)
			cv2[d]=0;
		if(cv2[d] > (scalar_grid.AxisSize(d)-1))
		{
			cv2[d] = scalar_grid.AxisSize(d)-1;
		}
	}
	v2 = scalar_grid.ComputeVertexIndex(cv2);


	if(debug){
		using namespace std;
		COORD_TYPE c[DIM3]={0.0,0.0,0.0};
		scalar_grid.ComputeCoord(v2,c);
		cout <<"neighbor of {"<<v1<<"}, is "<< v2 
			<<" coordv2 {"<<c[0]<<" "<<c[1]<<" "<<c[2]<<"} "
			<<" diff"<<cdiff[0]<<","<<cdiff[1]<<","<<cdiff[2]<<"\n";
	}
	//DEBUG end 
}
/*
* Tangent neighbor based computations
* tangent_neighbor_set
*/
bool Nt_based_computations_for_iv(
	const VERTEX_INDEX iv,
	float alpha,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	vector <VERTEX_INDEX> &tangent_neighbor_set
	)
{
	for (int v = 0; v < tangent_neighbor_set.size(); v++)
	{
		VERTEX_INDEX v1 = tangent_neighbor_set[v];
		VERTEX_INDEX v2 = 0;
		compute_v2(iv, v1, scalar_grid, v2);
		float deg =  lambda_computations
			(iv,v1,v2, scalar_grid, 
			gradient_grid, grad_mag_grid); 
		if( deg > alpha)
		{
			if(debug)
				cout <<" **fail** ";
			return false;
		}	
	}
	return true;
}


/*
* EXTENDED Tangent neighbor based computations
* Sets vertices to be true individually. 
* Returns a vector of correct gradients. 
*/
void Extended_Nt_based_computations_for_iv(
	const VERTEX_INDEX iv,
	float alpha,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	vector <VERTEX_INDEX> &tangent_neighbor_set,
	vector <VERTEX_INDEX> &extended_reliable_vertices
	)
{
	for (int v = 0; v < tangent_neighbor_set.size(); v++)
	{
		VERTEX_INDEX v1 = tangent_neighbor_set[v];
		VERTEX_INDEX v2 = 0;
		compute_v2(iv, v1, scalar_grid, v2);	
		// check if v1 and v2 are correct 
		if(reliable_grid.Scalar(v1) && reliable_grid.Scalar(v2))
			if(lambda_computations
				(iv,v1,v2, scalar_grid, 
				gradient_grid, grad_mag_grid) < alpha)
			{
				extended_reliable_vertices.push_back(iv);
			}	
	}

}



/*
* Orthogonal neighbor based computations
*/
bool No_based_computations_for_iv(
	const VERTEX_INDEX iv,
	float alpha,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	vector <VERTEX_INDEX> &ortho_neighbor_set)
{
	int m=0;
	//if ( ortho_neighbor_set.size() == 0 )
	//	return false; 

	for (int i = 0; i < ortho_neighbor_set.size(); i++)
	{
		VERTEX_INDEX v1 = ortho_neighbor_set[i];
		VERTEX_INDEX v2 = 0;
		compute_v2(iv,v1,scalar_grid,v2);	
		if(lambda_computations
			(iv,v1,v2, scalar_grid, gradient_grid, grad_mag_grid) > alpha)
		{
			if(debug)
				cout <<" **failNo** ";
			m++;
		}
	}
	if(m>1)
		return false;
	else
		return true;
}

/// OLD CODE NOT BEING CALLED 
void extended_curv_per_vertex
	(const VERTEX_INDEX iv,
	const bool increasing_iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const RELIGRADIENT::BOOL_GRID &boundary_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	vector<VERTEX_INDEX> & vertex_index_of_extended_correct_grads,
	INPUT_INFO & io_info
	)
{
	//set coord iv
	COORD_TYPE coord_iv[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(iv, coord_iv);
	float alpha = io_info.param_angle;

	GRADIENT_COORD_TYPE grad_iv[DIM3];  // unscaled gradient vector
	GRADIENT_COORD_TYPE normalized_grad_iv[DIM3]; // normalized grad_iv

	GRADIENT_COORD_TYPE mag_iv = grad_mag_grid.Scalar(iv);
	if (mag_iv > io_info.min_gradient_mag) {

		// grad_iv[] is an unscaled gradient
		for (int d = 0; d < DIM3; d++) 
		{
			COORD_TYPE spacing_d = scalar_grid.Spacing(d);
			grad_iv[d] = mag_iv*spacing_d*gradient_grid.Vector(iv,d);
		}

		IJK::normalize_vector
			(DIM3, grad_iv, io_info.min_gradient_mag, normalized_grad_iv);
		// Tangent neighbor set
		vector<VERTEX_INDEX> tangent_neighbor_set_iv;
		vector <VERTEX_INDEX> ortho_neighbor_set_iv;
		compute_Nt_and_No_iv (iv, grad_iv, scalar_grid, tangent_neighbor_set_iv,
			ortho_neighbor_set_iv, io_info);
		float alpha = 5;
		// alphs is half the angle parameter
		alpha = io_info.param_angle;
		Extended_Nt_based_computations_for_iv(iv, alpha, scalar_grid,
			gradient_grid, grad_mag_grid, reliable_grid, tangent_neighbor_set_iv, vertex_index_of_extended_correct_grads);
	}
}
/** 
* \Brief EXTENDED Curvature based computations for vertex iv, versionB.
* Check reliablity at iv using vertices distance two from v.
* Iv is not a boundary vertex.
* @param iv
*
**/
void extended_curvature_based_B_per_vertex(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> & extended_reliable,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{
	GRADIENT_COORD_TYPE magiv = grad_mag_grid.Scalar(iv); 
	GRADIENT_COORD_TYPE gradiv[DIM3], normgradiv[DIM3];
	if(magiv > io_info.min_gradient_mag){
		for (int d = 0; d < DIM3; d++) 
		{
			gradiv[d] = magiv*gradient_grid.Vector(iv,d);
			
		}
		IJK::normalize_vector
			(DIM3, gradiv, io_info.min_gradient_mag, normgradiv);
		vector <VERTEX_INDEX> tangent_neighbor_set_iv;
		vector <VERTEX_INDEX> ortho_neighbor_set_iv;
		compute_Nt_and_No_iv (iv, normgradiv, scalar_grid, tangent_neighbor_set_iv,
			ortho_neighbor_set_iv, io_info);

		float alpha = 5;
		alpha = io_info.param_angle;
		for (int i = 0; i < tangent_neighbor_set_iv.size(); i++)
		{
			VERTEX_INDEX v1 = tangent_neighbor_set_iv[i];
			COORD_TYPE diff[DIM3];
			COORD_TYPE v1c[DIM3], v2c[DIM3], v3c[DIM3], ivc[DIM3];
			scalar_grid.ComputeCoord(iv,ivc);
			scalar_grid.ComputeCoord(v1,v1c);
			IJK::subtract_coord_3D(v1c, ivc, diff);
			IJK::add_coord_3D(v1c, diff, v2c);
			IJK::add_coord_3D(v2c, diff, v3c);

			//boundary
			for (int d = 0; d < DIM3; d++)
			{
				if ( v3c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v3c[d] < 0)
					return;
				if(v2c[d] < 0)
					return;
			}
			int v2 = scalar_grid.ComputeVertexIndex(v2c); 
			int v3 = scalar_grid.ComputeVertexIndex(v3c);

			//extended version
			if(reliable_grid.Scalar(v1) && reliable_grid.Scalar(v2))
			{

				int k=1;
				float deg = lambda_computationsC(iv, v2, v1, k, scalar_grid, 
					gradient_grid, grad_mag_grid);
				if( deg < alpha)
				{
					if(io_info.cdist==2)
					{
						if(reliable_grid.Scalar(v3)){
							k=2;
							float deg = lambda_computationsC(iv, v3, v2, k, scalar_grid, 
								gradient_grid, grad_mag_grid);
							if( deg < alpha)
							{
								extended_reliable.push_back(iv);
							}
						}
					}
					else
					{
						extended_reliable.push_back(iv);
					}
				}
			}

		}
	}
}


/* Compute all neighbors using  boundary checking
*/
void compute_all_neighbors_iv
	(const VERTEX_INDEX iv, GRADIENT_COORD_TYPE *norm_grad_iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	vector <VERTEX_INDEX> &all_neighbor_set, INPUT_INFO & io_info)
{
	//insert_flag = update_tangent_vertex_list
	//(next_vertex, inn_pdt, degree_param, tangent_vertex_list);
	COORD_TYPE coord_iv[DIM3] = {0.0,0.0,0.0};
	scalar_grid.ComputeCoord(iv, coord_iv);
	for (int d = 0; d < DIM3; d++) 
	{
		int k = min(1, int(coord_iv[d]));
		VERTEX_INDEX prev_vertex = 
			iv - k * scalar_grid.AxisIncrement(d);
		all_neighbor_set.push_back (prev_vertex);
		k = min(1,(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));
		VERTEX_INDEX next_vertex = 
			iv + k * scalar_grid.AxisIncrement(d);
		all_neighbor_set.push_back(next_vertex);
	}
}
/*
* Run algorithim 1 or 2 based on cdist 
*/
void algo12_per_vertex(
	const VERTEX_INDEX iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{
	GRADIENT_COORD_TYPE magiv = grad_mag_grid.Scalar(iv); 
	GRADIENT_COORD_TYPE gradiv[DIM3], normgradiv[DIM3];
	if(magiv > io_info.min_gradient_mag){
		for (int d = 0; d < DIM3; d++) 
		{
			gradiv[d] = magiv*gradient_grid.Vector(iv,d);
		}
		IJK::normalize_vector
			(DIM3, gradiv, io_info.min_gradient_mag, normgradiv);

		vector <VERTEX_INDEX> neighbor_set_iv;
		// compute all neighbors (with  bondary checking)
		compute_all_neighbors_iv(iv, normgradiv, scalar_grid,neighbor_set_iv, io_info);
		for (int i = 0; i < neighbor_set_iv.size(); i++)
		{
			VERTEX_INDEX v1 = neighbor_set_iv[i];
			COORD_TYPE diff[DIM3];
			COORD_TYPE v1c[DIM3], v2c[DIM3], v3c[DIM3], ivc[DIM3];
			scalar_grid.ComputeCoord(iv,ivc);
			scalar_grid.ComputeCoord(v1,v1c);
			IJK::subtract_coord_3D(v1c, ivc, diff);
			IJK::add_coord_3D(v1c, diff, v2c);
			IJK::add_coord_3D(v2c, diff, v3c);

			//boundary
			for (int d = 0; d < DIM3; d++)
			{
				if ( v3c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v3c[d] < 0)
					return;
				if(v2c[d] < 0 )
					return;
			}
			int v2 = scalar_grid.ComputeVertexIndex(v2c); 
			int v3 = scalar_grid.ComputeVertexIndex(v3c);


			float alpha = 5;
			alpha = io_info.param_angle;
			if(io_info.cdist == 1) 
			{
				// algorithm 1
				float deg = lambda_computationsC(iv, v2, v1, 1, scalar_grid, 
					gradient_grid, grad_mag_grid);
				if( deg > alpha)
				{
					return;
				}
			}
			else {
				//algorithm 2
				float deg1 = lambda_computationsC(iv, v2, v1, 1, scalar_grid, 
					gradient_grid, grad_mag_grid);
				float deg2 = lambda_computationsC(iv, v3, v2, 2, scalar_grid, 
					gradient_grid, grad_mag_grid);
				if( deg1 >= alpha || deg2 >=alpha)
				{
					return;
				}
			}
		}
		reliable_grid.Set(iv, true);
	}
}



/** 
* \Brief Curvature based computations for vertex iv, versionB.
* Check reliablity at iv using vertices distance two from v.
* Iv is not a boundary vertex.
* @param iv
* This is the key function being called. 
**/
void curvature_based_B_per_vertex(
	const VERTEX_INDEX iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{
 	GRADIENT_COORD_TYPE magiv = grad_mag_grid.Scalar(iv); 
	GRADIENT_COORD_TYPE gradiv[DIM3], normgradiv[DIM3];
	if(magiv > io_info.min_gradient_mag){
		for (int d = 0; d < DIM3; d++) 
		{
			gradiv[d] = magiv*gradient_grid.Vector(iv,d);
		}
		IJK::normalize_vector
			(DIM3, gradiv, io_info.min_gradient_mag, normgradiv);

		vector <VERTEX_INDEX> tangent_neighbor_set_iv;
		vector <VERTEX_INDEX> ortho_neighbor_set_iv;
		compute_Nt_and_No_iv (iv, normgradiv, scalar_grid, tangent_neighbor_set_iv,
			ortho_neighbor_set_iv, io_info);

		for (int i = 0; i < tangent_neighbor_set_iv.size(); i++)
		{
			VERTEX_INDEX v1 = tangent_neighbor_set_iv[i];
			COORD_TYPE diff[DIM3];
			COORD_TYPE v1c[DIM3], v2c[DIM3], v3c[DIM3], ivc[DIM3];
			scalar_grid.ComputeCoord(iv,ivc);
			scalar_grid.ComputeCoord(v1,v1c);
			IJK::subtract_coord_3D(v1c, ivc, diff);
			IJK::add_coord_3D(v1c, diff, v2c);
			IJK::add_coord_3D(v2c, diff, v3c);

			//boundary
			for (int d = 0; d < DIM3; d++)
			{
				if ( v3c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
					return;
				if ( v3c[d] < 0)
					return;
				if(v2c[d] < 0 )
					return;
			}
			int v2 = scalar_grid.ComputeVertexIndex(v2c); 
			int v3 = scalar_grid.ComputeVertexIndex(v3c);


			float alpha = 5;
			alpha = io_info.param_angle;
			int k=1;

			float deg = lambda_computationsC(iv, v2, v1, k, scalar_grid, 
				gradient_grid, grad_mag_grid);
			if( deg > alpha)
			{
				return;
			}
			if(io_info.cdist==2){
				k=2;
				float deg2 = lambda_computationsC(iv, v3, v2, k, scalar_grid, 
					gradient_grid, grad_mag_grid);

				if( deg2 > alpha)
				{
					return;
				}
			}
		}

		//ortho match 
		if(ortho_neighbor_set_iv.size() > 0 )
		{
			// DEBUG
			
			bool oMatch = does_ortho_match(
				iv, ortho_neighbor_set_iv, scalar_grid,gradient_grid, 
				grad_mag_grid,io_info);
			
			if(oMatch)
			{
				reliable_grid.Set(iv, true);
			}
			
			/*
			if (does_ortho_matchA(
				iv, ortho_neighbor_set_iv, scalar_grid,gradient_grid, 
				grad_mag_grid,io_info))
				reliable_grid.Set(iv, true);
			else if(does_ortho_matchB(
				iv, ortho_neighbor_set_iv, scalar_grid,gradient_grid, 
				grad_mag_grid,io_info))
				reliable_grid.Set(iv, true);
			else
				return;
			*/
		}
		else
		{
				reliable_grid.Set(iv, true);
		}
		
	}
}

/**
* float angle_between (iv,v1)
**/
float angle_between(const VERTEX_INDEX v,
					const VERTEX_INDEX v2,
					const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
					const GRADIENT_GRID & gradient_grid,
					const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
					INPUT_INFO & io_info
					)
{
	GRADIENT_COORD_TYPE gv[DIM3], gv2[DIM3];

	for (int d = 0; d < DIM3; d++)
	{
		gv2[d] = grad_mag_grid.Scalar(v2)* gradient_grid.Vector(v2,d);
		gv[d] = grad_mag_grid.Scalar(v)*gradient_grid.Vector(v,d);
	}
	IJK::normalize_vector(DIM3, gv2, io_info.min_gradient_mag, gv2);
	IJK::normalize_vector(DIM3, gv, io_info.min_gradient_mag, gv);
	float innpdt = 0.0;
	IJK::compute_inner_product_3D(gv, gv2, innpdt);
	/*

	*/
	float deg = acos (innpdt) * 180.0 / M_PI;
	return  deg;

}
			
//debug
static int x11=0,x22=0,x33=0,x44=0,x55=0,x66=0;


bool does_ortho_matchA
	(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info)
{
	const int osize = ortho_neighbor_set.size();
	short numAgree=0;
	for (VERTEX_INDEX i = 0; i < osize; i++)
	{
		VERTEX_INDEX v1 = ortho_neighbor_set[i];
		COORD_TYPE diff[DIM3];
		COORD_TYPE v1c[DIM3], v2c[DIM3], v3c[DIM3], ivc[DIM3];

		scalar_grid.ComputeCoord(iv,ivc);
		scalar_grid.ComputeCoord(v1,v1c);
		IJK::subtract_coord_3D(v1c, ivc, diff);
		IJK::add_coord_3D(v1c, diff, v2c);
		IJK::add_coord_3D(v2c, diff, v3c);

		//boundary
		for (int d = 0; d < DIM3; d++)
		{
			if ( v3c[d] > (scalar_grid.AxisSize(d) - 1))
				return false;
			if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
				return  false;
			if ( v3c[d] < 0)
				return false;
			if(v2c[d] < 0 )
				return false;
		}
		VERTEX_INDEX v2 = scalar_grid.ComputeVertexIndex(v2c); 
		VERTEX_INDEX v3 = scalar_grid.ComputeVertexIndex(v3c);

		float alpha = 5, alpha1=5;
		alpha = io_info.param_angle;
		float angle_iv_v1 = angle_between(iv, v1, scalar_grid, gradient_grid,
			grad_mag_grid, io_info);
		float angle_iv_v2 = angle_between(iv, v2, scalar_grid, gradient_grid,
			grad_mag_grid, io_info); 
		if(angle_iv_v1 <= alpha1 && 0)
		{
			x44++;
			numAgree++;
		}
		else
		{
			if( (lambda_computationsC(iv, v2, v1, 1, scalar_grid, 
				gradient_grid, grad_mag_grid) <= alpha) && 
				(lambda_computationsC 
				(iv, v3, v2, 2, scalar_grid, gradient_grid, grad_mag_grid)
				<= alpha))
			{
				x55++;
				numAgree++;
			}
		}
	}
	if(numAgree >= 2)
	{
		x66++;
		return true;
	}
	else
		return false;
}

bool does_ortho_matchB
	(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info)
{
	const int osize = ortho_neighbor_set.size();
	for (VERTEX_INDEX i = 0; i < osize; i++)
	{
		VERTEX_INDEX v1 = ortho_neighbor_set[i];
		COORD_TYPE diff[DIM3];
		COORD_TYPE v1c[DIM3], v2c[DIM3], ivc[DIM3];

		scalar_grid.ComputeCoord(iv,ivc);
		scalar_grid.ComputeCoord(v1,v1c);
		IJK::subtract_coord_3D(v1c, ivc, diff);
		IJK::add_coord_3D(v1c, diff, v2c);

		//boundary
		for (int d = 0; d < DIM3; d++)
		{
			if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
				return  false;
			if(v2c[d] < 0 )
				return false;
		}
		VERTEX_INDEX v2 = scalar_grid.ComputeVertexIndex(v2c); 
		float alpha = 15, alpha1=15;
		alpha = io_info.param_angle;
			float angle_iv_v1 = angle_between(iv, v1, scalar_grid, gradient_grid,
			grad_mag_grid, io_info);
		float angle_iv_v2 = angle_between(iv, v2, scalar_grid, gradient_grid,
			grad_mag_grid, io_info); 
		if(angle_iv_v1 <= alpha1 && angle_iv_v2 <= alpha1)
		{
			COORD_TYPE v2diff[DIM3];
			IJK::subtract_coord_3D(v2c, ivc, v2diff);
			float angle_v2a,angle_v2b;
			GRADIENT_COORD_TYPE gv2[DIM3];
			for (int d = 0; d < DIM3; d++)
			{
				gv2[d] = grad_mag_grid.Scalar(v2)*gradient_grid.Vector(v2,d);
			}
			IJK::normalize_vector(DIM3, gv2, io_info.min_gradient_mag, gv2);
			IJK::normalize_vector(DIM3, v2diff, io_info.min_gradient_mag,v2diff);
			float innpdt=0.0f;
			IJK::compute_inner_product_3D(gv2, v2diff, innpdt);
			angle_v2a =  acos (innpdt) * 180.0 / M_PI;
			IJK::subtract_coord_3D(ivc, v2c, v2diff);
			IJK::normalize_vector(DIM3, v2diff, io_info.min_gradient_mag, v2diff);
			IJK::compute_inner_product_3D(gv2, v2diff, innpdt);
			angle_v2b =  acos (innpdt) * 180.0 / M_PI;

			if(angle_v2a < 20 || angle_v2b < 20)
			{
				return true;
			}
		}
	}
	return false;
}

/**
* Does orthogonal directions match
**/

bool does_ortho_match(
	const VERTEX_INDEX iv,
	vector <VERTEX_INDEX> &ortho_neighbor_set,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	INPUT_INFO & io_info)
{
	const int osize = ortho_neighbor_set.size();
	
	for (VERTEX_INDEX i = 0; i < osize; i++)
	{
		VERTEX_INDEX v1 = ortho_neighbor_set[i];
		COORD_TYPE diff[DIM3];
		COORD_TYPE v1c[DIM3], v2c[DIM3], v3c[DIM3], ivc[DIM3];
		scalar_grid.ComputeCoord(iv,ivc);
		scalar_grid.ComputeCoord(v1,v1c);
		IJK::subtract_coord_3D(v1c, ivc, diff);
		IJK::add_coord_3D(v1c, diff, v2c);
		IJK::add_coord_3D(v2c, diff, v3c);
		//boundary
		for (int d = 0; d < DIM3; d++)
		{
			if ( v3c[d] > (scalar_grid.AxisSize(d) - 1))
				return false;
			if ( v2c[d] > (scalar_grid.AxisSize(d) - 1))
				return  false;
			if ( v3c[d] < 0)
				return false;
			if(v2c[d] < 0 )
				return false;
		}
		int v2 = scalar_grid.ComputeVertexIndex(v2c); 
		int v3 = scalar_grid.ComputeVertexIndex(v3c);

		float alpha = 5;
		alpha = io_info.param_angle;
		int k=1;
		if(lambda_computationsC(iv, v2, v1, 1, scalar_grid, 
			gradient_grid, grad_mag_grid) < alpha)
		{
			if (io_info.cdist == 2)
			{
				k=2;
				if(lambda_computationsC \
					(iv, v3, v2, 2, scalar_grid, gradient_grid, grad_mag_grid)
					<= alpha)
				{
					x11++;
					return true;
				}
				float angle = angle_between(iv, v2, scalar_grid, gradient_grid,
					grad_mag_grid, io_info); 
				if( angle < 5)
				{
					x22++;
					return true;
				}
			}
			else
			{
				return true;
			}
		}
		/*
		debug
		float angle = angle_between(iv,v1,scalar_grid, gradient_grid,
			grad_mag_grid, io_info);
		if( angle < 5)
		{
			x33++;
			return true;
		}

		*/
		float anglev2 = angle_between(iv, v2, scalar_grid, gradient_grid,
					grad_mag_grid, io_info);
		float anglev1 = angle_between(iv,v1,scalar_grid, gradient_grid,
			grad_mag_grid, io_info);

		if( anglev1 <= 5)
		{
			x33++;
			return true;
		}
	}
	return false;
}



/** 
* Curvature based computations for vertex iv, 
* @param iv
*
* Called by compute_reliable_gradients_curvature_based.
**/
void curvature_based_per_vertex(
	const VERTEX_INDEX iv,
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{		
	GRADIENT_COORD_TYPE grad_iv[DIM3];  // unscaled gradient vector
	GRADIENT_COORD_TYPE normalized_grad_iv[DIM3]; // normalized grad_iv

	GRADIENT_COORD_TYPE mag_iv = grad_mag_grid.Scalar(iv);
	if (mag_iv > io_info.min_gradient_mag) {

		// grad_iv[] is an unscaled gradient
		// scaling
		for (int d = 0; d < DIM3; d++) 
		{
			COORD_TYPE spacing_d = scalar_grid.Spacing(d);
			grad_iv[d] = mag_iv*spacing_d*gradient_grid.Vector(iv,d);
		}

		IJK::normalize_vector
			(DIM3, grad_iv, io_info.min_gradient_mag, normalized_grad_iv);

		// Tangent neighbor set
		vector <VERTEX_INDEX> tangent_neighbor_set_iv;
		vector <VERTEX_INDEX> ortho_neighbor_set_iv;
		compute_Nt_and_No_iv (iv, grad_iv, scalar_grid, tangent_neighbor_set_iv,
			ortho_neighbor_set_iv, io_info);

		//DEBUG
		if(debug || iv == test_iv){
			cout <<"debug "<< debug ;
			using namespace std;
			COORD_TYPE c[DIM3]={0.0,0.0,0.0};
			scalar_grid.ComputeCoord(iv,c);
			cout <<"cube index{"<<iv<<"}"
				<<" coord{"<<c[0]<<" "<<c[1]<<" "<<c[2]<<"}\n"
				<<" normalized grad "<<normalized_grad_iv[0]<<","<< normalized_grad_iv[1]<<","<<normalized_grad_iv[2]
			<<" tangent neighbors {\n";

			for (int k = 0; k < tangent_neighbor_set_iv.size(); k++)
			{
				cout <<tangent_neighbor_set_iv[k]<<"[";
				scalar_grid.ComputeCoord(tangent_neighbor_set_iv[k],c);
				cout <<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
			}
			cout <<"}"<<endl;
			cout 	<<" ortho neighbors {\n";
			for (int k = 0; k < ortho_neighbor_set_iv.size(); k++)
			{
				cout <<ortho_neighbor_set_iv[k]<<"[";
				scalar_grid.ComputeCoord(ortho_neighbor_set_iv[k],c);
				cout <<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
			}
			cout <<"}"<<endl;
		}//DEBUG_END


		float alpha = 5;
		alpha = io_info.param_angle;

		bool flag_nt = Nt_based_computations_for_iv(iv, alpha, scalar_grid,
			gradient_grid, grad_mag_grid, tangent_neighbor_set_iv);
		bool flag_no = No_based_computations_for_iv(iv, alpha, scalar_grid,
			gradient_grid, grad_mag_grid, ortho_neighbor_set_iv);

		//DEBUG
		if(debug || iv == test_iv)
			cout <<"is reliable: Nt "<< flag_nt <<" No "<< flag_no<<endl;
		//DEBUG_END


		//bool is_rel = (flag_nt && flag_no);
		//reliable_grid.Set(iv, is_rel );

		if ( flag_nt ) 
			reliable_grid.Set(iv, true);
		if( ortho_neighbor_set_iv.size() > 0 ) 
		{
			if (flag_no)
				reliable_grid.Set(iv, true);
		}


		if(reliable_grid.Scalar(iv))
		{
			io_info.out_info.num_reliable++;
		}
		else
			io_info.out_info.num_unreliable++;

	}

}
/*
* Running algorithim 1 or 2 depending on cdist. 
*/
void compute_reliable_gradients_curvature_based_algo12(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const RELIGRADIENT::BOOL_GRID &boundary_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{
	//debug
	cerr <<"Running algorithm : " << io_info.cdist <<".\n";
	const int numVertices = scalar_grid.NumVertices();
	for (int iv = 0; iv < numVertices; iv++)
	{
		if(!boundary_grid.Scalar(iv))
		{
			 algo12_per_vertex(iv, scalar_grid, gradient_grid,
				grad_mag_grid, reliable_grid, io_info);
		}
	}
}
/*
*Curvature based reliable gradients computations B
*Check reliablity at vertex iv using vertices distance 2 from v.
*/
void compute_reliable_gradients_curvature_basedB(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const RELIGRADIENT::BOOL_GRID &boundary_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	)
{
	const int numVertices = scalar_grid.NumVertices();
	for (int iv = 0; iv < numVertices; iv++)
	{
		if(!boundary_grid.Scalar(iv))
		{
			curvature_based_B_per_vertex(iv, scalar_grid, gradient_grid,
				grad_mag_grid, reliable_grid, io_info);
		}
	}
	cerr <<"x11 "<< x11 <<", x22 "<< x22 <<", x33 "<<x33 <<endl;
	cerr <<"x44 "<<x44 <<",x55 "<< x55<<",x66 "<<x66 << endl;
	cerr <<"Starting Extended Correct Gradients." << endl;
	// check if extended
	if (io_info.extended_curv_based)
	{
		for (int i = 0; i < io_info.extend_max; i++)
		{
			vector<VERTEX_INDEX>  vertex_index_of_extended_correct_grads;
			//reset num unreliables. 
			io_info.out_info.num_unreliable = 0;
			io_info.out_info.num_reliable = 0;

			for (int iv = 0; iv < numVertices; iv++)
			{
				if(!boundary_grid.Scalar(iv))
				{
					extended_curvature_based_B_per_vertex(iv,  vertex_index_of_extended_correct_grads,
						scalar_grid, gradient_grid, grad_mag_grid, reliable_grid, io_info);
				}
			}

			for (int i = 0; i < vertex_index_of_extended_correct_grads.size(); i++)
			{
				reliable_grid.Set( vertex_index_of_extended_correct_grads[i], true);
			}
		}
	}
}

//
// Curvature based reliable gradients computations
//
void compute_reliable_gradients_curvature_based(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const RELIGRADIENT::BOOL_GRID &boundary_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	) 
{
	cout <<"Calling Function  " <<__FUNCTION__<< endl;
	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) 
	{
		if(!boundary_grid.Scalar(iv))
			curvature_based_per_vertex (iv, scalar_grid, gradient_grid,
			grad_mag_grid, reliable_grid, io_info);
	}
	cout <<"Calling Function  " <<__FUNCTION__<< endl;
}

// Extended version of curvature based reliable gradients computations
// Returns vertex indices of extended correct grads
void compute_reliable_gradients_extended_curvature_based(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const RELIGRADIENT::BOOL_GRID &boundary_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	vector<VERTEX_INDEX> & vertex_index_of_extended_correct_grads,
	INPUT_INFO & io_info)
{
	cout <<"Calling Function " <<__FUNCTION__<< endl;
	//increasing order
	bool increasing_indices = true;
	VERTEX_INDEX total = scalar_grid.NumVertices();
	for (VERTEX_INDEX iv = 0; iv < total; iv++)
	{
		if(!boundary_grid.Scalar(iv))
			extended_curv_per_vertex(iv, increasing_indices, scalar_grid, boundary_grid, gradient_grid,
			grad_mag_grid, reliable_grid, vertex_index_of_extended_correct_grads, io_info);
	}
}
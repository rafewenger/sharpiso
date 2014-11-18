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


// local type definition
namespace {

	typedef IJK::BOOL_GRID_BASE<RELIGRADIENT_GRID> BOOL_GRID_BASE;
	typedef IJK::BOOL_GRID<RELIGRADIENT_GRID> BOOL_GRID;

};

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
	if(inn_pdt > 0.0)
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

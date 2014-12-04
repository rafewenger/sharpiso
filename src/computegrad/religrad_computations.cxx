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






	//DEBUG
	const bool debug = false;
	/*
	* check if neighbor v of iv can be put into the
	* tangent neighbor set or orthogonal set.
	*/
	void check_v_for_insert2_Nt_and_No_iv(
		const VERTEX_INDEX iv, 
		GRADIENT_COORD_TYPE *grad_iv,
		const VERTEX_INDEX v, 
		const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
		vector <VERTEX_INDEX> &tangent_neighbor_set,
		vector <VERTEX_INDEX> &ortho_neighbor_set, INPUT_INFO & io_info)
	{
		COORD_TYPE edge_vec[DIM3] = {0, 0, 0};
		float degree_param = io_info.neighbor_angle_parameter*M_PI/180.0;
		compute_vector_between_grid_vertex
			(scalar_grid, iv, v, &(edge_vec[0]));
		IJK::normalize_vector(DIM3, &(edge_vec[0]), 0.0, &(edge_vec[0]));
		float inn_pdt = 0.0;
		IJK::compute_inner_product(DIM3, grad_iv, edge_vec, inn_pdt);
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
	*/
	void compute_Nt_and_No_iv
		(const VERTEX_INDEX iv, GRADIENT_COORD_TYPE *grad_iv,
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
			check_v_for_insert2_Nt_and_No_iv(iv, grad_iv, prev_vertex, scalar_grid, 
				tangent_neighbor_set, ortho_neighbor_set, io_info);

			k = min(1,(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));
			VERTEX_INDEX next_vertex = 
				iv + k * scalar_grid.AxisIncrement(d);
			check_v_for_insert2_Nt_and_No_iv(iv, grad_iv, next_vertex, scalar_grid, 
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


	/*
	* compute phi between v and v1
	*/
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
			phi_vector[d]=phi_vector[d]+2*(innpdt)*gradv1[d]-gradv[d];
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
	*/
	float lambda_computations(const VERTEX_INDEX iv, 
		const VERTEX_INDEX v1, const VERTEX_INDEX v2,
		const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID & gradient_grid,
		const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid)
	{
		float  phi_vector[DIM3]={0.0,0.0,0.0};
		compute_phi(iv, v1, scalar_grid, gradient_grid, grad_mag_grid,
			&(phi_vector[0]));

		GRADIENT_COORD_TYPE magv2 = grad_mag_grid.Scalar(v2);
		COORD_TYPE gradv2[DIM3]={0.0,0.0,0.0};
		for (int d = 0; d < DIM3; d++) 
		{
			COORD_TYPE spacing_d = scalar_grid.Spacing(d);
			gradv2[d] = magv2*spacing_d*gradient_grid.Vector(v1,d);
		}
		IJK::normalize_vector
			(DIM3, gradv2, 0.0, gradv2);
		
		//
		float inn_pdt = 0.0;
		IJK::compute_inner_product(DIM3, gradv2, phi_vector, inn_pdt);
		//DEBUG
		if(debug){
		cout <<"gradV2 "<< gradv2[0]<<","<<gradv2[1]<<","<<gradv2[2]
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
	/*
	* compute v2 from v1
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
			compute_v2(iv,v1,scalar_grid,v2);	
			if(lambda_computations
				(iv,v1,v2, scalar_grid, 
				gradient_grid, grad_mag_grid) > alpha)
			{
				if(debug)
					cout <<" **fail** ";
				return false;
			}	
		}
		return true;
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

// Curvature based computations for vertex iv, 
// @param iv , grid vertex
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
	// only run the test if gradient at vertex iv is reliable
	if (reliable_grid.Scalar(iv)) {
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

		
			// tangent neighbor set
			vector<VERTEX_INDEX> tangent_neighbor_set_iv;
			vector <VERTEX_INDEX> ortho_neighbor_set_iv;
			compute_Nt_and_No_iv(iv, grad_iv, scalar_grid, tangent_neighbor_set_iv,
				ortho_neighbor_set_iv, io_info);

			//DEBUG

			if(debug){
				cout <<"debug "<< debug<< endl;
			using namespace std;
			COORD_TYPE c[DIM3]={0.0,0.0,0.0};
			scalar_grid.ComputeCoord(iv,c);
			cout <<"cube index{"<<iv<<"}"
				<<" coord{"<<c[0]<<" "<<c[1]<<" "<<c[2]<<"}"
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
			cout <<"}"<<endl;}
			
			//DEBUG_END


			float  alpha = 5; //degrees
			bool flag_nt = Nt_based_computations_for_iv(iv, alpha, scalar_grid,
				gradient_grid, grad_mag_grid, tangent_neighbor_set_iv);
			bool flag_no = No_based_computations_for_iv(iv, alpha, scalar_grid,
				gradient_grid, grad_mag_grid, ortho_neighbor_set_iv);
			
			//DEBUG
			if(debug)
			cout <<"is reliable: Nt "<< flag_nt <<" No "<< flag_no<<endl;
			//DEBUG_END
			
			reliable_grid.Set(iv, (flag_nt && flag_no) );

		}
	}
}
//
// Curvature based reliable gradients computations
//
void compute_reliable_gradients_curvature_based(
	const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
	const GRADIENT_GRID & gradient_grid,
	const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
	IJK::BOOL_GRID<RELIGRADIENT_GRID> & reliable_grid,
	INPUT_INFO & io_info
	) 
{

	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) 
	{
		curvature_based_per_vertex(iv, scalar_grid, gradient_grid,
			grad_mag_grid, reliable_grid, io_info);
		
	}
}
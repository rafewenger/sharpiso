/// \file cgradient.cxx
/// compute gradients from scalar data
/// Version 0.0.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "cgradient.h"
#include "ijkscalar_grid.txx"
#include "ijkgrid.txx"
#include "ijkcoord.txx"
#include "isodual3D_datastruct.h"
#include "sharpiso_scalar.txx"
#include "ijkgrid_macros.h"
#include "sharpiso_types.h"
#include <algorithm>
#include <vector>
using namespace ISODUAL3D;
using namespace std;
using namespace SHARPISO;

// local type definition
namespace {

typedef IJK::BOOL_GRID_BASE<ISODUAL_GRID> BOOL_GRID_BASE;
typedef IJK::BOOL_GRID<ISODUAL_GRID> BOOL_GRID;

};

// local routines
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient, const SCALAR_TYPE & min_gradient_mag);

void compute_gradient_weighted_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient,
		const SCALAR_TYPE & min_gradient_mag,
		const float * weights);
void compute_boundary_gradient
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient);
void compute_cube_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & gradient_grid, const INPUT_INFO & io_info);
void compute_cube_grad
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient,
		const SCALAR_TYPE & min_gradient_mag, const INPUT_INFO & io_info);
void print_info (
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID & vertex_gradient_grid,
		const INPUT_INFO & io_info,
		const int direction,
		const VERTEX_INDEX iv,
		out e);

void compute_plane_point_dist (
		const GRADIENT_TYPE * normal, const SCALAR_TYPE * pt_on_plane,
		const SCALAR_TYPE * far_point, SCALAR_TYPE & dist);

SCALAR_TYPE zero_vector [3] = { 0.0, 0.0, 0.0};

bool gradients_agree
(const GRADIENT_GRID & vertex_gradient_grid,
		const GRADIENT_TYPE * gradient1,
		const VERTEX_INDEX vert2,
		const VERTEX_INDEX vert1,
		const INPUT_INFO & io_info)
{
	SCALAR_TYPE mag1=0.0, mag2=0.0;
	GRADIENT_TYPE   gradient2[DIM3]={0.0,0.0,0.0};
	std::copy(vertex_gradient_grid.VectorPtrConst(vert2),
			vertex_gradient_grid.VectorPtrConst(vert2)+DIM3,
			&(gradient2[0]));
	IJK::compute_magnitude_3D(gradient2, mag2);
	if (mag2 > 0.0)
	{
		// compute unit vertex gradient
		for (int l=0;l<DIM3;l++){
			gradient2[l]=gradient2[l]/mag2;
		}
		SCALAR_TYPE inn_pdt=0;
		IJK::compute_inner_product (DIM3, gradient1, gradient2, inn_pdt);
		/// DEBUG
		if(io_info.print_info && vert1==io_info.print_info_vertex){
			cout <<" (angle diff "<< (acos(inn_pdt)*180.0)/M_PI<<") ";
		}
		if ( inn_pdt > io_info.min_cos_of_angle){
			return true;
		}
	}
	else {
		if(io_info.print_info && vert1==io_info.print_info_vertex){
			cout <<" mag "<< mag2<<endl;
		}
	}
	return false;

}

/*
 * Compute reliable gradients by comparing how good the gradients
 * predict the scalar values of neighboring grids
 * OLD CODE : updating below to work only for certain neighbors
 */
void compute_reliable_gradients_SBP_2
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & vertex_gradient_grid,
		INPUT_INFO & io_info)
{
	using namespace SHARPISO;
	//compute the gradients using central difference
	compute_gradient_central_difference
	(scalar_grid, vertex_gradient_grid, io_info);

	// setup a secondary grid to keep track
	// of which gradients are reliable
	BOOL_GRID reliable_grid;
	reliable_grid.SetSize(scalar_grid);
	reliable_grid.SetAll(false);

	for (VERTEX_INDEX iv=0; iv < scalar_grid.NumVertices(); iv++)
	{
		int num_agree = 0;
		GRADIENT_TYPE  grad1[DIM3]={0.0,0.0,0.0};
		GRADIENT_TYPE  grad2[DIM3]={0.0,0.0,0.0};
		std::copy(vertex_gradient_grid.VectorPtrConst(iv),
				vertex_gradient_grid.VectorPtrConst(iv)+DIM3, grad1);
		SCALAR_TYPE mag1;
		IJK::compute_magnitude_3D(grad1, mag1);
		if (mag1 > io_info.min_gradient_mag)
		{
			IJK::multiply_coord_3D(-1.0, grad1, grad1);
			COORD_TYPE coord_iv[DIM3],coord_iv_prev[DIM3],coord_iv_next[DIM3];
			scalar_grid.ComputeCoord(iv, coord_iv);
			print_info (scalar_grid, vertex_gradient_grid, io_info, 0, iv, CURR_VERTEX);
			print_info (scalar_grid, vertex_gradient_grid, io_info, 0, iv, GRADIENT);
			for (unsigned int d = 0; d < DIM3; d++)
			{
				int k = min(io_info.scalar_prediction_dist, int(coord_iv[d]));
				VERTEX_INDEX prev_vertex = iv - k*scalar_grid.AxisIncrement(d);
				scalar_grid.ComputeCoord(prev_vertex, coord_iv_prev);

				print_info (scalar_grid, vertex_gradient_grid, io_info, d, iv, PREV_VERTEX);

				COORD_TYPE distance = 0;
				compute_signed_distance_to_gfield_plane
				(grad1, coord_iv, scalar_grid.Scalar(iv),
						coord_iv_prev, scalar_grid.Scalar(prev_vertex), distance);


				if(io_info.print_info && iv==io_info.print_info_vertex){
					cout <<"Distance: "<< distance <<endl;
				}
				k = min(io_info.scalar_prediction_dist,
						(scalar_grid.AxisSize(d)-int(coord_iv[d])-1));
				VERTEX_INDEX next_vertex = iv + k*scalar_grid.AxisIncrement(d);

				print_info (scalar_grid, vertex_gradient_grid, io_info, d, iv, NEXT_VERTEX);
				scalar_grid.ComputeCoord(next_vertex, coord_iv_next);
				compute_signed_distance_to_gfield_plane
				(grad1, coord_iv, scalar_grid.Scalar(iv),
						coord_iv_next, scalar_grid.Scalar(next_vertex), distance);
				if(io_info.print_info && iv==io_info.print_info_vertex){
					cout <<"Distance: "<< distance <<endl;
				}
			}


		}

	}

}

void compute_reliable_gradients_SBP(
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & vertex_gradient_grid,
		BOOL_GRID & reliable_grid,
		INPUT_INFO & io_info)
{
	using namespace SHARPISO;
	//compute the gradients using central difference
	compute_gradient_central_difference
	(scalar_grid, vertex_gradient_grid, io_info);

	for (VERTEX_INDEX iv=0; iv < scalar_grid.NumVertices(); iv++)
	{
		if (iv == int(scalar_grid.NumVertices()*0.25))
			cout <<"25% vertices done"<< endl;
		if (iv == int(scalar_grid.NumVertices()*0.50))
			cout <<"50% vertices done"<< endl;
		if (iv == int(scalar_grid.NumVertices()*0.75))
			cout <<"75% vertices done"<< endl;

		int num_agree = 0;
		// set up a vector to keep track of the distances
		vector <SCALAR_TYPE> vec_scalar_dists;

		GRADIENT_TYPE  grad1[DIM3]={0.0,0.0,0.0}, grad1_normalized[DIM3]={0.0,0.0,0.0};
		GRADIENT_TYPE  grad2[DIM3]={0.0,0.0,0.0};
		std::copy(vertex_gradient_grid.VectorPtrConst(iv),
				vertex_gradient_grid.VectorPtrConst(iv)+DIM3, grad1);
		SCALAR_TYPE mag1;
		IJK::compute_magnitude_3D(grad1, mag1);
		if (mag1 > io_info.min_gradient_mag)
		{
			// find the normalized gradient
			IJK::normalize_vector(DIM3, grad1, io_info.min_gradient_mag, grad1_normalized);
			// point on the plane
			COORD_TYPE coord_iv[DIM3];
			scalar_grid.ComputeCoord(iv, coord_iv);
			// find neighbor vertices
			vector <VERTEX_INDEX> near_vertices;
			scalar_grid.GetVertexNeighbors(iv, io_info.scalar_prediction_dist, near_vertices);

			if(io_info.draw && iv==io_info.draw_vert){
				// plot neighbor vertices
				cout <<"point "<< coord_iv[0] <<" "<<coord_iv[1]<<" "<<coord_iv[2]<<
						" 1 1 0"<<endl;
			}
			if(io_info.print_info && iv==io_info.print_info_vertex){
				print_info (scalar_grid, vertex_gradient_grid, io_info, 0, iv, CURR_VERTEX);
				print_info (scalar_grid, vertex_gradient_grid, io_info, 0, iv, GRADIENT);
				cout <<"num of neighbors: " << near_vertices.size() << endl;
			}


			// for all the neighboring points find the distance of points to plane
			// nv_ind = near_vertex_index

			// Number of vertices which are close to the gradient plane through the vertex iv
			VERTEX_INDEX close_vertices = 0;
			VERTEX_INDEX correct_prediction = 0;
			for (unsigned int nv_ind = 0; nv_ind < near_vertices.size(); nv_ind++)
			{
				// compute distance from plane
				COORD_TYPE coord_nv[DIM3]={0.0,0.0,0.0};
				scalar_grid.ComputeCoord(near_vertices[nv_ind], coord_nv);
				SCALAR_TYPE dist_to_plane = 0.0;
				compute_plane_point_dist (grad1_normalized, coord_iv, coord_nv, dist_to_plane);

				if(io_info.print_info && iv==io_info.print_info_vertex){
					cout <<"near vert "<< near_vertices[nv_ind] <<" coord (" << coord_nv[0] <<","<<coord_nv[1]<<","<<coord_nv[2]<<") ";
					cout <<"Distance to plane: "<< dist_to_plane <<endl;
					cout <<"Draw point "<<" "<< coord_nv[0] <<" "<<coord_nv[1]<<" "<<coord_nv[2]<<endl;
				}


				// if dist_to_plane is within the threshold
				if (abs (dist_to_plane) < 0.5){
					close_vertices++;

					// compute the distance between prediction and observed scalar value
					SCALAR_TYPE err_distance = 0.0;
					compute_signed_distance_to_gfield_plane
					(grad1, coord_iv, scalar_grid.Scalar(iv),
							coord_nv, scalar_grid.Scalar(near_vertices[nv_ind]), err_distance);
					// keep track of the error distances
					vec_scalar_dists.push_back (abs(err_distance));

					if(io_info.print_info && iv==io_info.print_info_vertex){
						cout <<"\n within threshold error in scalar pred is "<<
								err_distance<<"\n\n";
					}

					// Compare to error distance (set to .15)
					if (abs (err_distance) < io_info.scalar_prediction_err)
					{
						correct_prediction++;
						if(io_info.draw && iv==io_info.draw_vert){
							cout <<"point "<< coord_nv[0] <<" "<<coord_nv[1]<<" "<<coord_nv[2]<<
									" 0 1 0"<<endl;
						}
					}
					else
					{
						if(io_info.draw && iv==io_info.draw_vert){
							cout <<"point "<< coord_nv[0] <<" "<<coord_nv[1]<<" "<<coord_nv[2]<<
									" 1 0 0"<<endl;
						}
					}
				}
				else{
					if(io_info.draw && iv==io_info.draw_vert){
						cout <<"point "<< coord_nv[0] <<" "<<coord_nv[1]<<" "<<coord_nv[2]<<
								" 1 0 1"<<endl;
					}
				}
			}

			sort(vec_scalar_dists.begin(), vec_scalar_dists.end());

			if(io_info.print_info && iv==io_info.print_info_vertex)
			{
				cout << "io_info.scalar_prediction_err: "<<io_info.scalar_prediction_err<<endl;
				cout <<"correct prediction count: "<< correct_prediction
						<<" close vertices count: "<< close_vertices << endl;
				cout <<"prediction rate : " <<  correct_prediction/float(close_vertices) <<endl;

				cout <<"smallest scalar distance is "<< vec_scalar_dists[0]<<" largest "
						<< vec_scalar_dists[vec_scalar_dists.size()-1]<<endl;
			}

			if (vec_scalar_dists[vec_scalar_dists.size()-1] > io_info.scalar_prediction_err){
				if(io_info.print_info && iv==io_info.print_info_vertex){
					cout <<"Not Reliable" <<endl;
				}
				reliable_grid.Set(iv,false);
			}
		}

	}
}


/*
 * Compute reliable gradients by comparing with neighboring gradients
 * which are at a distance eliable_grad_far_dist
 */
void compute_reliable_gradients_far
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & vertex_gradient_grid,
		BOOL_GRID & reliable_grid,
		INPUT_INFO & io_info)
{
	using namespace SHARPISO;
	//Compute the vertex gradients using central difference
	compute_gradient_central_difference
	(scalar_grid, vertex_gradient_grid, io_info);
	int numAgree = 0;
	// DEBUG
	if(io_info.print_info)
	{
		cout <<"vertex Index " << io_info.print_info_vertex;
		COORD_TYPE coord0[DIM3];
		scalar_grid.ComputeCoord(io_info.print_info_vertex,coord0);
		cout <<" loc ["<<coord0[0]<<" "<<coord0[1]<<" "<<coord0[2]<<"]("
				<<scalar_grid.Scalar(io_info.print_info_vertex)<<")"<<endl;
	}
	for (VERTEX_INDEX iv=0; iv < scalar_grid.NumVertices(); iv++)
	{
		if (iv == int(scalar_grid.NumVertices()*0.25))
			cout <<"25% vertices done"<< endl;
		if (iv == int(scalar_grid.NumVertices()*0.50))
			cout <<"50% vertices done"<< endl;
		if (iv == int(scalar_grid.NumVertices()*0.75))
			cout <<"75% vertices done"<< endl;
		numAgree=0;
		GRADIENT_TYPE  gradient1[DIM3]={0.0,0.0,0.0}, gradient2[DIM3]={0.0,0.0,0.0};
		SCALAR_TYPE mag;
		std::copy(vertex_gradient_grid.VectorPtrConst(iv),
				vertex_gradient_grid.VectorPtrConst(iv)+DIM3,
				&(gradient1[0]));
		IJK::compute_magnitude_3D(gradient1, mag);
		if (mag > 0.0)
		{
			for (int l=0;l<DIM3;l++){
				gradient1[l]=gradient1[l]/mag;
			}

			if(io_info.print_info && iv==io_info.print_info_vertex){
				cout <<"vertex_gradient (cdiff) "<<gradient1[0]<<","<<gradient1[1]<<","<<gradient1[2];
				cout <<" mag "<<mag<<endl;
			}
			COORD_TYPE coord0[DIM3];
			scalar_grid.ComputeCoord(iv,coord0);

			for (int d=0; d<DIM3; d++)
			{
				int k = min(io_info.reliable_grad_far_dist, int(coord0[d]));
				VERTEX_INDEX prev_vertex = iv - k*scalar_grid.AxisIncrement(d);

				if (gradients_agree(vertex_gradient_grid, gradient1, prev_vertex, iv, io_info))
					numAgree++;

				if(io_info.print_info && iv==io_info.print_info_vertex){
					COORD_TYPE coord0[DIM3];
					scalar_grid.ComputeCoord(prev_vertex ,coord0);
					cout <<"prev vertex [" << coord0[0] <<" "<<coord0[1] <<" "<<coord0[2]<<"] ("
							<<scalar_grid.Scalar(prev_vertex)<<")"<<endl;
					cout <<"numAgree "<<numAgree <<endl;
				}

				k = min(io_info.reliable_grad_far_dist, (scalar_grid.AxisSize(d)-int(coord0[d])-1));
				VERTEX_INDEX next_vertex = iv + k*scalar_grid.AxisIncrement(d);

				if (gradients_agree(vertex_gradient_grid, gradient1, next_vertex, iv, io_info))
					numAgree++;

				/// DEBUG
				if(io_info.print_info && iv==io_info.print_info_vertex){
					COORD_TYPE coord0[DIM3];
					scalar_grid.ComputeCoord(next_vertex ,coord0);
					cout <<"next vertex [" << coord0[0] <<" "<<coord0[1] <<" "<<coord0[2]<<"] ("
							<<scalar_grid.Scalar(next_vertex)<<")"<<endl;
					cout <<"numAgree "<<numAgree <<endl;
				}

			}
			//DEBUG
			if (numAgree < io_info.min_num_agree){
				reliable_grid.Set(iv,false);
				io_info.out_info.num_unreliable++;
				io_info.out_info.un_reliable_grads_vert_info.push_back(iv);
			}
			else {
				io_info.out_info.num_reliable++;
			}
		}
	}
}

void compute_reliable_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & vertex_gradient_grid,
		INPUT_INFO & io_info)
{
	using namespace SHARPISO;
	//Compute the vertex gradients using central difference
	compute_gradient_central_difference
	(scalar_grid, vertex_gradient_grid, io_info);

	//Compute cube gradients
	GRADIENT_GRID  cube_gradient_grid;

	compute_cube_gradients
	(scalar_grid, cube_gradient_grid, io_info);

	BOOL_GRID boundary_grid;
	boundary_grid.SetSize(scalar_grid);
	compute_boundary_grid(boundary_grid);

	BOOL_GRID reliable_grid;
	reliable_grid.SetSize(scalar_grid);
	reliable_grid.SetAll(false);

	// DEBUG
	if(io_info.print_info)
	{
		cout <<"vertex Index " << io_info.print_info_vertex;
		COORD_TYPE coord0[DIM3];
		scalar_grid.ComputeCoord(io_info.print_info_vertex,coord0);
		cout <<" loc ["<<coord0[0]<<" "<<coord0[1]<<" "<<coord0[2]<<"]"<<endl;
	}

	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		if (!boundary_grid.Scalar(iv)) {
			int numAgree=0;
			GRADIENT_TYPE  gradient[DIM3]={0.0,0.0,0.0};
			SCALAR_TYPE mag,cube_grad_magnitude;

			std::copy(vertex_gradient_grid.VectorPtrConst(iv),
					vertex_gradient_grid.VectorPtrConst(iv)+DIM3,
					&(gradient[0]));
			IJK::compute_magnitude_3D(gradient, mag);
			/// DEBUG
			if(io_info.print_info && iv==io_info.print_info_vertex){
				cout <<"vertex_gradient (cdiff) "<<gradient[0]<<","<<gradient[1]<<","<<gradient[2];
				cout <<" mag "<<mag<<endl;
			}
			if (mag > 0.0)
			{
				// compute unit vertex gradient
				for (int l=0;l<DIM3;l++){
					gradient[l]=gradient[l]/mag;
				}
				GRADIENT_TYPE * cube_grad;
				VERTEX_INDEX cube_temp = iv-scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

				for (int k=0; k< NUM_CUBE_VERTICES3D; k++){
					VERTEX_INDEX in_cube = scalar_grid.CubeVertex(cube_temp,k);
					cube_grad=cube_gradient_grid.VectorPtr(in_cube);
					IJK::compute_magnitude_3D(cube_grad, cube_grad_magnitude);
					if (cube_grad_magnitude > 0.0)
						// compute unit cube gradient
						for (int l=0; l<DIM3; l++){
							cube_grad[l]=cube_grad[l]/cube_grad_magnitude;
						}
					/// DEBUG
					if(io_info.print_info && iv==io_info.print_info_vertex){
						cout << "cube ["<< in_cube<<"] ";
						COORD_TYPE coord0[DIM3];
						scalar_grid.ComputeCoord(in_cube,coord0);
						cout << "loc ["<<coord0[0]<<" "<<coord0[1]<<" "<<coord0[2]<<"] ";
						cout<<  "gradient ["<< cube_grad[0]<<","<<cube_grad[1]<<","<<cube_grad[2]<<"]";
					}

					SCALAR_TYPE inn_pdt=0;
					IJK::compute_inner_product (DIM3, cube_grad, gradient, inn_pdt);
					if ( inn_pdt > io_info.min_cos_of_angle){
						numAgree++;
						if(io_info.print_info && iv==io_info.print_info_vertex){
							cout <<" angle diff "<< (acos(inn_pdt)*180.0)/M_PI;
							cout <<" NumAgree "<< numAgree <<endl;
						}
					}
					else{
						if(io_info.print_info && iv==io_info.print_info_vertex)
							cout <<" angle diff "<< (acos(inn_pdt)*180.0)/M_PI <<endl;
					}
				}

				if (numAgree < io_info.min_num_agree){
					if(io_info.print_info && iv==io_info.print_info_vertex)
					{cout <<"Vertex "<<iv <<" not reliable, num agree " << numAgree<<endl;}

					//vertex_gradient_grid.Set(iv, zero_vector);
					io_info.out_info.num_unreliable++;
					io_info.out_info.un_reliable_grads_vert_info.push_back(iv);
				}
				else {
					reliable_grid.Set(iv,true);
					io_info.out_info.num_reliable++;
				}
			}
		}
		else
		{
			// boundary
			if(io_info.print_info && iv==io_info.print_info_vertex)
			{
				cout << "The vertex is in boundary."<<endl;
			}
			io_info.out_info.boundary_verts++;
		}

	}

	for (VERTEX_INDEX iv=0; iv < scalar_grid.NumVertices(); iv++){

		if (!reliable_grid.Scalar(iv)){
			vertex_gradient_grid.Set(iv, zero_vector);
		}
	}
}



/*
 * compute the central differnce
 */
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & gradient_grid, const INPUT_INFO & io_info)
{
	const int dimension = scalar_grid.Dimension();

	gradient_grid.SetSize(scalar_grid, dimension);

	BOOL_GRID boundary_grid;
	boundary_grid.SetSize(scalar_grid);
	compute_boundary_grid(boundary_grid);

	if (io_info.flag_weighted_cdiff)
		cout <<"Weighted central differnce being used.\n";

	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		if (boundary_grid.Scalar(iv)) {
			compute_boundary_gradient(scalar_grid, iv, gradient_grid.VectorPtr(iv));
		}
		else {
			if (io_info.flag_weighted_cdiff)
			{
				compute_gradient_weighted_central_difference
				(scalar_grid, iv, gradient_grid.VectorPtr(iv), io_info.min_gradient_mag, io_info.weights);
			}
			else
				compute_gradient_central_difference
				(scalar_grid, iv, gradient_grid.VectorPtr(iv), io_info.min_gradient_mag);

		}
	}
}

/// Compute central difference per vertex
void compute_gradient_weighted_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient,
		const SCALAR_TYPE & min_gradient_mag,
		const float * weights)
{
	const int dimension = scalar_grid.Dimension();
	for (int d = 0; d < dimension; d++) {
		VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
		VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
		gradient[d] = ((scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/(2*weights[d]));
	}
	// set min gradient to )
	SCALAR_TYPE mag=0.0;
	IJK::compute_magnitude_3D(gradient, mag);
	if (mag < min_gradient_mag ) {
		IJK::copy_coord_3D(zero_vector, gradient);
	}
}
/// Compute central difference per vertex
void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient,
		const SCALAR_TYPE & min_gradient_mag)
{
	const int dimension = scalar_grid.Dimension();
	for (int d = 0; d < dimension; d++) {
		VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
		VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
		gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2;
	}
	// set min gradient to )
	SCALAR_TYPE mag=0.0;
	IJK::compute_magnitude_3D(gradient, mag);
	if (mag < min_gradient_mag ) {
		IJK::copy_coord_3D(zero_vector, gradient);
	}
}

void compute_boundary_gradient
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient)
{
	const int dimension = scalar_grid.Dimension();
	GRID_COORD_TYPE coord[dimension];

	scalar_grid.ComputeCoord(iv1, coord);

	for (int d = 0; d < dimension; d++) {
		if (coord[d] > 0) {
			VERTEX_INDEX iv0 = scalar_grid.PrevVertex(iv1, d);
			if (coord[d]+1 < scalar_grid.AxisSize(d)) {
				VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
				// use central difference
				gradient[d] = (scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv0))/2;
			}
			else {
				gradient[d] = scalar_grid.Scalar(iv1) - scalar_grid.Scalar(iv0);
			}
		}
		else if (coord[d]+1 < scalar_grid.AxisSize(d)) {
			VERTEX_INDEX iv2 = scalar_grid.NextVertex(iv1, d);
			gradient[d] = scalar_grid.Scalar(iv2) - scalar_grid.Scalar(iv1);
		}
		else {
			gradient[d] = 0;
		}
	}
}

void compute_cube_gradients
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & gradient_grid,
		const INPUT_INFO & io_info)
{
	//compute the cube gradients
	const int dimension = scalar_grid.Dimension();
	gradient_grid.SetSize(scalar_grid, dimension);
	IJK_FOR_EACH_GRID_CUBE(iv, scalar_grid, VERTEX_INDEX){
		compute_cube_grad
		(scalar_grid, iv, gradient_grid.VectorPtr(iv), io_info.min_gradient_mag, io_info);
	}
}

/// Compute gradient per cube
void compute_cube_grad(
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const VERTEX_INDEX iv1, GRADIENT_TYPE * gradient,
		const SCALAR_TYPE & min_gradient_mag,
		const INPUT_INFO & io_info)
{
	const int num_facet_vertices = scalar_grid.NumFacetVertices();

	if(io_info.print_info && iv1==io_info.print_info_vertex)
	{cout <<" cube ID "<<iv1<<" "<<endl;}

	for (int d=0; d<DIM3; d++){
		SCALAR_TYPE diff =0;
		for (int k=0; k< num_facet_vertices; k++){
			VERTEX_INDEX iend0 = scalar_grid.FacetVertex(iv1,d,k);
			VERTEX_INDEX iend1 = scalar_grid.NextVertex(iend0,d);

			if(io_info.print_info && iv1==io_info.print_info_vertex){
				cout <<iend0<<"("<<scalar_grid.Scalar(iend0)<<")"<<"-"<<iend1
						<<"("<<scalar_grid.Scalar(iend1)<<")";
				cout <<" diff "<<(scalar_grid.Scalar(iend1)-scalar_grid.Scalar(iend0))<<endl;
			}

			diff = diff + ((scalar_grid.Scalar(iend1)-scalar_grid.Scalar(iend0)));
		}

		diff=diff/4.0;
		if(io_info.print_info && iv1==io_info.print_info_vertex){
			cout <<"Gradient direction ["<<d<<"] "<<diff<<endl;
		}
		gradient[d]=diff;
	}
	// set min gradient to )
	SCALAR_TYPE mag=0.0;
	IJK::compute_magnitude_3D(gradient, mag);
	if (mag < min_gradient_mag ) {
		IJK::copy_coord_3D(zero_vector, gradient);
	}
	if(io_info.print_info && iv1==io_info.print_info_vertex){
		cout <<"magnitude "<< mag << endl;
		if (mag > 0.0)
			cout <<"["<<gradient[0]/mag<<" "<<gradient[1]/mag<<" "<<gradient[2]/mag<<"]\n\n"<<endl;

	}
}
/*
 * Helper print function
 */
/*
 * Given a  plane ( point and a normal),
 *  find the distance of another point from this plane
 */
void compute_plane_point_dist (
		const GRADIENT_TYPE * normal,
		const SCALAR_TYPE * pt_on_plane,
		const SCALAR_TYPE * far_point,
		SCALAR_TYPE & dist
)
{
	SCALAR_TYPE d = 0; //  n.pt on plane
	SCALAR_TYPE temp = 0; // n.far_point
	IJK::compute_inner_product_3D(normal, pt_on_plane, d);
	IJK::compute_inner_product_3D(normal, far_point, temp);
	dist = temp - d;
}

void print_info (
		const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		const		GRADIENT_GRID & vertex_gradient_grid,
		const		INPUT_INFO & io_info,
		const 		int direction,
		const VERTEX_INDEX iv,
		out e)
{
	//CURR_VERTEX, PREV_VERTEX, NEXT_VERTEX
	if (io_info.print_info && io_info.print_info_vertex == iv )
	{
		COORD_TYPE coord_iv[DIM3], coord_iv_prev[DIM3], coord_iv_next[DIM3];
		scalar_grid.ComputeCoord(iv, coord_iv);
		int k;
		switch ( e )
		{
		case CURR_VERTEX:
		{
			cout <<"Curr Vertex: " << iv;
			cout << " ("<<coord_iv[0]<<","<<coord_iv[1]<<","<<coord_iv[2]<<")";
			cout <<" Scalar: "<< scalar_grid.Scalar(iv)<<endl;
			break;
		}
		case PREV_VERTEX:  // intentional fall-through
		{

			k = min(io_info.scalar_prediction_dist, int(coord_iv[direction ]));
			VERTEX_INDEX prev_vertex = iv - k*scalar_grid.AxisIncrement(direction);
			scalar_grid.ComputeCoord(prev_vertex, coord_iv_prev);
			cout <<"Prev Vertex: (direction "<<direction<<") "<< prev_vertex;
			cout << " ("<<coord_iv_prev[0]<<","<<coord_iv_prev[1]<<","<<coord_iv_prev[2]<<")";
			cout <<" Scalar: "<< scalar_grid.Scalar(prev_vertex)<<endl;
			break;
		}
		case NEXT_VERTEX:
		{
			k = min(io_info.scalar_prediction_dist,
					(scalar_grid.AxisSize(direction)-int(coord_iv[direction])-1));
			VERTEX_INDEX next_vertex = iv + k*scalar_grid.AxisIncrement(direction);
			cout <<"Next Vertex: (direction "<<direction<<") "<< next_vertex;
			scalar_grid.ComputeCoord(next_vertex, coord_iv_next);
			cout << " ("<<coord_iv_next[0]<<","<<coord_iv_next[1]<<","<<coord_iv_next[2]<<")";
			cout <<" Scalar: "<< scalar_grid.Scalar(next_vertex)<<endl;
			break;
		}
		case GRADIENT:
		{
			GRADIENT_TYPE  grad1[DIM3]={0.0,0.0,0.0};
			std::copy(vertex_gradient_grid.VectorPtrConst(iv),
					vertex_gradient_grid.VectorPtrConst(iv)+DIM3, grad1);
			SCALAR_TYPE mag1;
			IJK::compute_magnitude_3D(grad1, mag1);
			cout <<"Gradient: "<< grad1[0]<<" "<< grad1[1]<<" "<<grad1[2];
			IJK::normalize_vector
			(DIM3, grad1, io_info.min_gradient_mag, grad1);
			cout <<" Gradient Dir: "<< grad1[0]<<" "<< grad1[1]<<" "<<grad1[2];
			cout <<" Mag: "<<mag1<<endl;
			break;
		}
		default:
			cout <<" ERROR: "<< e << " is not defined."<<endl;
			break;
		}
	}
}

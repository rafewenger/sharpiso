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
#include "ijkgrid_macros.h"
#include "sharpiso_types.h"

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


SCALAR_TYPE zero_vector [3] = { 0.0, 0.0, 0.0};


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
			GRADIENT_TYPE * gradient;
			SCALAR_TYPE mag,cube_grad_magnitude;
			gradient=vertex_gradient_grid.VectorPtr(iv);
			IJK::compute_magnitude_3D(gradient, mag);

			/// DEBUG
			if(io_info.print_info && iv==io_info.print_info_vertex){
				cout <<"vertex_gradient (cdiff) "<<gradient[0]<<","<<gradient[1]<<","<<gradient[2];
				cout <<" mag "<<mag<<endl;
			}
			if (mag > 0.0){
				// compute unit vertex gradient
				for (int l=0;l<DIM3;l++){
					gradient[l]=gradient[l]/mag;
				}
				GRADIENT_TYPE * cube_grad;
				VERTEX_INDEX cube_temp=iv-scalar_grid.CubeVertexIncrement(NUM_CUBE_VERTICES3D-1);

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
							cout <<"angle < cos_angle NumAgree= " << numAgree<<endl;
						}
					}
					else{
						if(io_info.print_info && iv==io_info.print_info_vertex)
							cout <<" angle diff "<< (acos(inn_pdt)*180.0)/M_PI <<endl;
					}
				}
			}
			if (numAgree < 5){
				if(io_info.print_info && iv==io_info.print_info_vertex)
				{cout <<"Vertex "<<iv <<" not reliable, num agree " << numAgree<<endl;}

				vertex_gradient_grid.Set(iv, zero_vector);
				io_info.out_info.num_unreliable++;
				io_info.out_info.un_reliable_grads_vert_info.push_back(iv);
			}
			else {
				io_info.out_info.num_reliable++;
			}
		}
		else {
			// boundary
			if(io_info.print_info && iv==io_info.print_info_vertex){
				cout << "The vertex is in boundary."<<endl;
			}
			io_info.out_info.boundary_verts++;
		}
	}
}

void compute_gradient_central_difference
(const ISODUAL_SCALAR_GRID_BASE & scalar_grid,
		GRADIENT_GRID & gradient_grid, const INPUT_INFO & io_info)
{
	const int dimension = scalar_grid.Dimension();

	gradient_grid.SetSize(scalar_grid, dimension);

	BOOL_GRID boundary_grid;
	boundary_grid.SetSize(scalar_grid);
	compute_boundary_grid(boundary_grid);

	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		if (boundary_grid.Scalar(iv)) {
			compute_boundary_gradient(scalar_grid, iv, gradient_grid.VectorPtr(iv));
		}
		else {
			compute_gradient_central_difference
			(scalar_grid, iv, gradient_grid.VectorPtr(iv), io_info.min_gradient_mag);
		}
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
			VERTEX_INDEX iend1 = scalar_grid.NextVertex(iv1,d);

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


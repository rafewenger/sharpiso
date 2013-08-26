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
	GRID_COORD_TYPE coord[dimension];

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
// angle based 
// check if gradients agree

bool gradients_agree(const GRADIENT_GRID & gradient_grid,
		const GRADIENT_COORD_TYPE * gradient1, const VERTEX_INDEX vert2,
		const VERTEX_INDEX vert1, const INPUT_INFO & io_info) {
	GRADIENT_COORD_TYPE mag1 = 0.0, mag2 = 0.0;
	GRADIENT_COORD_TYPE gradient2[DIM3] = { 0.0, 0.0, 0.0 };

	std::copy(gradient_grid.VectorPtrConst(vert2),
			gradient_grid.VectorPtrConst(vert2) + DIM3, &(gradient2[0]));
	IJK::compute_magnitude_3D(gradient2, mag2);
	if (mag2 > io_info.min_gradient_mag) {
		// compute unit vertex gradient
		gradient2[0] = gradient2[0] / mag2;
		gradient2[1] = gradient2[1] / mag2;
		gradient2[2] = gradient2[2] / mag2;

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
		GRADIENT_GRID & gradient_grid, GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
		BOOL_GRID & reliable_grid, INPUT_INFO & io_info) {
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

			for (int d = 0; d < DIM3; d++) {
				int k = min(io_info.angle_based_dist, int(coord_iv[d]));
				VERTEX_INDEX prev_vertex = iv
						- k * scalar_grid.AxisIncrement(d);

				if (gradients_agree(gradient_grid, gradient_iv, prev_vertex, iv,
						io_info))
					numAgree++;

				k = min(io_info.angle_based_dist,
						(scalar_grid.AxisSize(d) - int(coord_iv[d]) - 1));
				VERTEX_INDEX next_vertex = iv
						+ k * scalar_grid.AxisIncrement(d);

				if (gradients_agree(gradient_grid, gradient_iv, next_vertex, iv,
						io_info))
					numAgree++;
			}

			if (numAgree < io_info.min_num_agree) {
				reliable_grid.Set(iv, false);
				io_info.out_info.num_unreliable++;
			} else {
				io_info.out_info.num_reliable++;
			}
		} else {
			io_info.out_info.grad_mag_zero++;
		}
	}
}

// scalar based 
void compute_plane_point_dist
(const GRADIENT_COORD_TYPE * normal,
 const COORD_TYPE * pt_on_plane, const COORD_TYPE * far_point,
 SCALAR_TYPE & dist) 
{
	COORD_TYPE d = 0; //  n.pt on plane
	COORD_TYPE temp = 0; // n.far_point
	IJK::compute_inner_product_3D(normal, pt_on_plane, d);
	IJK::compute_inner_product_3D(normal, far_point, temp);
	dist = temp - d;
}

// scalar based 
void compute_scaled_plane_point_dist
(const GRADIENT_COORD_TYPE normal[],
 const COORD_TYPE point0[], const COORD_TYPE point1[],
 const COORD_TYPE scale[], COORD_TYPE & dist) 
{
  dist = 0;
  for (int d = 0; d < DIM3; d++) {
    dist += (point1[d]-point0[d])*normal[d]/scale[d];
  }
}

// Scalar based prediction
void compute_reliable_gradients_SBP(
		const RELIGRADIENT_SCALAR_GRID_BASE & scalar_grid,
		const GRADIENT_GRID & gradient_grid,const  GRADIENT_MAGNITUDE_GRID & grad_mag_grid,
		IJK::BOOL_GRID<RELIGRADIENT_GRID> &reliable_grid,
		INPUT_INFO & io_info) {


	for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
		int num_agree = 0;
		// set up a vector to keep track of the distances
		vector<SCALAR_TYPE> vec_scalar_dists;

		GRADIENT_COORD_TYPE grad1_normalized[DIM3] = { 0.0, 0.0, 0.0 };
		GRADIENT_COORD_TYPE grad2[DIM3] = { 0.0, 0.0, 0.0 };
		std::copy(gradient_grid.VectorPtrConst(iv),
				gradient_grid.VectorPtrConst(iv) + DIM3, grad1_normalized);

		GRADIENT_COORD_TYPE mag1 = grad_mag_grid.Scalar(iv);

		if (mag1 > io_info.min_gradient_mag) {
			io_info.num_vertices_mag_grt_zero++;
			// find the normalized gradient
			// point on the plane
			COORD_TYPE coord_iv[DIM3];
			scalar_grid.ComputeScaledCoord(iv, coord_iv);
			// find neighbor vertices
			vector<VERTEX_INDEX> near_vertices;
			scalar_grid.GetVertexNeighbors(iv, io_info.scalar_prediction_dist,
					near_vertices);

			// for all the neighboring points find the distance of points to plane
			// nv_ind = near_vertex_index
			// Number of vertices which are close to the gradient plane through the vertex iv
			VERTEX_INDEX close_vertices = 0;
			VERTEX_INDEX correct_prediction = 0;

			bool flag_correct = true;
			for (int nv_ind = 0; nv_ind < near_vertices.size(); nv_ind++) {
				// compute distance from plane
				COORD_TYPE coord_nv[DIM3] = { 0.0, 0.0, 0.0 };
				scalar_grid.ComputeScaledCoord(near_vertices[nv_ind], coord_nv);

				COORD_TYPE dist_to_plane = 0.0;

				compute_scaled_plane_point_dist
          (grad1_normalized, coord_iv, coord_nv, scalar_grid.SpacingPtrConst(),
           dist_to_plane);

				// if dist_to_plane is within the threshold
				if (abs(dist_to_plane) < 0.5) {
					close_vertices++;
					// compute the distance between prediction and observed scalar value
					SCALAR_TYPE err_distance = 0.0;
					for (int l=0;l<DIM3;l++)
						grad2[l]=grad1_normalized[l]*mag1;
					compute_signed_distance_to_gfield_plane(grad2,
							coord_iv, scalar_grid.Scalar(iv), coord_nv,
							scalar_grid.Scalar(near_vertices[nv_ind]),
							err_distance);

					//keep track of the error distances
					//vec_scalar_dists.push_back (abs(err_distance));
					float abs_err = abs(err_distance);
					if (abs_err > io_info.scalar_prediction_err) {
						flag_correct = false;
						break;
					}

					// Compare to error distance (set to .15)
					if (abs(err_distance) < io_info.scalar_prediction_err) {
						correct_prediction++;
						io_info.out_info.num_reliable++;

					}
				}
			}
			if (!flag_correct) {
				reliable_grid.Set(iv, false);
				io_info.out_info.num_unreliable++;
			}
		}
	}
}


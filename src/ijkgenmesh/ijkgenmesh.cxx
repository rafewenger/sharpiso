/// \file ijkgenmesh.cxx
/// generate a mesh

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2014-2015 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define _USE_MATH_DEFINES
#include <math.h>

#include "ijkcube.txx"
#include "ijkmesh.txx"
#include "ijkgencoord.txx"

#include "ijkgenmesh.h"


using namespace IJK;
using namespace IJKGENMESH;
using namespace IJKGENCOORD;
using namespace IJKGENGEOM;
using namespace std;

// types
typedef IJK::CUBE<int, NUM_TYPE, COORD_TYPE, COORD_TYPE> CUBE_TYPE;
typedef IJK::CUBE_FACE_INFO<int, NUM_TYPE, NUM_TYPE> CUBE_FACE_INFO_TYPE;

// global variable
extern MESH_PARAM mesh_param;

// forward declarations
void add_vertices_on_circle
(const int dimension, const NUM_TYPE num_vert, const COORD_TYPE radius, 
 const COORD_TYPE zcoord, std::vector<COORD_TYPE> & coord);
void add_cube_vertices
(const CUBE_TYPE & cube, std::vector<COORD_TYPE> & coord);
void add_quad_between_circles
(const NUM_TYPE num_vert, const NUM_TYPE iv0_circleA, 
 const NUM_TYPE iv0_circleB, POLYMESH_TYPE & mesh);
void add_split_cube_facets
(const int icube, const CUBE_FACE_INFO_TYPE & cube_face_info, 
 const NUM_TYPE ivX, POLYMESH_TYPE & mesh);
void gen_single_cube
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);
void gen_two_cubes
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);
void find_vertex_in_interior
(const CUBE_TYPE & cubeA, const CUBE_TYPE & cubeB, 
 int & ivA, int & num_in_interior);
void compute_in_facet_vertex_coord
(const CUBE_TYPE & cube, const COORD_TYPE split_coord[],
 std::vector<COORD_TYPE> & in_facet_vertex_coord);
void compute_in_edge_vertex_coord
(const CUBE_TYPE & cube, const COORD_TYPE split_coord[],
 std::vector<COORD_TYPE> & in_edge_vertex_coord);
bool check_dimension
(const int dimension, const int dimension2, IJK::ERROR & error);
bool check_num_in_interior(const int num_in_interior, IJK::ERROR & error);
void get_center_diff
(const int dimension, const OBJECT_PROPERTIES & prop, 
 COORD_TYPE center_diff[], IJK::ERROR & error);

template <typename CTYPE0, typename CTYPE1>
void find_different_direction
(const int dimension, const CTYPE0 * dir0, CTYPE1 * dir1);
template <typename CTYPE0, typename CTYPE1>
void rotate_coord
(const int dimension, const CTYPE0 * dir0, const CTYPE1 * dir1,
 std::vector<COORD_TYPE> & coord);
template <typename CTYPE>
void translate_coord
(const int dimension, const CTYPE * translate_dir,
 std::vector<COORD_TYPE> & coord);
template<typename T>
void append_vector(const std::vector<T> & v0, std::vector<T> & v1);
template <typename CTYPE0, typename CTYPE1>
void compute_frame_coord3D
(const CTYPE0 * dir0, const CTYPE1 * dir1,
 const COORD_TYPE coord[], COORD_TYPE frame_coord[]);


// **************************************************
// GENERATE MESH
// **************************************************

void IJKGENMESH::gen_cube
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  IJK::PROCEDURE_ERROR error("gen_cube");

  if (mesh_param.NumObjects() <= 1) {
    gen_single_cube(dimension, prop, coord, mesh);
  }
  else if (mesh_param.NumObjects() == 2) {
    gen_two_cubes(dimension, prop, coord, mesh);
  }
  else {
    error.AddMessage
      ("Programming error. gen_cube not implemented for more than 2 cubes.");
    throw error;
  }
}

void gen_single_cube
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  const COORD_TYPE distance = prop.Distance();
  CUBE_TYPE cube(dimension);
  CUBE_FACE_INFO_TYPE cube_face_info(dimension);
  IJK::ARRAY<COORD_TYPE> vertex0_coord(dimension);
  IJK::PROCEDURE_ERROR error("gen_singel_cube");

  if (!check_dimension(DIM3, dimension, error)) { throw error; }

  set_coord(dimension, -distance, vertex0_coord.Ptr());

  cube.SetVertexCoord(vertex0_coord.PtrConst(), 2*distance);

  NUM_TYPE num_coord = dimension*cube.NumVertices();
  coord.resize(num_coord);
  COORD_TYPE * coord_ptr = &(coord[0]);

  std::copy(cube.VertexCoord(), cube.VertexCoord()+num_coord, coord_ptr);

  if (mesh_param.flag_tilt) { 
    rotate_coord
      (dimension, prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0), 
       coord); 
  }

  // translate by center
  translate_coord(dimension, prop.CenterPtrConst(0), coord);

  mesh.Clear();
  const NUM_TYPE num_vert_per_facet = cube_face_info.NumFacetVertices();
  for (NUM_TYPE ifacet = 0; ifacet < cube_face_info.NumFacets(); ifacet++) {
    mesh.num_poly_vert.push_back(num_vert_per_facet);
    mesh.first_poly_vert.push_back(num_vert_per_facet*ifacet);

    for (NUM_TYPE j = 0; j < num_vert_per_facet; j++) {
      NUM_TYPE iv = cube_face_info.FacetVertex(ifacet, j);
      mesh.poly_vert.push_back(iv);
    }
  }

  if (dimension == 3) {
    IJK::reorder_quad_vertices(mesh.poly_vert);
    IJK::reverse_quad_orientations(mesh.poly_vert);
  }

}

void gen_two_cubes
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  const COORD_TYPE distance = prop.Distance();
  CUBE<int, NUM_TYPE, COORD_TYPE, COORD_TYPE> cube[2];
  CUBE_FACE_INFO<int, int, int> cube_face_info(dimension);
  int ivX[2];
  int num_in_interior;
  vector<COORD_TYPE> in_facet_vertex_coord;
  vector<COORD_TYPE> in_edge_vertex_coord;
  IJK::ARRAY<COORD_TYPE> vertex0_coord(dimension);
  IJK::ARRAY<COORD_TYPE> center_diff(dimension);
  IJK::PROCEDURE_ERROR error("gen_two_cubes");

  if (!check_dimension(DIM3, dimension, error)) { throw error; }

  for (int i = 0; i < 2; i++) 
    { cube[i].SetDimension(DIM3); }

  set_coord(dimension, -distance, vertex0_coord.Ptr());

  cube[0].SetVertexCoord(vertex0_coord.PtrConst(), 2*distance);

  get_center_diff(dimension, prop, center_diff.Ptr(), error);

  add_coord(dimension, vertex0_coord.PtrConst(), 
            center_diff.PtrConst(), vertex0_coord.Ptr());
  cube[1].SetVertexCoord(vertex0_coord.PtrConst(), 2*distance);

  // Set coord[].
  coord.clear();
  for (int i = 0; i < 2; i++) 
    { add_cube_vertices(cube[i], coord); }

  for (int i = 0; i < 2; i++) {
    int i2 = (i+1)%2;
    find_vertex_in_interior(cube[i], cube[i2], ivX[i], num_in_interior);
    if (!check_num_in_interior(num_in_interior, error)) { throw error; }
  }

  for (int i = 0; i < 2; i++) {
    int i2 = (i+1)%2;
    compute_in_facet_vertex_coord
      (cube[i], cube[i2].VertexCoord(ivX[i2]), in_facet_vertex_coord);
    append_vector(in_facet_vertex_coord, coord);
  }

  for (int i = 0; i < 2; i++) {
    int i2 = (i+1)%2;
    compute_in_edge_vertex_coord
      (cube[i], cube[i2].VertexCoord(ivX[i2]), in_edge_vertex_coord);
    append_vector(in_edge_vertex_coord, coord);
  }

  if (mesh_param.flag_tilt) { 
    rotate_coord
      (dimension, prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0), 
       coord); 
  }

  // translate by center
  translate_coord(dimension, prop.CenterPtrConst(0), coord);

  // Set polymesh facets.
  mesh.Clear();
  for (int i = 0; i < 2; i++)
    { add_split_cube_facets(i, cube_face_info, ivX[i], mesh); }

  if (dimension == 3) {
    IJK::reverse_quad_orientations(mesh.poly_vert);
  }
}

void IJKGENMESH::gen_annulus
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  if (prop.flag_flange && 
      prop.flange_width.size() > 0 && prop.flange_height.size() > 0) {

    gen_flanged_annulus(dimension, prop, coord, mesh);
  }
  else {
    gen_annulus_no_flange(dimension, prop, coord, mesh);
  }
}


void IJKGENMESH::gen_annulus_no_flange
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  const COORD_TYPE distance = prop.Distance();
  COORD_TYPE radius0, radius1;
  COORD_TYPE half_height;
  NUM_TYPE num_angle;
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::PROCEDURE_ERROR error("gen_annulus");

  if (!check_dimension(DIM3, dimension, error)) { throw error; }

  radius0 = prop.radius[0]-distance;
  radius1 = prop.radius[0]+distance;

  if (radius0 < 1) {
    error.AddMessage("Error. Inner radius is too small.");
    error.AddMessage("  Increase distance or radius.");
    throw error;
  }

  half_height = distance+prop.length_difference[0];
  if (half_height < 1) {
    error.AddMessage("Error.  Height too small.");
    error.AddMessage("  Increase length_diff.");
    throw error;
  }
  
  num_angle = int(2*M_PI*radius1);
  if (num_angle < 4) { num_angle = 4; }

  coord.clear();
  add_vertices_on_circle
    (dimension, num_angle, radius0, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius0, half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, half_height, coord);

  find_different_direction
    (dimension, prop.DirectionPtrConst(0), dir1.Ptr());
  rotate_coord(dimension, prop.DirectionPtrConst(0), dir1.PtrConst(), coord);

  // translate by center
  translate_coord(dimension, prop.CenterPtrConst(0), coord);

  mesh.Clear();
  add_quad_between_circles(num_angle, num_angle, 0, mesh);
  add_quad_between_circles(num_angle, 2*num_angle, 3*num_angle, mesh);
  add_quad_between_circles(num_angle, 0, 2*num_angle, mesh);
  add_quad_between_circles(num_angle, 3*num_angle, num_angle, mesh);

  if (dimension == 3) {
    IJK::reverse_quad_orientations(mesh.poly_vert);
  }
}

void IJKGENMESH::gen_flanged_annulus
(const int dimension, const OBJECT_PROPERTIES & prop,
 std::vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  const COORD_TYPE distance = prop.Distance();
  COORD_TYPE radius0, radius1, radius2, radius3;
  COORD_TYPE flange_width, flange_height;
  COORD_TYPE half_height;
  NUM_TYPE num_angle;
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::PROCEDURE_ERROR error("gen_flanged_annulus");

  if (!check_dimension(DIM3, dimension, error)) { throw error; }

  if (!prop.flag_flange) {
    error.AddMessage("Programming error.  Annulus does not have a flange.");
    throw error;
  }

  if (prop.flange_width.size() < 1) {
    error.AddMessage("Programming error.  Missing flange width.");
    throw error;
  }

  if (prop.flange_height.size() < 1) {
    error.AddMessage("Programming error.  Missing flange height.");
    throw error;
  }

  flange_width = prop.flange_width[0];
  flange_height = prop.flange_height[0];
  radius0 = prop.radius[0]-distance-flange_width;
  radius1 = prop.radius[0]-distance;
  radius2 = prop.radius[0]+distance;
  radius3 = prop.radius[0]+distance+flange_width;

  if (radius0 < 1) {
    error.AddMessage("Error. Inner radius is too small.");
    error.AddMessage("  Increase distance or radius.");
    throw error;
  }

  half_height = distance+prop.length_difference[0];
  if (half_height < 1) {
    error.AddMessage("Error.  Height too small.");
    error.AddMessage("  Increase length_diff.");
    throw error;
  }
  
  num_angle = int(2*M_PI*radius3);
  if (num_angle < 4) { num_angle = 4; }

  coord.clear();
  add_vertices_on_circle
    (dimension, num_angle, radius0, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius0, half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, -half_height-flange_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius1, half_height+flange_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius2, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius2, half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius2, -half_height-flange_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius2, half_height+flange_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius3, -half_height, coord);
  add_vertices_on_circle
    (dimension, num_angle, radius3, half_height, coord);

  find_different_direction
    (dimension, prop.DirectionPtrConst(0), dir1.Ptr());
  rotate_coord(dimension, prop.DirectionPtrConst(0), dir1.PtrConst(), coord);

  // translate by center
  translate_coord(dimension, prop.CenterPtrConst(0), coord);

  mesh.Clear();
  add_quad_between_circles(num_angle, num_angle, 0, mesh);
  add_quad_between_circles(num_angle, 2*num_angle, 4*num_angle, mesh);
  add_quad_between_circles(num_angle, 5*num_angle, 3*num_angle, mesh);
  add_quad_between_circles(num_angle, 8*num_angle, 6*num_angle, mesh);
  add_quad_between_circles(num_angle, 7*num_angle, 9*num_angle, mesh);
  add_quad_between_circles(num_angle, 10*num_angle, 11*num_angle, mesh);

  add_quad_between_circles(num_angle, 0, 2*num_angle, mesh);
  add_quad_between_circles(num_angle, 3*num_angle, num_angle, mesh);
  add_quad_between_circles(num_angle, 4*num_angle, 8*num_angle, mesh);
  add_quad_between_circles(num_angle, 9*num_angle, 5*num_angle, mesh);
  add_quad_between_circles(num_angle, 6*num_angle, 10*num_angle, mesh);
  add_quad_between_circles(num_angle, 11*num_angle, 7*num_angle, mesh);

  if (dimension == 3) {
    IJK::reverse_quad_orientations(mesh.poly_vert);
  }
}


// **************************************************
// ROUTINES TO ADD VERTICES
// **************************************************

void add_vertices_on_circle
(const int dimension, const NUM_TYPE num_vert, const COORD_TYPE radius, 
 const COORD_TYPE zcoord, std::vector<COORD_TYPE> & coord)
{
  const ANGLE_TYPE angle0 = 2*M_PI/num_vert;
  NUM_TYPE s = coord.size();
  coord.resize(s+num_vert*dimension);

  for (NUM_TYPE iv = 0; iv < num_vert; iv++) {
    NUM_TYPE k = iv*dimension+s;
    coord[k] = radius*cos(iv*angle0);
    coord[k+1] = radius*sin(iv*angle0);
    coord[k+2] = zcoord;
  }
}

void add_cube_vertices
(const CUBE_TYPE & cube, std::vector<COORD_TYPE> & coord)
{
  const int dimension = cube.Dimension();
  const NUM_TYPE num_cube_coord = dimension*cube.NumVertices();
  NUM_TYPE s = coord.size();

  coord.resize(s+num_cube_coord);
  COORD_TYPE * coord_ptr = &(coord[0]) + s;

  std::copy
    (cube.VertexCoord(), cube.VertexCoord()+num_cube_coord, coord_ptr);
}

// Append vector v0 to v1.
template<typename T>
void append_vector(const std::vector<T> & v0, std::vector<T> & v1)
{
  const NUM_TYPE s0 = v0.size();
  const NUM_TYPE s1 = v1.size();

  v1.resize(s0+s1);
  T * v1_ptr = &(v1[0]) + s1;

  std::copy(v0.begin(), v0.end(), v1_ptr);
}

// **************************************************
// ROUTINES TO ADD POLYGONS
// **************************************************

void add_quad_between_circles
(const NUM_TYPE num_vert, const NUM_TYPE iv0_circleA, 
 const NUM_TYPE iv0_circleB, POLYMESH_TYPE & mesh)
{
  for (NUM_TYPE iv = 0; iv < num_vert; iv++) {
    ANGLE_TYPE iv2 = (iv+1)%num_vert;
    mesh.num_poly_vert.push_back(NUM_VERT_PER_QUAD);
    mesh.first_poly_vert.push_back(mesh.poly_vert.size());
    mesh.poly_vert.push_back(iv+iv0_circleA);
    mesh.poly_vert.push_back(iv+iv0_circleB);
    mesh.poly_vert.push_back(iv2+iv0_circleB);
    mesh.poly_vert.push_back(iv2+iv0_circleA);
  }
}


/// Add quads formed by splitting cube facets to polymesh.
/// Do not add quads incident on ivX.
/// @param icube = cube index.
void add_split_cube_facets
(const int icube, const CUBE_FACE_INFO_TYPE & cube_face_info, 
 const NUM_TYPE ivX, POLYMESH_TYPE & mesh)
{
  const int dimension = cube_face_info.Dimension();
  const int num_cube_vertices = cube_face_info.NumVertices();
  const int num_cube_facets = cube_face_info.NumFacets();
  const int num_cube_edges = cube_face_info.NumEdges();
  const int num_facet_vertices = cube_face_info.NumFacetVertices();
  //NUM_TYPE facet_vertex[num_facet_vertices];
  //NUM_TYPE facet_vertex[cube_face_info.NumFacetVertices()];
  IJK::ARRAY<NUM_TYPE> facet_vertex(num_facet_vertices);
  IJK::PROCEDURE_ERROR error("add_split_cube_facets");

  if (!check_dimension(DIM3, dimension, error)) { throw error; }


  for (NUM_TYPE ifacet = 0; ifacet < cube_face_info.NumFacets(); ifacet++) {

    int orth_dir = cube_face_info.FacetOrthDir(ifacet);
    int d1 = (orth_dir+1)%dimension;
    int d2 = (orth_dir+2)%dimension;

    // Copy facet vertices to facet_vert.
    for (NUM_TYPE k = 0; k < cube_face_info.NumFacetVertices(); k++) 
      { facet_vertex[k] = cube_face_info.FacetVertex(ifacet, k); }

    // Reorder vertices in cyclic order around facet.
    std::swap(facet_vertex[2], facet_vertex[3]);

    for (NUM_TYPE k = 0; k < cube_face_info.NumFacetVertices(); k++) {

      NUM_TYPE k1 = (k+1)%cube_face_info.NumFacetVertices();
      NUM_TYPE iv0 = facet_vertex[k];
      NUM_TYPE iv1 = facet_vertex[k1];

      if (iv0 != ivX) {
        // Add quadrilateral incident to iv0 on facet ifacet.

        NUM_TYPE ie1 = cube_face_info.IncidentEdge(iv0, d1);
        NUM_TYPE ie2 = cube_face_info.IncidentEdge(iv0, d2);

        if (cube_face_info.IsEdgeIncidentOnVertex(ie2, iv1))
          { std::swap(ie1, ie2); }

        if (!cube_face_info.IsEdgeIncidentOnVertex(ie1, iv1)) {
          error.AddMessage
            ("Programming error.  Vertex ", iv1, " on facet ", ifacet,
             " is not incident on edge ", ie1, ".");
          throw error;
        }

        mesh.num_poly_vert.push_back(NUM_VERT_PER_QUAD);
        mesh.first_poly_vert.push_back(mesh.poly_vert.size());
        mesh.poly_vert.push_back(iv0+icube*num_cube_vertices);

        NUM_TYPE jv = 2*(num_cube_vertices+num_cube_facets);
        jv += icube*(num_cube_edges) + ie1;
        mesh.poly_vert.push_back(jv);

        jv = 2*num_cube_vertices + icube*(num_cube_facets) + ifacet;
        mesh.poly_vert.push_back(jv);

        jv = 2*(num_cube_vertices+num_cube_facets);
        jv += icube*(num_cube_edges) + ie2;
        mesh.poly_vert.push_back(jv);
      }
    }
  }
}


// **************************************************
// ROUTINES TO ROTATE/TRANSLATE COORDINATES
// **************************************************

// Get vector from center 0 to center 1. 
void get_center_diff
(const int dimension, const OBJECT_PROPERTIES & prop, 
 COORD_TYPE center_diff[], IJK::ERROR & error)
{
  if (mesh_param.flag_stack) {
    IJK::copy_coord(dimension, mesh_param.TranslatePtrConst(0), center_diff);
  }
  else if (prop.NumCenters() > 1) {
    subtract_coord(dimension, prop.CenterPtrConst(1), 
                   prop.CenterPtrConst(0), center_diff);
  }
  else {
    error.AddMessage
      ("Programming error.  Missing location of second cube center.");
    throw error;
  }
}

// find arbitrary direction very different from dir0
template <typename CTYPE0, typename CTYPE1>
void find_different_direction
(const int dimension, const CTYPE0 * dir0, CTYPE1 * dir1)
{
  int d0 = 0;

  for (int d = 1; d < dimension; d++) {
    if (abs(dir0[d]) < abs(dir0[d0]))
      { d0 = d; }
  }

  set_coord(dimension, 0, dir1);
  dir1[d0] = 1;
}

// Rotate coord based to new coordinate frame.
// Compute new coordinate frame from dir0 and dir1.
// Rotate coord so that coordinates in new frame equal coordinates in old frame.
template <typename CTYPE0, typename CTYPE1>
void rotate_coord
(const int dimension, const CTYPE0 * dir0, const CTYPE1 * dir1,
 std::vector<COORD_TYPE> & coord)
{
  const NUM_TYPE numv = coord.size()/dimension;
  IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);
  IJK::ARRAY<COORD_TYPE> vertex_coord(dimension);

  compute_orthogonal_basis3D(dir0, dir1, zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  COORD_TYPE * coord_ptr = &(coord[0]);
  for (NUM_TYPE iv = 0; iv < numv; iv++) {

    multiply_coord(dimension, coord[iv*dimension], xdir.PtrConst(),
                   vertex_coord.Ptr());
    add_scaled_coord(dimension, coord[iv*dimension+1], ydir.PtrConst(),
                     vertex_coord.PtrConst(), vertex_coord.Ptr());
    add_scaled_coord(dimension, coord[iv*dimension+2], zdir.PtrConst(),
                     vertex_coord.PtrConst(), vertex_coord.Ptr());

    copy_coord
      (dimension, vertex_coord.PtrConst(), coord_ptr+iv*dimension);
  }

}

// translate coord
template <typename CTYPE>
void translate_coord
(const int dimension, const CTYPE * translate_dir,
 std::vector<COORD_TYPE> & coord)
{
  const NUM_TYPE numv = coord.size()/dimension;

  COORD_TYPE * coord_ptr = &(coord[0]);
  for (NUM_TYPE iv = 0; iv < numv; iv++) {
    add_coord(dimension, coord_ptr+iv*dimension, translate_dir,
              coord_ptr+iv*dimension);
  }
}

/// Compute coordinates in frame given by dir0 and dir1.
///  @pre coordinates are 3-dimensional.
template <typename CTYPE0, typename CTYPE1>
void compute_frame_coord3D
(const CTYPE0 * dir0, const CTYPE1 * dir1,
 const COORD_TYPE coord[], COORD_TYPE frame_coord[])
{
  IJK::ARRAY<COORD_TYPE> xdir(DIM3), ydir(DIM3), zdir(DIM3);
  COORD_TYPE new_coord[DIM3];

  compute_orthogonal_basis3D(dir0, dir1, zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  IJK::compute_inner_product(DIM3, xdir.PtrConst(), coord, new_coord[0]);
  IJK::compute_inner_product(DIM3, ydir.PtrConst(), coord, new_coord[1]);
  IJK::compute_inner_product(DIM3, zdir.PtrConst(), coord, new_coord[2]);
  IJK::copy_coord(DIM3, new_coord, frame_coord);
}


// **************************************************
// ROUTINES FOR VERTEX/EDGE/FACE/CUBE RELATIONSHIPS
// **************************************************

/// Find a vertex of cubeA contained in the interior of cubeB
/// @param ivA Some vertex in interior of cubeB (if one exists).
/// @param num_in_interior Number of vertices in interior of cubeB.
void find_vertex_in_interior
(const CUBE_TYPE & cubeA, const CUBE_TYPE & cubeB, 
 int & ivA, int & num_in_interior)
{
  ivA = 0;
  num_in_interior = 0;

  for (NUM_TYPE iv = 0; iv < cubeA.NumVertices(); iv++) {
    if (cubeB.InteriorContains(cubeA.VertexCoord(iv))) {
      ivA = iv;
      num_in_interior++;
    }
  }
}

/// Compute new vertex in interior of each facet of cube
/// @param split_coord[] Split cube facets using split_coord[].
void compute_in_facet_vertex_coord
(const CUBE_TYPE & cube, const COORD_TYPE split_coord[],
 std::vector<COORD_TYPE> & in_facet_vertex_coord)
{
  const int dimension = cube.Dimension();
  CUBE_FACE_INFO_TYPE cube_info(dimension);

  in_facet_vertex_coord.clear();
  for (NUM_TYPE ifacet = 0; ifacet < cube.NumFacets(); ifacet++) {

    int orth_dir = cube_info.FacetOrthDir(ifacet);
    int iv0 = cube_info.FacetVertex(ifacet, 0);
    COORD_TYPE c = cube.VertexCoord(iv0, orth_dir);

    for (int d = 0; d < dimension; d++) {
      if (d == orth_dir) 
        { in_facet_vertex_coord.push_back(c); }
      else 
        { in_facet_vertex_coord.push_back(split_coord[d]); }
    }
  }
}

/// Compute new vertex in interior of each edge of cube
/// @param split_coord[] Split cube facets using split_coord[].
void compute_in_edge_vertex_coord
(const CUBE_TYPE & cube, const COORD_TYPE split_coord[],
 std::vector<COORD_TYPE> & in_edge_vertex_coord)
{
  const int dimension = cube.Dimension();
  CUBE_FACE_INFO_TYPE cube_info(dimension);

  in_edge_vertex_coord.clear();
  for (NUM_TYPE ie = 0; ie < cube.NumEdges(); ie++) {

    int edge_dir = cube_info.EdgeDir(ie);
    int iv0 = cube_info.EdgeEndpoint(ie, 0);

    for (int d = 0; d < dimension; d++) {
      if (d == edge_dir)
        { in_edge_vertex_coord.push_back(split_coord[d]); }
      else
        { in_edge_vertex_coord.push_back(cube.VertexCoord(iv0, d)); }
    }
  }
}


// **************************************************
// CHECK ROUTINES
// **************************************************

bool check_dimension
(const int dimension, const int dimension2, IJK::ERROR & error)
{
  if (dimension2 != dimension) {
    error.AddMessage("Programming error. Illegal dimension ", dimension2, ".");
    error.AddMessage("  Dimension should be: ", dimension, ".");
    return(false);
  }

  return(true);
}

bool check_num_in_interior(const int num_in_interior, IJK::ERROR & error)
{
  if (num_in_interior < 1) {
    error.AddMessage("Error.  Cubes do not properly intersect.");
    error.AddMessage
      ("  Increase distance for larger cubes or change location of second cube.");
    return(false);
  }

  if (num_in_interior > 2) {
    error.AddMessage
      ("Programming error.  More than one cube vertex is in the interior of the other cube!?!.");
    return(false);
  }

  return(true);
}

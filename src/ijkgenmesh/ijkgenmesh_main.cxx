/// \file ijkgenmesh_main.cxx
/// generate a mesh
/// Version v0.1.0

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2014 Rephael Wenger

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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>


#include "ijk.txx"
#include "ijkIO.txx"
#include "ijkstring.txx"

#include "ijkgencoord.txx"

#include "ijkgenmesh.h"

using namespace IJKGENMESH;
using namespace IJKGENGEOM;
using namespace IJKGENCOORD;
using namespace std;


// global variables
char * ofilename = NULL;
SET_VALUE<int> grid_axis_size;
vector<MESH_INFO> mesh_info;
MESH_PARAM mesh_param;
bool output_gradients = false;
bool flag_silent = false;
double max_small_magnitude(0.0001);

// functions to read mesh parameters
void init_mesh_info
(vector<MESH_INFO> & mesh_info);
void read_num_vert_per_axis();
void read_mesh_distance(const MESH_INFO & mesh_info);
void read_mesh_name(const vector<MESH_INFO> & mesh_info);
void prompt_param(const std::string & mesh_name);
void set_mesh_param(const vector<MESH_INFO> & mesh_info);

// global functions
void set_random_param(const int dimension, const int num_objects);
void generate_mesh
(const int dimension, vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh);


// misc functions
void parse_command_line(const int argc, char **argv);
void usage_error(), help_msg();
void check_input_param(const int dimension);
void check_mesh_options();
void write_mesh
(const char * filename, const int dimension, const vector<COORD_TYPE> & coord, 
 const POLYMESH_TYPE & mesh);


int main(int argc, char **argv)
{
  IJK::PROCEDURE_ERROR error("ijkgenmesh");

  try {

    init_mesh_info(mesh_info);

    parse_command_line(argc, argv);

    if (ofilename == NULL) {
      error.AddMessage("Missing output file name");
      throw error;
    }

    if (!mesh_param.geom_info_index.IsSet()) 
      { read_mesh_name(mesh_info); }

    // Set parameters related to mesh.
    set_mesh_param(mesh_info);

    // Only dimension 3 is currently implemented.
    mesh_param.SetDimension(DIM3);

    int imesh = mesh_param.MeshIndex();
    std::string mesh_name = mesh_info[imesh].name;

    check_mesh_options();
    check_input_param(mesh_param.Dimension());

    if (!mesh_param.distance.IsSet()) 
      { read_mesh_distance(mesh_info[imesh]); }

    POLYMESH_TYPE mesh;
    vector<COORD_TYPE> coord;

    if (mesh_param.flag_grid_center && mesh_param.center.size() == 0) {

      if (!grid_axis_size.IsSet())
        { read_num_vert_per_axis(); }
      mesh_param.SetAxisSize(grid_axis_size.Value());

      for (int d = 0; d < mesh_param.Dimension(); d++) {
        COORD_TYPE x = (mesh_param.axis_size[d]-1)/2.0;
        mesh_param.center.push_back(x);
      }
    }

    // Prompt for any missing parameters.
    prompt_param(mesh_name);

    int dimension = mesh_param.Dimension();

    generate_mesh(dimension, coord, mesh);

    if (!flag_silent) {
      cout << "Writing mesh to " << ofilename << endl;
    }

    write_mesh(ofilename, dimension, coord, mesh);
  } 
  catch (IJK::ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// **************************************************
// Prompt for scalar mesh parameters
// **************************************************

template <typename CTYPE>
void prompt_coordinates
(const char * prompt, const int num_coord,
 std::vector<CTYPE> & coord)
{
  cout << prompt << " (" << num_coord << " values) : ";
  for (int ic = 0; ic < num_coord; ic++) {
    CTYPE c;
    cin >> c;
    coord.push_back(c);
  }
}

template <typename CTYPE>
void prompt_coordinates
(const std::ostringstream & prompt, const int num_coord,
 std::vector<CTYPE> & coord)
{
  prompt_coordinates(prompt.str().c_str(), num_coord, coord);
}

/// Prompt for a single scalar value and add to list
template <typename STYPE>
void prompt_scalar
(const char * prompt, std::vector<STYPE> & list)
{
  cout << prompt << " : ";
  STYPE s;
  cin >> s;
  list.push_back(s);
}

template <typename STYPE>
void prompt_scalar
(const std::ostringstream & prompt, std::vector<STYPE> & list)
{
  prompt_scalar(prompt.str().c_str(), list);
}

void prompt_num_objects(const int mesh_index)
{
  string mesh_name = mesh_info[mesh_index].name;

  int n(0);
  while (true) {

    if (mesh_name == "cube") {
      cout << "Enter number of cubes: ";
    }
    else if (mesh_name == "cylinder") {
      cout << "Enter number of cylinders: ";
    }
    else if (mesh_name == "annulus") {
      cout << "Enter number of annuli: ";
    }
    else if (mesh_name == "flange") {
      cout << "Enter number of flanges: ";
    }
    else if (mesh_name == "dist2planes") {
      cout << "Enter number of planes: ";
    }
    else {
      cout << "Enter number of objects: ";
    }

    cin >> n;
    if (n > 0) { break; };

    cout << "Illegal number.  Enter a positive integer." << endl;
  }

  mesh_param.num_objects.Set(n);
}

void prompt_center_translate(const int dimension)
{
  if (mesh_param.translate.size() == 0) {
    prompt_coordinates("Enter center translation vector",
                       dimension, mesh_param.translate);
  }
}

void prompt_center
(const int dimension, const int num_obj, const char * obj_name)
{
  int k = mesh_param.center.size()/dimension;

  // Clear any stray coordinates
  if (mesh_param.center.size()%dimension != 0) 
    { mesh_param.center.resize(k*dimension); }

  if (k >= num_obj) { return; };

  if (num_obj == 1) {
    ostringstream prompt;
    prompt << "Enter coordinates of " << obj_name << " center";
    prompt_coordinates(prompt, dimension, mesh_param.center);
  }
  else {

    if (mesh_param.flag_stack) {

      if (k < 1) {
        ostringstream prompt;
        prompt << "Enter coordinates of first " << obj_name << " center";
        prompt_coordinates(prompt, dimension, mesh_param.center);
      }

      prompt_center_translate(dimension);
    }
    else {
      for (int i = k; i < num_obj; i++) {
        if (mesh_param.center.size() <= i*dimension) {
          ostringstream prompt;
          prompt << "Enter coordinates of " << obj_name << " " 
                 << i << " center";
          prompt_coordinates(prompt, dimension, mesh_param.center);
        }
      }
    }
  }
}

void prompt_two_axis_directions(const int dimension)
{
  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter first axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.side_direction.size() == 0) {
    prompt_coordinates("Enter second axis direction",
                       dimension, mesh_param.side_direction);
  }
}

void prompt_wedge(const int dimension)
{
  const int mesh_index = mesh_param.MeshIndex();
  string mesh_name = mesh_info[mesh_index].name;

  if (mesh_param.wedge_angle.size() == 0) {
    prompt_scalar("Enter wedge angle", mesh_param.wedge_angle);
  }

  if (mesh_param.side_direction.size() == 0) {
    prompt_coordinates
      ("Enter coordinates of side direction",
       dimension, mesh_param.side_direction);
  }

  if (mesh_param.NumWedgeIsovalues() < 1) {
    prompt_scalar("Enter isovalue where wedge passes through the center", 
                  mesh_param.wedge_isovalue);
  }

}

void prompt_flange()
{
  if (mesh_param.flange_width.size() == 0) {
    prompt_scalar("Enter flange width", mesh_param.flange_width);
  }

  if (mesh_param.flange_height.size() == 0) {
    prompt_scalar("Enter flange height", mesh_param.flange_height);
  }
}

void prompt_cube
(const int dimension, const int num_cubes)
{
  prompt_center(dimension, num_cubes, "cube");

  if (mesh_param.flag_tilt) 
    { prompt_two_axis_directions(dimension); }
}

void prompt_cylinder
(const int dimension, const int num_cylinders, const bool flag_closed)
{
  prompt_center(dimension, num_cylinders, "cylinder");

  if (mesh_param.flag_tilt) {
    if (mesh_param.direction.size() == 0) {
      prompt_coordinates("Enter cylinder axis direction",
                         dimension, mesh_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    mesh_param.direction.clear();
    mesh_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { mesh_param.direction[0] = 1; }
  }

  if (flag_closed) {
    if (mesh_param.length_difference.size() == 0) {
      prompt_scalar("Enter difference between cylinder length and diameter", 
                    mesh_param.length_difference);
    }
  }

}

void prompt_cone(const int dimension, const bool flag_closed)
{
  const int num_cones = 1;

  prompt_center(dimension, num_cones, "cone");

  if (mesh_param.flag_tilt) {
    if (mesh_param.direction.size() == 0) {
      prompt_coordinates("Enter cone axis direction",
                         dimension, mesh_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    mesh_param.direction.clear();
    mesh_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { mesh_param.direction[0] = 1; }
  }

  if (mesh_param.angle.size() == 0) {
    prompt_scalar("Enter cone angle", mesh_param.angle);
  }

}

void prompt_frustrum(const int dimension)
{
  const int num_frustra = 1;

  prompt_center(dimension, num_frustra, "frustrum");

  if (mesh_param.flag_tilt) {
    if (mesh_param.direction.size() == 0) {
      prompt_coordinates("Enter frustrum axis direction",
                         dimension, mesh_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    mesh_param.direction.clear();
    mesh_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { mesh_param.direction[0] = 1; }
  }

  if (mesh_param.angle.size() == 0) {
    prompt_scalar("Enter frustrum angle", mesh_param.angle);
  }

  if (mesh_param.dist2near0.size() == 0) {
    prompt_scalar("Enter distance from apex to near plane", 
                  mesh_param.dist2near0);
  }

  if (mesh_param.dist2far0.size() == 0) {
    prompt_scalar("Enter distance from apex to far plane", 
                  mesh_param.dist2far0);
  }

}

void prompt_octahedron(const int dimension, const int num_objects)
{
  prompt_center(dimension, num_objects, "octahedron");
}

void prompt_dist2planes(const int dimension, const int num_planes)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point on all planes",
                       dimension, mesh_param.center);
  }

  int k = mesh_param.normal.size()/dimension;

  // Clear any stray coordinates
  if (mesh_param.normal.size()%dimension != 0) 
    { mesh_param.normal.resize(k*dimension); }

  if (k >= num_planes) { return; };

  for (int i = k; i < num_planes; i++) {
    ostringstream prompt;
    prompt << "Enter coordinates of normal to plane " << i;

    prompt_coordinates(prompt, dimension, mesh_param.normal);
  }
}

void prompt_edge(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point on edge",
                       dimension, mesh_param.center);
  }

  if (mesh_param.angle.size() == 0) {
    prompt_scalar
      ("Enter dihedral angle", mesh_param.angle);
  }

  if (mesh_param.normal.size() == 0) {
    prompt_coordinates
      ("Enter coordinates of normal to first plane", 
       dimension, mesh_param.normal);
  }

  if (mesh_param.rotation_direction.size() == 0) {
    prompt_coordinates
      ("Enter rotation direction", dimension, mesh_param.rotation_direction);
  }
}

void prompt_square_cylinder(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of cylinder center",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter cylinder axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.side_direction.size() == 0) {
    prompt_coordinates("Enter direction orthogonal to square side",
                       dimension, mesh_param.side_direction);
  }

  if (mesh_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  mesh_param.length_difference);
  }
}

void prompt_square_cylinder_rot45(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of cylinder center",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter cylinder axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.side_direction.size() == 0) {
    prompt_coordinates("Enter square corner direction",
                       dimension, mesh_param.side_direction);
  }

  if (mesh_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  mesh_param.length_difference);
  }
}

void prompt_crossing_square_cylinders(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of crossing center",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter direction of first cylinder",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.side_direction.size() == 0) {
    prompt_coordinates("Enter direction of second cylinder",
                       dimension, mesh_param.side_direction);
  }

  if (mesh_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  mesh_param.length_difference);
  }
}

void prompt_annulus(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of annulus center",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter annulus axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.radius.size() == 0) {
    prompt_scalar
      ("Enter annulus radius", mesh_param.radius);
  }

  if (mesh_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between annulus height and width", 
                  mesh_param.length_difference);
  }

}

void prompt_flange(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of flange center",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter flange axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.radius.size() == 0) {
    prompt_scalar
      ("Enter flange radius", mesh_param.radius);
  }

  if (mesh_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between flange height and width", 
                  mesh_param.length_difference);
  }

}

void prompt_torus(const int dimension, const int num_tori)
{
  prompt_center(dimension, num_tori, "torus");

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter torus axis direction",
                       dimension, mesh_param.direction);
  }

  if (mesh_param.radius.size() == 0) {
    prompt_scalar("Enter torus radius", mesh_param.radius);
  }
}

void prompt_constant_unit_gradient(const int dimension)
{
  if (mesh_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point with scalar value zero",
                       dimension, mesh_param.center);
  }

  if (mesh_param.direction.size() == 0) {
    prompt_coordinates("Enter gradient direction",
                       dimension, mesh_param.direction);
  }
}

template <typename T>
void prompt_maxval(IJKGENGEOM::SET_VALUE<T> & maxval)
{
  T x;

  cout << "Enter maximum scalar value: ";
  cin >> x;
  maxval.Set(x);
}

template <typename T>
void prompt_random_seed(IJKGENGEOM::SET_VALUE<T> & random_seed)
{
  T x;

  cout << "Enter random seed: ";
  cin >> x;
  random_seed.Set(x);
}

template <typename T>
void prompt_isotable_index
(const int dimension, IJKGENGEOM::SET_VALUE<T> & isotable_index)
{
  T num_table_entries;
  IJK::PROCEDURE_ERROR error("prompt_isotable_index");

  compute_num_table_entries(dimension, num_table_entries);

  if (num_table_entries < 1) {
    isotable_index.Set(0);
    return;
  }

  T x;
  do {
    cout << "Enter isosurface lookup table index [0," << num_table_entries-1
         << "]: ";
    cin >> x;

    if (x < 0) {
      cout << "Error.  Isosurface lookup table index must be non-negative."
           << endl;
    }
    else if (x >= num_table_entries) {
      cout << "Error.  Isosurface lookup table index must be less than "
           << num_table_entries << "." << endl;
    }

  } while (x < 0 || x >= num_table_entries);

  isotable_index.Set(x);
}

/// Prompt for any missing parameters.
void prompt_param(const std::string & mesh_name)
{
  const int dimension = mesh_param.Dimension();

  if (!mesh_param.num_objects.IsSet()) {
    if (mesh_param.flag_multi_centers || mesh_param.flag_multi_normals) {
      prompt_num_objects(mesh_param.MeshIndex());
    }
    else {
      mesh_param.num_objects.Set(1); 
    }
  }

  const int num_objects = mesh_param.NumObjects();

  set_random_param(dimension, num_objects);

  if (mesh_name == "cube") {
    prompt_cube(dimension, num_objects);
  }
  else if (mesh_name == "sphere") {
    prompt_center(dimension, num_objects, "sphere");
  }
  else if (mesh_name == "cylinder") {
    prompt_cylinder(dimension, num_objects, true);
  }
  else if (mesh_name == "open_cylinder") {
    prompt_cylinder(dimension, num_objects, false);
  }
  else if (mesh_name == "cone") {
    prompt_cone(dimension, true);
  }
  else if (mesh_name == "frustrum") {
    prompt_frustrum(dimension);
  }
  else if (mesh_name == "square") {
    prompt_square_cylinder(dimension);
  }
  else if (mesh_name == "square45") {
    prompt_square_cylinder_rot45(dimension);
  }
  else if (mesh_name == "squareX" ||
           mesh_name == "square45X") {
    prompt_crossing_square_cylinders(dimension);
  }
  else if (mesh_name == "annulus") {
    prompt_annulus(dimension);
  }
  else if (mesh_name == "flange") {
    prompt_flange(dimension);
  }
  else if (mesh_name == "torus") {
    prompt_torus(dimension, num_objects);
  }
  else if (mesh_name == "dist2planes") {
    prompt_dist2planes(dimension, num_objects);
  }
  else if (mesh_name == "edge") {
    prompt_edge(dimension);
  }
  else if (mesh_name == "octahedron") {
    prompt_octahedron(dimension, num_objects);
  }
  else if (mesh_name == "constant_unit_gradient") {
    prompt_constant_unit_gradient(dimension);
  }

  if (mesh_param.flag_flange) 
    { prompt_flange(); }

  if (mesh_param.flag_wedge)
    { prompt_wedge(dimension); }
}

// **************************************************
// Compute functions
// **************************************************

void compute_center_translation_vector
(const int dimension, COORD_TYPE translate[])
{
  const int imesh = mesh_param.MeshIndex();
  const std::string mesh_name = mesh_info[imesh].name;
  IJK::PROCEDURE_ERROR error("compute_center_translation_vector");

  if (mesh_param.flag_tilt) {

    if (mesh_param.NumDirections() <= 0) {
      error.AddMessage("Programming error. Missing tilt/axis direction.");
      throw error;
    };

    if (mesh_name == "cube") {
      IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

      if (mesh_param.NumSideDirections() <= 0) {
        error.AddMessage("Programming error. Missing side direction.");
        throw error;
      };

      compute_orthogonal_basis3D
        (mesh_param.DirectionPtrConst(0), mesh_param.SideDirectionPtrConst(0),
         zdir.Ptr(), xdir.Ptr(), ydir.Ptr());
      transform_basis3D
        (mesh_param.TranslatePtrConst(0), 
         xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), translate);
    }
    else {
      rotate_vector_ed
        (dimension, max_small_magnitude, mesh_param.DirectionPtrConst(0),
         mesh_param.TranslatePtrConst(0), translate);
    }
  }
  else {
    IJK::copy_coord
      (dimension, mesh_param.TranslatePtrConst(0), translate);
  }
}


// **************************************************
// Generate mesh
// **************************************************

/// Generate mesh
void generate_mesh
(const int dimension, vector<COORD_TYPE> & coord, POLYMESH_TYPE & mesh)
{
  const int num_objects = mesh_param.NumObjects();
  const int imesh = mesh_param.MeshIndex();
  const std::string mesh_name = mesh_info[imesh].name;
  OBJECT_PROPERTIES object_properties;
  std::vector<NUM_TYPE> tri_vert;
  IJK::PROCEDURE_ERROR error("generate_mesh");

  object_properties.Copy(mesh_param);
  mesh_info[imesh].function_ptr(dimension, object_properties, coord, mesh);

  if (mesh_param.flag_triangulate) {
    IJK::triangulate_polygon_list
      (mesh.num_poly_vert, mesh.poly_vert, mesh.first_poly_vert, tri_vert);

    mesh.Clear();
    mesh.AddPolytopes(tri_vert, NUM_VERT_PER_TRIANGLE);
  }
}

/// Compute a direction at a given angle from \a dir1.
/// @param dir1[] Direction.
/// @param dir2[] Second direction.
/// @param angle Angle in degrees.
/// @param[out] new_dir[] Direction at given with \a dir1 in plane
///                  determined by \a dir1 and \a dir2.
/// @pre dir1[] is a unit vector.
/// @pre dir2[] is a unit vector.
template <typename COORD1_TYPE, typename COORD2_TYPE, typename COORD3_TYPE>
void compute_dir_at_angle
(const int dimension, const COORD1_TYPE dir1[], 
 const COORD2_TYPE dir2[], const ANGLE_TYPE angle,
 COORD3_TYPE new_dir[])
{
  const double radians = M_PI*angle/180.0;

  IJK::multiply_coord(dimension, cos(radians), dir1, new_dir);
  IJK::add_scaled_coord
    (dimension, sin(radians), dir2, new_dir, new_dir);
}

template <typename T>
T ijk_mod(const T x, const int y)
{
  if (x >= 0) {
    return(fmod(x, y));
  }
  else {
    return(y-fmod(-x, y));
  }
}


/// Compute two normal vectors to wedge planes.
/// @param[out] normal[] Normals to wedge planes.
template <typename CTYPE>
void compute_wedge_plane_normals
(const OBJECT_PROPERTIES & prop, CTYPE normal[2*DIM3])
{
  IJK::ARRAY<COORD_TYPE> xdir(DIM3), ydir(DIM3), zdir(DIM3);
  IJK::PROCEDURE_ERROR error("compute_wedge_plane_normals");

  if (prop.NumWedgeAngles() < 1) {
    error.AddMessage("Programming error.  Missing wedge angle.");
    throw error;
  }

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0), 
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  IJK::copy_coord(DIM3, xdir.PtrConst(), normal);

  compute_dir_at_angle
    (DIM3, xdir.PtrConst(), ydir.PtrConst(),
     prop.wedge_angle[0], normal+DIM3);
  IJK::multiply_coord(DIM3, -1, normal+DIM3, normal+DIM3);
}


// **************************************************
// Read functions
// **************************************************

void read_num_vert_per_axis()
{
  int size;

  cout << "Enter # grid vertices per axis: ";
  cin >> size;

  grid_axis_size.Set(size);

  if (size < 2) {
    cerr << "Illegal grid axis size.  Must be at least 2." << endl;
    exit(30);
  };

}

void read_mesh_distance(const MESH_INFO & mesh_info)
{
  COORD_TYPE distance;

  if (mesh_info.flag_cylinder) 
    { cout << "Enter distance to circle: "; }
  else if (mesh_info.flag_center)
    { cout << "Enter distance to center: "; }
  else
    { cout << "Enter distance: "; }
  cin >> distance;

  mesh_param.distance.Set(distance);
}

// **************************************************
// MESH_INFO functions
// **************************************************

void init_mesh_info(vector<MESH_INFO> & mesh_info)
{
  MESH_INFO info;

  info.Set("cube", gen_cube);
  info.SetCubeFlags();
  info.SetDescription("Cube specified by distance from facets to cube center.");
  mesh_info.push_back(info);

  info.Set("annulus", gen_annulus);
  info.SetAnnulusFlags();
  info.SetDescription("Thickened annulus.");
  mesh_info.push_back(info);

}

/// Output mesh names.
/// @param with_gradients 
///           If true, output only meshes with gradients implemented.
void out_mesh_names(const vector<MESH_INFO> & mesh_info)
{
  bool newline_flag(false);

  int iout = 0;
  for (int i = 0; i < mesh_info.size(); i++) {
    if (iout > 0) { 
      cout << ","; 
      newline_flag = false;
      if (iout%5 == 0) { 
        newline_flag = true;
        cout << endl; 
      }
    }
    cout << mesh_info[i].name;
    iout++;
  }

  if (!newline_flag) { cout << endl; }
}

/// Output mesh descriptions.
void out_mesh_descriptions
(const vector<MESH_INFO> & mesh_info)
{
  const int COL1_WIDTH = 15;

  for (int i = 0; i < mesh_info.size(); i++) {
    int num_description_lines = mesh_info[i].description.size();
    string s = mesh_info[i].name + ":";

    cout << left << setw(COL1_WIDTH) << s;

    int j = 0;
    if ((num_description_lines > 0) && (s.size()+1 <= COL1_WIDTH)) {
      cout << mesh_info[i].description[0];
      j++;
    }
    cout << endl;

    while (j < num_description_lines) {
      cout << setw(COL1_WIDTH) << " " 
           << mesh_info[i].description[j] << endl;
      j++;
    }
  }

}

void read_mesh_name(const vector<MESH_INFO> & mesh_info)
{
  string s;
  int imesh;

  cout << "Enter scalar mesh name (h for help): ";
  cin >> s;

  while (true) {
    if (!find_geom_info_name(mesh_info, s, imesh)) {
      if (s != "h" && s != "help") {
        cout << "Illegal mesh name: " << s << endl;
      }
    }
    else {
      break;
    }

    cout << endl;
    cout << "Valid mesh names: ";
    out_mesh_names(mesh_info);
    cout << endl;
    cout << "Enter scalar mesh name (h for help): ";
    cin >> s;
  }

  mesh_param.SetMeshIndex(imesh);
}

// **************************************************
// MESH_PARAM functions
// **************************************************

void set_mesh_param(const vector<MESH_INFO> & mesh_info)
{
  int imesh = mesh_param.MeshIndex();

  if (mesh_info[imesh].flag_allow_multi_centers) {
    if (mesh_param.num_objects.IsSet()) {
      if (mesh_param.NumObjects() > 1) 
        { mesh_param.flag_multi_centers = true; }
    }
  }
  else
    { mesh_param.flag_multi_centers = false; }

  if (!mesh_info[imesh].flag_allow_tilt) 
    { mesh_param.flag_tilt = false; }

  if (!mesh_info[imesh].flag_allow_flange)
    { mesh_param.flag_flange = false; }

  if (!mesh_info[imesh].flag_allow_wedge)
    { mesh_param.flag_wedge = false; }
}

void set_random_param
(const int dimension, const int num_objects)
{
  const int MAX_GRAD_COORD = 5;
  const int MIN_GRAD_COORD = -MAX_GRAD_COORD;
  const int NUM_GRAD_COORD = (MAX_GRAD_COORD - MIN_GRAD_COORD) + 1;

  if (mesh_param.flag_random_centers) {

    if (mesh_param.randompos_seed.IsSet())
      { srand(mesh_param.randompos_seed.Value()); }

    for (int i = mesh_param.NumCenters(); i < num_objects; i++) {
      for (int d = 0; d < dimension; d++) {
        long int x = rand();
        x = x%(mesh_param.axis_size[d]);
        mesh_param.center.push_back(x);
      }
    }
  }

  if (mesh_param.flag_random_directions) {

    if (mesh_param.randomdir_seed.IsSet())
      { srand(mesh_param.randomdir_seed.Value()); }

    for (int i = mesh_param.NumDirections(); i < num_objects; i++) {
      for (int d = 0; d < dimension; d++) {
        long int x = rand();
        x = x%(NUM_GRAD_COORD);
        x = x + MIN_GRAD_COORD;
        mesh_param.direction.push_back(x);
      }
    }
  }

}


// **************************************************
// WRITE MESH
// **************************************************

void write_mesh
(const char * filename, const int dimension, const vector<COORD_TYPE> & coord, 
 const POLYMESH_TYPE & mesh)
{
  ofstream out(filename, ios::out);

  if (!out.good()) {
    cerr << "Unable to open output file " << filename << "." << endl;
    exit(65);
  };

  IJK::ijkoutPolytopeOFF(out, dimension, coord, mesh.num_poly_vert, 
                         mesh.poly_vert, mesh.first_poly_vert);

  out.close();

}


// **************************************************
// COMMAND LINE PARSER
// **************************************************

void check_input_param(const int dimension)
{
  int imesh = mesh_param.MeshIndex();
  std::string mesh_name = mesh_info[imesh].name;

  if (mesh_param.center.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -center does not match dimension "
         << dimension << ".";
    exit(50);
  }

  if (mesh_param.direction.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -dir does not match dimension "
         << dimension << ".";
    exit(50);
  }

  if (mesh_param.normal.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -dir does not match dimension "
         << dimension << ".";
    exit(50);
  }

}

void check_mesh_options()
{
  int imesh = mesh_param.MeshIndex();

  if (mesh_info[imesh].name == "cube") {
    if (mesh_param.flag_tilt && mesh_param.Dimension() != DIM3) {
      cerr << "Error."
           << " Illegal dimension " << mesh_param.Dimension() 
           << " used with mesh cube and option tilt."
           << endl;
      cerr << "  Mesh cube and option tilt can only be used with dimension " 
           << DIM3 << "." << endl;
      exit(50);
    }

    if (mesh_param.NumObjects() > 2) {
      cerr << "Error."
           << "  Number of cubes can be at most 2." << endl;
      exit(52);
    }
  }


  if (mesh_param.flag_wedge && mesh_param.Dimension() != DIM3) {
    cerr << "Usage error."
         << " Illegal dimension " << mesh_param.Dimension() 
         << " used with mesh " << mesh_info[imesh].name << "." << endl;
    cerr << "Wedge can only be used with dimension " << DIM3 << "." << endl;
    exit(55);
  }

}

int get_int(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  int x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

float get_float(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  float x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

template <typename T1, typename T2>
void get_arg(const int iarg, const int argc, char **argv,
             T1 & x1, T2 & x2)
{
  if (iarg+2 >= argc) { 
    cerr << "Usage error. Missing arguments for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  if (!IJK::string2val(argv[iarg+1], x1)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  if (!IJK::string2val(argv[iarg+2], x2)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+2] << endl;
    exit(50);
  }
}

/// Get argument for option argv[iarg].
/// Does not modify iarg.
template <typename ETYPE>
void get_float
(const int iarg, const int argc, char **argv, vector<ETYPE> & v)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!IJK::string2vector(argv[iarg+1], v)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }
}

int get_mesh_index
(const int iarg, const int argc, char **argv)
{
  int imesh(0);

  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!find_geom_info_name(mesh_info, argv[iarg+1], imesh)) {
    cerr << "Illegal mesh name: " << argv[iarg+1] << endl;
    exit(51);
  }

  return(imesh);
}

void parse_command_line(const int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];

    if (s == "-dim") {
      int dimension = get_int(iarg, argc, argv);
      mesh_param.SetDimension(dimension);
      iarg++;
    }
    else if (s == "-asize") {
      int size = get_int(iarg, argc, argv);
      grid_axis_size.Set(size);
      iarg++;
    }
    else if (s == "-mesh") {
      int imesh = get_mesh_index(iarg, argc, argv);
      mesh_param.SetMeshIndex(imesh);
      iarg++;
    }
    else if (s == "-distance") {
      COORD_TYPE distance = get_float(iarg, argc, argv);
      mesh_param.distance.Set(distance);
      iarg++;
    }
    else if (s == "-tilt") {
      mesh_param.flag_tilt = true;
    }
    else if (s == "-no_tilt") {
      mesh_param.flag_tilt = false;
    }
    else if (s == "-multi") {
      mesh_param.flag_multi_centers = true;
    }
    else if (s == "-single") {
      mesh_param.flag_multi_centers = false;
    }
    else if (s == "-two") {
      mesh_param.flag_multi_centers = true;
      mesh_param.num_objects.Set(2);
    }
    else if (s == "-n") {
      int num_obj = get_int(iarg, argc, argv);
      mesh_param.num_objects.Set(num_obj);
      iarg++;
    }
    else if (s == "-stack") {
      mesh_param.flag_stack = true;
      mesh_param.flag_multi_centers = true;
    }
    else if (s == "-flange") {
      mesh_param.flag_flange = true;
    }
    else if (s == "-flange_wh") {
      mesh_param.flag_flange = true;

      float w, h;
      get_arg(iarg, argc, argv, w, h);
      mesh_param.flange_width.push_back(w);
      mesh_param.flange_height.push_back(h);
      iarg = iarg+2;
    }
    else if (s == "-wedge") {
      mesh_param.flag_wedge = true;
    }
    else if (s == "-center") {
      get_float(iarg, argc, argv, mesh_param.center);
      iarg++;
    }
    else if (s == "-grid_center") {
      mesh_param.flag_grid_center = true;
    }
    else if (s == "-randompos") {
      mesh_param.flag_random_centers = true;
      int seed = get_int(iarg, argc, argv);
      mesh_param.randompos_seed.Set(seed);
      iarg++;
    }
    else if (s == "-randomdir") {
      mesh_param.flag_random_directions = true;
      int seed = get_int(iarg, argc, argv);
      mesh_param.randomdir_seed.Set(seed);
      iarg++;
    }
    else if (s == "-dir") {
      get_float(iarg, argc, argv, mesh_param.direction);
      iarg++;
    }
    else if (s == "-side_dir") {
      get_float(iarg, argc, argv, mesh_param.side_direction);
      iarg++;
    }
    else if (s == "-normal") {
      get_float(iarg, argc, argv, mesh_param.normal);
      iarg++;
    }
    else if (s == "-rot_dir") {
      get_float(iarg, argc, argv, mesh_param.rotation_direction);
      iarg++;
    }
    else if (s == "-spacing") {
      get_float(iarg, argc, argv, mesh_param.spacing);
      iarg++;
    }
    else if (s == "-translate") {
      get_float(iarg, argc, argv, mesh_param.translate);
      iarg++;
    }
    else if (s == "-radius") {
      float x = get_float(iarg, argc, argv);
      mesh_param.radius.push_back(x);
      iarg++;
    }
    else if (s == "-length_diff") {
      float x = get_float(iarg, argc, argv);
      mesh_param.length_difference.push_back(x);
      iarg++;
    }
    else if (s == "-wedge_angle") {
      ANGLE_TYPE angle = get_float(iarg, argc, argv);
      mesh_param.wedge_angle.push_back(angle);
      mesh_param.flag_wedge = true;
      iarg++;
    }
    else if (s == "-angle") {
      ANGLE_TYPE angle = get_float(iarg, argc, argv);
      mesh_param.angle.push_back(angle);
      iarg++;
    }
    else if (s == "-wedge_isovalue") {
      COORD_TYPE wedge_isovalue = get_float(iarg, argc, argv);
      mesh_param.wedge_isovalue.push_back(wedge_isovalue);
      iarg++;
    }
    else if (s == "-clip0") {
      float x1, x2;
      get_arg(iarg, argc, argv, x1, x2);
      mesh_param.dist2near0.push_back(x1);
      mesh_param.dist2far0.push_back(x2);
      iarg = iarg+2;
    }
    else if (s == "-seed") {
      int x = get_int(iarg, argc, argv);
      mesh_param.random_seed.Set(x);
      iarg++;
    }
    else if (s == "-triangulate") {
      mesh_param.flag_triangulate = true;
    }
    else if (s == "-poly") {
      mesh_param.flag_triangulate = false;
    }
    else if (s == "-list") {
      cout << "Valid mesh names: " << endl;
      out_mesh_descriptions(mesh_info);
      exit(20);
    }
    else if (s == "-s") {
      flag_silent = true;
    }
    else if (s == "-help") {
      help_msg();
    }
    else {
      cerr << "Usage error.  Illegal parameter: " << s << endl;
      usage_error();
    };
    iarg++;
  };

  if (iarg+1 > argc) { 
    cerr << "Error.  Missing file name." << endl << endl;
    usage_error(); 
  }
  else if (iarg+1 < argc) 
    { usage_error(); }

  ofilename = argv[iarg];
}

void usage_msg()
{
  cerr << "Usage: ijkgenscalar [OPTIONS] <filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-mesh <mesh name>] [-distance <d>]" << endl;
  cerr << "  [-asize <N>] [-tilt | -no_tilt]" << endl;
  cerr << "  [-single | -two | -multi | -n <num_obj>] [-stack]" << endl;
  cerr << "  [-center \"<coord>\" | -grid_center | -randompos <S>]" << endl;
  cerr << "  [-dir \"<coord>\" | -randomdir <S>]" << endl;
  cerr << "  [-side_dir \"<coord>\"] [-rot_dir \"<coord>\"]" << endl;
  cerr << "  [-radius <x>] [-length_diff <x>] [-angle <A>] [-clip0 <N> <F>]" 
       << endl;
  cerr << "  [-flange | -flange_wh <W> <H>] [-wedge | wedge_angle <A>]" 
       << endl;
  cerr << "  [-wedge_isovalue <S>]" << endl;
  cerr << "  [-normal \"<coord>\"] [-translate \"<coord>\"]" << endl;
  cerr << "  [-spacing \"<space_coord>\"]" << endl;
  cerr << "  [-seed <S>]" << endl;
  cerr << "  [-triangulate | -poly]" << endl;
  cerr << "  [-s] [-list] [-help]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cerr << "Usage: ijkgenscalar [OPTIONS] <filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -mesh <mesh name>: Set mesh."
       << "  (Use -list to list meshes.)" << endl;
  cerr << "  -distance <d>:  Distance to center/centerline/circle." << endl;
  cerr << "  -asize <N>: Set number of vertices along each grid axis to <N>."
       << endl;
  cerr << "  -tilt:      Use directions to tilt objects." << endl;
  cerr << "  -no_tilt:   Objects are not tilted." << endl; 
  cerr << "  -single:    Single object (default)." << endl;
  cerr << "  -two:       Two objects." << endl;
  cerr << "  -multi:     Multiple objects." << endl;
  cerr << "  -n <N>:     Number of objects is <N>." << endl;
  cerr << "  -stack:     Stack of translated objects." << endl;
  cerr << "  -center \"<coord>\": Set center to given coordinates." << endl;
  cerr << "  -grid_center: Set center to center of grid." << endl;
  cerr << "  -randompos <S>: Assign random positions to centers." << endl
       << "      Use <S> as random seed for positions." << endl;
  cerr << "  -dir \"<coord>\": Set direction to given coordinates." << endl;
  cerr << "  -randomdir <S>: Assign random directions." << endl
       << "      Use <S> as random seed for positions." << endl;
  cerr << "  -start_dir \"<coord>\": Set start direction to given coordinates." 
       << endl;
  cerr << "  -side_dir \"<coord>\": Set side direction to given coordinates." 
       << endl;
  cerr << "  -rot_dir \"<coord>\": Set rotation direction to given coordinates."
       << endl;
  cerr << "  -radius <x>: Set radius to <x>." << endl;
  cerr << "  -length_diff <x>:  Set (length-width) or (height-width) to <x>."
       << endl;
  cerr << "            Negative <x> sets width greater than length or height."
       << endl;
  cerr << "  -angle <A>:  Set cone, frustrum or edge angle to <A>." << endl;
  cerr << "  -clip0 <N> <F>: Set distance from isovalue 0 frustrum apex" << endl
       << "                    to near clipping plane to <N>." << endl
       << "                  Set distance to far clipping plane to <F>." 
       << endl;
  cerr << "  -flange: Add flange to object." << endl;
  cerr << "  -flange_wh <W> <H>: Set flange width to <W> and height to <H>."
       << endl;
  cerr << "  -wedge: Intersect object with wedge." << endl;
  cerr << "  -wedge_angle <A>: Set wedge angle to <A>." << endl;
  cerr << "  -wedge_isovalue <S>: Set scalar value where wedge passes" << endl
       << "                       through the center to <S>." << endl;
  cerr << "     Note: This determines how the wedge intersects objects."
       << endl;
  cerr << "  -normal \"<coord>\": Set plane normals to given coordinates." 
       << endl;
  cerr << "  -translate \"<coord>\": Set stack translation vector to given coordinates." 
       << endl;
  cerr << "  -spacing \"<space_coord>\": Set spacing between grid vertices"
       << endl
       << "                along each axis." << endl;
  cerr << "  -seed <S>: Set random seed to <S>." << endl;
  cerr << "  -triangulate:  Triangulate all polygons." << endl;
  cerr << "  -poly:  Output arbitrary polygons, no triangulation." << endl;
  cerr << "  -list:  List meshes." << endl;
  cerr << "  -s:     Silent mode." << endl;
  cerr << "  -help:  Print this help message." << endl;

  exit(20);
}


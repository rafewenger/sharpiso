/// \file ijkIO.txx
/// IO templates for reading/writing meshes.
/// - Input formats: Geomview .off.
/// - Output formats: Geomview .off, OpenInventor .iv (3D), Fig .fig (2D).
/// - Version 0.1.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011,2012 Rephael Wenger

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

#ifndef _IJKIO_
#define _IJKIO_

#include <iostream>
#include <vector>
#include <string>

#include "ijk.txx"

namespace IJK {

  // ******************************************
  // Miscellaneous local output routines
  // ******************************************

  namespace {

    /// Output vertex coordinates
    template <typename CTYPE> void ijkoutVertexCoord
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv);

    /// Output polygon vertices
    template <typename VTYPE> void ijkoutPolygonVertices
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump);

    /// Output polytope vertices
    template <typename NTYPE, typename VTYPE, typename ITYPE>
    void ijkoutPolytopeVertices
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly);

    /// Output vector endpoint coordinates.
    template <typename CTYPE0, typename CTYPE1, typename CTYPE2> 
    void ijkoutVectorEndpoints
    (std::ostream & out, const int dim, 
     const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
     const int num_vectors, const CTYPE2 scale);
  }

  // ******************************************
  // Write OpenInventor .iv file
  // ******************************************

  /// Output OpenInventor .iv file.
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // out = output stream
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    IJK::PROCEDURE_ERROR error("ijkoutIV");

    if (dim != 3) {
      throw error
        ("Illegal dimension.  OpenInventor files are only for dimension 3.");
    }

    out << "#Inventor V2.1 ascii" << std::endl;
    out << std::endl;

    out << "Separator {" << std::endl;

    // set vertex ordering to clockwise to turn on two-sided lighting
    out << "  ShapeHints {" << std::endl;
    out << "    vertexOrdering CLOCKWISE" << std::endl;
    out << "  }" << std::endl;
    out << std::endl;

    out << "  IndexedFaceSet {" << std::endl;

    out << "    vertexProperty VertexProperty {" << std::endl;
    out << "      vertex [" << std::endl;

    // output vertex coordinates
    out << std::endl << "# vertex coordinates" << std::endl;
    for (int i = 0; i < numv; i++) {
      for (int d = 0; d < dim; d++) {
        out << coord[dim*i+d];
        if (d < dim-1) { out << " "; }
        else {	
          if (i < numv-1) { out << "," << std::endl; };
        };
      };
    };

    out << " ]" << std::endl;
    out << "    }" << std::endl;

    out << "    coordIndex [" << std::endl;
    // output triangle vertices
    out << std::endl << "# triangle vertices" << std::endl;
    for (int it = 0; it < numt; it++) {
      for (int d = 0; d < dim; d++) {
        out << tri[dim*it+d] << ",";
      };
      out << "-1";
      if (it < numt-1) {
        out << "," << std::endl;
      };
    };
    out << " ]" << std::endl;

    out << "  }" << std::endl;

    out << "}" << std::endl;
  }

  /// Output OpenInventor .iv file to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // output OpenInventor .iv format to standard output
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    ijkoutIV(std::cout, dim, coord, numv, tri, numt); 
  }

  /// Output OpenInventor .iv file.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(out, dim, vector2pointer(coord), coord.size()/dim,
             vector2pointer(simplex_vert), simplex_vert.size()/dim);
  }

  /// Output OpenInventor .iv file to standard output.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutIV
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(dim, vector2pointer(coord), coord.size()/dim,
             vector2pointer(simplex_vert), simplex_vert.size()/dim);
  }	

  // ******************************************
  // Write Geomview OFF file
  // ******************************************

  /// Output Geomview .off file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices.
  ///        simplex_vert[numv_per_simplex*j+k] = k'th vertex of simplex j.
  /// @param nums = Number of simplices
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVertices(out, numv_per_simplex, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, numv, 
              simplex_vert, nums);
  }

  /// \brief Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // set numv_per_simplex to dim
  {
    ijkoutOFF(out, dim, dim, coord, numv, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF format to standard output
    // set numv_per_simplex to dim
  {
    ijkoutOFF(dim, dim, coord, numv, simplex_vert, nums);
  }


  /// Output Geomview .off file.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF files
    // vector input format
  {
    ijkoutOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex);
  }

  /// Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF files
    // vector input format
    // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(out, dim, numv_per_simplex, coord, simplex_vert);
  }


  /// Output Geomview .off file to standard output.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutOFF
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF format to standard output
    // vector input format
    // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file.
  /// Two types of polytopes.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param poly1_vlist = List of vertices of poly 1.
  /// @param numv_per_poly1 = Number of vertices per polygon in poly1_vlist.
  /// @param num_poly1 = Number of poly1 polytopes.
  /// @param poly2_vlist = List of vertices of poly 2.
  /// @param numv_per_poly2 = Number of vertices per polygon in poly2_vlist.
  /// @param num_poly2 = Number of poly2 polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVertices(out, numv_per_poly1, poly1_vlist, num_poly1);
    ijkoutPolygonVertices(out, numv_per_poly2, poly2_vlist, num_poly2);
  }

  /// Output Geomview .off file.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2> void ijkoutOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutOFF(out, dim, vector2pointer(coord), numv,
              vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
              vector2pointer(poly2_vlist), numv_per_poly2, num_poly2);
  }

  /// Output Geomview .off file.
  /// Multiple polytope types.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).

  template <typename CTYPE, typename NTYPE, typename VTYPE, typename ITYPE> 
  void ijkoutPolytopeOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const NTYPE * num_poly_vert, const VTYPE * poly_vert,
   const ITYPE * first_poly_vert, const int num_poly)
  {
    IJK::PROCEDURE_ERROR error("ijkoutPolytopeOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    ijkoutPolytopeVertices
      (out, num_poly_vert, poly_vert, first_poly_vert, num_poly);
  }

  /// Output Geomview .off file. Color vertices.
  template <typename T, typename colorT> void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
    // out = output stream
    // dim = dimension
    // numv_per_simplex = num vertices per simplex
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
    // front_color = array of front colors, 4 entries (RGBA) per vertex
    // back_color = array of backface colors, 4 entries (RGBA) per vertex
    //              May be NULL.
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorVertOFF");

    assert(front_color != NULL);

    out << "C";
    if (back_color != NULL) out << "C";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
        out << coord[iv*dim + d];
        if (d+1 < dim) out << " ";
      }
      out << "  ";
      for (int ic = 0; ic < 4; ic++) {
        out << front_color[4*iv+ic];
        if (ic < 3) out << " ";
      }
      out << "  ";
      if (back_color != NULL) {
        for (int ic = 0; ic < 4; ic++) {
          out << back_color[4*iv+ic];
          if (ic < 3) out << " ";
        }
      };
      out << std::endl;
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
        out << simplex_vert[is*numv_per_simplex + iv];
        if (iv < numv_per_simplex-1) { out << " "; }
        else { out << std::endl; };
      };
    };
  }

  /// Output Geomview .off file to standard output. Color vertices.
  template <typename T, typename colorT> void ijkoutColorVertOFF
  (const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorVertOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color vertices.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE, typename colorT> 
  void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorVertOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color simplices.
  template <typename T, typename colorT> void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
    // front_color = array of front colors, 4 entries (RGBA) per face
    // back_color = array of backface colors, 4 entries (RGBA) per face
    //              May be NULL.
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorFacesOFF");

    assert(front_color != NULL);

    out << "D";
    if (back_color != NULL) out << "D";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
        out << coord[iv*dim + d];
        if (d+1 < dim) out << " ";
      }
      out << std::endl;
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
        out << simplex_vert[is*numv_per_simplex + iv];
        if (iv < numv_per_simplex-1) { out << " "; }
      };
      out << "  ";
      for (int ic = 0; ic < 4; ic++) {
        out << front_color[4*is+ic];
        if (ic < 3) out << " ";
      }
      out << "  ";
      if (back_color != NULL) {
        for (int ic = 0; ic < 4; ic++) {
          out << back_color[4*is+ic];
          if (ic < 3) out << " ";
        }
      };
      out << std::endl;
    };

  }

  /// Output Geomview .off file. Color mesh faces.
  /// Two types of polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2,
            typename colorT> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE1 * poly1_vlist, const int numv_per_poly1, const int num_poly1,
   const colorT * poly1_front_color, const colorT * poly1_back_color,
   const VTYPE2 * poly2_vlist, const int numv_per_poly2, const int num_poly2,
   const colorT * poly2_front_color, const colorT * poly2_back_color)
  {
    const int num_poly = num_poly1 + num_poly2;
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    out << "D";
    if ((poly1_back_color != NULL) && (poly2_back_color != NULL))
      { out << "D"; }
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << num_poly << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutPolygonVerticesColor(out, numv_per_poly1, poly1_vlist, num_poly1, 
                               poly1_front_color, poly1_back_color);
    ijkoutPolygonVerticesColor(out, numv_per_poly2, poly2_vlist, num_poly2,
                               poly2_front_color, poly2_back_color);
  }

  /// Output Geomview .off file to standard output. Color simplices.
  template <typename T, typename colorT> void ijkoutColorFacesOFF
  (const int dim, const int numv_per_simplex, 
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorFacesOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color simplices.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename CTYPE, typename VTYPE, typename colorT> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorFacesOFF
      (out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
       vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  /// Output Geomview .off file. Color faces.
  /// Two types of polytopes.
  template <typename CTYPE, typename VTYPE1, typename VTYPE2, typename colorT> 
  void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE1> & poly1_vlist, const int numv_per_poly1,
   const colorT * poly1_front_color, const colorT * poly1_back_color,
   const std::vector<VTYPE2> & poly2_vlist, const int numv_per_poly2,
   const colorT * poly2_front_color, const colorT * poly2_back_color)
  {
    const int numv = coord.size()/dim;
    const int num_poly1 = poly1_vlist.size()/numv_per_poly1;
    const int num_poly2 = poly2_vlist.size()/numv_per_poly2;

    ijkoutColorFacesOFF
      (out, dim, vector2pointer(coord), numv,
       vector2pointer(poly1_vlist), numv_per_poly1, num_poly1,
       poly1_front_color, poly1_back_color,
       vector2pointer(poly2_vlist), numv_per_poly2, num_poly2,
       poly2_front_color, poly2_back_color);
  }

  /// Output Geomview .off file. Output vertex normals.
  template <typename CTYPE, typename NTYPE, typename VTYPE> 
  void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const NTYPE * normal, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF files with vertex normal information
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    IJK::PROCEDURE_ERROR error("ijkoutOFF");

    if (dim == 3) { out << "NOFF" << std::endl; }
    else if (dim == 4) { out << "N4OFF" << std::endl;}
    else {
      out << "NnOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
        out << coord[iv*dim + d] << " ";
      }
      for (int d = 0; d < dim; d++) {
        out << normal[iv*dim+d];
        if (d < dim-1) { out << " "; }
        else { out << std::endl; };
      }
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
        out << simplex_vert[is*numv_per_simplex + iv];
        if (iv < numv_per_simplex-1) { out << " "; }
        else { out << std::endl; };
      };

    };

  }

  /// Output Geomview .off file. Output vertex normals.
  /// C++ STL vector format for coord[], normal[] and simplex_vert[].
  template <typename CTYPE, typename NTYPE, typename VTYPE> 
  void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<NTYPE> & normal,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutNormalsOFF
      (out, dim, numv_per_simplex,
       vector2pointer(coord), vector2pointer(normal), coord.size()/dim,
       vector2pointer(simplex_vert),
       simplex_vert.size()/numv_per_simplex);
  }

  /// Output quadrilaterals to Geomview .off file.
  /// @param out = Output stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices.
  ///        simplex_vert[numv_per_simplex*j+k] = k'th vertex of simplex j.
  /// @param nums = Number of simplices
  /// @param flag_reorder_vertices = Flag for reordering vertices.
  ///        If true, output vertices in counter-clockwise order around quad.
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (std::ostream & out, const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    IJK::PROCEDURE_ERROR error("ijkoutQuadOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << numq << " " << 0 << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;
    ijkoutQuadVertices
      (out, quad_vert, numq, flag_reorder_vertices);
  }

  /// Output quadrilaterals in Geomview .off format to standard output.
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    ijkoutQuadOFF(std::cout, dim, coord, numv, quad_vert, numq);
  }

  /// Output quadrilaterals to Geomview .off
  /// C++ STL vector format for coord[] and quad_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & quad_vert,
   const bool flag_reorder_vertices)
  {
    const int NUMV_PER_QUAD = 4;

    ijkoutQuadOFF(out, dim, vector2pointer(coord), coord.size()/dim,
                  vector2pointer(quad_vert), quad_vert.size()/NUMV_PER_QUAD,
                  flag_reorder_vertices);
  }

  /// Output quadrilaterals in Geomview .off format to standard output.
  /// C++ STL vector format for coord[] and quad_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutQuadOFF
  (const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & quad_vert,
   const bool flag_reorder_vertices)
  {
    ijkoutQuadOFF(std::cout, dim, coord, quad_vert, flag_reorder_vertices);
  }

  // ******************************************
  // Write Geomview LINE file
  // ******************************************

  /// Output Geomview .line file.
  template <typename CTYPE, typename VTYPE> void ijkoutLINE
  (std::ostream & out, const int dim, 
   const CTYPE * coord, const int numv,
   const VTYPE * edge_vert, const int nume)
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // edge_vert[dim*j+k] = k'th vertex index of edge j
    // nums = number of simplices
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutLINE");

    if (dim == 3) { out << "LINE" << std::endl; }
    else if (dim == 4) { out << "4LINE" << std::endl;}
    else {
      error.AddMessage("Only dimensions 3 and 4 permitted for LINE output.");
      throw error;
    };

    out << numv << " " << nume << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }

  /// Output Geomview .line file.
  /// C++ STL vector format for coord[] and edge_vert[].
  template <typename CTYPE, typename VTYPE> void ijkoutLINE
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & edge_vert)
  // output Geomview LINE file
  // vector input format
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutLINE");

    if (coord.size()%dim != 0) {
      error.AddMessage("Illegal number of vertex coordinates: ",
                       edge_vert.size(), ".");
      error.AddMessage
        ("  Number of coordinates must be divisible by the dimension: ",
         dim, ".");
      throw error;
    }

    if (edge_vert.size()%NUM_EDGE_ENDPOINTS != 0) {
      error.AddMessage("Illegal number of edge endpoints: ",
                       edge_vert.size(), ".");
      error.AddMessage("  Number of edge endpoints must be even.");
      throw error;
    }

    ijkoutLINE(out, dim, vector2pointer(coord), coord.size()/dim,
               vector2pointer(edge_vert), edge_vert.size()/NUM_EDGE_ENDPOINTS);
  }

  /// Output Geomview .line file color header
  template <typename COLOR_TYPE> 
  void ijkoutColorLINEheader
  (std::ostream & out, const int dim, const int numv, const int nume,
   const COLOR_TYPE rgba[4])
  {
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    if (dim == 3) { out << "LINEC" << std::endl; }
    else if (dim == 4) { out << "4LINEC" << std::endl;}
    else {
      error.AddMessage("Only dimensions 3 and 4 permitted for LINE output.");
      throw error;
    };

    out << rgba[0] << " " << rgba[1] << " " 
        << rgba[2] << " " << rgba[3] << std::endl;
    out << numv << " " << nume << std::endl;
  }


  /// Output Geomview .line file with color
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorLINE
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * edge_vert, const int nume, const COLOR_TYPE rgba[4])
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // edge_vert[dim*j+k] = k'th vertex index of edge j
    // nums = number of simplices
    // rgba[] = array of R,G,B,A values
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    ijkoutColorLINEheader(out, dim, numv, nume, rgba);

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }

  /// Output Geomview .line file.
  /// C++ STL vector format for coord[] and edge_vert[].
  template <typename CTYPE, typename VTYPE, typename COLOR_TYPE> 
  void ijkoutColorLINE
  (std::ostream & out, const int dim,
   const std::vector<CTYPE> & coord, const std::vector<VTYPE> & edge_vert,
   const COLOR_TYPE rgba[4])
  {
    const int NUM_EDGE_ENDPOINTS = 2;
    IJK::PROCEDURE_ERROR error("ijkoutColorLINE");

    if (coord.size()%dim != 0) {
      error.AddMessage("Illegal number of vertex coordinates: ",
                       edge_vert.size(), ".");
      error.AddMessage
        ("  Number of coordinates must be divisible by the dimension: ",
         dim, ".");
      throw error;
    }

    if (edge_vert.size()%NUM_EDGE_ENDPOINTS != 0) {
      error.AddMessage("Illegal number of edge endpoints: ",
                       edge_vert.size(), ".");
      error.AddMessage("  Number of edge endpoints must be even.");
      throw error;
    }

    ijkoutColorLINE
      (out, dim, vector2pointer(coord), coord.size()/dim,
       vector2pointer(edge_vert), edge_vert.size()/NUM_EDGE_ENDPOINTS, rgba);
  }

  /// Output vectors to Geomview .line file.
  /// @param num_vectors Number of vectors
  /// @param scale Multiply all vectors by scale.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename COLOR_TYPE> 
  void ijkoutVectorsLINE
  (std::ostream & out, const int dim, 
   const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
   const int num_vectors, const CTYPE2 scale, const COLOR_TYPE rgba[4])
  {
    ijkoutColorLINEheader(out, dim, 2*num_vectors, num_vectors, rgba);

    ijkoutVectorEndpoints
      (out, dim, point_coord, dir_coord, num_vectors, scale);

    // Output line segments
    for (int iv = 0; iv < num_vectors; iv++) 
      { out << 2*iv << " " << 2*iv+1 << std::endl; }
  }

  /// Output vectors to Geomview .line file.
  /// @param num_vectors Number of vectors
  /// @param scale Multiply all vectors by scale.
  /// C++ STL vector format for point_coord[] and dir_coord[].
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename COLOR_TYPE> 
  void ijkoutVectorsLINE
  (std::ostream & out, const int dim, 
   const std::vector<CTYPE0> & point_coord, 
   const std::vector<CTYPE1> & dir_coord, 
   const CTYPE2 scale, const COLOR_TYPE rgba[4])
  {
    IJK::PROCEDURE_ERROR error("ijkoutVectorsLINE");

    if (dim <= 0) {
      error.AddMessage("Illegal dimension: ", dim, ".");
      throw error;
    }

    if (point_coord.size() != dir_coord.size()) {
      error.AddMessage
        ("Error.  Number of points does not equal number of directions.");
      error.AddMessage("Number of points: ", point_coord.size()/dim, ".");
      error.AddMessage("Number of directions: ", dir_coord.size()/dim, ".");
      throw error;
    }

    ijkoutVectorsLINE
      (out, dim, vector2pointer(point_coord), vector2pointer(dir_coord),
       point_coord.size()/dim, scale, rgba);
  }

  // ******************************************
  // Read Geomview OFF file
  // ******************************************

  // local namespace
  namespace {
    inline void gobble_line(std::istream & in)
    {
      while (in.good() && !in.eof() && in.get() != '\n') {};
    }
  }

  /// \brief Read header of Geomview .off file.
  /// @param in Input stream.
  /// @param dim Vertex dimension.
  /// @param numv Number of vertices.
  /// @param nums Number of simplices.
  /// @param nume Number of edges.
  /// @param flag_normals Normals flag.  True if file contains vertex normals.
  inline void ijkinOFFheader
  (std::istream & in, int & dim, int & numv, int & nums, int & nume,
   bool & flag_normals)
  {
    std::string header_keyword;
    IJK::PROCEDURE_ERROR error("ijkinOFF");

    // Initialize
    dim = 0;
    numv = 0;
    nums = 0;
    nume = 0;
    flag_normals = false;

    if (!in.good()) {
      error.AddMessage("Error reading from input stream in.");
      throw error;
    }

    // read input
    in >> header_keyword;

    if (!in.good()) {
      error.AddMessage
        ("Error reading from header keyword from input stream in.");
      throw error;
    }

    bool flag_valid_header = true;
    int hksize = header_keyword.size() ;
    if (hksize < 3) { flag_valid_header = false; }
    else if (header_keyword.substr(hksize-3, 3) != "OFF")
      { flag_valid_header = false; }
    else {
      dim = 3;  // default dimension

      if (hksize > 3) {
        if (header_keyword.substr(hksize-4, 4) == "4OFF")
          { dim = 4; }
        else if (header_keyword.substr(hksize-4, 4) == "nOFF") 
          { in >> dim; }
      }

      std::string prefix = header_keyword.substr(0, hksize-3);

      int prefix_size = prefix.size();
      if ( prefix_size >= 1) {

        if (prefix[prefix_size-1] == '3' ||
            prefix[prefix_size-1] == '4' ||
            prefix[prefix_size-1] == 'n') {
          prefix = prefix.substr(0, prefix_size-1); 
          prefix_size = prefix_size-1;
        }
      }

      if (prefix_size >= 1) {

        if (prefix[prefix_size-1] == 'N') 
          { flag_normals = true; }
      }
    }

    if (!flag_valid_header) {
      error.AddMessage("Illegal Geomview .off file header: "+
                       header_keyword);
      throw error;
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nums;
    in >> nume;

    if (!in.good()) {
      error.AddMessage("Error reading number of vertices and polyhedra from input stream in.");
      throw error;
    }
  }

  /// \brief Read header of Geomview .off file.
  /// @param in = Input stream.
  /// @param dim = Vertex dimension.
  /// @param numv = Number of vertices.
  /// @param nums = Number of simplices.
  /// @param nume = Number of edges.
  inline void ijkinOFFheader
  (std::istream & in, int & dim, int & numv, int & nums, int & nume)
  {
    bool flag_normals;

    ijkinOFFheader(in, dim, numv, nums, nume, flag_normals);
  }

  /// \brief Read coordinates from Geomview .off file.
  ///
  /// Ignores any normal, color information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @pre Array coord[] is preallocated with size at least \a numv * \a dim.
  template <typename T> void ijkinOFFcoord
  (std::istream & in, const int dim, const int numv, T * coord)
  {
    IJK::PROCEDURE_ERROR error("ijkinOFFcoord");

    if (!in.good()) {
      error.AddMessage("Error reading coordinates from input stream in.");
      throw error;
    }

    if (coord == NULL) 
      if (numv > 0 && dim > 0) {
        error.AddMessage("Programming error. Array coord[] not allocated.");
        throw error;
      }

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
        { in >> coord[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading coordinates of vertex ", iv,
                         " from input stream in.");
        throw error;
      }
    }
  }

  /// \brief Read coordinates and normals from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in Input stream.
  /// @param dim Dimension of vertices.
  /// @param numv Number of vertices.
  /// @param coord Array of coordinates. 
  ///        coord[dim*i+k] k'th coordinate of vertex i (k < dim).
  /// @pre Array coord[] is preallocated with size at least \a numv * \a dim.
  /// @param normal Array of normals.
  ///        normal[dim*i+k] k'th coordinate of normal of vertex i (k < dim).
  /// @pre Array normal[] is preallocated with size at least \a numv * \a dim.
  template <typename CTYPE, typename NTYPE> void ijkinNOFFcoord
  (std::istream & in, const int dim, const int numv, 
   CTYPE * coord, NTYPE * normal)
  {
    IJK::PROCEDURE_ERROR error("ijkinNOFFcoord");

    if (!in.good()) {
      error.AddMessage("Error reading coordinates from input stream in.");
      throw error;
    }

    if (coord == NULL) 
      if (numv > 0 && dim > 0) {
        error.AddMessage("Programming error. Array coord[] not allocated.");
        throw error;
      }

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
        { in >> coord[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading coordinates of vertex ", iv,
                         " from input stream in.");
        throw error;
      }

      for (int d = 0; d < dim; d++) 
        { in >> normal[iv*dim + d]; };

      if (!in.good()) {
        error.AddMessage("Error reading normals of vertex ", iv,
                         " from input stream in.");
        throw error;
      }
    }
  }

  /// \brief Read coordinates from Geomview .off file.
  template <typename T> void ijkinOFFcoord
  (std::istream & in, const int dim, const int numv, std::vector<T> & coord)
  {
    if (dim < 1) {
      coord.clear();
      return;
    }

    coord.resize(dim*numv);

    ijkinOFFcoord(in, dim, numv, &(coord[0]));
  }

  /// \brief Read coordinates from Geomview .off file.
  template <typename CTYPE, typename NTYPE> void ijkinNOFFcoord
  (std::istream & in, const int dim, const int numv, 
   std::vector<CTYPE> & coord, std::vector<NTYPE> & normal)
  {
    if (dim < 1) {
      coord.clear();
      normal.clear();
      return;
    }

    coord.resize(dim*numv);
    normal.resize(dim*numv);

    ijkinNOFFcoord(in, dim, numv, &(coord[0]), &(normal[0]));
  }

  /// \brief Read vertex indices from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param numv = Number of vertices.
  /// @param vert = Array of vertices.
  ///        vert[k] = k'th vertex.
  /// @pre Array vert[] has been preallocated
  ///      with size at least \a numv.
  template <typename T> void ijkinOFFvert
  (std::istream & in, const int numv, T * vert)
  {
    for (int k = 0; k < numv; k++) 
      { in >> vert[k]; }
    gobble_line(in);
  }

  /// \brief Read polytope vertices of polytope \a js from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param jpoly = Polytope index.
  /// @param numv_per_poly = Number of vertices per polytope.
  /// @param poly_vert = Array of polytope vertices. 
  ///        poly_vert[numv_per_poly*js+k] = 
  ///            k'th vertex index of polytope jp.
  /// @pre Array poly_vert[] has been preallocated
  ///      with size at least \a numv_per_poly * \a (jpoly+1).
  template <typename T> void ijkinOFFpolyVert
  (std::istream & in, const int jpoly, const int numv_per_poly,
   T * poly_vert)
  {
    for (int k = 0; k < numv_per_poly; k++) 
      { in >> poly_vert[jpoly*numv_per_poly + k]; }
    gobble_line(in);
  }

  /// \brief Read simplex vertices of polytope \a jpoly from Geomview .off file.
  /// C++ STL vector format for poly_vert[].
  /// @pre C++ vector poly_vert[] has size at least 
  ///        \a numv_per_poly * \a (jpoly+1).
  template <typename T> void ijkinOFFpolyVert
  (std::istream & in, const int jpoly, const int numv_per_poly,
   std::vector<T> & poly_vert)
  {
    ijkinOFFpolyVert(in, jpoly, numv_per_poly, &(poly_vert[0]));
  }

  /// \brief Read \a nums simplex vertices starting at simplex \a ifirst
  ///        from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param ifirst = Index of first simplex to read.
  /// @param nums = Number of simplices to read.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[numv_per_simplex*js+k] = 
  ///            k'th vertex index of simplex js.
  /// @pre Array simplex_vert[] has been preallocated
  ///      with size at least \a numv_per_simplex * \a (js+1).
  template <typename T> void ijkinOFFsimplex
  (std::istream & in, const int ifirst, const int nums, 
   const int numv_per_simplex, T * simplex_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinOFFsimplex");

    int nvert;
    for (int i = ifirst; i < ifirst + nums; i++) {

      in >> nvert;

      // Check that number of vertices is correct.
      if (nvert != numv_per_simplex) {
        error.AddMessage("Simplex ", i, " has ", nvert, " vertices.");
        error.AddMessage("All simplices are expected to have ",
                         numv_per_simplex, " vertices.");
        throw error;
      }

      ijkinOFFpolyVert(in, i, numv_per_simplex, simplex_vert);

      if (!in.good()) {
        error.AddMessage("Error reading vertices of simplex ", i, ".");
        throw error;
      }
    }
  }

  /// \brief Read \a nums simplex vertices starting at simplex \a ifirst
  ///        from Geomview .off file.
  /// C++ STL vector format for simplex_vert[].
  /// @pre C++ vector simplex_vert[] has been preallocated
  ///      with size at least \a numv_per_simplex * \a (js+1).
  template <typename T> void ijkinOFFsimplex
  (std::istream & in, const int ifirst, const int nums, 
   const int numv_per_simplex, std::vector<T> & simplex_vert)
  {
    ijkinOFFsimplex(in, ifirst, nums, numv_per_simplex,
                    &(simplex_vert.front()));
  }


  /// \brief Read Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param mesh_dim = Dimension of mesh.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[dim*j+k] = k'th vertex index of simplex j.
  /// @param nums = Number of simplices.
  /// @pre All simplices have the same dimension.
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   T * & coord, int & numv, int * & simplex_vert, int & nums)
  {
    int nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;
    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert = new int[nums*num_simplex_vert];

      // read in first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// \brief Read Geomview .off file.
  /// C++ STL vector format for coord[] and simplex_vert[].
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   std::vector<T> & coord, std::vector<int> & simplex_vert)
  {
    int numv, nums, nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord.clear();
    simplex_vert.clear();

    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert.resize(nums*num_simplex_vert);

      // read vertices of first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      // read remaining simplex vertices
      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// Read Geomview .off file from standard input.
  template <typename T> void ijkinOFF
  (int & dim, int & mesh_dim,
   T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF format to standard input
    // in = input stream
    // dim = dimension
    // mesh_dim = mesh dimension
    //   Precondition: All simplices have the same dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkinOFF(std::cin, dim, mesh_dim, coord, numv, simplex_vert, nums);
  }


  /// \brief Read Geomview .off file where each simplex has dimension \a dim-1.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[dim*j+k] = k'th vertex index of simplex j.
  /// @param nums = Number of simplices.
  /// @pre All simplices have the same dimension.
  template <typename T> void ijkinOFF
  (std::istream & in, int & dim, T * & coord, int & numv,
   int * & simplex_vert, int & nums)
  {
    const int numv_per_simplex = dim;
    int nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;

    if (!in.good()) {
      throw error("Error: Corrupted input stream. Unable to read input.");
    }

    ijkinOFFheader(in, dim, numv, nums, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    simplex_vert = new int[nums*dim];
    ijkinOFFsimplex(in, 0, nums, numv_per_simplex, simplex_vert);
  }

  /// \brief Read Geomview .off file where each simplex has dimension \a dim-1 
  ///   from standard input.
  template <typename T> void ijkinOFF
  (int & dim, T * & coord, int & numv, int * & simplex_vert, int & nums)
  {
    ijkinOFF(std::cin, dim, coord, numv, simplex_vert, nums);
  }

  /// \brief Read Geomview .off file.
  /// Read normal information from files with NOFF header.
  /// C++ STL vector format for coord[], normal[], and simplex_vert[].
  template <typename CTYPE, typename NTYPE> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
   std::vector<CTYPE> & coord, std::vector<NTYPE> & normal, 
   std::vector<int> & simplex_vert)
  {
    bool flag_normals;
    int numv, nums, nume;

    IJK::PROCEDURE_ERROR error("ijkinOFF");

    coord.clear();
    normal.clear();
    simplex_vert.clear();

    mesh_dim = 0;            // default mesh dimension

    ijkinOFFheader(in, dim, numv, nums, nume, flag_normals);

    if (flag_normals) {
      ijkinNOFFcoord(in, dim, numv, coord, normal);
    }
    else {
      ijkinOFFcoord(in, dim, numv, coord);
    }

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
        mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert.resize(nums*num_simplex_vert);

      // read vertices of first simplex
      ijkinOFFpolyVert(in, 0, num_simplex_vert, simplex_vert);

      // read remaining simplex vertices
      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// \brief Read triangles and quadrilaterals from Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param tri_vert = Array of triangle vertices. 
  ///        tri_vert[3*j+k] = k'th vertex index of triangle j.
  /// @param nums = Number of simplices.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        simplex_vert[4*j+k] = k'th vertex index of quad j.
  /// @param numq = Number of quadrilaterals.
  template <typename T> 
  void ijkinOFF_tri_quad
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & tri_vert, std::vector<int> & quad_vert)
  {
    const int NUM_TRIANGLE_VERTICES = 3;
    const int NUM_QUAD_VERTICES = 4;
    IJK::PROCEDURE_ERROR error("ijkinOFF_tri_quad");

    int nume;

    coord.clear();
    tri_vert.clear();
    quad_vert.clear();

    // nump: Total number of triangle and quads.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    int numt = 0;
    int numq = 0;
    for (int i = 0; i < nump; i++) {
      int nvert;

      in >> nvert;

      if (nvert == NUM_TRIANGLE_VERTICES) {
        tri_vert.resize(NUM_TRIANGLE_VERTICES*(numt+1));
        ijkinOFFpolyVert(in, numt, NUM_TRIANGLE_VERTICES, tri_vert);
        numt++;
      }
      else if (nvert == NUM_QUAD_VERTICES) {
        quad_vert.resize(NUM_QUAD_VERTICES*(numq+1));
        ijkinOFFpolyVert(in, numq, NUM_QUAD_VERTICES, quad_vert);
        numq++;
      }
      else {
        // Illegal number of vertices.
        error.AddMessage("Polygon ", i, " has ", nvert, " vertices.");
        error.AddMessage("All polygons are expected to have 3 or 4 vertices.");
        throw error;
      }

    }

  }

  /// \brief Read polygons from Geomview .off file.
  /// \brief Read triangles and quadrilaterals into separate arrays.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param tri_vert = Array of triangle vertices. 
  ///        tri_vert[3*j+k] = k'th vertex index of triangle j.
  /// @param quad_vert = Array of quadrilateral vertices.
  ///        quad_vert[4*j+k] = k'th vertex index of quad j.
  /// @param num_poly_vert = Array of number of polygon vertices.
  ///        num_poly_vert[i] = Number of vertices of polygon i.
  /// @param poly_vert = Array of polygon vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polygon.
  ///        Polygon i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolyOFF
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & tri_vert, std::vector<int> & quad_vert,
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    const int NUM_TRIANGLE_VERTICES = 3;
    const int NUM_QUAD_VERTICES = 4;
    IJK::PROCEDURE_ERROR error("ijkinPolyOFF");

    int nume;

    coord.clear();
    tri_vert.clear();
    quad_vert.clear();
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    // nump: Total number of polygons.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    int numt = 0;
    int numq = 0;
    int num_not_qt = 0;   // Number of polygons which are not tri or quad
    for (int i = 0; i < nump; i++) {
      int nvert;

      in >> nvert;

      if (nvert == NUM_TRIANGLE_VERTICES) {
        tri_vert.resize(NUM_TRIANGLE_VERTICES*(numt+1));
        ijkinOFFpolyVert(in, numt, NUM_TRIANGLE_VERTICES, tri_vert);
        numt++;
      }
      else if (nvert == NUM_QUAD_VERTICES) {
        quad_vert.resize(NUM_QUAD_VERTICES*(numq+1));
        ijkinOFFpolyVert(in, numq, NUM_QUAD_VERTICES, quad_vert);
        numq++;
      }
      else {
        num_poly_vert.push_back(nvert);
        int k = poly_vert.size();
        first_poly_vert.push_back(k);
        poly_vert.resize(k+nvert);
        ijkinOFFvert(in, nvert, &(poly_vert[k]));
        num_not_qt++;
      }
    }
  }

  /// \brief Read list of polytopes from Geomview .off file.
  /// @param in = Input stream.
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename NTYPE, typename VTYPE, typename ITYPE>
  void ijkinOFFpolyList
  (std::istream & in, const int nump, std::vector<NTYPE> & num_poly_vert,
   std::vector<VTYPE> & poly_vert, std::vector<ITYPE> & first_poly_vert)
  {
    for (int i = 0; i < nump; i++) {
      NTYPE nvert;

      in >> nvert;
      num_poly_vert.push_back(nvert);
      int k = poly_vert.size();
      first_poly_vert.push_back(k);
      poly_vert.resize(k+nvert);
      ijkinOFFvert(in, nvert, &(poly_vert[k]));
    }
  }

  /// \brief Read polytopes from Geomview .off file.
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolytopeOFF
  (std::istream & in, int & dim, T * & coord, int & numv, 
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinPolytopeOFF");

    int nume, nump;

    coord = NULL;
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    ijkinOFFheader(in, dim, numv, nump, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    ijkinOFFpolyList
      (in, nump, num_poly_vert, poly_vert, first_poly_vert);
  }

  /// \brief Read polytopes from Geomview .off file.
  /// C++ STL vector format for coord[].
  ///
  /// Ignores any color, normal information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param num_poly_vert = Array of number of polytope vertices.
  ///        num_poly_vert[i] = Number of vertices of polytope i.
  /// @param poly_vert = Array of polytope vertices.
  /// @param first_poly_vert = Array indexing first vertex of each polytope.
  ///        Polytope i has vertices poly_vert[first_poly_vert[i]] to
  ///        poly_vert[first_poly_vert[i]+num_poly_vert[i]-1].
  template <typename T> 
  void ijkinPolytopeOFF
  (std::istream & in, int & dim, std::vector<T> & coord, 
   std::vector<int> & num_poly_vert,
   std::vector<int> & poly_vert, std::vector<int> & first_poly_vert)
  {
    IJK::PROCEDURE_ERROR error("ijkinPolytopeOFF");

    int nume;

    coord.clear();
    num_poly_vert.clear();
    poly_vert.clear();
    first_poly_vert.clear();

    // nump: Total number of polytopes.
    int numv, nump;
    ijkinOFFheader(in, dim, numv, nump, nume);

    ijkinOFFcoord(in, dim, numv, coord);

    ijkinOFFpolyList
      (in, nump, num_poly_vert, poly_vert, first_poly_vert);
  }

  // ******************************************
  // Read Geomview LINE file
  // ******************************************

  /// \brief Read header of Geomview .line file.
  /// Ignores line color.
  /// @param in = Input stream.
  /// @param dim = Vertex dimension.
  /// @param numv = Number of vertices.
  /// @param nume = Number of edges.
  inline void ijkinLINEheader
  (std::istream & in, int & dim, int & numv, int & nume)
  {
    float r,g,b,a;
    std::string header_keyword;
    IJK::PROCEDURE_ERROR error("ijkinLINE");
    if (!in.good()) {
      error.AddMessage("Error reading from input stream in.");
      throw error;
    }

    // read input
    in >> header_keyword;

    if (!in.good()) {
      error.AddMessage
        ("Error reading from header keyword from input stream in.");
      throw error;
    }

    if (header_keyword == "LINEC") {
      dim = 3;
      in >> r >> g >> b >> a;

      if (!in.good()) {
        error.AddMessage("Error reading number of rgba edge color from input stream in.");
        throw error;
      }
    }
    else if (header_keyword == "LINE") {
      dim = 3;
    }
    else {
      error.AddMessage("Illegal Geomview .off file header: "+
                       header_keyword);
      throw error;
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nume;

    if (!in.good()) {
      error.AddMessage("Error reading number of vertices and edges from input stream in.");
      throw error;
    }
  }

  /// \brief Read \a nume edges from Geomview .off file.
  ///
  /// @param in = Input stream.
  /// @param nume = Number of edges to read.
  /// @param edge_endpoint = Array of edge endpoints. 
  ///        edge_enpdoint[2*je+k] = index of k'th endpoint of edge je.
  /// @pre Array edge_endpoint[] has been preallocated
  ///      with size at least 2 * \a nume.
  template <typename T> void ijkinLINEedge
  (std::istream & in, const int nume, T * & edge_endpoint)
  {
    IJK::PROCEDURE_ERROR error("ijkinLINEedge");

    for (int i = 0; i < nume; i++) {

      for (int k = 0; k < 2; k++) 
        { in >> edge_endpoint[2*i + k]; }
      gobble_line(in);

      if (!in.good()) {
        error.AddMessage("Error reading endpoints of edge ", i, ".");
        throw error;
      }
    }
  }

  /// \brief Read Geomview .line file.
  ///
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param coord = Array of coordinates. 
  ///                coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @param numv = Number of vertices.
  /// @param edge_endpoint = Array of edge endpoints. 
  ///        edge_endpoint[dim*j+k] = Vertex index of endpoint k of edge j.
  /// @param nume = Number of edges.
  template <typename T> void ijkinLINE
  (std::istream & in, int & dim, T * & coord, int & numv, 
   int * & edge_endpoint, int & nume)
  {
    IJK::PROCEDURE_ERROR error("ijkinLINE");

    coord = NULL;
    edge_endpoint = NULL;

    ijkinLINEheader(in, dim, numv, nume);

    coord = new T[numv*dim];
    ijkinOFFcoord(in, dim, numv, coord);

    int num_edge_endpoint = 0;
    edge_endpoint = new int[2*nume];

    // read in first edge
    ijkinLINEedge(in, nume, edge_endpoint);
  }

  // ******************************************
  // Fig file
  // ******************************************

  /// \brief Output Fig .fig file.
  /// @param out = output stream
  /// @param dim = dimension.  Must be 2.
  /// @param coord[2*i+k] = k'th coordinate of vertex j  (k < dim)
  /// @param numv = number of vertices 
  /// @param seg[2*j+k] = k'th vertex index of segment j (k < dim)
  /// @param nums = number of line segments
  /// @param flag_polyline: if true, join segments to form FIG polylines
  ///                if false, output segments individually
  /// @param flag_metric: if true, output "Metric" as units
  ///              if false, output "Inches"
  template <typename CTYPE, typename VTYPE, typename SCALE_TYPE> void ijkoutFIG
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord,
   const bool flag_polyline = false, const bool flag_metric = false)
  {
    IJK::PROCEDURE_ERROR error("ijkoutFIG");
    const int cap_style = 1;             // round cap

    if (dim != 2)
      throw error("Illegal dimension.  Fig files are only for dimension 2.");

    // FIG header
    out << "#FIG 3.2" << std::endl;

    out << "Landscape" << std::endl;
    out << "Center" << std::endl;
    if (flag_metric) { out << "Metric" << std::endl; }
    else { out << "Inches" << std::endl; }
    out << "Letter" << std::endl;
    out << "100.00" << std::endl;
    out << "Single" << std::endl;
    out << "-2" << std::endl;
    out << "1200 2" << std::endl;

    if (!flag_polyline) {
      // output each line segment separately line segments

      for (int i = 0; i < nums; i++) {
        out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
            << cap_style << " -1 0 0 2" << std::endl;
        out << "    ";
        for (int j = 0; j < 2; j++) {
          int iv = seg[2*i+j];
          for (int k = 0; k < 2; k++) {
            int x = int(scale_coord*coord[2*iv+k]);
            out << "  " << x;
          }
        }
        out << std::endl;
      }
    }
    else {
      ijkoutFIGpolyline(out, dim, coord, numv, seg, nums, scale_coord);
    }
  }

  // local namespace
  namespace {
    inline int index_other_endpoint(int k)
    {
      int kseg = k/2;
      int j = (k+1)%2;
      return(2*kseg+j);
    }
  }

  // Output Fig polygonal lines.
  template <typename CTYPE, typename VTYPE, typename SCALE_TYPE> 
  void ijkoutFIGpolyline
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord)
  {
    using std::vector;
    const int cap_style = 1;             // round cap
    IJK::ARRAY<bool> is_edge_processed(nums, false);

    // degree of vertex iv
    IJK::ARRAY<int> degree(numv, 0);

    // next[2*i] = next element (in next[]) around vertex 0 of segment i
    // next[2*i+1] = next element (in next[]) around vertex 1 of segment i
    IJK::ARRAY<int> next(2*nums);

    // first[iv] = first element in next[] incident on vertex iv
    // last[iv] = last element in next[] incident on vertex iv
    IJK::ARRAY<int> first(numv);
    IJK::ARRAY<int> last(numv);

    // set degree[], first[], last[], next[]
    for (int i = 0; i < nums; i++)
      for (int j = 0; j < 2; j++) {
        VTYPE iv = seg[2*i+j];
        if (degree[iv] == 0) {
          first[iv] = 2*i+j;
          last[iv] = 2*i+j;
          next[2*i+j] = 2*i+j;
        }
        else {
          next[2*i+j] = next[last[iv]];
          next[last[iv]] = 2*i+j;
          last[iv] = 2*i+j;
        }
        degree[iv]++;
      }

    vector<VTYPE> vlist;
    for (int i = 0; i < nums; i++) {
      if (is_edge_processed[i]) { continue; };

      vlist.clear();
      VTYPE iv0 = seg[2*i];
      VTYPE iv1 = seg[2*i+1];
      vlist.push_back(iv0);
      vlist.push_back(iv1);
      is_edge_processed[i] = true;

      iv0 = iv1;
      int k = next[2*i+1];
      k = index_other_endpoint(k);
      while (!is_edge_processed[k/2]) {
        iv1 = seg[k];
        vlist.push_back(iv1);
        is_edge_processed[k/2] = true;
        iv0 = iv1;
        k = next[k];
        k = index_other_endpoint(k);
      }

      out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
          << cap_style << " -1 0 0 " << vlist.size() << std::endl;

      out << "    ";
      for (int j = 0; j < vlist.size(); j++) {
        int iv = vlist[j];
        for (int k = 0; k < 2; k++) {
          int x = int(scale_coord*coord[2*iv+k]);
          out << "  " << x;
        }
      }
      out << std::endl;
    }

  }

  // ******************************************
  // Miscellaneous local output routines
  // ******************************************

  namespace {

    /// Output vertex coordinates
    /// @param out Output stream.
    /// @param dim Vertex dimension (number of vertex coordinates.)
    /// @param coord[] Array of vertex coordinates.
    ///                coord[dim*i+k] = k'th coordinate of vertex j  (k < dim).
    /// @param numv Number of vertices.
    template <typename CTYPE> void ijkoutVertexCoord
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv)
    {
      for (int iv = 0; iv < numv; iv++) {
        for (int d = 0; d < dim; d++) {
          out << coord[iv*dim + d];
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

    /// Output polygon vertices
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE> void ijkoutPolygonVertices
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv];
          if (iv < numv_per_polygon-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

    /// Output polytope vertices
    template <typename NTYPE, typename VTYPE, typename ITYPE>
    void ijkoutPolytopeVertices
    (std::ostream & out, 
     const NTYPE * num_poly_vert, const VTYPE * poly_vert,
     const ITYPE * first_poly_vert, const int num_poly)
    {
      for (int ipoly = 0; ipoly < num_poly; ipoly++) {
        const ITYPE * pvert = poly_vert + first_poly_vert[ipoly];
        NTYPE num_pvert = num_poly_vert[ipoly];

        out << num_pvert << " ";
        for (int i = 0; i < num_pvert; i++) {
          out << pvert[i];
          if (i+1 < num_pvert) { out << " "; }
          else { out << std::endl; };
        }
      }

    }

    /// Output rgba front and back color.
    template <typename ITYPE, typename colorT>
    void ijkoutColor(std::ostream & out, const ITYPE i, 
                     const colorT * front_color, const colorT * back_color)
    {
      for (int ic = 0; ic < 3; ic++) 
        { out << front_color[4*i+ic] << " "; }
      out << front_color[4*i+3] << "  ";
      for (int ic = 0; ic < 3; ic++)
        { out << back_color[4*i+ic] << " "; }
      out << back_color[4*i+3];
    }



    /// Output polygon vertices and polygon color.
    /// @param out Output stream.
    /// @param numv_per_polygon Number of vertices per polygon.
    /// @pre     All polygons have the same number of vertices.
    /// @param poly_vert[] Array of simplex vertices.
    ///        poly_vert[dim*j+k] = k'th vertex index of polygon j.
    /// @param nump Number of polygons.
    template <typename VTYPE, typename colorT> 
    void ijkoutPolygonVerticesColor
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump,
     const colorT * front_color, const colorT * back_color)
    {
      for (int is = 0; is < nump; is++) {
        out << numv_per_polygon << " ";
        for (int iv = 0; iv < numv_per_polygon; iv++) {
          out << poly_vert[is*numv_per_polygon + iv] << " ";
        }
        out << " ";  // Add another space
        ijkoutColor(out, is, front_color, back_color);
        out << std::endl;
      }
    }

    /// Output quad vertices.
    /// @param out Output stream.
    /// @param quad_vert[] Array of quad vertices.
    ///        quad_vert[4*j+k] = k'th vertex index of polygon j.
    /// @param numq Number of quadrilaterals.
    ///  @param flag_reorder_vertices = if true, change vertex order to be
    ///                         counter-clockwise around quad.
    template <typename VTYPE> void ijkoutQuadVertices
    (std::ostream & out, const VTYPE * quad_vert, const int numq,
     const bool flag_reorder_vertices)
    {
      const int NUMV_PER_QUAD = 4;

      if (flag_reorder_vertices) {
        for (int iq = 0; iq < numq; iq++) {
          out << NUMV_PER_QUAD << " ";
          const VTYPE * v = quad_vert+iq*NUMV_PER_QUAD;
          out << v[0] << " " << v[1] << " ";

          // Note change in order between v[2] and v[3]
          out << v[3] << " " << v[2] << std::endl;
        }
      }
      else {
        ijkoutPolygonVertices(out, NUMV_PER_QUAD, quad_vert, numq);
      }

    }

    /// Output vector endpoint coordinates.
    /// @param out Output stream.
    template <typename CTYPE0, typename CTYPE1, typename CTYPE2> 
    void ijkoutVectorEndpoints
    (std::ostream & out, const int dim, 
     const CTYPE0 * point_coord, const CTYPE1 * dir_coord, 
     const int num_vectors, const CTYPE2 scale)
    {
      for (int iv = 0; iv < num_vectors; iv++) {

        // Output point iv.
        for (int d = 0; d < dim; d++) {
          out << point_coord[iv*dim + d];
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }

        // Output coordinates of (point iv)+(direction iv)
        for (int d = 0; d < dim; d++) {
          int j = iv*dim + d;
          out << (point_coord[j] + dir_coord[j]*scale);
          if (d < dim-1) { out << " "; }
          else { out << std::endl; };
        }
      }
    }

  }

}

#endif

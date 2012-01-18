/// \file ijkIO.txx
/// IO templates for reading/writing meshes.
/// - Input formats: Geomview .off.
/// - Output formats: Geomview .off, OpenInventor .iv (3D), Fig .fig (2D).
/// - Version 0.1.1

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
    template <class CTYPE> void ijkoutVertexCoord
    (std::ostream & out, const int dim, const CTYPE * coord, const int numv);

    /// Output polygon vertices
    template <class VTYPE> void ijkoutPolygonVertices
    (std::ostream & out, const int numv_per_polygon,
     const VTYPE * poly_vert, const int nump);

  }

  // ******************************************
  // Write OpenInventor .iv file
  // ******************************************

  /// Output OpenInventor .iv file.
  template <class CTYPE, class VTYPE> void ijkoutIV
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

    if (dim != 3)
      throw error("Illegal dimension.  OpenInventor files are only for dimension 3.");

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
  template <class CTYPE, class VTYPE> void ijkoutIV
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
  template <class CTYPE, class VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(out, dim, vector2pointer(coord), coord.size()/dim,
             vector2pointer(simplex_vert), simplex_vert.size()/dim);
  }

  /// Output OpenInventor .iv file to standard output.
  template <class CTYPE, class VTYPE> void ijkoutIV
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
  template <class CTYPE, class VTYPE> void ijkoutOFF
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
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, numv, 
              simplex_vert, nums);
  }

  /// \brief Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // set numv_per_simplex to dim
  {
    ijkoutOFF(out, dim, dim, coord, numv, simplex_vert, nums);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF format to standard output
    // set numv_per_simplex to dim
  {
    ijkoutOFF(dim, dim, coord, numv, simplex_vert, nums);
  }


  /// Output Geomview .off file.
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF files
    // vector input format
  {
    ijkoutOFF(out, dim, numv_per_simplex, vector2pointer(coord), coord.size()/dim,
              vector2pointer(simplex_vert), simplex_vert.size()/numv_per_simplex);
  }

  /// Output Geomview .off file.
  /// Every simplex has \a dim vertices.
  template <class CTYPE, class VTYPE> void ijkoutOFF
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
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF format to standard output
    // vector input format
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file to standard output.
  /// Every simplex has \a dim vertices.
  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
    // output Geomview OFF format to standard output
    // vector input format
    // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  /// Output Geomview .off file. Color vertices.
  template <class T, class colorT> void ijkoutColorVertOFF
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
  template <class T, class colorT> void ijkoutColorVertOFF
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
  template <class CTYPE, class VTYPE, class colorT> void ijkoutColorVertOFF
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
  template <class T, class colorT> void ijkoutColorFacesOFF
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

  /// Output Geomview .off file to standard output. Color simplices.
  template <class T, class colorT> void ijkoutColorFacesOFF
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
  template <class CTYPE, class VTYPE, class colorT> void ijkoutColorFacesOFF
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

  /// Output Geomview .off file. Output vertex normals.
  template <class CTYPE, class NTYPE, class VTYPE> void ijkoutNormalsOFF
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
  template <class CTYPE, class NTYPE, class VTYPE> void ijkoutNormalsOFF
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
  template <class CTYPE, class VTYPE> void ijkoutQuadOFF
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
  template <class CTYPE, class VTYPE> void ijkoutQuadOFF
  (const int dim,
   const CTYPE * coord, const int numv,
   const VTYPE * quad_vert, const int numq,
   const bool flag_reorder_vertices)
  {
    ijkoutQuadOFF(std::cout, dim, coord, numv, quad_vert, numq);
  }

  /// Output quadrilaterals to Geomview .off
  /// Vector input format
  template <class CTYPE, class VTYPE> void ijkoutQuadOFF
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
  /// Vector input format
  template <class CTYPE, class VTYPE> void ijkoutQuadOFF
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
  template <class CTYPE, class VTYPE> void ijkoutLINE
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
  template <class CTYPE, class VTYPE> void ijkoutLINE
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

  /// Output Geomview .line file with color
  template <class CTYPE, class VTYPE, class COLOR_TYPE> 
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

    if (dim == 3) { out << "LINEC" << std::endl; }
    else if (dim == 4) { out << "4LINEC" << std::endl;}
    else {
      error.AddMessage("Only dimensions 3 and 4 permitted for LINE output.");
      throw error;
    };

    out << rgba[0] << " " << rgba[1] << " " 
        << rgba[2] << " " << rgba[3] << std::endl;
    out << numv << " " << nume << std::endl;

    ijkoutVertexCoord(out, dim, coord, numv);
    out << std::endl;

    for (int ie = 0; ie < nume; ie++) {
      out << edge_vert[ie*NUM_EDGE_ENDPOINTS]
          << " " << edge_vert[ie*NUM_EDGE_ENDPOINTS+1] << std::endl;
    };

  }

  /// Output Geomview .line file.
  /// Vector input format
  template <class CTYPE, class VTYPE, class COLOR_TYPE> 
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
  /// @param in = Input stream.
  /// @param dim = Vertex dimension.
  /// @param numv = Number of vertices.
  /// @param nums = Number of simplices.
  /// @param nume = Number of edges.
  inline void ijkinOFFheader
  (std::istream & in, int & dim, int & numv, int & nums, int & nume)
  {
    std::string header_keyword;
    IJK::PROCEDURE_ERROR error("ijkinOFF");

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

  /// \brief Read coordinates from Geomview .off file.
  ///
  /// Ignores any normal, color information.
  /// @param in = Input stream.
  /// @param dim = Dimension of vertices.
  /// @param numv = Number of vertices.
  /// @param coord = Array of coordinates. 
  ///        coord[dim*i+k] = k'th coordinate of vertex i (k < dim).
  /// @pre Array coord[] is preallocated with size at least \a numv * \a dim.
  template <class T> void ijkinOFFcoord
  (std::istream & in, const int dim, const int numv, T * & coord)
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

  /// \brief Read simplex vertices of simplex \a js from Geomview .off file.
  ///
  /// Ignores any color information.
  /// @param in = Input stream.
  /// @param js = Simplex index.
  /// @param numv_per_simplex = Number of vertices per simplex.
  /// @param simplex_vert = Array of simplex vertices. 
  ///        simplex_vert[numv_per_simplex*js+k] = 
  ///            k'th vertex index of simplex js.
  /// @pre Array simplex_vert[] has been preallocated
  ///      with size at least \a numv_per_simplex * \a (js+1).
  template <class T> void ijkinOFFsimplexVert
  (std::istream & in, const int js, const int numv_per_simplex, 
   T * & simplex_vert)
  {
    for (int k = 0; k < numv_per_simplex; k++) 
      { in >> simplex_vert[js*numv_per_simplex + k]; }
    gobble_line(in);
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
  template <class T> void ijkinOFFsimplex
  (std::istream & in, const int ifirst, const int nums, 
   const int numv_per_simplex, T * & simplex_vert)
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

      ijkinOFFsimplexVert(in, i, numv_per_simplex, simplex_vert);

      if (!in.good()) {
        error.AddMessage("Error reading vertices of simplex ", i, ".");
        throw error;
      }
    }
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
  template <class T> void ijkinOFF
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
      ijkinOFFsimplexVert(in, 0, num_simplex_vert, simplex_vert);

      ijkinOFFsimplex(in, 1, nums-1, num_simplex_vert, simplex_vert);
    };
  }

  /// Read Geomview .off file from standard input.
  template <class T> void ijkinOFF
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
  template <class T> void ijkinOFF
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
  template <class T> void ijkinOFF
  (int & dim, T * & coord, int & numv, int * & simplex_vert, int & nums)
  {
    ijkinOFF(std::cin, dim, coord, numv, simplex_vert, nums);
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
  template <class CTYPE, class VTYPE, class SCALE_TYPE> void ijkoutFIG
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
  template <class CTYPE, class VTYPE, class SCALE_TYPE> void ijkoutFIGpolyline
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
    template <class CTYPE> void ijkoutVertexCoord
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
    template <class VTYPE> void ijkoutPolygonVertices
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

    /// Output quad vertices.
    /// @param out Output stream.
    /// @param quad_vert[] Array of quad vertices.
    ///        quad_vert[4*j+k] = k'th vertex index of polygon j.
    /// @param numq Number of quadrilaterals.
    ///  @param flag_reorder_vertices = if true, change vertex order to be
    ///                         counter-clockwise around quad.
    template <class VTYPE> void ijkoutQuadVertices
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

  }

}

#endif

/// \file ijktable_poly.h
/// Class containing isosurface table polyhedron.
/// Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012 Rephael Wenger

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

#ifndef _IJKTABLE_POLY_
#define _IJKTABLE_POLY_

#include "ijk.txx"

/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKTABLE {

  typedef unsigned char 
  ISOSURFACE_VERTEX_INDEX;  ///< Index of isosurface vertex.
  typedef unsigned char EDGE_INDEX;    ///< Index of edge.
  typedef unsigned char FACET_INDEX;   ///< Index of facet.
  typedef int FACET;          ///< Bits representing vertices in facet.
  typedef int FACET_SET;      ///< Bits representing set of facets.

  const int NO_VERTEX = -1;

  // **************************************************
  // ISOSURFACE TABLE POLYHEDRON
  // **************************************************

  /// Isosurface table polyhedron.
  class ISOSURFACE_TABLE_POLYHEDRON {

  protected:
    int dimension;         ///< Polyhedron dimension.
    int num_vertices;      ///< Number of polyhedron vertices.
    int num_edges;         ///< Number of polyhedron edges.
    int num_facets;        ///< Number of polyhedron facets.
    int * vertex_coord;    ///< Polyhedron vertex coordinates.
    int * edge_endpoint;   ///< Polyhedron edge endpoints.
    int * num_facet_vertices; ///< Number of vertices of each facet.
    int ** facet_vertex_list; ///< List of vertices in each facet.
    FACET * facet;         ///< Polyhedron facets.
    void Init();           ///< Initialize.
    void FreeFacets();     ///< Free all facet arrays.

  public:
    ISOSURFACE_TABLE_POLYHEDRON(const int d);  ///< Constructor
    ~ISOSURFACE_TABLE_POLYHEDRON();            ///< Destructor
    ISOSURFACE_TABLE_POLYHEDRON
    (const ISOSURFACE_TABLE_POLYHEDRON & init);  ///< Copy constructor.
    const ISOSURFACE_TABLE_POLYHEDRON & operator = 
    (const ISOSURFACE_TABLE_POLYHEDRON &);  ///< Assignment.

    /// @name Get Functions
    //@{
    int Dimension() const   ///< Polyhedron dimension.
    { return(dimension); };
    int NumVertices() const ///< Number of polyhedron vertices.
    { return(num_vertices); };
    int NumEdges() const    ///< Number of polyhedron edges.
    { return(num_edges); };
    int NumFacets() const   ///< Number of polyhedron facets.
    { return(num_facets); };
    int NumFacetVertices    ///< Number of facet vertices of facet \a jf.
    (const FACET_INDEX jf) const
    { return(num_facet_vertices[jf]); };
    int VertexCoord         ///< \a ic'th vertex coordinate of vertex \a iv.
    (const int iv, const int ic) const
    { return(vertex_coord[iv*dimension + ic]); };
    int EdgeEndpoint        ///< \a j'th endpoint of edge \a ie. \a j = 0 or 1.
    (const EDGE_INDEX ie, const int j) const
    { return(edge_endpoint[int(ie)*2 + j]); };
    int MidpointCoord       ///< \a ic'th coordinate of midpoint of edge \a ie.
    (const EDGE_INDEX ie, const int ic) const;
    FACET Facet             ///< Bits representing vertices in facet \a jf.
    (const FACET_INDEX jf) const
    { return(facet[jf]); };
    bool IsVertexInFacet    ///< Return true if vertex \a iv is in facet \a jf.
    (const FACET_INDEX jf, const int iv) const
    { return((facet[jf] & ((1L) << iv)) != 0); };
    int FacetVertex         ///< Return \a k'th vertex in facet \a jf.
    (const FACET_INDEX jf, const int k) const
    { return(facet_vertex_list[jf][k]); };
    //@}

    /// @name Set Functions
    //@{
    void SetDimension(const int d);      ///< Set polyhedron dimension.
    void SetNumVertices(const int numv); ///< Set number of polyhedron vertices.
    void SetNumEdges(const int nume);    ///< Set number of polyhedron edges.
    void SetNumFacets(const int numf);   ///< Set number of polyhedron facets.

    /// Set number of polyhedron vertices, edges and facets.
    void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };

    /// Set \a ic'th coordinate of vertex \a iv.
    /// @pre SetNumVertices or SetSize must be called before SetVertexCoord.
    void SetVertexCoord 
    (const int iv, const int ic, const int coord);

    /// Set endpoints of edge \a ie.
    /// @pre SetNumEdges or SetSize must be called before SetEdge.
    void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1);

    /// Set number of vertices in facet \a jf.
    /// @pre SetNumFacets or SetSize must be called before SetNumFacetVertices.
    void SetNumFacetVertices 
    (const FACET_INDEX jf, const int numv);

    /// Set \a k'th facet vertex of facet \a jf to vertex \a iv.
    /// @pre SetNumFacetVertices(jf) must be called before SetFacetVertex.
    void SetFacetVertex(const FACET_INDEX jf, const int k, const int iv);
    //@}

    /// @name Memory Management Functions
    //@{
    void FreeAll();             ///< Free all memory.
    //@}

    /// @name Check Functions
    //@{
    bool CheckDimension() const;
    bool Check(IJK::ERROR & error_msg) const;
    //@}

    /// @name Generate Polyhedron
    //@{

    /// Generate a square, cube or hypercube.
    void GenCube(const int cube_dimension); 

    /// Generate a square, cube or hypercube.
    ///   Use order of facets and edges given by template class CUBE_FACE_INFO.
    void GenCubeOrderA(const int cube_dimension);

    /// Generate a triangle, tetrahedron or simplex.      
    void GenSimplex(const int simplex_dimension);

    /// Generate a pyramid over a square, cube or hypercube base.
    void GenPyramid(const int pyramid_dimension);

    //@}

  };

  typedef ISOSURFACE_TABLE_POLYHEDRON * ISOSURFACE_TABLE_POLYHEDRON_PTR;


  // **************************************************
  // ROUTINES FOR GENERATING POLYHEDRA
  // **************************************************

  void generate_prism(const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
                      ISOSURFACE_TABLE_POLYHEDRON & prism);
}

#endif

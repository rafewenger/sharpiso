/// \file ijkcube.txx
/// ijk templates defining cube (hypercube) classes and functions.
/// Version 0.1.0

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

#ifndef _IJKCUBE_
#define _IJKCUBE_

#include <algorithm>
#include <limits>

#include "ijk.txx"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS CUBE_INFO
  // **************************************************

  /// Base cube class
  template <typename DTYPE, typename NTYPE> 
  class CUBE_INFO {

  protected:
    DTYPE dimension;           ///< Cube dimension
    NTYPE num_vertices;        ///< Number of cube vertices
    NTYPE num_edges;           ///< Number of cube edges.
    NTYPE num_facets;          ///< Number of cube facets.
    NTYPE num_facet_vertices;  ///< Number of cube facet vertices.
    NTYPE num_ridge_vertices;  ///< Number of cube ridge vertices.

    void Init                  /// Initialize cube.
    (const DTYPE dimension);
    void FreeAll();            ///< Free all allocated memory.


  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    // Constructors, destructor.
    CUBE_INFO(const DTYPE dimension);     ///< Constructor.
    CUBE_INFO();                          ///< Constructor.
    ~CUBE_INFO();                         ///< Destructor.

    // set functions
    template <typename DTYPE2>
    void SetDimension                     ///< Set cube dimension
    (const DTYPE2 dimension);

    // get functions
    DTYPE Dimension() const               ///< Dimension.
    { return(dimension); }
    NTYPE NumVertices() const             ///< Number of cube vertices. 
    { return(num_vertices); }
    NTYPE NumEdges() const                ///< Number of cube edges.
    { return(num_edges); }
    NTYPE NumFacets() const               ///< Number of cube facets. 
    { return(num_facets); }
    NTYPE NumFacetVertices() const        ///< Number of cube facet vertices. 
    { return(num_facet_vertices); }
    NTYPE NumRidgeVertices() const        ///< Number of cube ridge vertices. 
    { return(num_ridge_vertices); }
    NTYPE NumDiagonals() const            ///< Number of cube diagonals.
    { return(NumFacetVertices()); }

    /// Maximum cube vertex index. Undefined if dimension < 1.
    NTYPE MaxVertexIndex() const
    { return(NumVertices()-1); }

    NTYPE OppositeVertex                  ///< Index of vertex opposite iv
    (const NTYPE iv) const
    { return(MaxVertexIndex()-iv); }

    /// Return neighbor of vertex \a iv0.
    /// @param dir Edge direction.
    /// @param iv0 Vertex.
    NTYPE VertexNeighbor
    (const NTYPE iv0, const DTYPE dir)
    {
      int mask = (1L << dir);
      return((iv0^mask));
    };

  };

  // **************************************************
  // TEMPLATE CLASS CUBE_FACE_INFO
  // **************************************************

  /// Information about cube faces
  template <typename DTYPE, typename NTYPE, typename VTYPE>
  class CUBE_FACE_INFO:public CUBE_INFO<DTYPE,NTYPE> {

  protected:
    VTYPE * facet_vertex;
    VTYPE * edge_endpoint;

    void Init();
    void Init(const int dimension);
    void FreeAll();

    /// Set facet vertices and edge endpoints.
    void SetFaces();

    /// Set facet vertices.
    void SetFacetVertices();

    /// Set edge endoints.
    /// @pre Array facet_vertex[] should be set before edge_endpoint[].
    void SetEdgeEndpoints();

  public:
    typedef VTYPE VERTEX_INDEX;           ///< Vertex index type.

  public:
    CUBE_FACE_INFO();
    CUBE_FACE_INFO(const DTYPE dimension);
    ~CUBE_FACE_INFO() { FreeAll(); };

    // set functions
    template <typename DTYPE2>
    void SetDimension                     ///< Set cube dimension
    (const DTYPE2 dimension);

    /// Return pointer to edge endpoints
    const VTYPE * EdgeEndpoint() const
    { return(edge_endpoint); }

    /// Return edge endpoint.
    /// @param ie Edge index.
    /// @param iend Endpoint (0 or 1).
    VTYPE EdgeEndpoint(const VTYPE ie, const NTYPE iend) const
    { return(edge_endpoint[2*ie+iend]); };

    /// Return pointer to facet vertices
    const VTYPE * FacetVertex() const
    { return(facet_vertex); }

    /// Return facet vertex.
    /// @param ifacet Facet index. 
    ///    Facet is orthogonal to axis (ifacet%dimension).
    ///    If ifacet < dimension, facet contains vertex 0.
    //     If ifacet >= dimension, facet contains largest vertex.
    VTYPE FacetVertex(const NTYPE ifacet, const NTYPE k) const
    { return(facet_vertex[ifacet*this->NumFacetVertices()+k]); };

    /// Return facet index.
    VTYPE FacetIndex(const VTYPE facet_vertex_index) const
    { return(int(facet_vertex_index)/int(this->NumFacetVertices())); }
  };

  // **************************************************
  // TEMPLATE CLASS CUBE
  // **************************************************

  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  class CUBE:public CUBE_INFO<DTYPE,NTYPE> {

  protected:
    /// vertex_coord[dimension*k+j] = 
    ///   j'th coordinate of k'th vertex of cube
    CTYPE * vertex_coord;

    /// max_vertex_index.  Stored for faster processing.
    NTYPE max_vertex_index;

    /// cube edge length
    LTYPE edge_length;

    void ZeroLocal();
    void InitLocal();
    template <typename CTYPE2, typename LTYPE2>
    void InitLocal
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);    
    void FreeLocal();
    void AllocateLocal();           ///< Allocate data structures in CUBE.
    void SetToUnitCube();           ///< Set coordinates to unit cube.
    void SetMaxVertexIndex();       ///< Set max vertex index.

  public:
    typedef CTYPE COORD_TYPE;       ///< Coordinate type.

  public:
    // Constructors and destructors.
    template <typename CTYPE2>
    CUBE(const DTYPE dimension, const CTYPE2 * vertex0_coord);
    template <typename CTYPE2, typename LTYPE2>
    CUBE(const DTYPE dimension, 
         const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
    CUBE(const DTYPE dimension);
    CUBE();
    ~CUBE();

    // set functions
    template <typename DTYPE2>
    void SetDimension               ///< Set cube dimension.
    (const DTYPE2 dimension);
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

    // *** get functions ***

    NTYPE MaxVertexIndex() const    ///< Maximum cube vertex index.
    { return(max_vertex_index); }

    NTYPE OppositeVertex            ///< Index of vertex opposite iv
    (const NTYPE iv) const
    { return(MaxVertexIndex()-iv); }

    LTYPE EdgeLength() const
    { return(edge_length); }        ///< Cube edge length.

    /// Return pointer to vertex coordinates
    const CTYPE * VertexCoord() const
    { return(vertex_coord); }

    /// Return pointer to coordinates of k'th cube vertex
    const CTYPE * VertexCoord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return j'th coordinate of k'th vertex
    const CTYPE VertexCoord         
    (const NTYPE k, const NTYPE j) const
    { return(vertex_coord[this->Dimension()*k+j]); }

    /// Return true if cube contains point.
    /// @param point_coord[] Point coordinates.
    /// @pre point_coord[] contains Dimension() point coordinates.
    template <typename CTYPE2>
    bool Contains(const CTYPE2 * point_coord) const;

    // *** get diagonal functions ***

    /// Return pointer to coordinates of endpoint 0 of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalEnd0Coord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return pointer to coordinates of endpoint 1 of diagonal k.
    /// Corner opposite cube corner k is endpoint 1 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalEnd1Coord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*OppositeVertex(k)); }

    /// Return d'th coordinate of endpoint 0 of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalEnd0Coord
    (const NTYPE k, const DTYPE d) const
    { return(vertex_coord[this->Dimension()*k+d]); }

    /// Return d'th coordinate of endpoint 1 of diagonal k.
    /// Cube corner k is endpoint 1 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalEnd1Coord
    (const NTYPE k, const DTYPE d) const
    { return(vertex_coord[this->Dimension()*OppositeVertex(k)+d]); }

    /// Return pointer to coordinates of endpoint iend of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalCoord
    (const NTYPE k, const NTYPE iend) const
    { if (iend == 0) { return(DiagonalEnd0Coord(k)); }
      else { return(DiagonalEnd1Coord(k)); }
    }

    /// Return d'th coordinate of endpoint iend of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalCoord
    (const NTYPE k, const NTYPE iend, const DTYPE d) const
    { if (iend == 0) { return(DiagonalEnd0Coord(k, d)); }
      else { return(DiagonalEnd1Coord(k, d)); }
    }

  };

  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class UNIT_CUBE:public CUBE<DTYPE,NTYPE,CTYPE,CTYPE> {

  public:
    // Constructors and destructors.
    UNIT_CUBE(const DTYPE dimension);
    UNIT_CUBE();
    ~UNIT_CUBE() {};

    // get functions
    CTYPE EdgeLength() const
    { return(1); }                  ///< Unit cube edge length.

    // Undefine set functions which do not apply to UNIT_CUBE.

    /// Undefine SetVertexCoord().
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

    /// Undefine Set().
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

  };

  // **************************************************
  // TEMPLATE CLASS CUBE3D
  // **************************************************

  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  class CUBE3D:public CUBE<DTYPE,NTYPE,CTYPE,LTYPE> {

  public:
    CUBE3D():CUBE<DTYPE,NTYPE,CTYPE,LTYPE>(3) {};
    template <typename CTYPE2>
    CUBE3D(const CTYPE2 * vertex0_coord):
      CUBE<DTYPE,NTYPE,CTYPE,LTYPE> (3, vertex0_coord) {};
    template <typename CTYPE2, typename LTYPE2>
    CUBE3D(const CTYPE2 * vertex0_coord, const LTYPE2 edge_length):
      CUBE<DTYPE,NTYPE,CTYPE,LTYPE> (vertex0_coord, edge_length) {};

    // Undefine set functions which do not apply to CUBE3D.

    /// Undefine SetDimension().
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 dimension);

    /// Undefine Set().
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

  };

  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE3D
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class UNIT_CUBE3D:public CUBE3D<DTYPE,NTYPE,CTYPE,CTYPE> {

  public:
    // Constructors and destructors.
    UNIT_CUBE3D() {};

    // Undefine set functions which do not apply to UNIT_CUBE3D.

    /// Undefine SetVertexCoord().
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
  };

  // **************************************************
  // TEMPLATE FUNCTIONS: COUNTING
  // **************************************************

  /// Return number of cube vertices
  template <typename DTYPE> 
  long compute_num_cube_vertices(const DTYPE dimension)
  { return(1L << dimension); }

  /// Return number of cube facet vertices
  /// @pre dimension > 0
  template <typename DTYPE> 
  long compute_num_cube_facet_vertices(const DTYPE dimension)
  { return(1L << (dimension-1)); }

  /// Return number of cube ridge vertices
  /// @pre dimension > 0
  template <typename DTYPE> 
  long compute_num_cube_ridge_vertices(const DTYPE dimension)
  { return(1L << (dimension-2)); }

  /// Return number of cube facets
  template <typename DTYPE> 
  long compute_num_cube_facets(const DTYPE dimension)
  { return(2*dimension); }

  /// Return number of cube edges
  template <typename DTYPE> 
  long compute_num_cube_edges(const DTYPE dimension)
  { 
    long numv = compute_num_cube_vertices(dimension);
    return((numv*dimension)/2);
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTE COORDINATES
  // **************************************************

  /// Compute cube vertex coordinates.
  /// @param dimension  Dimension of grid.
  /// @param vertex0_coord Coordinates of lowest/leftmost cube vertex.
  /// @param cube_edge_length Cube edge length.
  /// @param[out] coord[] Cube vertex coordinates.
  /// @pre Array coord[] is allocated with size at least 
  ///      (number of cube vertices)*dimension
  template <typename DTYPE, typename CTYPE0, typename LTYPE0, typename CTYPE1>
  void compute_cube_vertex_coord
  (const DTYPE dimension, const CTYPE0 * vertex0_coord,
   const LTYPE0 cube_edge_length, CTYPE1 * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_vertex_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(vertex0_coord, "vertex0_coord", error)) 
      { throw error; }
    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        int iend = j0 % 2;
        coord[j*dimension+d] = vertex0_coord[d] + iend*cube_edge_length;
        j0 = j0/2;
      }
    }
  }

  /// Compute unit cube vertex coordinates (0,0,...,0) to (1,1,...,1)
  /// @param dimension  Dimension of grid.
  /// @param[out] coord[] = Unit cube vertex coordinates.
  /// @pre Array coord[] is allocated with size at least 
  ///      (number of cube vertices)*dimension
  template <typename DTYPE, typename CTYPE>
  void compute_unit_cube_vertex_coord(const DTYPE dimension, CTYPE * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_unit_cube_vertex_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        coord[j*dimension+d] = j0 % 2;
        j0 = j0/2;
      }
    }
  }

  /// Compute coordinates of cube diagonal endpoints.
  /// @param dimension  Dimension of grid.
  /// @param vertex0_coord Coordinates of lowest/leftmost cube vertex.
  /// @param cube_edge_length Cube edge length.
  /// @param[out] coord[] = Unit cube diagonal coordinates.
  ///   coord[(2*k+i)*dimension+j] = 
  ///     j'th coordinate of endpoint i of diagonal k.
  /// @pre Array coord[] is allocated with size at least 
  ///      2*(number of cube vertices)*dimension
  template <typename DTYPE, typename CTYPE0, 
            typename LTYPE0, typename CTYPE1>
  void compute_cube_diagonal_coord
  (const DTYPE dimension, const CTYPE0 * vertex0_coord,
   const LTYPE0 cube_edge_length, CTYPE1 * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_diagonal_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(vertex0_coord, "vertex0_coord", error)) 
      { throw error; }
    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long k = 0; k < num_cube_vertices; k++) {
      long k0 = k;
      for (DTYPE d = 0; d < dimension; d++) {
        int iend = k0 % 2;
        // Diagonal k is pointing away from cube vertex k.
        coord[(2*k)*dimension+d] = vertex0_coord[d] + iend*cube_edge_length;
        // Diagonal k is pointing toward the opposite vertex.
        coord[(2*k+1)*dimension+d] = 
          vertex0_coord[d] + (1-iend)*cube_edge_length;
        k0 = k0/2;
      }
    }

  }

  /// Compute coordinates of unit cube diagonal endpoints.
  /// @param dimension  Dimension of grid.
  /// @param[out] coord[] = Cube diagonal coordinates.
  ///   coord[(2*k+i)*dimension+j] = 
  ///     j'th coordinate of endpoint i of diagonal k.
  /// @pre Array coord[] is allocated with size at least 
  ///      2*(number of cube vertices)*dimension
  template <typename DTYPE, typename CTYPE>
  void compute_unit_cube_diagonal_coord
  (const DTYPE dimension, CTYPE * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_unit_cube_diagonal_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long k = 0; k < num_cube_vertices; k++) {
      long k0 = k;
      for (DTYPE d = 0; d < dimension; d++) {
        CTYPE c = k0 % 2;
        // Diagonal k is pointing away from cube vertex k.
        coord[(2*k)*dimension+d] = c;
        // Diagonal k is pointing toward the opposite vertex.
        coord[(2*k+1)*dimension+d] = 1-c;
        k0 = k0/2;
      }
    }

  }

  // **************************************************
  // TEMPLATE FUNCTION: VERTEX NEIGHBOR
  // **************************************************

  /// Compute the neighbor of cube vertex iv in given direction.
  template <typename VTYPE1, typename DTYPE, typename VTYPE2>
  void compute_cube_vertex_neighbor
  (const VTYPE1 iv1, const DTYPE d, VTYPE2 & iv2)
    {
      VTYPE1 mask = (VTYPE1(1) << d);
      iv2 = (iv1^mask);
    }

  template <typename NTYPE, typename FTYPE1, typename FTYPE2>
  void compute_opposite_cube_facet
  (const NTYPE num_facets, const FTYPE1 kf, FTYPE2 & kf2)
  { 
    kf2 = (kf+(num_facets/2))%num_facets;
  }

  // **************************************************
  // TEMPLATE CUBE FACET FUNCTIONS
  // **************************************************

  /// Return orthogonal direction to facet.
  template <typename DTYPE, typename FTYPE>
  inline FTYPE cube_facet_orth_dir(const DTYPE dimension, const FTYPE kf)
  {
    return(kf%dimension);
  }

  /// Return side containing facet.
  /// Return 0 if facet contains vertex 0.
  /// Return 1 if facet does not contain vertex 0.
  template <typename DTYPE, typename FTYPE>
  inline FTYPE cube_facet_side(const DTYPE dimension, const FTYPE kf)
  {
    return(kf/dimension);
  }
  
  /// Return true if facet jf contains vertex i.
  template <typename DTYPE, typename FTYPE, typename VTYPE>
  bool cube_facet_contains
  (const DTYPE dimension, const FTYPE jf, const VTYPE iv)
  {
    FTYPE side = jf/dimension;
    FTYPE orth_dir = jf%dimension;
    FTYPE mask = (FTYPE(1) << orth_dir);

    bool flag_contains = ((side << orth_dir) == (mask & iv));
    return(flag_contains);
  }

  // **************************************************
  // TEMPLATE BOUNDARY BIT FUNCTIONS
  // **************************************************


  /// Return index of boundary bit.
  template <typename DTYPE, typename STYPE, typename ITYPE>
  inline void compute_boundary_bit_index
  (const DTYPE orth_dir, const STYPE side, ITYPE & bit_index)
  {
    bit_index = 2*orth_dir + DTYPE(side);
  }


  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTE CORNER IN DIRECTION
  // **************************************************

  /// Compute cube corner nearest to given direction.
  template <typename DTYPE, typename CTYPE, typename NTYPE>
  void compute_corner_nearest_direction
  (const DTYPE dimension, const CTYPE * dir, NTYPE & icorner)
  {
    icorner = 0;
    NTYPE bit_flag = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (dir[d] > 0) 
        { icorner = (icorner | bit_flag); }
      bit_flag = (bit_flag << 1);
    }
  }

  // **************************************************
  // TEMPLATES TO CHECK VALUES
  // **************************************************

  /// Check dimension
  /// Return true if dimension is non-negative
  template <typename DTYPE>
  bool check_cube_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension < 0) {
      error.AddMessage("Illegal dimension ", dimension, ".");
      error.AddMessage("Dimension must be non-negative.");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS CUBE_INFO MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::CUBE_INFO(const DTYPE dimension)
  {
    Init(dimension);
  }

  /// Constructor
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::CUBE_INFO()
  {
    Init(0);
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::~CUBE_INFO()
  {
    dimension = 0;
    num_vertices = 0;
    num_edges = 0;
    num_facets = 0;
    num_facet_vertices = 0;
    num_ridge_vertices = 0;
  }

  /// Initialize
  template <typename DTYPE, typename NTYPE> 
  void CUBE_INFO<DTYPE,NTYPE>::Init(const DTYPE dimension)
  {
    SetDimension(dimension);
  }

  /// Set cube dimension
  template <typename DTYPE, typename NTYPE> 
  template <typename DTYPE2>
  void CUBE_INFO<DTYPE,NTYPE>::SetDimension(const DTYPE2 dimension)
  {
    this->dimension = dimension;
    this->num_vertices = compute_num_cube_vertices(dimension);
    this->num_facets = compute_num_cube_facets(dimension);
    this->num_edges = compute_num_cube_edges(dimension);
    this->num_facet_vertices = num_vertices/2;
    this->num_ridge_vertices = num_facet_vertices/2;
  }

  // **************************************************
  // TEMPLATE CLASS CUBE_FACE_INFO MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  CUBE_FACE_INFO(const DTYPE dimension):CUBE_INFO<DTYPE,NTYPE>(dimension)
  {
    Init(dimension);
  }

  /// Constructor
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::CUBE_FACE_INFO()
  {
    Init();
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::Init()
  {
    facet_vertex = NULL;
    edge_endpoint = NULL;
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::Init(const int dimension)
  {
    facet_vertex = NULL;
    edge_endpoint = NULL;

    SetDimension(dimension);
  }

  // Free all arrays
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::FreeAll()
  {
    if (facet_vertex != NULL) {
      delete [] facet_vertex;
      facet_vertex = NULL;
    }

    if (edge_endpoint != NULL) {
      delete [] edge_endpoint;
      edge_endpoint = NULL;
    }
  }


  // Set dimension.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  template <typename DTYPE2>
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    FreeAll();
    CUBE_INFO<DTYPE,NTYPE>::SetDimension(dimension);
    facet_vertex = new VTYPE[this->NumFacetVertices()*this->NumFacets()];
    edge_endpoint = new VTYPE[this->NumEdges()*2];

    SetFacetVertices();
    SetEdgeEndpoints();
  }

  // Set facet vertices.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetFacetVertices()
  {
    const DTYPE dimension = this->Dimension();
    const NTYPE num_cube_vertices = this->NumVertices();
    const NTYPE num_facet_vertices = this->NumFacetVertices();
    IJK::PROCEDURE_ERROR error("SetFacetVertices");

    if (dimension < 1) { return; };

    if (facet_vertex == NULL) {
      error.AddMessage("Programming error.  Array facet_vertex[] not allocated.");
      throw error;
    }

    for (NTYPE ifacet = 0; ifacet < this->NumFacets(); ifacet++) {

      DTYPE orth_dir = ifacet%dimension;

      // Multiply by s mod (num_cube_vertices-1) to compute vertex index.
      NTYPE s = 2;
      for (DTYPE d = 0; d < orth_dir; d++)
        { s = s*2; };
      s = s%(num_cube_vertices-1);

      NTYPE i0 = ifacet*num_facet_vertices;

      if (ifacet < dimension) {

        for (NTYPE i = 0; i < num_facet_vertices; i++) 
          { facet_vertex[i0+i] = (i*s)%(num_cube_vertices-1); }
      }
      else {

        // Translate by t mod num_cube_vertices to compute vertex index.
        NTYPE t = 1;
        for (DTYPE d = 0; d < orth_dir; d++)
          { t = t*2; };

        for (NTYPE i = 0; i < num_facet_vertices; i++) 
          { facet_vertex[i0+i] = ((i*s)%(num_cube_vertices-1))+t; }

        // Swap subfacets to get consistent facet orientation.
        if (num_facet_vertices > 1) {
          NTYPE num_subfacet_vertices = num_facet_vertices/2;
          for (NTYPE i = 0; i < num_subfacet_vertices; i++) {
            std::swap(facet_vertex[i0+i], 
                      facet_vertex[i0+i+num_subfacet_vertices]);
          }
        }
      }
    }

  }

  // Set edge endpoints.
  // @pre Array facet_vertex[] should be set before edge_endpoints[].
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetEdgeEndpoints()
  {
    const DTYPE dimension = this->Dimension();
    IJK::PROCEDURE_ERROR error("SetEdgeEndpoints");

    if (dimension < 1) { return; };

    if (edge_endpoint == NULL) {
      error.AddMessage("Programming error.  Array edge_endpoint[] not allocated.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      for (NTYPE i = 0; i < this->NumFacetVertices(); i++) {
        NTYPE j = d*this->NumFacetVertices()+i;
        VTYPE iv0 = FacetVertex(d, i);
        edge_endpoint[2*j] = iv0;
        VTYPE iv1 = this->VertexNeighbor(iv0, d);
        edge_endpoint[2*j+1] = iv1;
      }
    }

  }

  // **************************************************
  // TEMPLATE CLASS CUBE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::CUBE(const DTYPE dimension):
    CUBE_INFO<DTYPE,NTYPE>(dimension)
  {
    InitLocal();
  }

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::CUBE():CUBE_INFO<DTYPE,NTYPE>()
  {
    InitLocal();
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::~CUBE()
  {
    FreeLocal();
  }

  /// Initialize.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::InitLocal()
  {
    ZeroLocal();
    AllocateLocal();
    SetMaxVertexIndex();
    SetToUnitCube();
  }

  /// Allocate local data structures in CUBE
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::AllocateLocal()
  {
    IJK::PROCEDURE_ERROR error("CUBE::AllocateLocal");

    if (!check_cube_dimension(this->Dimension(), error)) 
      { throw error; }
    if (!IJK::check_is_NULL(vertex_coord, "vertex_coord", error))
      { throw error; }

    NTYPE numv = this->NumVertices();
    this->vertex_coord = new CTYPE[numv*this->Dimension()];
  }

  // Return true if cube contains coordinate.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2>
  bool CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::
  Contains(const CTYPE2 * point_coord) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      COORD_TYPE c = VertexCoord(0,d);
      if (point_coord[d] < c)
        { return(false); }

      if (point_coord[d] > c + EdgeLength())
        { return(false); }
    }

    return(true);
  }

  // Set max vertex index.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::SetMaxVertexIndex()
  {
    if (this->NumVertices() > 0) 
      { max_vertex_index = this->NumVertices()-1; }
    else
      { max_vertex_index = 0; }
  }

  // Set vertex coordinates to unit cube vertex coordinates.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::SetToUnitCube()
  {
    IJK::PROCEDURE_ERROR error("CUBE::SetToUnitCube");

    if (!check_cube_dimension(this->Dimension(), error)) 
      { throw error; }
    if (!check_array_allocated(vertex_coord, "vertex_coord", error)) 
      { throw error; }

    edge_length = 1;
    compute_unit_cube_vertex_coord
      (this->Dimension(), this->vertex_coord);
  }

  /// Set pointers defined in CUBE to NULL.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::ZeroLocal()
  {
    max_vertex_index = 0;
    this->vertex_coord = NULL;
  }


  /// Free memory allocated in CUBE.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::FreeLocal()
  {
    if (vertex_coord != NULL) { delete [] vertex_coord; };
    ZeroLocal();
  }

  /// Set cube dimension.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename DTYPE2>
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    FreeLocal();

    CUBE_INFO<DTYPE,NTYPE>::SetDimension(dimension);
    AllocateLocal();
    SetMaxVertexIndex();
    SetToUnitCube();
  }

  /// Set cube vertex coordinates.
  /// @param vertex0_coord Coordinates of lowest/leftmost vertex.
  /// @param edge_length Length of cube edge.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2, typename LTYPE2>
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::
  SetVertexCoord             ///< Set cube vertex coordinates.
  (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length)
  {
    this->edge_length = edge_length;
    compute_cube_vertex_coord
      (this->Dimension(), vertex0_coord, edge_length, vertex_coord);
  }

  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  UNIT_CUBE<DTYPE,NTYPE,CTYPE>::UNIT_CUBE(const DTYPE dimension):
    CUBE<DTYPE,NTYPE,CTYPE,CTYPE>(dimension)
  {}

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  UNIT_CUBE<DTYPE,NTYPE,CTYPE>::UNIT_CUBE():
    CUBE<DTYPE,NTYPE,CTYPE,CTYPE>()
  {}

}

#endif

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
    { return(num_vertices); }
    NTYPE NumFacetVertices() const        ///< Number of cube facet vertices. 
    { return(num_facet_vertices); }
    NTYPE NumDiagonals() const            ///< Number of cube diagonals.
    { return(NumFacetVertices()); }

    /// Maximum cube vertex index. Undefined if dimension < 1.
    NTYPE MaxVertexIndex() const
    { return(NumVertices()-1); }

    NTYPE OppositeVertex                  ///< Index of vertex opposite iv
    (const NTYPE iv) const
    { return(MaxVertexIndex()-iv); }
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

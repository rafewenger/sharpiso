/// \file ijkgrid.txx
/// ijk templates defining regular grid classes and functions.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008-2013 Rephael Wenger

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

#ifndef _IJKGRID_
#define _IJKGRID_

#include <algorithm>
#include <limits>
#include <vector>

#include "ijk.txx"
#include "ijkcube.txx"

namespace IJK {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  /// Default type for grid size
  typedef long GRID_SIZE_TYPE;

  // **************************************************
  // TEMPLATE CLASS GRID
  // **************************************************

  /// Base grid class
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  class GRID {

  protected:
    DTYPE dimension;     ///< grid dimension
    ATYPE * axis_size;   ///< axis_size[i] = # grid points along axis i
    NTYPE num_vertices;  ///< number of grid vertices

    void Init            /// Initialize grid.
    (const DTYPE dimension, const ATYPE * axis_size);
    void FreeAll();      ///< Free all allocated memory.


  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef ATYPE AXIS_SIZE_TYPE;         ///< Axis size type.
    typedef VTYPE VERTEX_INDEX_TYPE;      ///< Vertex index type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    // Constructors, destructors, assignment.
    GRID(const DTYPE dimension, const ATYPE * axis_size);
    GRID();
    ~GRID();
    GRID(const GRID & grid);
    const GRID & operator = (const GRID & right);

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    DTYPE Dimension() const         /// Dimension.
    { return(dimension); }
    const ATYPE * AxisSize() const  /// Array axis size.
    { return(axis_size); }
    ATYPE AxisSize(const DTYPE i) const /// Axis size[i].
    { return(axis_size[i]); }
    NTYPE NumVertices() const       /// Number of grid vertices. 
    { return(num_vertices); }

    // compute functions
    NTYPE ComputeNumCubes() const;
    NTYPE ComputeNumEdges() const;
    NTYPE ComputeNumInteriorCubes() const;
    NTYPE ComputeNumBoundaryCubes() const;
    template <typename DTYPE2>
    NTYPE ComputeNumVerticesInFacet(const DTYPE2 orth_dir) const;
    template <typename DTYPE2, typename WTYPE>
    NTYPE ComputeNumVerticesInFacet
    (const DTYPE2 orth_dir, const WTYPE boundary_width) const;
    template <typename DTYPE2>
    NTYPE ComputeNumCubesInFacet(const DTYPE2 orth_dir) const;
    template <typename GTYPE>
    VTYPE ComputeVertexIndex(const GTYPE * coord) const;
    template <typename GTYPE>
    VTYPE ComputeVertexIndex(const std::vector<GTYPE> coord) const;
    template <typename GTYPE>
    void ComputeCoord(const VTYPE iv, GTYPE * coord) const;
    template <typename GTYPE>
    void ComputeCubeCenterCoord(const VTYPE iv, GTYPE * coord) const;
    template <typename BTYPE>
    void ComputeBoundaryBits(const VTYPE iv, BTYPE & boundary_bits) const;
    template <typename BTYPE>
    void ComputeBoundaryCubeBits
    (const VTYPE icube, BTYPE & boundary_bits) const;
    template <typename PTYPE, typename ATYPE2>
    void ComputeSubsampledAxisSizes
    (const PTYPE subsample_period, ATYPE2 subsampled_axis_size[]) const;

    // compare
    template <typename DTYPE2, typename ATYPE2>
    bool CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const;
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid) const;

    /// Return true if grid contains specified point.
    template <typename CTYPE>
    bool ContainsPoint(const CTYPE * coord) const;
    template <typename CTYPE>
    bool ContainsPoint(const std::vector<CTYPE> & coord) const;

    /// Return true if cube contains specified point.
    template <typename CTYPE, typename VTYPE2>
    bool CubeContainsPoint
    (const VTYPE2 icube, const CTYPE * coord) const;
    template <typename CTYPE, typename VTYPE2>
    bool CubeContainsPoint
    (const VTYPE2 icube, const std::vector<CTYPE> & coord) const;

    /// Return true if grid contains specified region.
    /// @pre Box.MinCoord(d) <= Box.MaxCoord(d) for every d < Box.Dimension().
    template <typename BOX_TYPE>
    bool ContainsRegion(const BOX_TYPE & region) const;

    /// Return true if cube facet is on grid boundary.
    template <typename CTYPE, typename DTYPE2>
    bool IsCubeFacetOnGridBoundary
    (const CTYPE cube_index, const DTYPE2 facet_orth_dir, 
     const bool facet_side) const;

    // check function
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool CheckDimension
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid,
     const char * grid_label1, const char * grid_label2,
     IJK::ERROR & error) const;
    bool Check(const DTYPE dimension, const ATYPE * axis_size,
               IJK::ERROR & error) const;
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid,
               IJK::ERROR & error) const;

    /// Return true if current grid size (grid1) matches grid2 size.
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
               const char * grid_label1, const char * grid_label2,
               IJK::ERROR & error) const;

    template <typename GTYPE>
    bool CheckCoord(const GTYPE * coord, IJK::ERROR & error) const;
    template <typename GTYPE>
    bool CheckCoord(const std::vector<GTYPE> & coord, 
                    IJK::ERROR & error) const;
    template <typename GTYPE>
    bool CheckCubeCoord(const GTYPE * coord, IJK::ERROR & error) const;
    template <typename GTYPE>
    bool CheckCubeCoord(const std::vector<GTYPE> & coord, 
                    IJK::ERROR & error) const;
    template <typename ITYPE>
    bool CheckVertexIndex(const ITYPE vertex_index,
                          IJK::ERROR & error) const;
    template <typename ITYPE>
    bool CheckCubeIndex(const ITYPE cube_index,
                        IJK::ERROR & error) const;

    /// Check that grid contains specified region.
    template <typename VTYPE2, typename ATYPE2>
    bool CheckContainsRegion
    (const VTYPE2 region_v0, const ATYPE2 * region_axis_size,
     IJK::ERROR & error) const;
  };

  // **************************************************
  // TEMPLATE CLASS GRID_PLUS
  // **************************************************

  /// Inherits class GRID and adds other indexes and operators 
  ///   for fast accessing of grid vertices.
  /// @tparam DTYPE  Dimension data type.
  /// @tparam ATYPE  Axis size type.
  /// @tparam VTYPE  Vertex index type.
  /// @tparam NTYPE  Number type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  class GRID_PLUS:public GRID<DTYPE,ATYPE,VTYPE,NTYPE> {

  protected:

    /// iv+axis_increment[d] is vertex next to iv
    VTYPE * axis_increment;    

    /// iv0+cube_vertex_increment[k] = k'th vertex of cube with primary vertex iv0
    VTYPE * cube_vertex_increment;

    /// \brief Increment for computing facet vertices.
    /// iv0+facet_vertex_increment[k+num_facet_vertices*d] = 
    ///   k'th vertex of facet orthogonal to d with primary vertex iv0
    VTYPE * facet_vertex_increment;

    /// \brief Increment for computing ridge vertices.
    /// iv0+ridge_vertex_increment[k+nv*d0+nv*nv*d1] = 
    ///   k'th vertex of ridge orthogonal to d0 and d1 with primary vertex iv0
    ///   nv = number of ridge vertices.
    VTYPE * ridge_vertex_increment;

    /// unit_cube_coord[dimension*k+j] = j'th coordinate of k'th vertex of unit cube
    NTYPE * unit_cube_coord;

    NTYPE num_cube_vertices;   ///< Number of cube vertices.
    NTYPE num_facet_vertices;  ///< Number of facet vertices.
    NTYPE num_cube_facets;     ///< Number of cube facets.
    NTYPE num_cube_edges;      ///< Number of cube edges.

    /// Number of cube ridge vertices.
    NTYPE num_cube_ridge_vertices;

    void ZeroLocal();
    void InitLocal();
    void FreeLocal();
    void Create();            ///< Allocate and set data in GRID_PLUS.

  public:
    /// Constructors.
    GRID_PLUS(const DTYPE dimension, const ATYPE * axis_size);
    GRID_PLUS();
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_PLUS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_PLUS(const GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> & grid2);

    ~GRID_PLUS();                       ///< Destructor

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    NTYPE NumCubeVertices() const       /// Return number of cube vertices.
    { return(num_cube_vertices); };
    NTYPE NumCubeFacets() const         /// Return number of cube facets.
    { return(num_cube_facets); };
    NTYPE NumCubeEdges() const          /// Return number of cube edges.
    { return(num_cube_edges); };
    NTYPE NumDiagonals() const          /// Return number of cube diagonals.
    { return(num_facet_vertices); };
    const VTYPE * AxisIncrement() const /// Return array axis_increment[]
    { return(axis_increment); }
    const VTYPE AxisIncrement(const DTYPE d) const /// Return axis_increment[d]
    { return(axis_increment[d]); }
    const VTYPE * CubeVertexIncrement() const /// Return cube_vertex_increment[]
    { return(cube_vertex_increment); }
    const VTYPE CubeVertexIncrement     /// Return cube_vertex_increment[k]
    (const VTYPE k) const              
    { return(cube_vertex_increment[k]); }

    /// Return number of cube facet vertices.
    NTYPE NumFacetVertices() const
    { return(num_facet_vertices); };

    /// Return number of cube ridge vertices.
    NTYPE NumCubeRidgeVertices() const
    { return(num_cube_ridge_vertices); };

    /// Return vertex increment of k'th vertex of facet ifacet
    const VTYPE FacetVertexIncrement
    (const DTYPE ifacet, const VTYPE k) const              
    { return(facet_vertex_increment[k+ifacet*num_facet_vertices]); }

    /// Return vertex increment of k'th vertex of ridge orthogonal
    ///   to orth_dir0 and orth_dir1.
    const VTYPE RidgeVertexIncrement    /// Return ridge_vertex_increment[k]
    (const DTYPE orth_dir0, const DTYPE orth_dir1, const VTYPE k) const
    { 
      const DTYPE dimension = this->Dimension();
      NTYPE j = 
        k+this->num_cube_ridge_vertices*(orth_dir0+dimension*orth_dir1);
      return(ridge_vertex_increment[j]);
    }

    // *** DEPRECATED. REPLACE BY class UNIT_CUBE. ***
    /// Return pointer to coordinates of k'th cube vertex
    const NTYPE * UnitCubeCoord(const NTYPE k) const
    { return(unit_cube_coord+this->Dimension()*k); }

    // *** DEPRECATED. REPLACE BY class UNIT_CUBE. ***
    /// Return j'th coordinate of k'th vertex
    const NTYPE UnitCubeCoord         
    (const NTYPE k, const NTYPE j) const
    { return(unit_cube_coord[this->Dimension()*k+j]); }

    /// \brief Return next vertex in direction d.
    /// @pre iv is not the last vertex in direction d.
    VTYPE NextVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv+axis_increment[d]); }

    /// \brief Return previous vertex in direction d.
    /// @pre iv is not the first vertex in direction d.
    VTYPE PrevVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv-axis_increment[d]); }

    /// \brief Return adjacent vertex in direction d.
    /// @param side = 0 (or false) or 1 (or true).
    template <typename STYPE>
    VTYPE AdjacentVertex(const VTYPE iv, const DTYPE d, const STYPE side) const
    {
      if (DTYPE(side) == 0) { return(PrevVertex(iv, d)); }
      else { return(NextVertex(iv, d)); };
    }

    /// \brief Return k'th cube vertex.
    /// @param iv0 is a primary cube vertex.
    /// @param k k'th cube vertex.
    /// @pre k is less than the number of unit cube vertices.
    VTYPE CubeVertex(const VTYPE iv0, const int k) const  
    { return(iv0+cube_vertex_increment[k]); }

    /// \brief Return k'th facet vertex.
    /// @param iv0 is a primary facet vertex.
    /// @param ifacet Facet index. In range [0..(Nf-1)]
    ///        where Nf is the number of cube facets.
    ///        If \a ifacet < (Nf/2), then facet \a ifacet is the
    ///          lower facet orthogonal to direction \a ifacet.
    ///        If \a ifacet >= (Nf/2), then facet \a ifacet is the
    ///          upper facet orthogonal to direction (\a ifacet/2).
    /// @param k \a k'th facet vertex.
    /// @pre \a ifacet is less than the number of cube facets.
    /// @pre \a k is less than the number of facet vertices.
    VTYPE FacetVertex(const VTYPE iv0, const DTYPE ifacet, const int k) const
    { return(iv0+facet_vertex_increment[k+ifacet*num_facet_vertices]); }

    /// \brief Return k'th ridge vertex.
    VTYPE RidgeVertex(const VTYPE iv0, 
                      const DTYPE orth_dir0, const DTYPE orth_dir1, 
                      const int k) const
    { return(iv0+RidgeVertexIncrement(orth_dir0, orth_dir1, k)); }

    /// \brief Compute vertex neighbors.
    template <typename VTYPE0, typename VTYPE1, typename DIST_TYPE>
    void GetVertexNeighbors
    (const VTYPE0 iv0, const DIST_TYPE distance, std::vector<VTYPE1> & vlist)
      const;
  };

  // **************************************************
  // TEMPLATE CLASS GRID_NEIGHBORS
  // **************************************************

  /// Inherits class GRID_PLUS and adds other indexes and operators 
  ///   for fast accessing neighbors of grid vertices.
  /// @tparam DTYPE  Dimension data type.
  /// @tparam ATYPE  Axis size type.
  /// @tparam VTYPE  Vertex index type.
  /// @tparam DIFFTYPE  Index difference type.  Must be signed.
  /// @tparam NTYPE  Number type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename DIFFTYPE, 
            typename NTYPE> 
  class GRID_NEIGHBORS:public GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> {

  protected:

    /// Number of neighbors of a vertex in cubes containing the vertex.
    NTYPE num_vertex_neighborsC; 

    /// iv + vertex_neighborC[k] = k'th neighbor of vertex iv in cubes containing iv.
    DIFFTYPE * vertex_neighborC;

    /// Number of vertices which share an edge with a vertex.
    NTYPE num_vertex_neighborsE;

    /// iv + vertex_neighborE[k] = k'th neighbor sharing an edge with vertex iv.
    DIFFTYPE * vertex_neighborE;

    /// cube_index + cube_neighborE[k] = 
    ///   k'th neighbor sharing an edgex with cube cube_index.
    DIFFTYPE * cube_neighborE;

    /// cube_index + cube_neighborV[k] = 
    ///   k'th neighbor sharing a vertex with cube cube_index.
    DIFFTYPE * cube_neighborV;

    /// \brief Number of vertices in cubes containing a facet, 
    ///   not including facet vertices.
    NTYPE num_facet_neighborsC;

    /// \brief iv + facet_neighborC[d*num_facet_neighborsC+k] = 
    ///   k'th neighbor of facet d, primary vertex iv
    DIFFTYPE * facet_neighborC;

    /// \brief Number of vertices in 2-faces containing an edge,
    ///   not including edge vertices.
    NTYPE num_edge_neighborsF2;

    /// \brief iv + edge_neighborC[d*num_edge_neighborsF2 + K] =
    ///   k'th neighbor of edge d, primary vertex iv
    DIFFTYPE * edge_neighborF2;

    void ZeroLocal();
    void InitLocal();
    void FreeLocal();
    void CreateLocal();       ///< Allocate and set data in GRID_NEIGHBORS.

  public:
    /// Constructors.
    GRID_NEIGHBORS(const DTYPE dimension, const ATYPE * axis_size);
    GRID_NEIGHBORS();
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_NEIGHBORS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_NEIGHBORS
    (const GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE, NTYPE> & grid2);

    ~GRID_NEIGHBORS();        ///< Desctructor.

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions

    /// \brief Return number of neighbors of a vertex.
    /// \details All vertices in a cube containing the vertex 
    ///    are counted as neighbors.
    /// The vertex itself is not counted in this number.
    NTYPE NumVertexNeighborsC() const   
    { return(num_vertex_neighborsC); };

    /// \brief Return number of vertices which share an edge with a vertex.
    /// \details The vertex itself is not counted in this number.
    NTYPE NumVertexNeighborsE() const   
    { return(num_vertex_neighborsE); };

    /// \brief Return number of cubes which share a facet with a cube.
    NTYPE NumCubeNeighborsF() const   
    { return(this->NumCubeFacets()); };

    /// \brief Return number of cubes which share an edge with a cube.
    NTYPE NumCubeNeighborsE() const   
    { return(this->NumCubeEdges()); };

    /// \brief Return number of cubes which share a vertex with a cube.
    NTYPE NumCubeNeighborsV() const   
    { return(this->NumCubeVertices()); };

    /// \brief Return number of neighbors of a facet.
    /// \details All vertices in a cube containing the facet are counted as neighbors,
    ///          not including vertices lying on the facet.
    NTYPE NumFacetNeighborsC() const   
    { return(num_facet_neighborsC); };

    /// \brief Return number of neighbors of an edge.
    /// \details All vertices in a 2-face containing the fedge are counted as neighbors,
    ///          not including the edge endpoints.
    NTYPE NumEdgeNeighborsF2() const   
    { return(num_edge_neighborsF2); };

    /// \brief Return k'th neighbor of vertex \a iv in cubes containing \a iv.
    /// \details A grid vertex is not a neighbor of itself.
    /// @pre Vertex \a iv must be an internal grid vertex, 
    ///    i.e., not on the grid boundary.
    VTYPE VertexNeighborC               
    (const VTYPE iv, const NTYPE k) const
    { return(iv+vertex_neighborC[k]); };

    /// \brief Return k'th vertex which shares an edge with the vertex.
    /// \details A grid vertex is not a neighbor of itself.
    /// @pre Vertex \a iv must be an internal grid vertex, 
    ///    i.e., not on the grid boundary.
    VTYPE VertexNeighborE
    (const VTYPE iv, const NTYPE k) const
    { return(iv+vertex_neighborE[k]); };

    /// \brief Return k'th cube which shares a facet with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborF
    (const VTYPE cube_index, const NTYPE k) const
    { return(this->VertexNeighborE(cube_index, k)); };

    /// \brief Return k'th cube which shares an edge with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborE
    (const VTYPE cube_index, const NTYPE k) const
    { return(cube_index+cube_neighborE[k]); };

    /// \brief Return k'th cube which shares a vertex with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborV
    (const VTYPE cube_index, const NTYPE k) const
    { return(cube_index+ cube_neighborV[k]); };

    /// \brief Return k'th neighbor of facet.
    /// @param iv  Primary facet vertex (facet vertex with lowest coordinates.)
    /// @param orth_dir  Direction orthogonal to facet.
    /// @param k   Return \a k'th neighbor.
    /// @pre Facet must NOT be contained in the grid boundary.
    VTYPE FacetNeighborC
    (const VTYPE iv, const DTYPE orth_dir, const NTYPE k) const
    { return(iv+facet_neighborC[orth_dir*num_facet_neighborsC+k]); };

    /// \brief Return k'th neighbor of edge.
    /// @param iv  Primary edge vertex (edge endpoint with lowest coordinates.)
    /// @param dir  Direction of edge.
    /// @param k   Return \a k'th neighbor.
    /// @pre Edge must NOT be contained in the grid boundary.
    VTYPE EdgeNeighborF2
    (const VTYPE iv, const DTYPE dir, const NTYPE k) const
    { return(iv+edge_neighborF2[dir*num_edge_neighborsF2+k]); };

  };

  // **************************************************
  // TEMPLATE CLASS GRID_SPACING
  // **************************************************

  /// Template to add grid spacing to grid.
  template <typename STYPE, typename GRID_TYPE>
  class GRID_SPACING:public GRID_TYPE
  {
  protected:
    STYPE * spacing;   ///< spacing[i] = grid spacing along axis i.

    void InitLocal();                    ///< Initialize GRID_SPACING.
    void CreateLocal();                  ///< Create local data structure.
    void FreeLocal();                    ///< Free memory.

  public:
    GRID_SPACING();                      ///< Constructructor.
    template <typename DTYPE2, typename ATYPE2>
    GRID_SPACING(const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_SPACING(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_SPACING(const GRID_SPACING<STYPE, GRID_TYPE> & grid2);
    ~GRID_SPACING();                     ///< Destructor.

    template <typename DTYPE2, typename STYPE2>    ///< Set spacing[d] to \a c.
    void SetSpacing(const DTYPE2 d, const STYPE2 c)
    { spacing[d] = c; };

    template <typename STYPE2>           ///< Set all grid spacing to \a c.
    void SetAllSpacing(const STYPE2 c);

    template <typename STYPE2>           ///< Set spacing[d] to spacing2[d].
    void SetSpacing(const STYPE2 * spacing2);

    /// Set spacing[d] to (c*spacing2[d]).
    template <typename SCALE_TYPE, typename STYPE2>
    void SetSpacing(const SCALE_TYPE c, const STYPE2 * spacing2);

    /// Set spacing[d] to grid2.Spacing(d).
    template <typename STYPE2, typename GRID_TYPE2>
    void SetSpacing(const GRID_SPACING<STYPE2, GRID_TYPE2> & grid2);

    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    /// Get spacing.
    template <typename DTYPE2>
    STYPE Spacing(const DTYPE2 d) const
    { return(spacing[d]); };

    /// Get const pointer to spacing.
    const STYPE * SpacingPtrConst() const
    { return(spacing); }

    /// Compute scaled coord.
    template <typename VTYPE2, typename CTYPE>
    void ComputeScaledCoord
    (const VTYPE2 iv, CTYPE * coord) const;

    /// Compute scaled coordinates of cube center.
    template <typename VTYPE2, typename CTYPE>
    void ComputeCubeCenterScaledCoord
    (const VTYPE2 iv, CTYPE * coord) const;
  };

  // **************************************************
  // TEMPLATE CLASS FACET_LIST
  // **************************************************

  /// Base class for list of facets.
  template <typename DTYPE>
  class FACET_LIST {

  protected:
    FACET_LIST(){};        ///< Constructor.

  public:
    inline bool Contains(const DTYPE d) const;
  };

  /// list containing zero facets
  template <typename DTYPE>
  class FACET_LIST0:public FACET_LIST<DTYPE> {

  public:
    FACET_LIST0() {};      ///< Constructor.

    /// Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { return(false); };
  };

  /// List containing one facet.
  template <typename DTYPE>
  class FACET_LIST1:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;       ///< Direction orthogonal to facet0.

  public:

    /// Set facet1.
    inline void Set(const DTYPE facet0) { this->facet0 = facet0; };

    /// Constructor.
    FACET_LIST1(const DTYPE facet0) { Set(facet0); };

    /// Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0) { return(true); };
      return(false);
    };
  };

  /// List containing two facets.
  template <typename DTYPE>
  class FACET_LIST2:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;       ///< Direction orthogonal to facet0.
    DTYPE facet1;       ///< Direction orthogonal to facet1.

  public:

    /// Set facet0 and facet1.
    inline void Set(const DTYPE facet0, const DTYPE facet1) 
    { this->facet0 = facet0; this->facet1 = facet1; };

    /// Constructor.
    FACET_LIST2(const DTYPE facet0, const DTYPE facet1) 
    { Set(facet0, facet1); };

    /// Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0 || d == facet1) { return(true); };
      return(false);
    };
  };

  // **************************************************
  // TEMPLATE CLASS GRID_VERTEX_LIST
  // **************************************************


  /// Class containing list of grid vertices.
  template <typename VTYPE>
  class GRID_VERTEX_LIST {

  protected:
    VTYPE * vertex_list;         ///< Vertex list.
    VTYPE list_length;           ///< Length of array vertex_list[].
    VTYPE num_vertices;          ///< Number of vertices in vertex_list[].

    void Init();                 ///< Initialize class.
    void FreeAll();              ///< Free all allocated memory.
    void AllocateList(const VTYPE n); ///< Allocate vertex list.

  public:
    GRID_VERTEX_LIST() 
    { Init(); };
    ~GRID_VERTEX_LIST()
    { FreeAll(); }

    // get functions
    VTYPE VertexIndex(const VTYPE i) const
    { return(vertex_list[i]); }

    VTYPE NumVertices() const
    { return(num_vertices); }

    VTYPE ListLength() const
    { return(list_length); }
  };

  /// Class containing list of grid boundary vertices.
  template <typename VTYPE>
  class GRID_BOUNDARY_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:
    /// GRID_BOUNDARY_VERTEX_LIST constructor.
    /// @param grid Grid.
    template<typename GCLASS>
    GRID_BOUNDARY_VERTEX_LIST(const GCLASS & grid);
  };


  /// Class containing list of grid cubes.
  template <typename VTYPE>
  class GRID_CUBE_LIST:protected GRID_VERTEX_LIST<VTYPE> {

  public:
    GRID_CUBE_LIST() {};
    ~GRID_CUBE_LIST() {};

    // get functions
    VTYPE CubeIndex(const VTYPE i) const
    { return(this->VertexIndex(i)); }

    VTYPE NumCubes() const
    { return(this->NumVertices()); }

    VTYPE ListLength() const
    { return(GRID_VERTEX_LIST<VTYPE>::ListLength()); }
  };

  /// Class containing list of facet 0 cubes.
  template <typename VTYPE>
  class FACET0_CUBE_LIST:public GRID_CUBE_LIST<VTYPE> {

  public:
    template <typename GCLASS>
    FACET0_CUBE_LIST(const GCLASS & grid)
    { GetCubes(grid); }

    /// Get cubes from grid and store in list.
    /// Reallocates list if (current list length < num cubes in facet0)
    template <typename GCLASS>
    void GetCubes(const GCLASS & grid);
  };

  /// Class containing list of facet vertices.
  /// List length is max size
  template <typename VTYPE>
  class FACET_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:
    /// FACET_VERTEX_LIST constructor.
    /// Get and store vertices in facet orthogonal to \a orth_dir.
    /// @param grid Grid.
    /// @param orth_dir Directional orthogonal to facet.
    /// @param allocate_max If true, allocate the list length
    ///            to be max number of interior vertices over all facets.
    template <typename GCLASS>
    FACET_VERTEX_LIST
    (const GCLASS & grid, const VTYPE orth_dir,
     const bool allocate_max=true);

    /// Get vertices from grid and store in list.
    /// Reallocates list if (current list length < num vertices in facet)
    template <typename GCLASS>
    void GetVertices(const GCLASS & grid, const VTYPE orth_dir);
  };

  /// Class containing list of facet interior vertices.
  /// List length is max size
  template <typename VTYPE>
  class FACET_INTERIOR_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:
    /// FACET_INTERIOR_VERTEX_LIST constructor.
    /// Get and store interior vertices in facet orthogonal to \a orth_dir.
    /// @param grid Grid.
    /// @param orth_dir Directional orthogonal to facet.
    /// @param allocate_max If true, allocate the list length
    ///            to be max number of interior vertices over all facets.
    template <typename GCLASS>
    FACET_INTERIOR_VERTEX_LIST
    (const GCLASS & grid, const VTYPE orth_dir,
     const bool allocate_max=true);

    /// Get vertices from grid and store in list.
    /// Reallocates list if (current list length < num vertices in facet)
    template <typename GCLASS>
    void GetVertices(const GCLASS & grid, const VTYPE orth_dir);
  };

  // **************************************************
  // inline UTILITY FUNCTIONS
  // **************************************************

  /// Integer divide.
  template <typename ATYPE, typename BTYPE>
  inline long integer_divide(const ATYPE a, const BTYPE b)
  { return(long(a)/long(b)); }

  /// Integer divide.
  inline long integer_divide(const long a, const long b)
  { return(a/b); }

  /// Integer divide.
  inline int integer_divide(const int a, const int b)
  { return(a/b); }
  
  /// Integer divide.
  inline unsigned long 
  integer_divide(const unsigned long a, const unsigned long b)
  { return(a/b); }

  /// Integer divide.
  inline unsigned int
  integer_divide(const unsigned int a, const unsigned int b)
  { return(a/b); }

  // **************************************************
  // THROW ERROR ROUTINES
  // **************************************************

  /// Throw subsample period error.
  template <typename STRING_TYPE, typename PTYPE>
  void throw_subsample_period_error
  (const STRING_TYPE proc_name, const PTYPE subsample_period)
  {
    IJK::PROCEDURE_ERROR error(proc_name);
    error.AddMessage("Subsample period must be a positive integer.");
    error.AddMessage("  Subsampling period = ", subsample_period, ".");
    throw error;
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: COUNTING AND INDEXING
  // **************************************************

  ///
  /// \defgroup counting Counting and indexing
  /* \ingroup counting */
  /* \{ */

  /// Return coordinate of vertex \a iv.
  /// @param iv = Vertex index.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param[out] coord Array. <em>coord[k]</em> = \a k'th coordinate of vertex \a iv.
  /// @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
  /// @pre \li Array coord[] is pre-allocated to size at least \a dimension.
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename GTYPE>
  void compute_coord(const VTYPE iv, const DTYPE dimension,
                     const ATYPE * axis_size, GTYPE * coord)
  {
    VTYPE k = iv;
    for (DTYPE d = 0; d < dimension; d++) {
      coord[d] = GTYPE(k % axis_size[d]);
      k = k / axis_size[d];
    };
  }

  /// Return coordinate of vertex \a iv.
  /// @param iv = Vertex index.
  /// @param dimension  Dimension of grid.
  /// @param axis_size[]
  ///    <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param scale[]
  ///    <em>scale[d]</em> = Scale along axis \a d.
  /// @param[out] scaled_coord[] 
  ///    <em>scaled_coord[k]</em> = \a k'th coordinate of vertex \a iv.
  /// @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
  /// @pre \li scale[d] > 0 for all \a d = 0,..., \a dimension-1.
  /// @pre \li Array coord[] is pre-allocated to size at least \a dimension.
  template <typename VTYPE, typename DTYPE, typename ATYPE, 
            typename STYPE, typename CTYPE>
  void compute_scaled_coord
  (const VTYPE iv, const DTYPE dimension, const ATYPE * axis_size, 
   const STYPE * scale, CTYPE * scaled_coord)
  {
    compute_coord(iv, dimension, axis_size, scaled_coord);
    for (DTYPE d = 0; d < dimension; d++)
      { scaled_coord[d] *= scale[d]; }
  }

  /// Return index of vertex with specified coord.
  template <typename VTYPE, typename GTYPE, typename DTYPE, typename ATYPE>
  VTYPE compute_vertex_index
  (const GTYPE * coord, const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE iv = 0;
    VTYPE inc = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      iv += inc*coord[d];
      inc = inc*axis_size[d];
    }

    return(iv);
  }

  /// Return number of vertices in grid or subgrid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param facet_list = List of facets determining subgrid.
  ///   If empty, return number of vertices in entire grid.
  ///   Otherwise, return number of vertices contained in the intersection of all facets in facet_list.
  /// @param[out] num_vertices Number of vertices.
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const FTYPE & facet_list,
   NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) 
        { num_vertices = num_vertices * axis_size[d]; }
    }
  }

  /// Return number of subsampled vertices along axis
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param subsample_period = only count every k'th vertex along axis where k = subsample_period.
  /// @pre subsample_period is a positive integer.
  template <typename ATYPE, typename PTYPE>
  inline ATYPE compute_num_subsampled_vertices_along_axis
  (const ATYPE axis_size, const PTYPE subsample_period)
  { 
    return(integer_divide(axis_size+subsample_period-1, subsample_period)); 
  }

  /// Return number of vertices in subsampled grid or grid subspace.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename FTYPE, 
            typename NTYPE>
  void compute_num_subsampled_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const PTYPE subsample_period, const FTYPE & facet_list,
   NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
    // subsample_period[d] = only count every k'th vertex along axis d
    //   where k = subsample_period[d]
    // facet_list = ignore facets in facet_list
    // num_vertices = number of vertices.
  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {

        if (subsample_period[d] < 1) 
          { throw_subsample_period_error
              ("compute_num_subsampled_vertices", subsample_period[d]); }

        ATYPE num_subsampled_vertices_along_axis =
          compute_num_subsampled_vertices_along_axis
          (axis_size[d], subsample_period[d]);
        num_vertices = num_vertices * num_subsampled_vertices_along_axis;
      }
    }
  }

  /// Return number of grid vertices
  /// NOTE: REPLACE BY FASTER VERSION WHICH DOES NOT USE FACET_LIST0.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  { 
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_vertices(dimension, axis_size, facet_list0, num_vertices);
  }

  /// Return number of grid vertices between two vertices.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param iv0 = Vertex index.
  /// @param iv1 = Vertex index.
  /// @param[out] num_vertices = Number of grid vertices between \a iv0 and \a iv1.
  /// @pre 0 <= iv0 <= iv1 < (total number of grid vertices)
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, NTYPE & num_vertices)
  {
    VTYPE coord0[dimension];
    VTYPE coord1[dimension];
    IJK::PROCEDURE_ERROR error("compute_num_grid_vertices");

    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] > coord1[d]) {
        error.AddMessage("Programming error in calculating ", d,
                         "'th coordinate.");
        error.AddMessage("  coord0 (", coord0[d], 
                         ") > coord1 (", coord1[d], ").");
        throw error;
      }

      num_vertices = num_vertices * (coord1[d]-coord0[d]+1);
    }
  }

  /// Return number of vertices in grid interior.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    ATYPE interior_axis_size[dimension];

    num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 2*boundary_width) { return; };
      interior_axis_size[d] = axis_size[d]-2*boundary_width;
    }
    compute_num_grid_vertices(dimension, interior_axis_size, num_vertices);
  }

  /// Return number of vertices in grid interior for boundary width 1.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  {
    compute_num_interior_grid_vertices
      (dimension, axis_size, 1, num_vertices);
  }

  /// Return number of vertices in grid boundary.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, NTYPE & num_boundary_vertices)
  {
    NTYPE num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    NTYPE num_interior_vertices;
    compute_num_interior_grid_vertices
      (dimension, axis_size, boundary_width, num_interior_vertices);
    num_boundary_vertices = num_grid_vertices - num_interior_vertices;
  }

  /// Return number of vertices in grid boundary.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_boundary_vertices)
  {
    compute_num_boundary_grid_vertices
      (dimension, axis_size, 1, num_boundary_vertices);
  }

  /// Return number of cubes in grid or grid subspace
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   const FTYPE & facet_list, NTYPE & num_cubes)
    // facet_list = ignore facets in facet_list
  {
    if (dimension < 1) {
      num_cubes = 0;
      return;
    }

    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
        if (axis_size[d] < 2) { 
          num_cubes = 0;
          return; 
        };
        num_cubes = num_cubes * (axis_size[d]-1);
      }
    }
  }

  /// Return number of grid cubes
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_cubes(dimension, axis_size, facet_list0, num_cubes);
  }

  /// Return number of grid edges
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_edges
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_edges)
  {
    num_edges = 0;
    for (DTYPE orth_dir = 0; orth_dir < dimension; orth_dir++) {

      if (axis_size[orth_dir] >= 2) {
        NTYPE num_vertices_in_facet;
        compute_num_vertices_in_grid_facet
          (dimension, axis_size, orth_dir, num_vertices_in_facet);

        num_edges += (num_vertices_in_facet*(axis_size[orth_dir]-1));
      }
    }
  }

  /// Return number of cubes in grid interior
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_interior_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_interior_cubes)
  {
    if (dimension <= 0) {
      num_interior_cubes = 0;
      return;
    }

    num_interior_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 3) { 
        num_interior_cubes = 0;
        return; 
      };
      num_interior_cubes = num_interior_cubes*(axis_size[d]-3);
    }
  }

  /// Return number of cubes in grid boundary
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_boundary_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_boundary_cubes)
  {
    NTYPE num_cubes, num_interior_cubes;

    compute_num_grid_cubes(dimension, axis_size, num_cubes);
    compute_num_interior_grid_cubes
      (dimension, axis_size, num_interior_cubes);

    num_boundary_cubes = num_cubes - num_interior_cubes;
  }

  /// Return number of inner vertices in grid or grid subspace.
  /// Inner vertices are vertices which are not on an outer facet of the grid.
  /// Outer facets are grid facets which do not contain the origin.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_inner_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    if (dimension < 1) { 
      num_vertices = 0;
      return; 
    };

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE asize = axis_size[d];
      if (asize < 2) { 
        num_vertices = 0;
        return; 
      }
      else { num_vertices = num_vertices * (asize-1);  }
    }
  }

  /// Return number of outer vertices in grid or grid subspace.
  /// Outer vertices are vertices which are are on an outer facet of the grid.
  /// Outer facets are grid facets which do not contain the origin.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_outer_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_outer_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    NTYPE num_grid_vertices;
    NTYPE num_inner_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    compute_num_inner_vertices(dimension, axis_size, num_inner_vertices);
    num_outer_vertices = num_grid_vertices-num_inner_vertices;
  }

  /// Return number of grid vertices in a region whose edges all have the same length.
  /// @param dimension  Dimension of grid.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @param[out] num_region_vertices  Number of grid vertices in a region.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_vertices)
  {
    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_vertices = num_region_vertices * (region_edge_length+1);
    }
  }

  /// Return number of grid cubes in a region whose edges all have the same length.
  /// @param dimension  Dimension of grid.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @param[out] num_region_cubes  Number of grid cubes in a region.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_cubes)
  {
    if (region_edge_length < 1) { 
      num_region_cubes = 0; 
      return;
    };

    num_region_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_cubes = num_region_cubes * region_edge_length;
    }
  }

  /// \brief Return number of vertices in a region.  
  ///        Regions completely contained in the interior of the grid all have the
  ///        same number of vertices.
  ///        Regions bounded by the grid boundary will have fewer vertices.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param iv0  Primary vertex (lowest x,y,z,... coordinates) in region.
  /// @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
  /// @param[out] num_region_vertices  Number of vertices in region.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const ATYPE max_region_edge_length,
   NTYPE & num_region_vertices)
  {
    ATYPE coord[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE numv_along_axis = max_region_edge_length + 1;
      if (coord[d] + max_region_edge_length >= axis_size[d]) {
        if (coord[d] < axis_size[d]) 
          { numv_along_axis = axis_size[d] - coord[d]; }
        else
          { numv_along_axis = 0; };
      }
      num_region_vertices = num_region_vertices * numv_along_axis;
    }
  }

  /// Return number of regions along a single axis.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <typename ATYPE>
  ATYPE compute_num_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    ATYPE num_regions_along_axis = 
      long(axis_size+region_edge_length-2)/long(region_edge_length);
    return(num_regions_along_axis);
  }

  /// Return total number of regions in grid or subgrid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @param facet_list = List of facets determining subgrid.
  ///   If empty, return number of regions in entire grid.
  ///   Otherwise, return number of regions in subgrid formed by intersection of all facets in facet_list.
  /// @param[out] num_regions Number of regions.
  /// @pre \a region_edge_length > 0.
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (!facet_list.Contains(d)) {
        if (axis_size[d] <= 1) { 
          num_regions = 0;
          return; 
        };

        ATYPE num_regions_along_axis = 
          compute_num_regions_along_axis(axis_size[d], region_edge_length);
        num_regions = num_regions * num_regions_along_axis; 
      }
    }
  }

  /// Return total number of regions.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @param[out] num_regions Number of regions.
  /// @pre \a region_edge_length > 0.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_regions
      (dimension, axis_size, region_edge_length, facet_list0, num_regions);
  }

  /// Return number of full regions along a single axis.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <typename ATYPE>
  ATYPE compute_num_full_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    ATYPE num_full_regions_along_axis = 
      long(axis_size-1)/long(region_edge_length);
    return(num_full_regions_along_axis);
  }

  /// Return number of full regions in grid or subgrid.
  /// @param dimension Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length Number of grid edges contained in each region edge.
  /// @param facet_list List of facets determining subgrid.
  ///   If empty, return number of vertices in entire grid.
  ///   Otherwise, return number of vertices in subgrid. formed by intersection of all facets in facet_list.
  /// @param[out] num_full_regions Number of full regions in grid or subgrid.
  /// @pre \a region_edge_length > 0.
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_full_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_full_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_full_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
        NTYPE k = compute_num_full_regions_along_axis
          (axis_size[d], region_edge_length);
        num_full_regions = num_full_regions*k;
      }
    }
  }

  /// Return number of full regions in grid.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_full_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_full_regions);
  }

  /// Return number of partial regions along a single axis (0 or 1.)
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <typename ATYPE>
  ATYPE compute_num_partial_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    if (long(axis_size-1)%long(region_edge_length) == 0) { return(0); };
    return(1);
  }

  /// Return number of partial regions in grid or subgrid.
  template <typename DTYPE, typename ATYPE, 
            typename FTYPE, typename NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_partial_regions)
  {
    NTYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_regions);
    NTYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_full_regions);
    num_partial_regions = num_regions-num_full_regions;
  }

  /// Return number of partial regions in grid.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, 
   NTYPE & num_partial_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_partial_regions);
  }

  /// Return total number of regions in grid facet.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param orth_dir = Direction orthogonal to facet.
  template <typename NTYPE, typename DTYPE, typename ATYPE>
  NTYPE compute_num_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE * region_edge_length, const DTYPE orth_dir)
  {
    NTYPE num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (d != orth_dir) {
        if (axis_size[d] <= 1) { return(0); };

        ATYPE num_regions_along_axis = 
          compute_num_regions_along_axis(axis_size[d], region_edge_length[d]);
        num_regions = num_regions * num_regions_along_axis; 
      }
    }
    return(num_regions);
  }

  /// Return total number of full regions in grid facet.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length  Number of grid edges contained in each region edge.
  /// @param orth_dir  Direction orthogonal to facet.
  /// @param[out] num_full_regions  Number of full regions in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_full_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_full_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_full_regions);
  }

  /// Return total number of partial regions in grid facet.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_partial_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_partial_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_partial_regions);
  }

  /// Return number of vertices along subsampled axis.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param subsample_period = Only count every k'th vertex along axis where k = \a subsample period.
  /// @pre \a subsample_period > 0.
  template <typename ATYPE, typename PTYPE>
  ATYPE compute_subsample_size
  (const ATYPE axis_size, const PTYPE subsample_period)
  {
    ATYPE subsample_size =
      integer_divide(axis_size+subsample_period-1, subsample_period);
    return(subsample_size);
  }

  /// Compute axis size for each axis in subsampled grid.
  /// @param dimension Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = 
  ///                     Number of vertices along axis \a d.
  /// @param[out] subsampled_axis_size Array: 
  ///     <em>subsampled_axis_size[d]</em> = 
  ///        Number of vertices along subsampled axis \a d.
  /// @pre Array subsampled_axis_size[] is preallocated 
  ///        to size at least \a dimension.
  template <typename DTYPE, typename ATYPE, typename PTYPE, 
            typename ATYPE2>
  void compute_subsample_axis_sizes
  (const DTYPE dimension, const ATYPE axis_size[], 
   const PTYPE subsample_period, ATYPE2 subsampled_axis_size[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      subsampled_axis_size[d] = 
        compute_subsample_size(axis_size[d], subsample_period); 
    }
  }

  /// Return number of vertices in subsampled grid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param subsample_period  
  ///        Array: <em>subsample_period[d]</em> = subsample period along axis \a d.
  /// @param[out] num_vertices Number of vertices.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename NTYPE>
  void compute_subsample_size
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, NTYPE & num_vertices)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_subsampled_vertices
      (dimension, axis_size, subsample_period, facet_list0, num_vertices);
  }

  /// Return number of vertices along supersampled axis.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param supersample_period = Only count every k'th vertex along axis where k = \a supersample period.
  /// @pre \a supersample_period > 0.
  template <typename ATYPE, typename PTYPE>
  ATYPE compute_supersample_size
  (const ATYPE axis_size, const PTYPE supersample_period)
  {
    if (axis_size < 1) { return(0); };
    ATYPE supersample_size = (axis_size-1)*supersample_period+1;
    return(supersample_size);
  }

  /* \} */

  // **************************************************
  // COMPUTE BOUNDARY BITS
  // **************************************************

  /// Compute boundary bits for vertex \a iv.
  /// @param iv Vertex index.
  /// @param dimension Dimension of grid.
  /// @param axis_size  Array. <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param [out] boundary_bits Bits flagging boundaries containing vertex \a iv.
  ///       If bit \a 2d is true, then <em>d</em>'th coordinate of vertex \a iv is zero.
  ///       If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate of vertex \a iv equals
  ///            <em>axis_size[d]-1</em>.
  /// @pre \li Variable \a boundary_bits has at least <em>(2*dimension)</em> bits.
  /// @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename BTYPE>
  void compute_boundary_bits
  (const VTYPE iv, const DTYPE dimension,
   const ATYPE * axis_size, BTYPE & boundary_bits)
  {
    VTYPE k = iv;
    BTYPE flag = 1;
    boundary_bits = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
      if (c+1 >= axis_size[d]) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
    };
  }

  /// Compute boundary bits for cube \a icube.
  /// @param icube Cube index.
  /// @param dimension Dimension of grid.
  /// @param axis_size  Array. 
  ///       <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param [out] boundary_bits Bits flagging boundaries 
  ///              containing vertex \a iv.
  ///       If bit \a 2d is true, then <em>d</em>'th coordinate 
  ///              of cube \a icube is zero.
  ///       If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate 
  ///              of cube \a icube equals <em>axis_size[d]-2</em>.
  /// @pre \li Variable \a boundary_bits has at least <em>(2*dimension)</em> bits.
  /// @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename BTYPE>
  void compute_boundary_cube_bits
  (const VTYPE iv, const DTYPE dimension,
   const ATYPE * axis_size, BTYPE & boundary_bits)
  {
    VTYPE k = iv;
    BTYPE flag = 1;
    boundary_bits = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
      if (c+2 >= axis_size[d]) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
    };
  }

  /// Return true if edge is on grid boundary.
  template <typename VTYPE, typename DIR_TYPE,
            typename DTYPE, typename ATYPE>
  bool is_edge_on_grid_boundary
  (const VTYPE iend0, const DIR_TYPE edge_dir, 
   const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE k = iend0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) {
        if (edge_dir != d)
          { return(true); }
      }
      if (c+1 >= axis_size[d]) { 
        if (edge_dir != d)
          { return(true); }
      }
    }

    return(false);
  }
              
  // **************************************************
  // COMPUTE REGION
  // **************************************************

  /// Compute region within boundary around given cube.
  template <typename VTYPE0, typename DTYPE, typename ATYPE0,
            typename DIST_TYPE, typename VTYPE1, typename ATYPE1>
  void compute_region_around_cube
  (const VTYPE0 icube, const DTYPE dimension, const ATYPE0 * axis_size,
   const DIST_TYPE dist2cube, VTYPE1 & region_iv0, ATYPE1 * region_axis_size)
  {
    IJK::ARRAY<ATYPE0> cube_coord(dimension);
    IJK::ARRAY<ATYPE0> region_iv0_coord(dimension);

    compute_coord(icube, dimension, axis_size, cube_coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (cube_coord[d] > dist2cube) 
        { region_iv0_coord[d] = cube_coord[d]-dist2cube; }
      else
        { region_iv0_coord[d] = 0; }

      ATYPE0 c = cube_coord[d]+dist2cube+2;
      if (c <= axis_size[d]) 
        { region_axis_size[d] = c-region_iv0_coord[d]; }
      else 
        { region_axis_size[d] = axis_size[d]-region_iv0_coord[d]; }
    }

    region_iv0 = 
      compute_vertex_index<VTYPE1>
      (region_iv0_coord.PtrConst(), dimension, axis_size);
  }

  /// Compute region within boundary around given vertex coord
  template <typename CTYPE0, typename DTYPE, typename ATYPE0, 
            typename DIST_TYPE, typename CTYPE1, typename ATYPE1>
  void compute_region_around_vertex_coord
  (const CTYPE0 * vertex_coord, 
   const DTYPE dimension, const ATYPE0 * axis_size, 
   const DIST_TYPE dist2vertex, 
   CTYPE1 * region_iv0_coord, ATYPE1 * region_axis_size)
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (vertex_coord[d] > dist2vertex) 
        { region_iv0_coord[d] = vertex_coord[d]-dist2vertex; }
      else
        { region_iv0_coord[d] = 0; }

      ATYPE0 c = vertex_coord[d]+dist2vertex;
      if (c < axis_size[d]) 
        { region_axis_size[d] = c-region_iv0_coord[d]+1; }
      else 
        { region_axis_size[d] = axis_size[d]-region_iv0_coord[d]; }
    }

  }

  // **************************************************
  // TEMPLATES TO CHECK VALUES AND ARRAYS.
  // **************************************************

  /// Check dimension
  /// Return true if dimension is non-negative
  template <typename DTYPE>
  bool check_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension < 0) {
      error.AddMessage("Illegal dimension ", dimension, ".");
      error.AddMessage("Dimension must be non-negative.");
      return(false);
    }

    return(true);
  }

  /// Check range of vertices
  /// Return true if 0 <= iv0 <= iv1 < total_num_grid_vertices
  template <typename VTYPE, typename DTYPE, typename ATYPE>
  bool check_range
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, IJK::ERROR & error)
  {
    GRID_SIZE_TYPE total_num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, total_num_grid_vertices);

    if (iv0 > iv1) {
      error.AddMessage("Illegal vertex range. Vertex index ", iv0, 
                       " is greater than vertex index ", iv1, ".");
      return(false);
    }

    if (0 > iv0 || iv1 >= total_num_grid_vertices) {
      error.AddMessage("Illegal vertex indices: ", iv0, ", ", iv1, ".");
      error.AddMessage("Vertex indices should be in range: [0,",
                       total_num_grid_vertices, "].");
      return(false);
    }

    return(true);
  }

  /// Check region coordinates.
  template <typename DTYPE, typename VTYPE, typename CTYPE>
  bool check_region_coordinates
  (const DTYPE dimension, const VTYPE iv0, const CTYPE * coord0, 
   const VTYPE iv1, const CTYPE * coord1, IJK::ERROR & error)
    // return true if (coord0[d] <= coord1[d]) for all d
  {
    for (DTYPE d = 0; d < dimension; ++d) {
      if (coord0[d] > coord1[d]) {
        error.AddMessage("Illegal coordinates.  Coordinates of vertex 0 should be less than or equal to coordinates of vertex 1.");
        error.AddMessage(" Vertex 0 = ", iv0, 
                         ".  Coordinate ", d , " = ", coord0[d], ".");
        error.AddMessage(" Vertex 1 = ", iv1, 
                         ".  Coordinate ", d, " = ", coord1[d], ".");
        return(false);
      }
    }

    return(true);
  }

  /* INCORRECT
  /// Check that axis size is positive.
  template <typename DTYPE, typename ATYPE>
  bool check_positive_axis_size
  (const DTYPE dimension, ATYPE * axis_size, IJK::ERROR & error)
  {
    for (DTYPE d = 0; d < dimension; d++) {
      error.AddMessage("Illegal axis size.  Axis size must be positive.");
      error.AddMessage("  Axis ", d, " has size ", axis_size[d], ".");
      return(false);
    }
    return(true);
  }
  */

  /// Check that array vertex_list[] is not NULL.
  template <typename VTYPE>
  bool check_vertex_list(const VTYPE * vertex_list, IJK::ERROR & error)
  {
    if (vertex_list == NULL) {
      error.AddMessage("Vertex list is NULL");
      return(false);
    }
    else { return(true); }
  }

  /// Check that region edge length is positive.
  template <typename LTYPE>
  bool check_region_edge_length(const LTYPE length, IJK::ERROR & error)
  {
    if (length <= 0) {
      error.AddMessage("Region edge length must be positive.");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices added to vertex list.
  template <typename NTYPE0, typename NTYPE1>
  bool check_num_vertices_added
  (const NTYPE0 num_added, const NTYPE1 numv, IJK::ERROR & error)
  {
    if (num_added != numv) {
      error.AddMessage("Added ", num_added, " vertices to vertex list.");
      error.AddMessage("Number of vertices in list should be ", numv, ".");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices equals number of grid vertices
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  bool check_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const NTYPE numv,
   IJK::ERROR & error)
  {
    NTYPE numv2;
    compute_num_grid_vertices(dimension, axis_size, numv2);
    if (numv != numv2) {
      error.AddMessage("Programming error. Incorrect number of vertices.");
      return(false);
    }
    return(true);
  }

  /// Check that DIFFTYPE is signed and has range [-n:n]
  template <typename DIFFTYPE, typename NTYPE>
  bool check_difftype(const NTYPE n, IJK::ERROR & error)
  {
    if (std::numeric_limits<DIFFTYPE>::min() >= 0) {
      error.AddMessage("Programming error. Template parameter DIFFTYPE must be signed.");
      return(false);
    }
    else if (n > std::numeric_limits<DIFFTYPE>::max()) {
      error.AddMessage("Error. Template parameter DIFFTYPE is not large enough to store integer ",
                       n, ".");
      return(false);
    }
    else if (-n < std::numeric_limits<DIFFTYPE>::min()) {
      error.AddMessage("Error. Template parameter DIFFTYPE is not large enough to store integer ",
                       -n, ".");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTING INCREMENTS
  // **************************************************

  /// Compute increment to add to index of a vertex to get next vertices along the axes.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param[out] increment[] = Axis increment. iv+increment[d] is the index of the vertex after iv along axis \a d.
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_axis_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_increment");

    if (dimension <= 0) { return; };

    if (axis_size == NULL || increment == NULL) {
      error.AddMessage("Programming error. axis_size == NULL or increment == NULL.");
      throw error;
    }

    increment[0] = 1;
    for (DTYPE d = 1; d < dimension; d++)
      { increment[d] = increment[d-1]*axis_size[d-1]; }
  }

  /// DEPRECATED
  /// OLD function name.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  {
    compute_axis_increment(dimension, axis_size, increment);
  }

  /// Compute increment to add to index of current vertex to get
  ///   next vertex along each axis
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE, typename ITYPE>
  void compute_increment
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, ITYPE * increment)
  {
    compute_increment(grid.Dimension(), grid.AxisSize(), increment);
  }

  /// Compute increment to add to index of vertex to get next subsampled vertex along each axis.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param subsample_period  Array: <em>subsample_period[d]</em> = Subsample period along axis \a d.
  /// @param[out] increment  Array: <em>increment[d]</em> = Increment to add to index of a vertex to get next subsampled vertex along axis \a d.
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename ITYPE>
  void compute_subsample_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, ITYPE * increment)
  {
    compute_increment(dimension, axis_size, increment);

    for (DTYPE d = 0; d < dimension; ++d)
      { increment[d] *= subsample_period[d]; };
  }

  /// Compute increment to add to index of primary vertex to get
  ///   k'th corner vertex of region.
  /// @param dimension  Dimension of grid.
  /// @param grid_axis_size  Array: <em>grid_axis_size[d]</em> = Number of vertices along grid axis \a d.
  /// @param region_axis_size  Array: <em>region_axis_size[d]</em> = Number of vertices along region axis \a d.
  /// @param[out] increment  Array: Region corner increment. iv0+increment[k] is k'th corner vertex of region.
  /// @pre Array increment[] is allocated with size at least number of corner regions.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE * region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_region_corner_increment");
    IJK::ARRAY<ITYPE> axis_increment(dimension);

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (region_axis_size == NULL) {
      error.AddMessage("Programming error. region_axis_size == NULL.");
      throw error;
    }    
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size[d] < 1) {
        error.AddMessage("Programming error.  Region axis size must be at least 1.");
        throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment.Ptr());
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          increment[j] = increment[j] + (region_axis_size[d]-1)*axis_increment[d];
        };
        j0 = j0/2;
      };
    }
  }

  /// Compute increment to add to index of primary vertex to get
  ///   k'th corner vertex of cubic region (all axes have same size.)
  /// @param dimension  Dimension of grid.
  /// @param grid_axis_size[] = Grid axis size. grid_axis_size[i] is the number of grid vertices along axis i.
  /// @param region_axis_size = Region axis size. Number of region vertices along any axis.
  /// @param[out] increment[] = Region corner increment. iv0+increment[k] is the k'th corner vertex of the region.
  /// @pre Array increment[] is allocated with size at least number of corner regions.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_cubic_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_cubic_region_corner_increment");
    IJK::ARRAY<ITYPE> axis_increment(dimension);

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size < 1) {
        error.AddMessage("Programming error.  Region axis size must be at least 1.");
        throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment.Ptr());
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          increment[j] = increment[j] + (region_axis_size-1)*axis_increment[d];
        };
        j0 = j0/2;
      };
    }
  }


  /// Compute increment to add to index of primary vertex to get k'th vertex in region.
  /// @pre Array increment is allocated with size at least number of region vertices
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_vertex_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const ATYPE scale, ITYPE * increment)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(scale);

    for (DTYPE d = 0; d < dimension; d++) {
      subgrid_axis_size[d] = region_edge_length*scale+1; 

      if (subgrid_axis_size[d] > axis_size[d]) {
        IJK::PROCEDURE_ERROR error("compute_region_vertex_increment");
        error.AddMessage("Region is larger than grid.");
        error.AddMessage("  Grid axis size[", d, "] = ", axis_size[d], ".");
        error.AddMessage("  Region size[", d, "] = ", subgrid_axis_size[d],
                         ".");
        throw error;
      }
    }

    subsample_subgrid_vertices
      (dimension, axis_size, ITYPE(0), subgrid_axis_size.PtrConst(), 
       subsample_period, increment);
  }

  /// Compute increment to add to index of primary vertex of region to get
  ///   primary vertex of k'th grid cube in cubic region.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges along every edge of the cubic region.
  /// @param increment = Region increment. iv0+increment[k] is the k'th cube of the region.
  /// @pre Array increment[] is allocated with size at least number of grid cubes contained in the region.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_cube_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, ITYPE * increment)
  {
    IJK::ARRAY<ATYPE> subgrid_size(dimension);

    if (region_edge_length < 1)  // No cubes.
      { return; }

    for (DTYPE d = 0; d < dimension; d++) 
      { subgrid_size[d] = region_edge_length; }

    get_subgrid_vertices
      (dimension, axis_size, 0, subgrid_size.PtrConst(), increment);
  }

  /// Compute increment to add to vertex 0 to compute vertex i of hypercube.
  /// @param dimension  Dimension of grid.
  /// @param axis_increment[] = Axis increment. iv+axis_increment[i] 
  ///        is the next vertex after vertex iv along axis i.
  /// @param[out] cube_vertex_increment[] = Cube vertex increment. iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube with primary vertex iv0.
  /// @pre Array cube_vertex_increment[] is allocated with size at least number of cube vertices
  template <typename DTYPE, typename ITYPE1, typename ITYPE2>
  void compute_cube_vertex_increment
  (const DTYPE dimension, const ITYPE1 * axis_increment, 
   ITYPE2 * cube_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_vertex_increment");

    if (dimension <= 0) { return; };

    if (axis_increment == NULL) {
      error.AddMessage("Programming error. Array axis_increment[] must be allocated and set before calling compute_cube_vertex_increment.");
      throw error;
    }

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated before calling compute_cube_vertex_increment.");
      throw error;
    }
    
    ITYPE2 num_cube_vertices = compute_num_cube_vertices(dimension);

    for (ITYPE2 j = 0; j < num_cube_vertices; j++) {
      cube_vertex_increment[j] = 0;
      ITYPE2 j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          cube_vertex_increment[j] += axis_increment[d];
        };
        j0 = j0/2;
      }
    }
  }

  /// Compute increment to add to vertex 0 to compute vertex i of hypercube.
  /// @param grid Grid.
  /// @param[out] cube_vertex_increment[] = Cube vertex increment. iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube with primary vertex iv0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE, typename ITYPE>
  void compute_cube_vertex_increment
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, ITYPE * cube_vertex_increment)
  {
    const DTYPE dimension = grid.Dimension();
    VTYPE axis_increment[dimension];

    compute_increment(grid, axis_increment);
    compute_cube_vertex_increment(dimension, axis_increment, cube_vertex_increment);
  }

  /// Compute increment to add to vertex 0 to compute vertex i of facet.
  /// @param dimension  Dimension of grid.
  /// @param cube_vertex_increment[] = Cube vertex increment. iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube with primary vertex iv0.
  /// @param[out] facet_vertex_increment[] = Facet vertex increment. iv0+facet_vertex_increment[i] is the i'th vertex of the facet with primary vertex iv0.
  /// @pre Array facet_vertex_increment[] is allocated with size at least number of facet vertices
  template <typename DTYPE, typename ITYPE>
  void compute_facet_vertex_increment
  (const DTYPE dimension, const DTYPE ifacet,
   const ITYPE * cube_vertex_increment, ITYPE * facet_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_facet_vertex_increment");

    if (dimension <= 0) { return; };

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated and set before calling compute_facet_vertex_increment.");
      throw error;
    }

    if (facet_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array facet_vertex_increment[] must be allocated before calling compute_facet_vertex_increment.");
      throw error;
    }
    
    DTYPE orth_dir = ifacet%dimension;
    ITYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    ITYPE num_facet_vertices = compute_num_cube_facet_vertices(dimension);

    // Multiply by s mod (num_cube_vertices-1) to compute vertex index.
    ITYPE s = 2;
    for (ITYPE d = 0; d < orth_dir; d++)
      { s = s*2; };
    s = s%(num_cube_vertices-1);

    if (ifacet < dimension) {

      for (ITYPE i = 0; i < num_facet_vertices; i++) {
        ITYPE k = (i*s)%(num_cube_vertices-1);
        facet_vertex_increment[i] = cube_vertex_increment[k];
      }
    }
    else {

      // Translate by t mod num_cube_vertices to compute vertex index.
      ITYPE t = 1;
      for (ITYPE d = 0; d < orth_dir; d++)
        { t = t*2; };

      for (ITYPE i = 0; i < num_facet_vertices; i++) {
        ITYPE k = (i*s)%(num_cube_vertices-1);
        facet_vertex_increment[i] = cube_vertex_increment[k+t];
      }

      // Swap subfacets to get consistent facet orientation.
      if (num_facet_vertices > 1) {
        ITYPE num_subfacet_vertices = num_facet_vertices/2;
        for (ITYPE i = 0; i < num_subfacet_vertices; i++) {
          std::swap(facet_vertex_increment[i], 
                    facet_vertex_increment[i+num_subfacet_vertices]);
        }
      }
    }
  }

  /// Compute increment to add to vertex 0 to compute vertex i of ridge.
  /// @param dimension  Dimension of grid.
  /// @param cube_vertex_increment[] = Cube vertex increment. iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube with primary vertex iv0.
  /// @param[out] ridge_vertex_increment[] = Ridge vertex increment. iv0+ridge_vertex_increment[i] is the i'th vertex of the ridge with primary vertex iv0.
  /// @pre Array ridge_vertex_increment[] is allocated with size at least number of ridge vertices
  template <typename DTYPE, typename ITYPE>
  void compute_ridge_vertex_increment
  (const DTYPE dimension, const DTYPE orth_dir0, const DTYPE orth_dir1,
   const ITYPE * cube_vertex_increment, ITYPE * ridge_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_ridge_vertex_increment");

    if (dimension <= 0) { return; };

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated and set before calling compute_ridge_vertex_increment.");
      throw error;
    }

    if (ridge_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array ridge_vertex_increment[] must be allocated before calling compute_ridge_vertex_increment.");
      throw error;
    }
    
    ITYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    ITYPE num_cube_ridge_vertices = compute_num_cube_ridge_vertices(dimension);
    ITYPE mask0, mask1;

    mask0 = (ITYPE(1) << orth_dir0);
    mask1 = (ITYPE(1) << orth_dir1);

    ITYPE i = 0;
    for (ITYPE j = 0; j < num_cube_vertices; j++) {

      if (((j & mask0) == 0) && ((j & mask1) == 0)) {
        ridge_vertex_increment[i] = cube_vertex_increment[j];
        i++;
      }
    }

    if (i != num_cube_ridge_vertices) {
      error.AddMessage("Programming error.  Added ", i,
                       " values to ridge_vertex_increment[].");
      error.AddMessage("  Should have added ",
                       num_cube_ridge_vertices, " values.");
      throw error;
    }
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: GETTING VERTICES
  // **************************************************

  /// Subsample vertices in subgrid.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param  subgrid_origin  Subgrid origin.
  /// @param subgrid_axis_size 
  ///        Array: <em>subgrid_axis_size[d]</em> = Number of vertices along subgrid axis d.
  /// @param subsample_period 
  ///             Array: <em>subsample_period[d]</em> = Only report every k'th vertex along subgrid axis \a d where k = \a subsample_period[d].
  /// @param[out] vlist[]  List of vertices.
  /// @pre \li Subgrid is contained in grid, i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li \a subsample_period[d] is a positive integer.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE>
  void subsample_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("subsample_subgrid_vertices");

    GRID_SIZE_TYPE num_vertices;
    compute_subsample_size
      (dimension, subgrid_axis_size, subsample_period, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: subgrid_axis_size[d] >= 1 for all d

    if (!check_vertex_list(vlist, error)) { throw error; };

    IJK::ARRAY<VTYPE> subsample_increment(dimension);
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment.Ptr());

    // initialize prev_num_vertices
    vlist[0] = subgrid_origin;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_subsampled_along_axis =
        compute_num_subsampled_vertices_along_axis
        (subgrid_axis_size[d], subsample_period[d]);

      VTYPE iv0 = subsample_increment[d];

      for (VTYPE i = 1; i < num_subsampled_along_axis; i++) {
        for (VTYPE * vprev_ptr = vlist; 
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        };
        iv0 = iv0 + subsample_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_subsampled_along_axis;
    }

    if (!check_num_vertices_added(prev_num_vertices, num_vertices, error))
      throw error;

  }

  /// Subsample vertices in subgrid.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE, typename NTYPE>
  void subsample_subgrid_vertices
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices
      (grid.Dimension(), grid.AxisSize(), subgrid_origin, subgrid_axis_size,
       subsample_period, vlist);
  }

  /// Subsample vertices in grid.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param subsample_period Array: <em>subsample_period[d]</em> = Only report every k'th vertex along subgrid axis \a d where k = \a subsample_period[d].
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li \a subsample_period is a positive integer.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE>
  void subsample_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices(dimension, axis_size, VTYPE(0), axis_size, 
                               subsample_period, vlist);
  }

  /// Get vertices in subgrid.
  /// Returns vertex indices sorted in increasing order.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = 
  ///        Number of vertices along axis \a d.
  /// @param axis_increment[] = Axis increment. iv+axis_increment[i] 
  ///        is the next vertex after vertex iv along axis i.
  /// @param  subgrid_origin  Subgrid origin.
  /// @param subgrid_axis_size
  ///        Array: <em>subgrid_axis_size[d]</em> = Number of vertices along subgrid axis d.
  /// @param[out] vlist[]  List of vertices.
  /// @pre \li Subgrid is contained in grid, i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  /// @pre dimension >= 1.
  template <typename DTYPE, typename ATYPE, typename ITYPE, typename VTYPE,
            typename NTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ITYPE * axis_increment,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist, NTYPE & num_subgrid_vertices)
  {
    IJK::PROCEDURE_ERROR error("get_subgrid_vertices");

    compute_num_grid_vertices
      (dimension, subgrid_axis_size, num_subgrid_vertices);
    if (num_subgrid_vertices < 1) { return; };
    // Note: subgrid_axis_size[d] >= 1 for all d

    if (!check_vertex_list(vlist, error)) { throw error; };

    // add vertices along dimension 0 to vlist.
    vlist[0] = subgrid_origin;
    for (VTYPE i = 1; i < subgrid_axis_size[0]; i++) 
      { vlist[i] = subgrid_origin + i; }

    // initialize prev_num_vertices
    VTYPE prev_num_vertices = subgrid_axis_size[0];

    // Process axes 1,2,...,dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 1; d < dimension; d++) {

      VTYPE iv0 = axis_increment[d];

      for (VTYPE i = 1; i < subgrid_axis_size[d]; i++) {
        for (VTYPE * vprev_ptr = vlist;
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        }
        iv0 = iv0+axis_increment[d];
      }

      prev_num_vertices = prev_num_vertices*subgrid_axis_size[d];
    }

    if (!check_num_vertices_added
        (prev_num_vertices, num_subgrid_vertices, error))
      throw error;
  }

  /// Get vertices in subgrid.
  /// Returns vertex indices sorted in increasing order.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = 
  ///        Number of vertices along axis \a d.
  /// @param  subgrid_origin  Subgrid origin.
  /// @param subgrid_axis_size
  ///        Array: <em>subgrid_axis_size[d]</em> = 
  ///        Number of vertices along subgrid axis d.
  /// @param[out] vlist[]  List of vertices.
  /// @pre \li Subgrid is contained in grid, 
  ///     i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li Array vlist[] is preallocated to length at least 
  ///          number of vertices in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist)
  {
    VTYPE num_subgrid_vertices;

    IJK::ARRAY<VTYPE> axis_increment(dimension);
    compute_axis_increment(dimension, axis_size, axis_increment.Ptr());

    get_subgrid_vertices(dimension, axis_size, axis_increment.Ptr(),
                         subgrid_origin, subgrid_axis_size, 
                         vlist, num_subgrid_vertices);
  }

  /// Get cubes in subgrid.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param  subgrid_origin  Subgrid origin.
  /// @param subgrid_axis_size
  ///        Array: <em>subgrid_axis_size[d]</em> = Number of vertices along subgrid axis d.
  /// @param[out] vlist[]  List of vertices.
  /// @pre \li Subgrid is contained in grid, i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li Array cube_list[] is preallocated to length at least number of cubes in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_subgrid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * cube_list)
  {
    IJK::ARRAY<ATYPE> subgrid2_axis_size(dimension);

    for (DTYPE d = 0; d < dimension; d++) {
      if (subgrid_axis_size[d] < 2) { return; }
      subgrid2_axis_size[d] = subgrid_axis_size[d]-1;
    }

    get_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid2_axis_size.PtrConst(), 
       cube_list);
  }

  /// Get vertices in grid.
  /// @param dimension = Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   VTYPE * vlist)
  {
    subsample_grid_vertices(dimension, axis_size, 1, vlist);
  }

  /// Get grid vertices in region between two grid vertices (inclusive).
  /// @param dimension = Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param iv0 = Lower grid vertex.
  /// @param iv1 = Upper grid vertex.
  /// @param[out] vlist[] = List of vertices between \a iv0 and \a iv1.
  /// @pre \li 0 <= iv0 <= iv1 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices_between
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const VTYPE iv1, VTYPE * vlist)
  {
    ATYPE coord0[dimension];
    ATYPE coord1[dimension];
    ATYPE region_size[dimension];
    IJK::PROCEDURE_ERROR error("get_grid_vertices_between");

    if (dimension < 0) { return; };
    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    for (DTYPE d = 0; d < dimension; ++d) {
      if (!check_region_coordinates
          (dimension, iv0, coord0, iv1, coord1, error)) { throw error; };

      region_size[d] = coord1[d]-coord0[d]+1;
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }

  /// Get grid vertices in region.
  /// @param dimension = Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param iv0 = Primary vertex of region.
  /// @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
  /// @param[out] vlist = List of vertices in region.
  /// @pre \li 0 <= iv0 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist)
  {
    ATYPE coord[dimension];
    ATYPE region_size[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
        { region_size[d] = max_region_edge_length + 1; }
      else if (coord[d] < axis_size[d])
        { region_size[d] = axis_size[d] - coord[d]; }
      else {
        // Vertex iv0 does not lie inside grid.
        return;
      }
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }

  /// Get grid cubes in region.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param iv0  Primary vertex of region.
  /// @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
  /// @param[out] vlist[] List of primary vertices of cubes in region.
  /// @param num_cubes Number of cubes in region (= number of vertices in vlist[].)
  /// @pre \li 0 <= iv0 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void get_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist,
   NTYPE & num_cubes)
  {
    IJK::ARRAY<ATYPE> coord(dimension);
    IJK::ARRAY<ATYPE> subgrid_size(dimension);

    num_cubes = 0;
    if (max_region_edge_length < 1) { return; };

    compute_coord(iv0, dimension, axis_size, coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
        { subgrid_size[d] = max_region_edge_length; }
      else if (coord[d]+1 < axis_size[d])
        { subgrid_size[d] = axis_size[d] - coord[d] - 1; }
      else {
        // Region contains no cubes
        return;
      }
    }

    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) 
      { num_cubes *= subgrid_size[d]; }

    get_subgrid_vertices
      (dimension, axis_size, iv0, subgrid_size.PtrConst(), vlist);
  }

  /// Get grid vertices in neighborhood around vertex \a iv.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param axis_increment[] = Axis increment. iv+axis_increment[i] 
  ///        is the next vertex after vertex iv along axis i.
  /// @param iv  Neighborhood around vertex \a iv.
  /// @param distance Distance to \a iv.
  /// @param[out] vlist[] List of vertices in neighborhood around \a iv.
  /// @pre dimension > 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE0, typename VTYPE1,
            typename DIST_TYPE>
  void get_grid_vertices_in_neighborhood
  (const DTYPE dimension, const ATYPE * axis_size, const ATYPE * axis_increment,
   const VTYPE0 iv, const DIST_TYPE distance,
   std::vector<VTYPE1> & vlist)
  {
    VTYPE0 region_iv0;
    IJK::ARRAY<ATYPE> region_iv0_coord(dimension);
    IJK::ARRAY<ATYPE> region_axis_size(dimension);
    IJK::ARRAY<ATYPE> vertex_coord(dimension);
    VTYPE0 num_subgrid_vertices;

    compute_coord(iv, dimension, axis_size, vertex_coord.Ptr());

    compute_region_around_vertex_coord
      (vertex_coord.PtrConst(), dimension, axis_size, 
       distance, region_iv0_coord.Ptr(), region_axis_size.Ptr());
    region_iv0 =
      compute_vertex_index<VTYPE0>
      (region_iv0_coord.PtrConst(), dimension, axis_size);

    compute_num_grid_vertices
      (dimension, region_axis_size.PtrConst(), num_subgrid_vertices);

    if (num_subgrid_vertices == 0) {
      vlist.clear();
      return;
    }

    // get all vertices including iv, then remove iv.
    vlist.resize(num_subgrid_vertices);
    get_subgrid_vertices(dimension, axis_size, axis_increment,
                         region_iv0, region_axis_size.PtrConst(),
                         &(vlist[0]), num_subgrid_vertices);

    // remove iv from list
    std::remove(vlist.begin(), vlist.end(), iv);
    vlist.pop_back();
  }

  /// Get primary cube vertices in subgrid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_primary_cube_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size,
   VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size2(dimension);

    for (DTYPE d = 0; d < dimension; d++) {
      if (subgrid_axis_size[d] < 2) 
        { return; }                      // zero cubes
      subgrid_axis_size2[d] = subgrid_axis_size[d]-1;
    }
    get_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size2.PtrConst(), 
       vlist);
  }

  // ********************************************************
  // TEMPLATE FUNCTIONS: FACET VERTICES, CUBES AND REGIONS
  // ********************************************************

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param orth_dir  Direction orthogonal to the facet.
  /// @param boundary_width  Width of boundary, (Number of vertices. Must be non-negative.)
  /// @param[out] num_vertices Number of vertices.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE2 orth_dir, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    if (dimension < 1) { num_vertices = 0; }
    else if (dimension == 1 && boundary_width > 0)
      { num_vertices = 0; }
    else { num_vertices = 1; }

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
        if (axis_size[d] > 2*boundary_width) 
          { num_vertices *= (axis_size[d]-2*boundary_width); }
        else
          { num_vertices = 0; };
      }
      else if (axis_size[d] < 1)
        { num_vertices = 0; }
    }
  }

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, 0, num_vertices);
  }

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename VTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, const DTYPE2 orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, num_vertices);
  }

  /// Return number of vertices in grid facet interior.
  /// Same as compute_num_vertices_in_grid_facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE2 orth_dir, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, boundary_width, num_vertices);
  }

  /// Return number of vertices in specified grid ridge
  /// Specify grid facet by the directions orthogonal to ridge
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_vertices_in_grid_ridge
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir0, const DTYPE orth_dir1,
   NTYPE & num_vertices)
  {
    FACET_LIST2<DTYPE> facet_list2(orth_dir0, orth_dir1);
    compute_num_vertices(dimension, axis_size, facet_list2, num_vertices);
  }

  /// Return maximum number of vertices over all grid facets
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_max_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   NTYPE & max_num_vertices)
  {
    max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      NTYPE num_face_vertices; 
      compute_num_vertices_in_grid_facet
        (dimension, axis_size, d, num_face_vertices);
      if (num_face_vertices > max_num_vertices)
        max_num_vertices = num_face_vertices;
    };
  }

  /// Return maximum number of vertices over all grid facet interiors.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_max_num_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size, const WTYPE boundary_width,
   NTYPE & max_num_vertices)
  {
    max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      NTYPE num_face_vertices; 
      compute_num_vertices_in_grid_facet_interior
        (dimension, axis_size, d, boundary_width, num_face_vertices);
      if (num_face_vertices > max_num_vertices)
        max_num_vertices = num_face_vertices;
    };
  }


  /// Return number of cubes in specified grid facet.
  /// Facet is lower facet orthogonal to the specified direction
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir,
   NTYPE & num_cubes)
  {
    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) { 
        num_cubes = 0;
        return; 
      };
      if (d != orth_dir) {
        num_cubes = num_cubes*(axis_size[d]-1);
      };
    }
  }

  /// Return number of cubes in grid facet orthogonal to axis 0
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    if (dimension < 1) 
      { num_cubes = 0; }
    else {
      compute_num_cubes_in_grid_facet
        (dimension, axis_size, DTYPE(0), num_cubes);
    }
  }

  /// Get vertices in specified grid facet.
  /// Does not return any vertices if \a axis_size[orth_dir] == 0.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param orth_dir = Direction orthogonal to facet.
  /// @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
  /// @param boundary_width = Width of boundary, (Number of vertices. Must be non-negative.)
  /// @param[out] vlist[] = List of primary vertices.
  /// @pre Array vlist[] must be pre-allocated to size at least number of vertices in facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, const BTYPE boundary_width, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::ARRAY<ATYPE> coord(dimension);

    if (dimension < 1) { return; };
    if (axis_size[orth_dir] < 1) { return; };
    if (dimension == 1 && boundary_width > 0) { return; };
    VTYPE subgrid_origin = 0;

    if (boundary_width == 0) {
      std::copy(axis_size, axis_size+dimension, 
                subgrid_axis_size.Ptr());
    
      if (side) {
        for (DTYPE d = 0; d < dimension; d++) { coord[d] = 0; };
        coord[orth_dir] = axis_size[orth_dir]-1;
        subgrid_origin = 
          compute_vertex_index<VTYPE>(coord.PtrConst(), dimension, axis_size);
      }
    }
    else {
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != orth_dir) {
          if (axis_size[d] < 2*boundary_width) { return; };

          coord[d] = boundary_width; 
          subgrid_axis_size[d] = axis_size[d]-2*boundary_width;
        }
      }

      if (side) { coord[orth_dir] = axis_size[orth_dir]-1; }
      else { coord[orth_dir] = 0; }

      subgrid_origin = 
        compute_vertex_index<VTYPE>(coord.PtrConst(), dimension, axis_size);
    }
    subgrid_axis_size[orth_dir] = 1;

    get_subgrid_vertices(dimension, axis_size, subgrid_origin, 
                         subgrid_axis_size.PtrConst(), vlist);
  }

  /// Get vertices in specified grid facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, side, 0, vlist);
  }

  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, 
   const BTYPE boundary_width, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, DTYPE(0), false, boundary_width, vlist);
  }

  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0(dimension, axis_size, 0, vlist);
  }

  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void get_vertices_in_grid_facet0
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0
      (grid.Dimension(), grid.AxisSize(), vlist);
  }

  /// Get vertices in grid facet interior.
  /// Same as get_vertices_in_grid_facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, const BTYPE boundary_width, VTYPE * vlist)
  {
    get_vertices_in_grid_facet(dimension, axis_size, orth_dir, side,
                               boundary_width, vlist);
  }

  /// Get vertices in specified grid ridge
  /// Ridge is lower ridge orthogonal to the specified directions
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param orth_dir0 is first direction orthogonal to the ridge.
  /// @param orth_dir1 is second direction orthogonal to the ridge.
  /// @param[out] vlist[] = List of vertices in grid ridge.
  /// @pre Array vlist[] is preallocated to length at least number of vertices in grid ridge.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_ridge
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir0, const DTYPE orth_dir1, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    if (axis_size[orth_dir0] < 1 || axis_size[orth_dir1] < 1) { return; }
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size.Ptr());
    subgrid_axis_size[orth_dir0] = 1;
    subgrid_axis_size[orth_dir1] = 1;
    
    get_subgrid_vertices
      (dimension, axis_size, 0, subgrid_axis_size.PtrConst(), vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in specified grid facet where d = \a dimension.
  /// Does not return any vertices if \a axis_size[orth_dir] == 0.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param orth_dir = Direction orthogonal to facet.
  /// @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
  /// @param[out] vlist[] = List of primary vertices.
  /// @pre Array vlist[] must be pre-allocated to size at least number of primary vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, const bool side, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    if (axis_size[orth_dir] < 2) { return; };
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size.Ptr());
    subgrid_axis_size[orth_dir] = 2;

    VTYPE subgrid_origin = 0;

    if (side) {
      IJK::ARRAY<VTYPE> axis_increment(dimension);
      compute_increment(dimension, axis_size, axis_increment.Ptr());

      subgrid_origin = axis_increment[orth_dir]*(axis_size[orth_dir]-2);
    }

    get_primary_cube_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size.PtrConst(), 
       vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in specified grid facet where d = \a dimension.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param orth_dir  Direction orthogonal to facet.  Facet contains (0,0,...,0).
  /// @param[out] vlist[] = List of primary vertices.
  /// @pre Array vlist[] must be pre-allocated to size at least number of primary vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, orth_dir, false, vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in grid facet 0 where d = \a dimension.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, DTYPE(0), vlist);
  }

  /// Get outer grid vertices
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_outer_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
    // Precondition: vlist[] is preallocated to size 
    //   at least num_outer_vertices
  {
    VTYPE axis_increment[dimension];

    if (dimension < 1) { return; };

    DTYPE d_last = dimension - 1;

    if (dimension == 1) {
      if (axis_size[0] < 2) { return; }
      else {
        vlist[0] = axis_size[0]-1;
        return;
      }
    }
    else {
      if (axis_size[d_last] < 2) { return; }

      get_outer_grid_vertices(d_last, axis_size, vlist);
      VTYPE num_vertices; 
      compute_num_outer_vertices(d_last, axis_size, num_vertices);

      compute_increment(dimension, axis_size, axis_increment);

      for (VTYPE i = 1; i < axis_size[d_last]-1; i++) {
        VTYPE k = i*num_vertices;
        VTYPE k_increment = i*axis_increment[d_last];
        for (VTYPE j = 0; j < num_vertices; j++) {
          vlist[k+j] = vlist[j]+k_increment;
        } 
      }

      num_vertices = num_vertices * (axis_size[d_last]-1);

      VTYPE * vlist2 = vlist + num_vertices;
      FACET_LIST1<DTYPE> facet_list1(d_last);

      get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist2);
    }
  }

  /// Get primary vertices of regions in grid or subgrid.
  /// @param dimension = Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of full regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of full regions.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);
    
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < 2) { return; };

      subgrid_axis_size[d] = axis_size[d]-1;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }

  /// Get primary vertices of regions in grid or subgrid.
  /// @param dimension  Grid dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length Number of grid edges contained in each region edge.
  /// @param[out] vlist[] Array of primary vertices of full regions.
  /// @param[out] is_full Array: <em>is_full[k]</em> = True if region \a k is a full \a LxLx...xL region where \a L = \a region_edge_length.
  /// @pre Array vlist[] must be pre-allocated to size at least number of regions.
  /// @pre Array is_full[] must be pre-allocated to size at least number of regions.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist, bool * is_full)
  {
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    GRID_SIZE_TYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, num_regions);

    if (num_regions < 1) { return; };
    // Note: axis_size[d] >= 2 for all d

    GRID_SIZE_TYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, num_full_regions);

    if (!check_vertex_list(vlist, error)) { throw error; };

    IJK::ARRAY<VTYPE> subsample_increment(dimension);
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment.Ptr());

    // set vlist[0], is_full[0] and initialize prev_num_vertices
    vlist[0] = 0;
    if (num_full_regions > 0) { is_full[0] = true; }
    else { is_full[0] = false; };
    VTYPE prev_num_regions = 1;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_regions_along_axis =
        compute_num_regions_along_axis(axis_size[d], region_edge_length);

      ATYPE num_full_along_axis =
        compute_num_full_regions_along_axis(axis_size[d], region_edge_length);

      VTYPE iv0 = subsample_increment[d];
      for (VTYPE i = 1; i < num_full_along_axis; i++) {
        VTYPE i2 = i * prev_num_regions;
        for (VTYPE j = 0; j < prev_num_regions; j++) {
          VTYPE k = j + i2;
          vlist[k] = vlist[j] + iv0;
          is_full[k] = is_full[j];
        }
        iv0 = iv0 + subsample_increment[d];
      }

      if (num_regions_along_axis != num_full_along_axis &&
          num_regions_along_axis > 1) {
        VTYPE i2 = num_full_along_axis * prev_num_regions;
        for (VTYPE j = 0; j < prev_num_regions; j++) {
          VTYPE k = j + i2;
          vlist[k] = vlist[j] + iv0;
          is_full[k] = false;
        }
      }

      prev_num_regions = prev_num_regions*num_regions_along_axis;
    }

    if (!check_num_vertices_added(prev_num_regions, num_regions, error))
      throw error;

  }

  /// Get primary vertices of full regions in grid.
  /// @param dimension = Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of full regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of full regions.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_full_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    IJK::PROCEDURE_ERROR error("get_full_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= region_edge_length) { return; };

      subgrid_axis_size[d] = axis_size[d]-region_edge_length;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }

  /// Get primary vertices of partial regions in grid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of partial regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of partial regions.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_partial_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("get_partial_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    GRID_SIZE_TYPE num_vertices;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: axis_size[d] >= 1 for all d not in facet_list

    if (!check_vertex_list(vlist, error)) { throw error; };

    VTYPE region_axis_increment[dimension];

    compute_subsample_increment(dimension, axis_size, subsample_period,
                                region_axis_increment);

    // initialize prev_num_vertices
    VTYPE prev_num_vertices = 0;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {
      VTYPE * vcur_ptr = vlist + prev_num_vertices;

      ATYPE num_full_regions_along_axis = 
        compute_num_full_regions_along_axis
        (axis_size[d], region_edge_length);

      VTYPE iv0 = region_axis_increment[d];
      for (VTYPE i = 1; i < num_full_regions_along_axis; i++) {
        for (VTYPE * vprev_ptr = vlist; 
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        };
        iv0 = iv0 + region_axis_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_full_regions_along_axis;

      ATYPE num_partial_regions_along_axis =
        compute_num_partial_regions_along_axis
        (axis_size[d], region_edge_length);

      if (num_partial_regions_along_axis > 0) {

        VTYPE inc = 
          region_axis_increment[d]*num_full_regions_along_axis;

        if (d == 0) {
          vlist[prev_num_vertices] = inc;
          prev_num_vertices++;
        }
        else {

          get_region_primary_vertices(d, axis_size, region_edge_length, 
                                      vlist + prev_num_vertices);

          ATYPE k;
          compute_num_regions(d, axis_size, region_edge_length, k);
          for (VTYPE i = prev_num_vertices; i < prev_num_vertices+k; i++)
            { vlist[i] += inc; }

          prev_num_vertices += k;
        }
      }
    }

    if (prev_num_vertices != num_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
                       " vertices to vertex list.");
      error.AddMessage("Number of vertices in list should be ", 
                       num_vertices, ".");
      throw error;
    }

  }

  /// Get vertices on boundary of grid.
  /// Allows boundary_width to be greater than 1.
  /// @param boundary_width Width of boundary, in vertices.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename WTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, VTYPE * vlist)
  {
    if (dimension < 1) { return; }
    if (boundary_width < 1) { return; }

    const DTYPE d_last = dimension - 1;

    if (axis_size[d_last] <= 2*boundary_width) {
      // all vertices are on the boundary
      VTYPE num_grid_vertices;
      compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
      for (VTYPE j = 0; j < num_grid_vertices; j++) 
        { vlist[j] = j; }
      return;
    }

    if (dimension == 1) {
      for (VTYPE j = 0; j < boundary_width; j++) 
        { vlist[j] = j; }

      for (VTYPE j = 0; j < boundary_width; j++) {
        VTYPE iv = axis_size[0]-boundary_width + j;
        vlist[j+boundary_width] = iv;
      };
      return;
    }

    ATYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // get vertices in lower facet
    VTYPE num_vertices_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d_last, num_vertices_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d_last, false, vlist);

    // get remaining vertices in lower boundary
    for (VTYPE i = 1; i < boundary_width; i++) {
      VTYPE inc = i*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
        vlist[i*num_vertices_in_grid_facet + j] = vlist[j] + inc;
      }
    }

    VTYPE * vlist2 = vlist+boundary_width*num_vertices_in_grid_facet;
    get_boundary_grid_vertices(dimension-1, axis_size, 
                               boundary_width, vlist2);

    VTYPE num_boundary_grid_vertices;
    compute_num_boundary_grid_vertices
      (dimension-1, axis_size, boundary_width, num_boundary_grid_vertices);
    for (VTYPE * vcur_ptr = vlist2; 
         vcur_ptr != vlist2+num_boundary_grid_vertices; vcur_ptr++)
      { *vcur_ptr += boundary_width*axis_increment[d_last]; }

    VTYPE * vlist3 = vlist2+num_boundary_grid_vertices;
    for (ATYPE j = boundary_width+1; j+boundary_width < axis_size[d_last]; 
         j++) {
      VTYPE inc = axis_increment[d_last]*(j-boundary_width);

      for (VTYPE i = 0; i < num_boundary_grid_vertices; i++)  
        { vlist3[i] = vlist2[i] + inc; }

      vlist3 += num_boundary_grid_vertices;
    }

    VTYPE * vlist4 = vlist3 + (boundary_width-1)*num_vertices_in_grid_facet;

    // get vertices in upper facet
    get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist4);

    // get remaining vertices in upper boundary
    for (VTYPE i = 0; i+1 < boundary_width; i++) {
      VTYPE inc = (boundary_width-i-1)*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
        vlist3[i*num_vertices_in_grid_facet + j] = vlist4[j] - inc;
      }
    }
  
  }

  /// Get vertices on boundary of grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_boundary_grid_vertices(dimension, axis_size, 1, vlist);
  }

  /// Get boundary grid cubes
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_boundary_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * cube_list)
  {
    if (dimension < 1) { return; }
    if (dimension == 1) {
      if (axis_size[0] > 0) { cube_list[0] = 0; }
      if (axis_size[0] > 2) { cube_list[1] = axis_size[0]-2; };
      return;
    }

    DTYPE d_last = dimension - 1;
    if (axis_size[d_last] < 2) { return; };

    // get vertices in lower facet
    VTYPE num_cubes_in_grid_facet;
    compute_num_cubes_in_grid_facet(dimension, axis_size, d_last,
                                    num_cubes_in_grid_facet);
    get_cubes_in_grid_facet(dimension, axis_size, d_last, false, cube_list);

    VTYPE * cube_list2 = cube_list+num_cubes_in_grid_facet;
    VTYPE * cube_list3 = cube_list2;
    if (axis_size[d_last] > 3) {
      ATYPE axis_increment[dimension];

      compute_increment(dimension, axis_size, axis_increment);
      get_boundary_grid_cubes(dimension-1, axis_size, cube_list2);

      VTYPE n;
      compute_num_boundary_grid_cubes(dimension-1, axis_size, n);
      for (VTYPE * vcur_ptr = cube_list2; vcur_ptr != cube_list2+n; vcur_ptr++)
        { *vcur_ptr += axis_increment[d_last]; }

      cube_list3 = cube_list2 + n;
      for (ATYPE j = 2; j+2 < axis_size[d_last]; j++) {
        VTYPE inc = axis_increment[d_last]*(j-1);

        for (VTYPE i = 0; i < n; i++)  { cube_list3[i] = cube_list2[i] + inc; }

        cube_list3 += n;
      }
    }

    get_cubes_in_grid_facet(dimension, axis_size, d_last, true, cube_list3);
  }

  // ********************************************************
  // TEMPLATE FUNCTIONS: GRID BOUNDARIES
  // ********************************************************

  /// Return true if cube_facet is on grid boundary.
  template <typename DTYPE, typename ATYPE, 
            typename DTYPE2, typename CTYPE>
  bool is_cube_facet_on_grid_boundary
  (const DTYPE dimension, const ATYPE * axis_size,
   const CTYPE cube_index, const DTYPE2 facet_orth_dir,
   const bool facet_side)
  {
    long boundary_bits;
    compute_boundary_cube_bits
      (cube_index, dimension, axis_size, boundary_bits);

    // Convert kf to index into boundary bits.
    long bit_index;
    IJK::compute_boundary_bit_index(facet_orth_dir, facet_side, bit_index);
    long mask = (long(1) << bit_index);
    if ((boundary_bits & mask) == 0) 
      { return(false); }
    else 
      { return(true); }
  }

  // ********************************************************
  // TEMPLATE FUNCTIONS: COMPUTING NEIGHBORS
  // ********************************************************

  /// \brief Compute number of neighbors of a vertex in all cubes containing the vertex.
  /// Does not count the vertex itself.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_vertex_neighborsC
  (const DTYPE dimension, NTYPE & num_neighbors)
  { 
    num_neighbors = 1;
    for (DTYPE d = 0; d < dimension; d++) 
      { num_neighbors = num_neighbors*3; }
    num_neighbors = num_neighbors-1;
  }

  /// \brief Compute number of vertices which share an edge with a vertex.
  /// Does not count the vertex itself.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_vertex_neighborsE
  (const DTYPE dimension, NTYPE & num_neighbors)
  { 
    num_neighbors = 2*dimension;
  }

  /// \brief Compute number of vertices in cubes containing a facet,
  ///        not including facet vertices.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_facet_neighborsC
  (const DTYPE dimension, NTYPE & num_neighbors)
  {
    NTYPE num_facet_vertices = compute_num_cube_facet_vertices(dimension);
    num_neighbors = 2*num_facet_vertices;
  }

  /// \brief Compute number of vertices in 2-faces containing an edge,
  ///        not including edge vertices.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_edge_neighborsF2
  (const DTYPE dimension, NTYPE & num_neighbors)
  {
    num_neighbors = 0;
    if (dimension > 0) 
      { num_neighbors = 4*(dimension-1); }
  }

  /// \brief Compute integer to add to vertex index to compute vertex neighbors.
  ///        Use only for vertex neighbors of internal vertices.
  /// @param[out] vertex_neighborC Array. 
  ///        iv + vertex_neighborC[k] = index of k'th vertex neighbor of iv
  ///        where iv is the index of an internal grid vertex.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_vertex_neighborC
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * vertex_neighborC)
  {
    IJK::PROCEDURE_ERROR error("compute_vertex_neighborC");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (dimension == 0) { return; }

    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    compute_increment(dimension, axis_size, axis_increment.Ptr());

    // iv0 = index of vertex (1,1,...,1).
    DIFFTYPE iv0 = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { iv0 += axis_increment[d]; };

    vertex_neighborC[0] = -iv0;
    DIFFTYPE k = 1;

    // Process axes 0,1,2,..., dimension-1
    DIFFTYPE * vcur_ptr = vertex_neighborC + k;
    DIFFTYPE * vcenter_ptr = vertex_neighborC;
    for (DTYPE d = 0; d < dimension; d++) {

      for (ATYPE j = 1; j < 3; j++) {
        for (DIFFTYPE * vprev_ptr = vertex_neighborC; 
             vprev_ptr != vertex_neighborC+k; vprev_ptr++) {
          DIFFTYPE iv = j*axis_increment[d] + *(vprev_ptr);
          if ((d+1 < dimension) || (j != 1) || 
              (vprev_ptr != vcenter_ptr)) { 
            *(vcur_ptr) = iv; 
            vcur_ptr++;
          }
        }
      }
      vcenter_ptr += k;

      k = vcur_ptr - vertex_neighborC;
    }

    DIFFTYPE num_neighbors; 
    compute_num_vertex_neighborsC(dimension, num_neighbors);

    if (!check_num_vertices_added(k, num_neighbors, error))
      throw error;
  }

  /// Compute integer to add to vertex index to compute vertex neighbors 
  ///   across edges.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_vertex_neighborE
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * vertex_neighborE)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_vertex_neighborE");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    DIFFTYPE num_neighbors; 
    compute_num_vertex_neighborsE(dimension, num_neighbors);

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    for (DIFFTYPE d = 0; d < dimension; d++) 
      { vertex_neighborE[d] = -axis_increment[d]; };

    for (DIFFTYPE d = 0; d < dimension; d++) 
      { vertex_neighborE[d + dimension] = axis_increment[d]; };
  }

  /// Compute integer to add to cube_index to compute cubes sharing
  ///   an edge with cube cube_index.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_cube_neighborE
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * cube_neighborE)
  {
    const DIFFTYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    const DIFFTYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    const DIFFTYPE num_cube_edges = compute_num_cube_edges(dimension);
    const DIFFTYPE facet_vlast = num_facet_vertices-1;
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);
    IJK::ARRAY<DIFFTYPE> facet_vertex_increment(num_facet_vertices);
    IJK::PROCEDURE_ERROR error("compute_cube_neighborE");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (num_cube_vertices <= 0) { return; }
    if (num_cube_edges <= 0) { return; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    facet_vertex_increment[0] = -cube_vertex_increment[0];

    DIFFTYPE k = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      compute_facet_vertex_increment
        (dimension, d, cube_vertex_increment.PtrConst(), 
         facet_vertex_increment.Ptr());

      for (int i = 0; i < num_facet_vertices; i++) {
        cube_neighborE[k] = 
          2*facet_vertex_increment[i] - facet_vertex_increment[facet_vlast];
        k++;
      }
    }

    if (k != num_cube_edges) {
      error.AddMessage
        ("Programming error.  Wrong number of cubes added to cube_neighborE.");
      throw error;
    }

  }

  /// Compute integer to add to cube_index to compute cubes sharing
  ///   a vertex with cube cube_index.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_cube_neighborV
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * cube_neighborV)
  {
    const DIFFTYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);
    IJK::PROCEDURE_ERROR error("compute_cube_neighborV");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (num_cube_vertices <= 0) { return; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    cube_neighborV[0] = 0;
    for (DTYPE d = 0; d < dimension; d++)
      { cube_neighborV[0] -= axis_increment[d]; }

    for (DIFFTYPE k = 0; k < num_cube_vertices; k++)
      { cube_neighborV[k] = cube_neighborV[0] + 2*cube_vertex_increment[k]; }
  }

  /// Compute integer to add to vertex index to compute facet neighbors.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_facet_neighborC
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * facet_neighborC)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_facet_neighborC");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    DIFFTYPE num_neighbors;
    compute_num_facet_neighborsC(dimension, num_neighbors);

    const DIFFTYPE num_cube_vertices = 
      compute_num_cube_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);

    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    const DIFFTYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> facet_vertex_increment(dimension*num_facet_vertices);
    
    // for each facet defined by a different orthogonal direction
    for (DTYPE orth_dir = 0; orth_dir < dimension; orth_dir++) {

      DIFFTYPE * facet_ptr = 
        facet_vertex_increment.Ptr()+num_facet_vertices*orth_dir;

      compute_facet_vertex_increment
        (dimension, orth_dir, cube_vertex_increment.PtrConst(), facet_ptr);
      
      // for each vertex in the prvious facet
      for (DIFFTYPE k = 0; k < num_facet_vertices; k++) {
        facet_neighborC[k+num_neighbors*orth_dir] =
          facet_ptr[k] - axis_increment[orth_dir];
      }

      // for each vertex in the next facet
      for (DIFFTYPE k = 0; k < num_facet_vertices; k++) {
        facet_neighborC[k+num_facet_vertices+num_neighbors*orth_dir] =
          facet_ptr[k] + axis_increment[orth_dir];
      }
    }
  }

  /// Compute integer to add to vertex index to compute edge neighbors.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_edge_neighborF2
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * edge_neighborF2)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_edge_neighborF2");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    DIFFTYPE num_neighbors;
    compute_num_edge_neighborsF2(dimension, num_neighbors);

    DIFFTYPE k = 0;
    // for each edge determined by a direction
    for (DTYPE dir = 0; dir < dimension; dir++) {

      // for each edge incident on the lower vertex
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = -axis_increment[d];
          k++;
        }
      }

      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[d];
          k++;
        }
      }

      // for each edge incident on the upper vertex
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[dir]-axis_increment[d];
          k++;
        }
      }

      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[dir]+axis_increment[d];
          k++;
        }
      }
    }

    if (k != num_neighbors*dimension) {
      error.AddMessage("Programming error.  Wrong number of edge neighbors added to edge_neighborF2[].");
      throw error;
    }

  }

  // **************************************************
  // TEMPLATE CLASS GRID MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::GRID
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    Init(dimension, axis_size);
  }

  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::GRID()
  {
    Init(0, NULL);
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID()
  {
    FreeAll();
  }

  /// Initialize grid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Init
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    this->axis_size = NULL;
    this->dimension = 0;
    this->num_vertices = 1;
    if (dimension > 0) 
      { SetSize(dimension, axis_size); };
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeAll();

    if (dimension < 0) {
      IJK::PROCEDURE_ERROR error("Grid::SetSize");
      error.AddMessage("Programming error.  Illegal dimension ",
                       dimension, ".");
      error.AddMessage("Dimension should be non-negative.");
      throw error;
    }

    this->dimension = dimension;
    this->axis_size = NULL;
    if (dimension > 0)
      { this->axis_size = new ATYPE[dimension]; }
    else {
      // allocate axis_size even if dimension equals 0
      this->axis_size = new ATYPE[1];
      this->axis_size[0] = 0;
    }
      
    for (DTYPE d = 0; d < dimension; d++)
      { this->axis_size[d] = axis_size[d]; }

    compute_num_grid_vertices(dimension, axis_size, this->num_vertices);
  }

  /// Set size of \a grid to size of \a grid2.
  /// @param grid2  Grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }


  /// Free all allocated memory in typename GRID.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::FreeAll()
  {
    if (axis_size != NULL) { delete [] axis_size; };
    axis_size = NULL;
    dimension = 0;
    num_vertices = 0;
  }

  /// Copy constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID(const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid)
  {
    Init(grid.Dimension(), grid.AxisSize());
  }

  /// Copy assignment.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::operator = 
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & right)
  {
    if (&right != this) {         // avoid self-assignment
      SetSize(right.Dimension(), right.AxisSize());
    }
  }

  /// Compute and return number of grid cubes.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumCubes() const
  {
    NTYPE num_grid_cubes;
    compute_num_grid_cubes(Dimension(), AxisSize(), num_grid_cubes);
    return(num_grid_cubes);
  }

  /// Compute and return number of grid edges
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumEdges() const
  {
    NTYPE num_grid_edges;
    compute_num_grid_edges(Dimension(), AxisSize(), num_grid_edges);
    return(num_grid_edges);
  }

  /// Compute and return number of cubes in grid interior.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumInteriorCubes() const
  {
    NTYPE num_interior_cubes;
    compute_num_interior_grid_cubes
      (Dimension(), AxisSize(), num_interior_cubes);
    return(num_interior_cubes);
  }

  /// Compute and return number of cubes in grid boundary.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumBoundaryCubes() const
  {
    NTYPE num_boundary_cubes;
    compute_num_boundary_grid_cubes
      (Dimension(), AxisSize(), num_boundary_cubes);
    return(num_boundary_cubes);
  }

  /// Compute and return number of vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumVerticesInFacet(const DTYPE2 orth_dir) const
  {
    NTYPE num_vertices_in_facet;
    compute_num_vertices_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, num_vertices_in_facet);
    return(num_vertices_in_facet);
  }

  /// Compute and return number of vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename WTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumVerticesInFacet
  (const DTYPE2 orth_dir, const WTYPE boundary_width) const
  {
    NTYPE num_vertices_in_facet;
    compute_num_vertices_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, boundary_width,
       num_vertices_in_facet);
    return(num_vertices_in_facet);
  }

  /// Compute and return number of cubes in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumCubesInFacet(const DTYPE2 orth_dir) const
  {
    NTYPE num_cubes_in_facet;
    compute_num_cubes_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, num_cubes_in_facet);
    return(num_cubes_in_facet);
  }

  /// Compute index of vertex with given coordinates.
  /// @param coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeVertexIndex(const GTYPE * coord) const
  {
    return(compute_vertex_index<VTYPE>(coord, Dimension(), AxisSize()));
  }

  /// Compute index of vertex with given coordinates.
  /// @param coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeVertexIndex(const std::vector<GTYPE> coord) const
  {
    return(compute_vertex_index<VTYPE>(&(coord[0]), Dimension(), AxisSize()));
  }

  /// Compute coordinates of given vertex.
  /// @param iv  Vertex index.
  /// @param[out] coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  /// @pre Array <em>coord[]</em> is preallocated to length 
  ///      at least \a Dimension().
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCoord(const VTYPE iv, GTYPE * coord) const
  {
    compute_coord(iv, Dimension(), AxisSize(), coord);
  }

  /// Compute coordinates of given cube center.
  ///   Cube center coord is (vertex coord) + (0.5,0.5,...,0.5).
  /// @param iv  Index of primary vertex of cube.
  /// @param[out] coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  /// @pre Array <em>coord[]</em> is preallocated to length 
  ///      at least \a Dimension().
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCubeCenterCoord(const VTYPE iv, GTYPE * coord) const
  {
    ComputeCoord(iv, coord);
    for (DTYPE d = 0; d < Dimension(); d++)
      { coord[d] += 0.5; }
  }

  /// Compute bits identifying which boundary contains vertex \a iv.
  /// @param iv  Vertex index.
  /// @param [out] boundary_bits Bits flagging boundaries containing 
  ///         vertex \a iv.
  ///       If bit \a 2d is true, then <em>d</em>'th coordinate of 
  ///         vertex \a iv is zero.
  ///       If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate 
  ///         of vertex \a iv equals <em>axis_size[d]-1</em>.
  /// @pre \li Variable \a boundary_bits has at least 
  ///          <em>(2*dimension)</em> bits.
  /// @pre \li AxisSize(d) > 0 for all \a d = 0,..., \a dimension-1.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeBoundaryBits(const VTYPE iv, BTYPE & boundary_bits) const
  {
    compute_boundary_bits(iv, Dimension(), AxisSize(), boundary_bits);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename PTYPE, typename ATYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeSubsampledAxisSizes
  (const PTYPE subsample_period, ATYPE2 subsampled_axis_size[]) const
  {
    compute_subsample_axis_sizes
      (Dimension(), AxisSize(), subsample_period, subsampled_axis_size);
  }

  /// Compute bits identifying which boundary contains cube \a icube.
  /// @param icube  Cube index.
  /// @param [out] boundary_bits Bits flagging boundaries containing cube \a icube.
  ///       If bit \a 2d is true, then <em>d</em>'th coordinate 
  ///              of cube \a icube is zero.
  ///       If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate 
  ///              of vertex \a iv equals <em>axis_size[d]-2</em>.
  /// @pre \li Variable \a boundary_bits has at least <em>(2*dimension)</em> bits.
  /// @pre \li AxisSize(d) > 0 for all \a d = 0,..., \a dimension-1.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeBoundaryCubeBits(const VTYPE icube, BTYPE & boundary_bits) const
  {
    compute_boundary_cube_bits(icube, Dimension(), AxisSize(), boundary_bits);
  }

  /// Return true if grid dimension and axis size match parameters.
  /// @param dimension  Dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const
  {
    if (dimension != this->Dimension()) { return(false); };
    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->AxisSize(d)) { return(false); };
    }

    return(true);
  }

  /// Return true if dimensions and axis size match dimensions and axis size of \a grid2.
  /// @param grid2  Grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2) const
  {
    return(this->CompareSize(grid2.Dimension(), grid2.AxisSize()));
  }

  /// Return true if grid contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsPoint(const CTYPE * coord) const
  {
    const DTYPE dimension = this->Dimension();

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] < 0 || coord[d]+1 > this->AxisSize(d))
        { return(false); }
    }

    return(true);
  }

  /// Return true if grid contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsPoint(const std::vector<CTYPE> & coord) const
  {
    return(this->ContainsPoint(&(coord[0])));
  }

  /// Return true if cube contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename VTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CubeContainsPoint(const VTYPE2 icube, const CTYPE * coord) const
  {
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<VTYPE> cube_coord(dimension);

    this->ComputeCoord(icube, cube_coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] < cube_coord[d] || coord[d] > cube_coord[d]+1)
        { return(false); }
    }

    return(true);
  }

  /// Return true if cube contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename VTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CubeContainsPoint
  (const VTYPE2 icube, const std::vector<CTYPE> & coord) const
  {
    return(this->CubeContainsPoint(icube, &(coord[0])));
  }

  /// Return true if grid contains specified region.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BOX_TYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsRegion(const BOX_TYPE & region) const
  {
    DTYPE dimension = this->Dimension();
    if (region.Dimension() < dimension)
      { dimension = region.Dimension(); };

    for (DTYPE d = 0; d < dimension; d++) {
      if (region.MinCoord(d) < 0 ||
          this->AxisSize(d) <= region.MaxCoord(d))
        { return(false); }
    }

    return(true);
  }

  /// Return true if cube facet is on grid boundary.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename DTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  IsCubeFacetOnGridBoundary
  (const CTYPE cube_index, const DTYPE2 facet_orth_dir, 
   const bool facet_side) const
  {
    return(is_cube_facet_on_grid_boundary
           (Dimension(), AxisSize(), cube_index, facet_orth_dir, facet_side));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckDimension(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
                 const char * grid_label1, const char * grid_label2,
                 IJK::ERROR & error) const
  {
    if (this->dimension != grid2.Dimension()) {
      error.AddMessage("Mismatch of volume dimensions.");
      error.AddMessage("  ", grid_label1, " has dimension ",
                       this->dimension, ".");
      error.AddMessage("  ", grid_label2, " has dimension ",
                       grid2.Dimension(), ".");
      return(false);
    }

    return(true);
  }

  /// Return true if grid dimension and axis size match parameters.
  /// @param dimension  Dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param[out] error Error message if grid dimension or axis_size do not match.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Check
  (const DTYPE dimension, const ATYPE * axis_size, IJK::ERROR & error) const
  {
    if (dimension != this->dimension) {
      error.AddMessage("Incorrect grid dimension ", this->dimension, ".");
      error.AddMessage("  Dimension should be ", dimension, ".");
      return(false);
    }

    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->axis_size[d]) {
        error.AddMessage("Illegal axis size[", d, "] = ", 
                         this->axis_size[d], ".");
        error.AddMessage("  Axis size[", d, "] should be ", axis_size[d], ".");
        return(false);
      }
    }

    NTYPE num_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_vertices);

    if (num_vertices != this->num_vertices) {
      error.AddMessage("Incorrect number of grid vertices ", 
                       this->num_vertices, ".");
      error.AddMessage("  Number of grid vertices should be ", 
                       num_vertices, ".");
      return(false);
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
        IJK::ERROR & error) const
  {
    return(Check(grid2.Dimension(), grid2.AxisSize(), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
        const char * grid_label1, const char * grid_label2,
        IJK::ERROR & error) const
  {
    if (!CheckDimension(grid2, grid_label1, grid_label2, error))
      { return(false); }

    for (int d = 0; d < dimension; d++) {
      if (this->AxisSize(d) != grid2.AxisSize(d)) {
        error.AddMessage("Mismatch of axis size ", d, ".");
        error.AddMessage("  ", grid_label1, " axis_size[", d,
                         "] = ", this->AxisSize(d), ".");
        error.AddMessage("  ", grid_label2, " axis_size[", d,
                         "] = ", grid2.AxisSize(d), ".");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCoord(const GTYPE * coord, IJK::ERROR & error) const
  {
    for (DTYPE d = 0; d < Dimension(); d++) {
      if (coord[d] < 0) {
        error.AddMessage("Coordinate should be non-negative.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        return(false);
      }
      else if (AxisSize(d) < 1) {
        error.AddMessage("Coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        error.AddMessage
          ("  axis size[", d, "] = ", AxisSize(d), 
           " so all coordinates on axis ", d, " are out of bounds.");
        return(false);

      }
      else if (coord[d] >= AxisSize(d)) {
        error.AddMessage("Coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], 
                         ".  axis size[", d, "] = ", AxisSize(d), ".");
        error.AddMessage("  coord[", d, 
                         "] should be less than axis size[", d, "].");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCoord(const std::vector<GTYPE> & coord, IJK::ERROR & error) const
  {
    return(this->CheckCoord(&(coord[0]), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeCoord(const GTYPE * coord, IJK::ERROR & error) const
  {
    for (DTYPE d = 0; d < Dimension(); d++) {
      if (coord[d] < 0) {
        error.AddMessage("Cube coordinates should be non-negative.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        return(false);
      }
      else if (AxisSize(d) < 2) {
        error.AddMessage("Cube coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        error.AddMessage
          ("  axis size[", d, "] = ", AxisSize(d), 
           " so all cube coordinates on axis ", d, " are out of bounds.");
        return(false);

      }
      else if (coord[d]+1 >= AxisSize(d)) {
        error.AddMessage("Cube coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], 
                         ".  axis size[", d, "] = ", AxisSize(d), ".");
        error.AddMessage("  coord[", d, 
                         "] should be less than (axis size[", d, "]-1).");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeCoord(const std::vector<GTYPE> & coord, IJK::ERROR & error) const
  {
    return(this->CheckCubeCoord(&(coord[0]), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename ITYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckVertexIndex(const ITYPE vertex_index, IJK::ERROR & error) const
  {
    if (vertex_index < 0) {
      error.AddMessage("Illegal vertex index ", vertex_index, ".");
      error.AddMessage("  Vertex index should be non-negative.");
      return(false);
    }

    if (vertex_index >= NumVertices()) {
      error.AddMessage("Illegal vertex index ", vertex_index, ".");
      error.AddMessage
        ("  Vertex index should be less than number of grid vertices.");
      error.AddMessage("  Number of grid vertices = ", 
                       NumVertices(), ".");
      return(false);
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename ITYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeIndex(const ITYPE cube_index, IJK::ERROR & error) const
  {
    if (cube_index < 0) {
      error.AddMessage("Illegal cube index ", cube_index, ".");
      error.AddMessage("  Cube index should be non-negative.");
      return(false);
    }

    for (DTYPE d = 0; d < Dimension(); d++) {
      if (AxisSize(d) <= 0) {
        // Never reached since NumVertices() == 0, but just in case.
        error.AddMessage("Illegal cube index ", cube_index, ".");
        if (AxisSize(d) <= 0) {
          error.AddMessage("  Grid axis size[", d, "] = ",
                           AxisSize(d), " so grid has no cubes.");
          return(false);
        }
      }
    }

    if (cube_index >= NumVertices()) {
      error.AddMessage("Illegal cube index ", cube_index, ".");
      error.AddMessage
        ("  Cube index should be less than number of grid vertices.");
      error.AddMessage("  Number of grid vertices = ", 
                       NumVertices(), ".");
      return(false);
    }

    VTYPE iv = cube_index;
    for (DTYPE d = 0; d < Dimension(); d++) {
      // Note: Already checked that AxisSize(d) > 0.
      VTYPE c = iv%AxisSize(d);
      iv = iv/AxisSize(d);
      if (c+1 >= AxisSize(d)) {
        error.AddMessage("  Vertices on right/top of cube with index ",
                         cube_index, " would have coordinate[", d,
                         "] = ", c+1, ".");
        error.AddMessage
          ("  Maximum coordinate[", d, "] in grid is ", AxisSize(d)-1, ".");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename VTYPE2, typename ATYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::CheckContainsRegion
  (const VTYPE2 region_v0, const ATYPE2 * region_axis_size,
   IJK::ERROR & error) const
  {
    const DTYPE dimension = this->Dimension();

    if (region_v0 < 0 || region_v0 >= this->NumVertices()) {
      error.AddMessage("Illegal lower/leftmost region vertex ",
                       region_v0, ".");
      return(false);
    }

    IJK::ARRAY<VTYPE> coord0(dimension);
    ComputeCoord(region_v0, coord0.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size[d] < 0) {
        error.AddMessage("Illegal region_axis_size[", d, 
                         "] = ", region_axis_size[d], ".");
        return(false); 
      }

      if (coord0[d] + region_axis_size[d] > this->AxisSize(d)) {
        error.AddMessage("Error.  Region extends beyond grid.");
        error.AddMessage("  lower/leftmost coord[", d, 
                         "] = ", coord0[d], ".");
        error.AddMessage("  region_axis_size[", d, "] = ",
                         region_axis_size[d], ".");
        error.AddMessage("  grid axis_size[", d, "] = ",
                         this->AxisSize(d), ".");
        return(false); 
      }
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS GRID_PLUS MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GRID_PLUS
  (const DTYPE dimension, const ATYPE * axis_size):
    GRID<DTYPE,ATYPE,VTYPE,NTYPE> (dimension,axis_size)
  {
    InitLocal();
  }

  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GRID_PLUS()
  {
    InitLocal();
  }

  /// Constructor from another GRID_PLUS grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID_PLUS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID<DTYPE, ATYPE, VTYPE, NTYPE>(grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID_PLUS(const GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> & grid2):
    GRID<DTYPE,ATYPE,VTYPE,NTYPE>(grid2)
  {
    InitLocal();
  }

  /// \brief Set all local (not inherited) arrays to NULL.
  /// Set all local variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::ZeroLocal()
  {
    this->axis_increment = NULL;
    this->cube_vertex_increment = NULL;
    this->facet_vertex_increment = NULL;
    this->unit_cube_coord = NULL;
    this->num_cube_vertices = 0;
    this->num_facet_vertices = 0;
    this->num_cube_ridge_vertices = 0;
    this->num_cube_facets = 0;
    this->num_cube_edges = 0;
  }

  /// \brief Initialize data structures in GRID_PLUS.
  /// @pre \a dimension and \a axis_size are already set.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::InitLocal()
  {
    ZeroLocal();
    if (this->Dimension() > 0) { Create(); };
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID_PLUS()
  {
    FreeLocal();
  }

  /// Allocate arrays and compute data in GRID_PLUS.
  /// @pre \a dimension and \a axis_size[] are already set.
  /// @pre All other arrays are set to NULL.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::Create()
  {
    const DTYPE dimension = this->Dimension();
    IJK::PROCEDURE_ERROR error("GRID_PLUS::Create");

    if (!check_dimension(dimension, error)) { throw error; };
    if (!IJK::check_is_NULL(axis_increment, "axis_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (cube_vertex_increment, "cube_vertex_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (facet_vertex_increment, "facet_vertex_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (unit_cube_coord, "unit_cube_coord", error))
      { throw error; }

    this->num_cube_vertices = compute_num_cube_vertices(dimension);
    this->num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    this->num_cube_ridge_vertices = 
      compute_num_cube_ridge_vertices(dimension);
    this->num_cube_facets = compute_num_cube_facets(dimension);
    this->num_cube_edges = compute_num_cube_edges(dimension);
    this->axis_increment = new VTYPE[dimension];
    this->cube_vertex_increment = new VTYPE[num_cube_vertices];
    this->facet_vertex_increment = 
      new VTYPE[num_facet_vertices*this->num_cube_facets];
    this->ridge_vertex_increment =
      new VTYPE[num_cube_ridge_vertices*dimension*dimension];
    this->unit_cube_coord = new NTYPE[num_cube_vertices*dimension];

    compute_increment
      (dimension, this->AxisSize(), this->axis_increment);
    compute_cube_vertex_increment
      (dimension, this->AxisIncrement(), this->cube_vertex_increment);

    for (DTYPE ifacet = 0; ifacet < this->NumCubeFacets(); ifacet++) {
      compute_facet_vertex_increment
        (dimension, ifacet, this->cube_vertex_increment, 
         this->facet_vertex_increment+this->num_facet_vertices*ifacet);
    }

    // Initialize ridge_vertex_increment to zero.
    for (NTYPE i = 0; i < this->num_cube_ridge_vertices*dimension*dimension; 
         i++) 
      {this->ridge_vertex_increment[i] = 0; }

    // Compute ridge_vertex_increment.
    for (DTYPE orth_dir1 = 0; orth_dir1 < dimension; orth_dir1++)
      for (DTYPE orth_dir0 = 0; orth_dir0 < dimension; orth_dir0++)
        if (orth_dir0 != orth_dir1) {
          compute_ridge_vertex_increment
            (dimension, orth_dir0, orth_dir1,
             this->cube_vertex_increment, 
             this->ridge_vertex_increment+
             this->num_cube_ridge_vertices*(orth_dir0+dimension*orth_dir1));
        }

    compute_unit_cube_vertex_coord
      (dimension, this->unit_cube_coord);
  }

  /// Free memory in the derived typename GRID_PLUS.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::FreeLocal()
  {
    if (axis_increment != NULL) { delete [] axis_increment; }
    if (cube_vertex_increment != NULL) 
      { delete [] cube_vertex_increment; };
    if (facet_vertex_increment != NULL)
      { delete [] facet_vertex_increment; };
    if (unit_cube_coord != NULL) { delete [] unit_cube_coord; };
    ZeroLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    IJK::PROCEDURE_ERROR error("GRID_PLUS::SetSize");

    FreeLocal();

    GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize(dimension, axis_size);
    Create();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename VTYPE0, typename VTYPE1, typename DIST_TYPE>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GetVertexNeighbors
  (const VTYPE0 iv0, const DIST_TYPE distance, std::vector<VTYPE1> & vlist) 
    const
  {
    get_grid_vertices_in_neighborhood
      (this->Dimension(), this->AxisSize(), this->AxisIncrement(), 
       iv0, distance, vlist);
  }

  // **************************************************
  // TEMPLATE CLASS GRID_NEIGHBORS MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::GRID_NEIGHBORS
  (const DTYPE dimension, const ATYPE * axis_size):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (dimension,axis_size)
  {
    InitLocal();
  }

  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::GRID_NEIGHBORS()
  {
    InitLocal();
  }

  /// Constructor from another GRID_NEIGHBORS grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::
  GRID_NEIGHBORS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::
  GRID_NEIGHBORS
  (const GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE, NTYPE> & grid2):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (grid2)
  {
    InitLocal();
  }

  /// Set all local (not inherited) arrays to NULL.
  /// Set all local variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::ZeroLocal()
  {
    this->vertex_neighborC = NULL;
    this->vertex_neighborE = NULL;
    this->cube_neighborV = NULL;
    this->facet_neighborC = NULL;
    this->edge_neighborF2 = NULL;
    this->num_vertex_neighborsC = 0;
    this->num_vertex_neighborsE = 0;
    this->num_facet_neighborsC = 0;
    this->num_edge_neighborsF2 = 0;
  }

  /// \briefInitialize data structures in GRID_NEIGHBORS
  /// @pre \a dimension and \a axis_size are already set
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::InitLocal()
  {
    ZeroLocal();
    if (this->Dimension() > 0) { CreateLocal(); };
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::~GRID_NEIGHBORS()
  {
    FreeLocal();
  }

  /// Allocate arrays and compute data in GRID_NEIGHBORS.
  /// @pre \a dimension and \a axis_size[] are already set.
  /// @pre All other arrays are set to NULL.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::CreateLocal()
  {
    const DTYPE dimension = this->Dimension();
    const ATYPE * axis_size = this->AxisSize();
    IJK::PROCEDURE_ERROR error("GRID_NEIGHBORS::CreateLocal");

    FreeLocal();

    compute_num_vertex_neighborsC(dimension, num_vertex_neighborsC);
    vertex_neighborC = new DIFFTYPE[num_vertex_neighborsC];
    compute_vertex_neighborC(dimension, axis_size, vertex_neighborC);

    compute_num_vertex_neighborsE
      (dimension, num_vertex_neighborsE);
    vertex_neighborE = new DIFFTYPE[num_vertex_neighborsE];
    compute_vertex_neighborE(dimension, axis_size, vertex_neighborE);

    cube_neighborV = new DIFFTYPE[this->NumCubeVertices()];
    compute_cube_neighborV(dimension, axis_size, cube_neighborV);

    cube_neighborE = new DIFFTYPE[this->NumCubeEdges()];
    compute_cube_neighborE(dimension, axis_size, cube_neighborE);

    compute_num_facet_neighborsC
      (dimension, num_facet_neighborsC);
    facet_neighborC = 
      new DIFFTYPE[num_facet_neighborsC*(dimension)];
    compute_facet_neighborC(dimension, axis_size, facet_neighborC);

    compute_num_edge_neighborsF2
      (dimension, num_edge_neighborsF2);
    this->edge_neighborF2 = 
      new DIFFTYPE[num_edge_neighborsF2*(dimension)];
    compute_edge_neighborF2(dimension, axis_size, edge_neighborF2);
  }

  // \brief Free all local (not inherited) arrays.
  // Set all local arrays to NULL and variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::FreeLocal()
  {
    if (vertex_neighborC != NULL) { delete [] vertex_neighborC; };
    if (vertex_neighborE != NULL) { delete [] vertex_neighborE; };
    if (facet_neighborC != NULL) { delete [] facet_neighborC; };
    if (edge_neighborF2 != NULL) { delete [] edge_neighborF2; };
    ZeroLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    IJK::PROCEDURE_ERROR error("GRID_NEIGHBORS::SetSize");

    FreeLocal();

    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize(dimension, axis_size);
    CreateLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  // **************************************************
  // TEMPLATE CLASS GRID_SPACING MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename STYPE, typename GRID_TYPE>
  GRID_SPACING<STYPE, GRID_TYPE>::GRID_SPACING():GRID_TYPE()
  {
    InitLocal();
  }

  /// Constructor.
  template <typename STYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  GRID_SPACING<STYPE, GRID_TYPE>::
  GRID_SPACING(const DTYPE2 dimension, const ATYPE2 * axis_size):
    GRID_TYPE(dimension, axis_size)
  {
    InitLocal();
  }

  /// Constructor from another grid.
  template <typename STYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  GRID_SPACING<STYPE, GRID_TYPE>::
  GRID_SPACING(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID_TYPE(grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename STYPE, typename GRID_TYPE>
  GRID_SPACING<STYPE, GRID_TYPE>::
  GRID_SPACING(const GRID_SPACING<STYPE, GRID_TYPE> & grid2):
    GRID_TYPE(grid2)
  {
    InitLocal();
    SetSpacing(grid2.SpacingPtrConst());
  }

  /// Initialize.
  template <typename STYPE, typename GRID_TYPE>
  void GRID_SPACING<STYPE, GRID_TYPE>::InitLocal()
  {
    spacing = NULL;
    CreateLocal();
    SetAllSpacing(1);
  }

  /// Create data structure.
  template <typename STYPE, typename GRID_TYPE>
  void GRID_SPACING<STYPE, GRID_TYPE>::CreateLocal()
  {
    FreeLocal();
    spacing = new STYPE[this->dimension];
  }

  /// Free memory.
  template <typename STYPE, typename GRID_TYPE>
  void GRID_SPACING<STYPE, GRID_TYPE>::FreeLocal()
  {
    if (spacing != NULL) {
      delete [] spacing;
      spacing = NULL;
    };
  }

  /// Set spacing along all axes to c.
  template <typename STYPE, typename GRID_TYPE>
  template <typename STYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::SetAllSpacing(const STYPE2 c)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, c); }
  }

  /// Set spacing[d] to spacing2[d].
  template <typename STYPE, typename GRID_TYPE>
  template <typename STYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::SetSpacing(const STYPE2 * spacing2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, spacing2[d]); }
  }

  /// Set spacing[d] to (c*spacing2[d]).
  template <typename STYPE, typename GRID_TYPE>
  template <typename SCALE_TYPE, typename STYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::
  SetSpacing(const SCALE_TYPE c, const STYPE2 * spacing2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, c*spacing2[d]); }
  }

  /// Set spacing[d] to grid2.Spacing(d).
  template <typename STYPE, typename GRID_TYPE>
  template <typename STYPE2, typename GRID_TYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::
  SetSpacing(const GRID_SPACING<STYPE2, GRID_TYPE2> & grid2)
  {
    SetSpacing(grid2.SpacingPtrConst());
  }

  template <typename STYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeLocal();
    GRID_TYPE::SetSize(dimension, axis_size);
    CreateLocal();
    SetAllSpacing(1);
  }

  template <typename STYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  void GRID_SPACING<STYPE, GRID_TYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  template <typename STYPE, typename GRID_TYPE>
  template <typename VTYPE2, typename CTYPE>
  void GRID_SPACING<STYPE, GRID_TYPE>::
  ComputeScaledCoord(const VTYPE2 iv, CTYPE * coord) const
  {
    compute_scaled_coord
      (iv, this->Dimension(), this->AxisSize(), SpacingPtrConst(), coord);
  }

  template <typename STYPE, typename GRID_TYPE>
  template <typename VTYPE2, typename CTYPE>
  void GRID_SPACING<STYPE, GRID_TYPE>::
  ComputeCubeCenterScaledCoord(const VTYPE2 iv, CTYPE * coord) const
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    ComputeScaledCoord(iv, coord);
    for (DTYPE d = 0; d < this->Dimension(); d++) 
      { coord[d] += Spacing(d)/2.0; }
  }

  /// Destructor.
  template <typename STYPE, typename GRID_TYPE>
  GRID_SPACING<STYPE, GRID_TYPE>::~GRID_SPACING()
  {
    FreeLocal();
  }

  // **************************************************
  // TEMPLATE CLASS GRID_VERTEX_LIST MEMBER FUNCTIONS
  // **************************************************
  
  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::Init()
  {
    this->vertex_list = NULL;
    this->list_length = 0;
    this->num_vertices = 0;
  }

  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::FreeAll()
  {
    if (vertex_list != NULL) { delete [] vertex_list; };
    vertex_list = NULL;
    this->list_length = 0;
    num_vertices = 0;
  }

  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::AllocateList
  (const VTYPE list_length)
  {
    FreeAll();

    if (list_length > 0) {
      vertex_list = new VTYPE[list_length];
      this->list_length = list_length;
      num_vertices = 0;
    }
  }

  template <typename VTYPE>
  template<typename GCLASS>
  void FACET0_CUBE_LIST<VTYPE>::GetCubes(const GCLASS & grid)
  {
    VTYPE num_cubes;
    compute_num_cubes_in_grid_facet0
      (grid.Dimension(), grid.AxisSize(), num_cubes);

    if (num_cubes > this->ListLength()) 
      { this->AllocateList(num_cubes); }

    if (num_cubes > 0) {
      get_cubes_in_grid_facet0
        (grid.Dimension(), grid.AxisSize(), this->vertex_list);
    }

    this->num_vertices = num_cubes;
  }

  /// GRID_BOUNDARY_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  GRID_BOUNDARY_VERTEX_LIST<VTYPE>::GRID_BOUNDARY_VERTEX_LIST
  (const GCLASS & grid)
  {
    VTYPE numv = 0;

    compute_num_boundary_grid_vertices
      (grid.Dimension(), grid.AxisSize(), numv);

    AllocateList(numv);

    if (numv > 0) {
      get_boundary_grid_vertices
        (grid.Dimension(), grid.AxisSize(), this->vertex_list);
    }

    this->num_vertices = numv;
  }

  /// FACET_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  FACET_VERTEX_LIST<VTYPE>::FACET_VERTEX_LIST
  (const GCLASS & grid, const VTYPE orth_dir,
   const bool allocate_max)
  {
    VTYPE numv = 0;

    if (allocate_max) {
      compute_max_num_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), numv);
    }
    else {
      compute_num_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), orth_dir, numv);
    }

    this->AllocateList(numv);
    GetVertices(grid, orth_dir);
  }

  /// Get vertices in grid facet
  template <typename VTYPE>
  template<typename GCLASS>
  void FACET_VERTEX_LIST<VTYPE>::GetVertices
  (const GCLASS & grid, const VTYPE orth_dir)
  {
    VTYPE numv;
    const bool side = false;

    compute_num_vertices_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, numv);

    if (numv > this->ListLength()) 
      { this->AllocateList(numv); }

    if (numv > 0) {
      get_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), orth_dir, side,
         this->vertex_list);
    }

    this->num_vertices = numv;
  }

  /// FACET_INTERIOR_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  FACET_INTERIOR_VERTEX_LIST<VTYPE>::FACET_INTERIOR_VERTEX_LIST
  (const GCLASS & grid, const VTYPE orth_dir,
   const bool allocate_max)
  {
    const VTYPE boundary_width = 1;
    VTYPE numv = 0;

    if (allocate_max) {
      compute_max_num_vertices_in_grid_facet_interior
        (grid.Dimension(), grid.AxisSize(), boundary_width, numv);
    }
    else {
      compute_num_vertices_in_grid_facet_interior
        (grid.Dimension(), grid.AxisSize(), orth_dir, boundary_width,
         numv);
    }

    this->AllocateList(numv);
    GetVertices(grid, orth_dir);
  }

  /// Get vertices in grid facet interior
  template <typename VTYPE>
  template<typename GCLASS>
  void FACET_INTERIOR_VERTEX_LIST<VTYPE>::GetVertices
  (const GCLASS & grid, const VTYPE orth_dir)
  {
    VTYPE numv;
    const VTYPE boundary_width = 1;
    const bool side = false;

    compute_num_vertices_in_grid_facet_interior
      (grid.Dimension(), grid.AxisSize(), orth_dir, boundary_width,
       numv);

    if (numv > this->ListLength()) 
      { this->AllocateList(numv); }

    if (numv > 0) {
      get_vertices_in_grid_facet_interior
        (grid.Dimension(), grid.AxisSize(), orth_dir, side,
         boundary_width, this->vertex_list);
    }

    this->num_vertices = numv;
  }

  // **************************************************
  // TEMPLATE OUTPUT FUNCTIONS
  // **************************************************

  /// Output coord (for debugging purposes)
  template <typename DTYPE, typename CTYPE>
  void ijkgrid_output_coord
  (std::ostream & out,
   const DTYPE dimension, const CTYPE * coord)
  {
    out << "(";
    for (DTYPE d = 0; d < dimension; d++) {
      out << coord[d];
      if (d+1 < dimension) 
        { out << ","; }
    }
    out << ")";
  }

  /// Output vertex coord (for debugging purposes)
  template <typename GTYPE, typename VTYPE>
  void ijkgrid_output_vertex_coord
  (std::ostream & out, 
   const GTYPE & grid, const VTYPE iv)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<VTYPE> coord(dimension);

    grid.ComputeCoord(iv, coord.Ptr());
    ijkgrid_output_coord(out, dimension, coord.PtrConst());
  }

}


#endif

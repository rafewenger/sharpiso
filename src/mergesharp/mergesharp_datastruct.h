/// \file mergesharp_datastruct.h
/// Data structure definitions for mergesharp.

/*
  Copyright (C) 2011-2013 Rephael Wenger

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

#ifndef _MERGESHARP_DATASTRUCT_
#define _MERGESHARP_DATASTRUCT_

#include <string>
#include <unordered_map>
#include <vector>

#include "ijk.txx"
#include "ijkobject_grid.txx"
#include "ijkscalar_grid.txx"
#include "ijkvector_grid.txx"
#include "ijkmerge.txx"

#include "ijkdualtable.h"

#include "sharpiso_eigen.h"
#include "sharpiso_grids.h"
#include "sharpiso_feature.h"

#include "mergesharp_types.h"
#include "mergesharp_isovert.h"


namespace MERGESHARP {

  class MERGE_DATA;
  class MULTIRES_GRID;

  // **************************************************
  // DATA STRUCTURES: ARRAY, HASH TABLE
  // **************************************************

  /// Array of vertices.
  typedef IJK::ARRAY<VERTEX_INDEX> VERTEX_ARRAY;

  /// Vertex hash table.
  typedef std::unordered_map<VERTEX_INDEX, NUM_TYPE>
    VERTEX_HASH_TABLE;

  /// Cube hash table.
  typedef VERTEX_HASH_TABLE CUBE_HASH_TABLE;

  // **************************************************
  // DUAL ISOSURFACE VERTICES, CUBE DATA
  // **************************************************

  typedef IJK::DUAL_ISOVERT
    <ISO_VERTEX_INDEX, FACET_VERTEX_INDEX, IJKDUALTABLE::TABLE_INDEX>
    DUAL_ISOVERT;

  // **************************************************
  // DUAL CONTOURING ISOSURFACE CLASS
  // **************************************************

  /// Dual contouring isosurface.
  /// Representation of isosurface returned by Dual Contouring algorithm.
  class DUAL_ISOSURFACE {

  public:
    /// List of vertices of each quadrilateral in isosurface mesh.
    VERTEX_INDEX_ARRAY quad_vert;

    /// List of vertices of each triangle in isosurface mesh.
    /// Not always set.
    VERTEX_INDEX_ARRAY tri_vert;

    /// List of vertex coordinates.
    COORD_ARRAY vertex_coord;

  public:
    DUAL_ISOSURFACE() {};

    NUM_TYPE NumVerticesPerTri() const
    { return(NUM_VERT_PER_TRI); }

    NUM_TYPE NumVerticesPerQuad() const
    { return(NUM_VERT_PER_QUAD); }

    NUM_TYPE NumTri() const
    { return(tri_vert.size()/NumVerticesPerTri()); };

    NUM_TYPE NumQuad() const
    { return(quad_vert.size()/NumVerticesPerQuad()); };

    NUM_TYPE NumVertices() const
    { return(vertex_coord.size()/DIM3); }

    void Clear();
  };

  // **************************************************
  // DUAL CONTOURING INPUT DATA AND DATA STRUCTURES
  // **************************************************

  /// dual contouring parameters
  class MERGESHARP_PARAM:public SHARPISO::SHARP_ISOVERT_PARAM {

  protected:
    void Init();

  public:

    /// Interpolation type. Linear or multi-linear interpolation.
    INTERPOLATION_TYPE interpolation_type;

    /// Isosurface vertex position method.
    /// Position at cube center or at of centroid of grid edge-isosurface 
    ///   intersections or use vertex gradients to position,
    ///   or use isosurface normals computed using interpolation
    ///   or use isosurface normals computed using gradients.
    VERTEX_POSITION_METHOD vertex_position_method; 

    /// If true, convert isosurface quadrilaterals to isosurface triangles.
    bool flag_convert_quad_to_tri;

    /// Quadrilateral triangulation method.
    /// Uniform triangulation or split max quadrilateral angle.
    QUAD_TRI_METHOD quad_tri_method;

    /// If true, merge grid cube near sharp vertices.
    bool flag_merge_sharp;

    /// If true, allow multiple isosurface vertices in a single cube
    ///   to represent multiple isosurface patches.
    bool allow_multiple_iso_vertices;

    /// If true, isosurfaced patches separate negative vertices.
    bool flag_separate_neg;

    /// If true, resolve ambiguous facets based on location of sharp features.
    bool flag_resolve_ambiguous_facets;

    /// If true, split isosurface vertices to avoid non-manifold edges.
    bool flag_split_non_manifold;

    /// If true, select which cube has configuration of split isosurface 
    /// vertices where adjacent cubes share an ambiguous facet and one will 
    /// have one isosurface vertex while the other has two isosurface vertices.
    bool flag_select_split;

    /// If true, remove isolated isosurface vertices.  Isolated vertices
    ///   do not lie in any isosurface polygon.
    bool flag_delete_isolated_vertices;

    /// If true, store isosurface vertex information.
    bool flag_store_isovert_info;

    /// If true, convert gradient to hermite data.
    bool flag_grad2hermite;

    /// If true, convert gradient to hermite data using linear interpolation.
    bool flag_grad2hermiteI;

  public:
    MERGESHARP_PARAM() { Init(); };
    ~MERGESHARP_PARAM() { Init(); };

    // Set functions
    void SetInterpolationType       ///< Set type of interpolation.
      (const INTERPOLATION_TYPE interpolation_type);
    void SetVertexPositionMethod    ///< Set isosurface vertex position method.
      (const VERTEX_POSITION_METHOD vertex_position_method);
    void SetUseSelectedGradients    ///< Set flag for using selected gradients.
      (const bool flag);
    void SetUseOnlyCubeGradients    ///< Set flag for using only cube gradients.
      (const bool flag);

    /// Set flag for using gradients at endpoints of intersected edges.
    void SetUseIntersectedEdgeEndpointGradients    
      (const bool flag);

    /// Set
    void Set(const MERGESHARP_PARAM & param);

    /// Return interpolation type.
    INTERPOLATION_TYPE InterpolationType() const
      { return(interpolation_type); };

    /// Return isosurface vertex position method.
    VERTEX_POSITION_METHOD VertexPositionMethod() const
      { return(vertex_position_method); };

    /// Return true if gradient data required.
    bool GradientsRequired() const;

    /// Return true if edgeI normal data required.
    bool NormalsRequired() const;
  };

  /// Input data to Dual Contouring and related algorithms
  class MERGESHARP_DATA:public MERGESHARP_PARAM {

  protected:
    SHARPISO_SCALAR_GRID scalar_grid;  ///< Regular grid of scalar values.
    GRADIENT_GRID gradient_grid;       ///< Regular grid of vertex gradients.

    /// Coordinate of edge-isosurface intersections
    std::vector<COORD_TYPE> edgeI_coord;  

    /// Coordinates of normal vectors at edge-isosurface intersections.
    std::vector<GRADIENT_COORD_TYPE> edgeI_normal_coord;
    

    // flags
    bool is_scalar_grid_set;
    bool is_gradient_grid_set;
    bool are_edgeI_set;

    void Init();
    void FreeAll();

  public:
    MERGESHARP_DATA() { Init(); };
    ~MERGESHARP_DATA() { FreeAll(); };

    // Set functions
    void CopyScalarGrid             /// Copy scalar_grid to MERGESHARP_DATA
      (const SHARPISO_SCALAR_GRID_BASE & scalar_grid2);
    void CopyGradientGrid             /// Copy gradient_grid to MERGESHARP_DATA
      (const GRADIENT_GRID_BASE & gradient_grid2);
    void SubsampleScalarGrid        /// Subsample scalar_grid.
      (const SHARPISO_SCALAR_GRID_BASE & scalar_grid2, 
       const int subsample_resolution);
    void SubsampleGradientGrid
      (const GRADIENT_GRID_BASE & gradient_grid2, 
       const int subsample_resolution);

    void SupersampleScalarGrid      /// Supersample scalar_grid.
      (const SHARPISO_SCALAR_GRID_BASE & scalar_grid2, 
       const int supersample_resolution);

    /// Copy, subsample or supersample scalar grid.
    /// Precondition: flag_subsample and flag_supersample are not both true.
    void SetScalarGrid
      (const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution);

    /// Copy, subsample or supersample scalar and gradient grids.
    /// Precondition: flag_subsample and flag_supersample are not both true.
    void SetGrids
      (const SHARPISO_SCALAR_GRID_BASE & full_scalar_grid,
       const GRADIENT_GRID_BASE & full_gradient_grid,
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution);

    /// Set edge-isosurface intersections and normals.
    void SetEdgeI(const std::vector<COORD_TYPE> & edgeI_coord,
                  const std::vector<GRADIENT_COORD_TYPE> & edgeI_normal_coord);

    // Get functions
    bool IsScalarGridSet() const     /// Return true if scalar grid is set.
      { return(is_scalar_grid_set); };
    bool IsGradientGridSet() const   /// Return true if gradient grid is set.
      { return(is_gradient_grid_set); };

    /// Return true if edge-intersections (and normals) are set.
    bool AreEdgeISet() const         
      { return(are_edgeI_set); };

    /// Return scalar_grid.
    const SHARPISO_SCALAR_GRID_BASE & ScalarGrid() const
      { return(scalar_grid); };

    /// Return gradient_grid.
    const GRADIENT_GRID_BASE & GradientGrid() const     
      { return(gradient_grid); };

    /// Return edgeI coordinates.
    const std::vector<COORD_TYPE> & EdgeICoord() const
      { return(edgeI_coord); }

    /// Return edgeI normal coordinates.
    const std::vector<GRADIENT_COORD_TYPE> & EdgeINormalCoord() const  
      { return(edgeI_normal_coord); }

    /// Check data structure
    bool Check(IJK::ERROR & error) const;
  };


  // **************************************************
  // MERGESHARP TIME
  // **************************************************

  /// dual contouring time.
  /// Uses system clock() function to determine time.  
  /// System clock() function should return cpu time, 
  ///   but may return wall time.
  class MERGESHARP_TIME {

  public:
    // all times are in seconds
    float preprocessing;  
    // time to create data structure for faster isosurface extraction
    float extract;          // time to extract isosurface mesh
    float merge_identical;  // time to merge identical vertices
    float position;         // time to position isosurface vertices
    float merge_sharp;      // time to merge sharp isosurface vertices
    float total;            // extract_time+merge_time+position_time

    MERGESHARP_TIME();
    void Clear();
    void Add(const MERGESHARP_TIME & isodual_time);
  };

  // **************************************************
  // GRID INFO
  // **************************************************

  /// Regular grid information.
  class GRID_INFO {

  public:
    GRID_INFO();                    ///< Constructor.

    VERTEX_INDEX num_cubes;         ///< Number of grid cubes.

    void Clear();                   ///< Clear all data.
  };

  // **************************************************
  // SCALAR INFO
  // **************************************************

  /// Scalar grid information.
  class SCALAR_INFO {

  protected:
    int dimension;

    void Copy(const SCALAR_INFO & info);  ///< Copy scalar info.
    void Init(const int dimension);       ///< Initialize scalar info.
    void FreeAll();                       ///< Free all memory.

  public:
    SCALAR_INFO() { Init(3); };
    SCALAR_INFO(const int dimension) { Init(dimension); };
    ~SCALAR_INFO();
    SCALAR_INFO(const SCALAR_INFO & info) ///< Copy constructor. 
      { Copy(info); };       
    const SCALAR_INFO & operator =        ///< Copy assignment.
      (const SCALAR_INFO & right);

    VERTEX_INDEX num_non_empty_cubes;     ///< Number of cubes containing an isosurface simplex.

    VERTEX_INDEX num_bipolar_edges; ///< Number of bipolar edges.
    //   (number of edges with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)


    // Set functions.
    void SetDimension(const int dimension); ///< Set dimension.

    // Get functions.
    int Dimension() const { return(dimension); };

    void Clear();     // clear all data
  };

  // **************************************************
  // SHARPISO INFO
  // **************************************************

  class DUAL_ISOVERT_INFO:public DUAL_ISOVERT {

  public:
    unsigned char num_eigenvalues;  ///< Number of eigenvalues.

    /// If true, location is centroid of (grid edge)-isosurface intersections.
    bool flag_centroid_location;

    /// Set DUAL_ISOVERT_INFO from DUAL_ISOVERT.
    void Set(const MERGESHARP::DUAL_ISOVERT & isov);
  };

  /// Sharpiso information.
  class SHARPISO_INFO:public ISOVERT_INFO {

  public:
    int num_edge_collapses;
    int num_repositioned_vertices;

    // Ambiguity information.
    int num_cube_not_ambiguous;
    int num_cube_separate_pos;
    int num_cube_separate_neg;
    int num_cube_unresolved_ambiguity;

    // Multiple isosurface vertices information.
    int num_cube_single_isov;  ///< Number of cubes with single iso vertices.
    int num_cube_multi_isov;   ///< Number of cubes with multiple iso vertices.

    /// Number of merges blocked by non-disk isosurface patches.
    int num_non_disk_isopatches;

    /// Number of splits to avoid non-manifold edges.
    int num_non_manifold_split;

    /// Number of cube configuration changes from 1 to 2 isosurface vertices
    ///   or vice versa.
    int num_1_2_change;

    SHARPISO_INFO();  ///< Constructor.
    void Clear();     ///< Clear all data.

    /// Increment num_sharp_corners or num_sharp_edges or num_smooth_vertices
    ///   depending upon the number of large eigenvalues.
    void IncrementIsoVertexNum(const int num_large_eigenvalues);

    /// Isosurface vertex information.
    std::vector<DUAL_ISOVERT_INFO> vertex_info;
  };

  // **************************************************
  // MERGESHARP INFO
  // **************************************************

  /// dual contouring information.
  /// Statistical and timing information from the Dual Contouring algorithm.
  class MERGESHARP_INFO {

  public:
    GRID_INFO grid;
    SCALAR_INFO scalar;
    MERGESHARP_TIME time;
    SHARPISO_INFO sharpiso;

    MERGESHARP_INFO();
    MERGESHARP_INFO(const int dimension);

    void Clear();     // clear all data
  };

  // **************************************************
  // MERGE DATA
  // **************************************************

  /// Internal data structure for merge_identical_vertices 
  class MERGE_DATA: public IJK::INTEGER_LIST<MERGE_INDEX, MERGE_INDEX> {

  protected:
    MERGE_INDEX num_edges;             ///< Number of edges.
    MERGE_INDEX num_vertices;          ///< Number of vertices.
    MERGE_INDEX num_obj_per_vertex;    ///< Number of objects per vertex.
    MERGE_INDEX num_obj_per_edge;      ///< Number of objects per edge.
    MERGE_INDEX num_obj_per_grid_vertex; ///< Number of objects per grid vertex.
    MERGE_INDEX vertex_id0;            ///< First vertex identifier.

    /// Initialize.
    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size,
              const MERGE_INDEX num_obj_per_vertex,
              const MERGE_INDEX num_obj_per_edge);

  public:
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size)
      { Init(dimension, axis_size, 0, 1); };
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size,
               const MERGE_INDEX num_obj_per_vertex, 
               const MERGE_INDEX num_obj_per_edge)
      { Init(dimension, axis_size, num_obj_per_vertex, num_obj_per_edge); };

    // get functions
    MERGE_INDEX NumEdges() const        /// Number of edges.
      { return(num_edges); };
    MERGE_INDEX NumVertices() const     /// Number of vertices.
      { return(num_vertices); };
    MERGE_INDEX NumObjPerVertex() const { return(num_obj_per_vertex); };
    MERGE_INDEX NumObjPerEdge() const { return(num_obj_per_edge); };
    MERGE_INDEX NumObjPerGridVertex() const
      { return(num_obj_per_grid_vertex); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv) const { return(vertex_id0 + iv); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv, const MERGE_INDEX j) const 
      { return(vertex_id0 + iv + j * num_vertices); };
    MERGE_INDEX EdgeIdentifier         /// Edge identifier.
      (const MERGE_INDEX ie) const { return(ie); };
    MERGE_INDEX EdgeIdentifier         /// edge identifier
      (const MERGE_INDEX ie, const MERGE_INDEX j) const 
      { return(ie + j * num_edges); };

    /// Get first endpoint of edge containing isosurface vertex isov.
    inline VERTEX_INDEX GetFirstEndpoint(const MERGE_INDEX isov) const
      { return(isov/NumObjPerGridVertex()); };

    /// Get direction of edge containing isosurface vertex isov.
    inline MERGE_INDEX GetEdgeDir(const MERGE_INDEX isov) const
      { return(isov%NumObjPerGridVertex()); };

    bool Check(IJK::ERROR & error) const;     ///< Check allocated memory.
  };

  /// Merge data structure for isosurface vertices.
  /// Vertices are identified by a single integer.
  class ISO_MERGE_DATA: public MERGE_DATA {
  public:
    ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 0) {};
  };

  // **************************************************
  // MESH VERTEX LIST
  // **************************************************

  /// List of mesh vertices.
  template <class DTYPE, class CTYPE, class STYPE>
    class MESH_VERTEX_LIST {

    protected:
    DTYPE dimension;        ///< Dimension.

    /// Coordinates of mesh vertices.
    /// coord[i*dimension+j] = j'th coordinate of i'th mesh vertex.
    std::vector<CTYPE> coord;

    /// scalar[i] = scalar value of i'th mesh vertex.
    std::vector<STYPE> scalar;   

    public:
    MESH_VERTEX_LIST(const DTYPE dim) { this->dimension = dim; };

    // Get functions.
    DTYPE Dimension() const { return(dimension); } ///< Return dimension.

    /// Return j'th coordinate of i'th mesh vertex.
    CTYPE Coord(const int i, const int j) const
      { return(coord[i*Dimension()+j]); };

    /// Copy coordinates of \a i'th mesh vertex into array \a coord[].
    /// @pre coord[] is preallocated to length at least \a dimension.
    template <class CTYPE2>
      void CopyCoord(const int i, CTYPE2 * coord2) const
      {
        const CTYPE * coord_begin = &(coord[i*Dimension()]);
        std::copy(coord_begin, coord_begin+Dimension(), coord2);
      }

    /// Return scalar value of i'th mesh vertex.
    STYPE Scalar(const int i) const { return(scalar[i]); };

    int NumVertices() const     ///< Return number of vertices in list.
      { return(scalar.size()); };

    // Set functions.

    void Add()                 ///< Add a vertex to the list.
      {
        const int k = coord.size();
        coord.resize(k+Dimension());
        scalar.push_back(0);
      }

    void SetScalar(const int i, const STYPE s)  ///< Set scalar of vertex i.
      { scalar[i] = s; };

    /// Set j'th coordinate of vertex i.
    void SetCoord(const int i, const int j, const CTYPE c)
      { coord[i*Dimension()+j] = c; };

    /// Set coordinates of vertex i
    template <class CTYPE2>
      void SetCoord(const int i, CTYPE2 * coord2)
      {
        for (DTYPE d = 0; d < Dimension(); d++)
          { SetCoord(i, d, coord2[d]); };
      }

  };

  typedef MESH_VERTEX_LIST<int, COORD_TYPE, SCALAR_TYPE> 
    MERGESHARP_MESH_VERTEX_LIST;    ///< Dual Contouring mesh vertex list.


  // **************************************************
  // SUBROUTINES
  // **************************************************

  /// Store isosurface vertex information in isovert_info.
  void set_isovert_info
    (const std::vector<MERGESHARP::DUAL_ISOVERT> & iso_vlist,
     const std::vector<GRID_CUBE> & gcube_list,
     std::vector<DUAL_ISOVERT_INFO> & isovert_info);

  /// Delete vertices i where flag_keep[i] = false.
  void delete_vertices
  (const std::vector<MERGESHARP::DUAL_ISOVERT> & iso_vlist,
   const std::vector<NUM_TYPE> & new_isovert_index,
   const std::vector<bool> & flag_keep,
   std::vector<MERGESHARP::DUAL_ISOVERT> & new_iso_vlist);

}

#endif

/// \file ijkcontourtree.txx
/// ijk templates defining contour tree routines and data structures
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009,2012 Rephael Wenger

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

#ifndef _IJKCONTOURTREE_
#define _IJKCONTOURTREE_

#include <algorithm>
#include <vector>

#include "ijk.txx"
#include "ijkscalar_grid.txx"

namespace IJK {

  /// Contour tree classes and functions
  namespace CONTOUR_TREE {

    // **************************************************
    // TYPE DEFINITIONS
    // **************************************************

    typedef enum { AUGMENT_NONE, AUGMENT_JS, AUGMENT_ALL } AUGMENT_FLAG;
    typedef enum { MERGE_IDENT_NONE, MERGE_IDENT_CTREE, 
                   MERGE_IDENT_ALL } MERGE_IDENT_FLAG;
    typedef enum { JOIN_EDGE, SPLIT_EDGE } EDGE_TYPE;
    typedef enum { EDGE_CONNECT, CUBE_CONNECT, FACET_CONNECT,
                   SIMPLEX_CONNECT } CONNECTIVITY;

    // **************************************************
    // TEMPLATE CLASS TREE_NODE
    // **************************************************

    /// Join/split tree node base class.
    class TREE_NODE {

    protected:
      TREE_NODE * parent;

    public:
      TREE_NODE() { parent = NULL; };

      // set functions
      void SetParent(TREE_NODE * p) { parent = p; };

      // get functions
      TREE_NODE * Parent() const { return(parent); };
      bool IsRoot() const { return(parent == NULL); };

      // IsDeleted() always returns false
      // Allows use of TREE_NODE by templates calling IsDeleted()
      bool IsDeleted() const { return(false); };
    };

    template <typename TREE_NODE_TYPE, typename WTYPE>
    class WEIGHTED_TREE_NODE:public TREE_NODE_TYPE {

    protected:
      WTYPE weight;

    public:
      typedef WTYPE WEIGHT_TYPE;         ///< Weight type.

    public:
      WEIGHTED_TREE_NODE() { weight = 1; };

      // set functions
      void SetWeight(const WTYPE w)
      { weight = w; };
      void AddWeight(const WTYPE w)
      { weight += w; }

      // get functions
      WTYPE Weight() const { return(weight); };
    };

    template <typename TREE_NODE_TYPE, typename STYPE>
    class SCALAR_TREE_NODE:public TREE_NODE_TYPE {

    protected:
      STYPE scalar;

    public:
      typedef STYPE SCALAR_TYPE;         ///< Scalar type.

    public:
      SCALAR_TREE_NODE() { scalar = 1; };

      // set functions
      void SetScalar(const STYPE s)
      { scalar = s; };

      // get functions
      STYPE Scalar() const { return(scalar); };

      // Copy functions
      // Note: Copying tree nodes corrupts the Parent() pointer.
      void Copy(const TREE_NODE_TYPE & x)
      {
        if (&x != this) {
          this->TREE_NODE_TYPE::Copy(x);
        }
      }

    };

    /// Tree node with delete flag.
    template <typename TREE_NODE_TYPE>
    class TREE_NODE_DEL:public TREE_NODE_TYPE {

    protected:
      bool is_deleted;           /// true, if node is deleted.

    public:
      TREE_NODE_DEL() { is_deleted = false; };

      /// Redefine set parent
      void SetParent(TREE_NODE_DEL * p) 
      { TREE_NODE_TYPE::SetParent(static_cast<TREE_NODE_TYPE *>(p)); };

      // get functions
      bool IsDeleted() const     /// Return true if node is deleted.
      { return(is_deleted); };
      TREE_NODE_DEL * Parent() const;  /// Get parent, skipping deleted nodes.

      // set functions
      void Delete()              /// delete function
      { is_deleted = true; };

      TREE_NODE_DEL * ParentCompress();     /// Compress tree path to parent
    };

    /// Weighted tree node with delete flag.
    template <typename TREE_NODE_TYPE, typename WTYPE>
    class WEIGHTED_TREE_NODE_DEL:
      public TREE_NODE_DEL<WEIGHTED_TREE_NODE<TREE_NODE_TYPE,WTYPE> > {

    public:
      WEIGHTED_TREE_NODE_DEL() {};

      /// Redefine get parent
      WEIGHTED_TREE_NODE_DEL * Parent() const
      { return(static_cast<WEIGHTED_TREE_NODE_DEL *>
               (TREE_NODE_DEL<WEIGHTED_TREE_NODE<TREE_NODE_TYPE,WTYPE> >::
                Parent())); }

      /// Redefine ParentCompress.
      /// Compress tree path to parent.
      WEIGHTED_TREE_NODE_DEL * ParentCompress();     
    };

    /// Extended tree node class.
    template <typename NTYPE>
    class TREE_NODE_EXT:public TREE_NODE_DEL<TREE_NODE> {

    protected:
      NTYPE num_children;

      void Init();

    public:
      TREE_NODE_EXT() { Init(); };

      // set functions
      void SetNumChildren(const NTYPE n) { num_children = n; };
      void IncrementNumChildren() { num_children++; };
      void DecrementNumChildren() { num_children--; };

      /// Redefine set parent
      void SetParent(TREE_NODE_EXT * p) 
      { parent = static_cast<TREE_NODE *>(p); };

      // get functions
      NTYPE NumChildren() const { return(num_children); };
      bool IsDeleted() const { return(is_deleted); };

      /// Redefine get parent
      TREE_NODE_EXT * Parent() const
      { return(static_cast<TREE_NODE_EXT *>
	       (TREE_NODE_DEL<TREE_NODE>::Parent())); }

      /// Redefine ParentCompress()
      TREE_NODE_EXT * ParentCompress()
      { return(static_cast<TREE_NODE_EXT *>
               (TREE_NODE_DEL<TREE_NODE>::ParentCompress())); }

      // Redefine delete function
      void Delete();
    };

    /// Tree node with node identifier
    template <typename ITYPE, typename TREE_NODE_TYPE>
    class TREE_NODE_ID:public TREE_NODE_TYPE {

    protected:
      ITYPE ident;               /// node identifier

    public:
      typedef ITYPE IDENT_TYPE;

    public:
      TREE_NODE_ID() {};

      // set functions
      void SetIdent(const ITYPE ident) { this->ident = ident; };

      /// Redefine Parent()
      TREE_NODE_ID * Parent() const
      { return(static_cast<TREE_NODE_ID *>(TREE_NODE_TYPE::Parent())); }

      /// Redefine set parent
      void SetParent(TREE_NODE_ID * p)
      { TREE_NODE_TYPE::SetParent(static_cast<TREE_NODE_TYPE *>(p)); };

      // get functions
      ITYPE Ident() const { return(ident); };
    };


    /// Contour tree node template.
    template <typename TREE_NODE_TYPE>
    class CONTOUR_TREE_NODE:public TREE_NODE_TYPE {

    protected:
      EDGE_TYPE edge_type;

    public:
      CONTOUR_TREE_NODE() { edge_type = JOIN_EDGE; };

      /// Redefine set parent
      void SetParent(CONTOUR_TREE_NODE<TREE_NODE_TYPE> * p) 
      { TREE_NODE_TYPE::SetParent(static_cast<TREE_NODE_TYPE *>(p)); };

      // get functions
      EDGE_TYPE EdgeType() const { return(edge_type); };

      /// Redefine get parent
      CONTOUR_TREE_NODE<TREE_NODE_TYPE> * Parent() const
      { return(static_cast<CONTOUR_TREE_NODE<TREE_NODE_TYPE> *>
	       (TREE_NODE_TYPE::Parent())); }

      // set functions
      void SetEdgeType(const EDGE_TYPE edge_type)
      { this->edge_type = edge_type; };
      void SetJoinEdge() { this->edge_type = JOIN_EDGE; };
      void SetSplitEdge() { this->edge_type = SPLIT_EDGE; };

      /// Copy
      // Note: Copying tree nodes corrupts the Parent() pointer.
      void Copy(const CONTOUR_TREE_NODE & x)
      {
        *this = x;
      };

    };

    /// Join tree node template.
    template <typename TREE_NODE_TYPE>
    class JOIN_TREE_NODE:public TREE_NODE_TYPE {

    protected:
      EDGE_TYPE edge_type;

    public:
      JOIN_TREE_NODE() {};

      /// Redefine set parent
      void SetParent(JOIN_TREE_NODE<TREE_NODE_TYPE> * p) 
      { TREE_NODE_TYPE::SetParent(static_cast<TREE_NODE_TYPE *>(p)); };

      // Always return edge type JOIN_EDGE
      EDGE_TYPE EdgeType() const { return(JOIN_EDGE); };

      /// Redefine get parent
      JOIN_TREE_NODE<TREE_NODE_TYPE> * Parent() const
      { return(static_cast<JOIN_TREE_NODE<TREE_NODE_TYPE> *>
	       (TREE_NODE_TYPE::Parent())); }

      /// Copy
      // Note: Copying tree nodes corrupts the Parent() pointer.
      void Copy(const JOIN_TREE_NODE & x)
      {
        *this = x;
      }
    };

    /// Split tree node template.
    template <typename TREE_NODE_TYPE>
    class SPLIT_TREE_NODE:public TREE_NODE_TYPE {

    protected:
      EDGE_TYPE edge_type;

    public:
      SPLIT_TREE_NODE() {};

      /// Redefine set parent
      void SetParent(SPLIT_TREE_NODE<TREE_NODE_TYPE> * p) 
      { TREE_NODE_TYPE::SetParent(static_cast<TREE_NODE_TYPE *>(p)); };

      // Always return edge type SPLIT_EDGE
      EDGE_TYPE EdgeType() const { return(SPLIT_EDGE); };

      /// Redefine get parent
      SPLIT_TREE_NODE<TREE_NODE_TYPE> * Parent() const
      { return(static_cast<SPLIT_TREE_NODE<TREE_NODE_TYPE> *>
	       (TREE_NODE_TYPE::Parent())); }

      /// Copy
      void Copy(const SPLIT_TREE_NODE & x)
      {
        *this = x;
      }
    };

    /// Return index of parent node in array node_array.
    /// Precondition: Parent() is not NULL.
    template <typename ITYPE, typename NODE_TYPE>
    ITYPE ParentIndex(NODE_TYPE * node_array,
                      const ITYPE i)
    { return(node_array[i].Parent() - node_array); };

    /// Return index of parent node in array node_array.
    /// Compress tree path to parent.
    /// Precondition: Parent() is not NULL.
    template <typename ITYPE, typename NODE_TYPE>
    ITYPE ParentCompressIndex
    (NODE_TYPE * node_array, const ITYPE i)
    { return(static_cast<NODE_TYPE *>(node_array[i].ParentCompress()) 
             - node_array); };

    template <typename ITYPE, typename NODE_TYPE>
    ITYPE ParentIndex
    (IJK::ARRAY<NODE_TYPE> & node_array,const ITYPE i)
    {
      return(ParentIndex(node_array.Ptr(), i));
    }

    // **************************************************
    // TEMPLATE CLASS UNION_FIND
    // **************************************************

    /// Data structure supporting union-find operations.
    template <typename NTYPE>
    class UNION_FIND_BASE {
    protected:
      NTYPE * parent;
      NTYPE num_elements;

      void Init(const NTYPE num_elements);
      void FreeAll();
      NTYPE FindRoot(const NTYPE e);

    public:
      UNION_FIND_BASE(const NTYPE num_elements) { Init(num_elements); };
      ~UNION_FIND_BASE() { FreeAll(); };

      // get functions
      NTYPE NumElements() const { return(num_elements); };

      // Union-Find functions
      NTYPE Find(const NTYPE e);
      void Union(const NTYPE e1, const NTYPE e2);
    };

    /// Union-find data structure with data associated with each set.
    template <typename NTYPE, typename DTYPE>
    class UNION_FIND:public UNION_FIND_BASE<NTYPE> {
    protected:
      DTYPE * data;

      void Init(const NTYPE num_elements);
      void FreeAll();

    public:
      UNION_FIND(const NTYPE num_elements):
        UNION_FIND_BASE<NTYPE>(num_elements)
      { Init(num_elements); };
      ~UNION_FIND() { FreeAll(); };

      // get functions
      DTYPE Data(const NTYPE e);

      // set functions
      void SetData(const NTYPE e, const DTYPE d) { data[e] = d; };

      // Redeclare union function to set data
      void Union(const NTYPE e1, const NTYPE e2, const DTYPE d);
    };

    // **************************************************
    // TEMPLATE CLASS VERTEX_NEIGHBORS
    // **************************************************

    template <typename DTYPE, typename ATYPE, typename VTYPE, 
	      typename NTYPE, typename BTYPE>
    class VERTEX_NEIGHBORS:public GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> {

    protected:
      NTYPE num_neighbors;
      VTYPE origin_decrement;
      VTYPE origin;
      VTYPE vertex_index;
      VTYPE * neighbor_increment;
      BTYPE * boundary_bits;

      void SetAllToZero();

    public:
      VERTEX_NEIGHBORS(const DTYPE dimension, const ATYPE * axis_size);
      ~VERTEX_NEIGHBORS();

      // Set function
      void SetVertex(const VTYPE vertex_index);
      void SetOrigin(const VTYPE origin);

      // Get functions
      NTYPE NumNeighbors() const { return(num_neighbors); };
      VTYPE Neighbor(const NTYPE i)
      { return(origin+neighbor_increment[i]); };
      VTYPE VertexIndex() const
      { return(vertex_index); };
      VTYPE Origin() const
      { return(origin); };
      VTYPE OriginDecrement() const { return(origin_decrement); };
      BTYPE BoundaryBits(const NTYPE i) const
      { return(boundary_bits[i]); };
    };

    /// Set the vertex index. Neighbors surround vertex.
    /// @pre Vertex is not on grid boundary.
    /// @pre vertex_index >= origin_decrement.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE,
	      typename BTYPE>
    void VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NTYPE,BTYPE>::
    SetVertex(const VTYPE vertex_index)
    {
      this->vertex_index = vertex_index;
      this->origin = vertex_index - origin_decrement;
    };

    /// Set the origin. vertex_index = origin + origin_decrement.
    /// @pre Origin >= 0.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE,
	      typename BTYPE>
    void VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NTYPE,BTYPE>::
    SetOrigin(const VTYPE origin)
    {
      this->origin = origin;
      this->vertex_index = origin + origin_decrement;
    };

    /// Set all values to zero.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE,
	      typename BTYPE>
    void VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NTYPE,BTYPE>::SetAllToZero()
    {
      origin = 0;
      vertex_index = 0;
      origin_decrement = 0;
      num_neighbors = 0;
      neighbor_increment = NULL;
      boundary_bits = NULL;
    }
    
    /// Constructor.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE,
	      typename BTYPE>
    VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NTYPE,BTYPE>::
    VERTEX_NEIGHBORS(const DTYPE dimension, const ATYPE * axis_size):
      GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>(dimension, axis_size) 
    {
      IJK::PROCEDURE_ERROR error("VERTEX_NEIGHBORS");
      const ATYPE REGION_EDGE_LENGTH = 2;
      //ATYPE coord2[dimension];
	  IJK::ARRAY < ATYPE > coord2 ( dimension );
      SetAllToZero();
      VTYPE * vlist = NULL;

      try {

        const VTYPE ilast = this->NumCubeVertices()-1;
        origin_decrement = this->CubeVertexIncrement(ilast);

        VTYPE numv_in_region;
        compute_num_grid_vertices_in_region
          (dimension, REGION_EDGE_LENGTH, numv_in_region);
        num_neighbors = numv_in_region-1;

        neighbor_increment = new VTYPE[num_neighbors];

        VTYPE * vlist = new VTYPE[numv_in_region];
        get_grid_vertices_in_region(this->dimension, this->axis_size, VTYPE(0),
                                    REGION_EDGE_LENGTH, vlist);

        NTYPE j = 0;
        for (NTYPE i = 0; i < numv_in_region; i++) {
          VTYPE iv = vlist[i];
          if (iv != origin_decrement) {
            if (j >= num_neighbors) {
              error.AddMessage("Programming error. Too many neighbors.");
              throw error;
            }
            neighbor_increment[j] = iv;
            j++;
          }
        }

        if (j != num_neighbors) {
          error.AddMessage("Programming error. Too few neighbors.");
          throw error;
        }

        // Set boundary bits.
        boundary_bits = new BTYPE[num_neighbors];
        for (NTYPE i = 0; i < num_neighbors; i++) {
          boundary_bits[i] = 0;
          VTYPE iv = neighbor_increment[i];
          ComputeCoord(iv, &(coord2[0]));
          BTYPE flag = 1;
          for (DTYPE d = 0; d < dimension; d++) {
            if (coord2[d] == 0) 
              { boundary_bits[i] = boundary_bits[i] | flag; };
            flag = (flag << 1);
            if (coord2[d] == 2) 
              { boundary_bits[i] = boundary_bits[i] | flag; };
            flag = (flag << 1);
          }
        }

      }
      catch(...) {
        delete [] vlist;
        vlist = NULL;
        delete [] neighbor_increment;
        delete [] boundary_bits;
        SetAllToZero();
        throw;
      }

      delete [] vlist;
      vlist = NULL;
    };

    /// Destructor.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE,
              typename BTYPE>
    VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NTYPE,BTYPE>::
    ~VERTEX_NEIGHBORS()
    {
      delete [] neighbor_increment;
      delete [] boundary_bits;

      SetAllToZero();
    }


    // **************************************************
    // Construct join/split tree from contour tree
    // **************************************************

    // local namespace
    namespace {

      inline bool join_compare(const int i, const int j)
      { return(j < i); };

      inline bool split_compare(const int i, const int j)
      { return(j > i); };

      /// Add contour tree vertex to a join/split tree.
      template<typename NUM_TYPE, typename NODE_TYPE, typename VTYPE, 
               typename LTYPE, typename JS_NODE_TYPE, typename COMPARE_FUNC>
      void add_vertex
      (const NUM_TYPE i, const NODE_TYPE * tree,
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       JS_NODE_TYPE * js_tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        VTYPE iv = sorted_vertices[i];
        set.SetData(i, i);
        js_tree[iv].SetParent(NULL);

        if (tree[iv].IsDeleted()) 
          { 
            js_tree[iv].Delete(); 
          }
        else {
          if (!tree[iv].IsRoot()) {
            NUM_TYPE iv2 = ParentIndex(tree, iv);
            NUM_TYPE j = vertex_loc[iv2];

            if (compare(j,i)) {
              js_tree[iv].SetParent(js_tree+iv2);
            }
            else {
              // among vertices less than i, get root of subtree containing j
              VTYPE k = set.Data(j);
              NUM_TYPE iv2 = sorted_vertices[k];
              JS_NODE_TYPE * node = js_tree[iv2].Parent();
              while (node != NULL) {
                iv2 = node - js_tree;
                k = vertex_loc[iv2];
                if (compare(i,k)) 
                  { set.Union(k, j, k); }
                else 
                  { break; }

                node = node->Parent();
              }
              k = set.Data(j);
              NUM_TYPE iv3 = sorted_vertices[k];

              // insert iv above iv3 in tree
              js_tree[iv].SetParent(js_tree[iv3].Parent());
              js_tree[iv3].SetParent(js_tree+iv);
              set.Union(k, i, i);
            }
          }
        }
      }

    }

    /// Construct a join tree from another tree.
    /// @param sorted_vertices[] = Array of grid vertex indices sorted by increasing scalar value.
    /// @param vertex_loc[iv] = Location of vertex \a iv in array \a sorted_vertics[].
    /// @param num_vertices = Number of grid vertices.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename NODE_TYPE, typename VTYPE, typename LTYPE, 
              typename NUM_TYPE, typename JNODE_TYPE>
    void construct_join_tree
    (const NODE_TYPE * tree,
     const VTYPE * sorted_nodes, const LTYPE * node_loc,
     const NUM_TYPE num_nodes, JNODE_TYPE * join_tree)
    {
      IJK::PROCEDURE_ERROR error("construct_join_tree");

      UNION_FIND<NUM_TYPE,NUM_TYPE> set(num_nodes);

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        add_vertex(i, tree, sorted_nodes, node_loc,
                   join_tree, set, join_compare);
      }
    }


    /// Construct a split tree from another tree.
    /// @param sorted_vertices[] = Array of grid vertex indices sorted by increasing scalar value.
    /// @param vertex_loc[iv] = Location of vertex \a iv in array \a sorted_vertics[].
    /// @param num_vertices = Number of grid vertices.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename NODE_TYPE, typename VTYPE, typename LTYPE, 
              typename NUM_TYPE, typename SNODE_TYPE>
    void construct_split_tree
    (const NODE_TYPE * tree,
     const VTYPE * sorted_nodes, const LTYPE * node_loc,
     const NUM_TYPE num_nodes, SNODE_TYPE * split_tree)
    {
      IJK::PROCEDURE_ERROR error("construct_split_tree");

      UNION_FIND<NUM_TYPE,NUM_TYPE> set(num_nodes);

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        add_vertex(num_nodes-1-i, tree, sorted_nodes, node_loc,
                   split_tree, set, split_compare);
      }
    }

    // **************************************************
    // Contour tree main functions
    // **************************************************

    // local namespace
    namespace {

      /// Union sets containing vertices iv1 and iv2 if compare is true.
      template <typename NUM_TYPE, typename VTYPE, typename LTYPE,
                typename NODE_TYPE, typename COMPARE_FUNC>
      void union_vertices
      (const VTYPE iv1, const VTYPE iv2, const VTYPE * sorted_vertices,
       const LTYPE * vertex_loc, NODE_TYPE * tree,
       UNION_FIND<NUM_TYPE,NUM_TYPE> & set, COMPARE_FUNC compare)
      {
        NUM_TYPE i1 = vertex_loc[iv1];
        NUM_TYPE i2 = vertex_loc[iv2];
        if (compare(i1,i2) && set.Find(i1) != set.Find(i2)) {
          VTYPE i3 = set.Data(i2);
          VTYPE iv3 = sorted_vertices[i3];
          tree[iv3].SetParent(tree+iv1);
          set.Union(i1, i2, i1);
        }
      }

      /// Add boundary grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices by grid edges.
      template<typename NUM_TYPE, typename DTYPE, typename ATYPE, typename VTYPE, 
	       typename LTYPE, typename ITYPE, typename NODE_TYPE, 
	       typename COMPARE_FUNC>
      void add_boundary_vertex_edge_connect
      (const NUM_TYPE i, const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       const ITYPE * axis_increment,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        //ATYPE coord[dimension];
		IJK::ARRAY < ATYPE > coord ( dimension );

        VTYPE iv = sorted_vertices[i];

        set.SetData(i, i);
        tree[iv].SetParent(NULL);
        compute_coord(iv, dimension, axis_size, &(coord[0]));

        for (DTYPE d = 0; d < dimension; d++) {
          for (DTYPE dir = 0; dir < 2; dir++) {
            VTYPE iv2 = iv;
            if (dir == 0) {
              if (coord[d] <= 0) { continue; };
              iv2 -= axis_increment[d];
            }
            else {
              if (coord[d]+1 >= axis_size[d]) { continue; };
              iv2 += axis_increment[d];
            }

            union_vertices
              (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
          }
        }
      }

      /// Add grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices by grid edges.
      template<typename NUM_TYPE, typename VTYPE, typename LTYPE, 
               typename ITYPE, typename BGRID_TYPE, typename NODE_TYPE, 
               typename COMPARE_FUNC>
      void add_vertex_edge_connect
      (const NUM_TYPE i, 
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       const BGRID_TYPE & boundary_grid,
       const ITYPE * axis_increment,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        typedef typename BGRID_TYPE::DIMENSION_TYPE DTYPE;
        const DTYPE dimension = boundary_grid.Dimension();

        VTYPE iv = sorted_vertices[i];

        if (boundary_grid.Scalar(iv)) {
          add_boundary_vertex_edge_connect
            (i, dimension, boundary_grid.AxisSize(), 
             sorted_vertices, vertex_loc,
             axis_increment, tree, set, compare);
        }
        else {
          set.SetData(i, i);
          tree[iv].SetParent(NULL);

          for (DTYPE d = 0; d < dimension; d++) {
            for (DTYPE dir = 0; dir < 2; dir++) {

              VTYPE iv2 = iv;
              if (dir == 0) { iv2 -= axis_increment[d]; }
              else { iv2 += axis_increment[d]; }

              union_vertices
                (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
            }
          }
        }

      }

      /// Get grid vertices sharing a grid cube with vertex iv
      /// Note: iv is included in the list.
      /// Precondition: Array vlist[] is preallocated to length
      ///   at least num_vertices.
      template<typename DTYPE, typename ATYPE, typename VTYPE, typename VTYPE2,
               typename NUM_TYPE>
      void get_cube_connected_vertices
      (const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE iv, VTYPE2 * vlist, NUM_TYPE & num_vertices)
      {
        //ATYPE coord[dimension];
        //ATYPE subgrid_origin_coord[dimension];
        //ATYPE subgrid_axis_size[dimension];

		IJK::ARRAY <ATYPE> coord ( dimension );
		IJK::ARRAY <ATYPE> subgrid_origin_coord ( dimension );
		IJK::ARRAY <ATYPE> subgrid_axis_size ( dimension );


        compute_coord(iv, dimension, axis_size, &(coord[0]));

        // compute coordinates of subgrid origin
        for (DTYPE d = 0; d < dimension; d++) {
          if (coord[d] <= 0 )
            { subgrid_origin_coord[d] = coord[d]; }
          else
            { subgrid_origin_coord[d] = coord[d] - 1; };

          if (coord[d]+1 >= axis_size[d]) 
            { subgrid_axis_size[d] = coord[d]+1 - subgrid_origin_coord[d]; }
          else
            { subgrid_axis_size[d] = coord[d]+2 - subgrid_origin_coord[d]; };
        }

        VTYPE subgrid_origin = 
          IJK::compute_vertex_index<VTYPE>
          (&(subgrid_origin_coord[0]), dimension, axis_size);

        IJK::compute_num_grid_vertices
          (dimension, &(subgrid_axis_size[0]), num_vertices);

        IJK::get_subgrid_vertices(dimension, axis_size, subgrid_origin,
                                  &(subgrid_axis_size[0]), vlist);
      }

      /// Add boundary grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices sharing a grid cube.
      template<typename NUM_TYPE, typename DTYPE, typename ATYPE, typename VTYPE,
               typename BTYPE, typename LTYPE, typename NODE_TYPE, 
               typename COMPARE_FUNC>
      void add_boundary_vertex_cube_connect
      (const NUM_TYPE i, 
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NUM_TYPE,BTYPE> & neighbors,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        long boundary_bits;

        VTYPE iv = sorted_vertices[i];
        set.SetData(i, i);
        tree[iv].SetParent(NULL);

        neighbors.SetOrigin(0);
        neighbors.ComputeBoundaryBits(iv, boundary_bits);

        for (NUM_TYPE k = 0; k < neighbors.NumNeighbors(); k++) {
	  
          if (boundary_bits & neighbors.BoundaryBits(k)) { continue; }

          VTYPE iv2 = 
            (iv + neighbors.Neighbor(k)) - neighbors.OriginDecrement();

          union_vertices
            (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
        }
      }

      /// Add grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices sharing a grid cube.
      template<typename NUM_TYPE, typename DTYPE, typename ATYPE, typename VTYPE, 
               typename BGRID_TYPE, typename BBITS_TYPE,
               typename LTYPE, typename NODE_TYPE, 
               typename COMPARE_FUNC>
      void add_vertex_cube_connect
      (const NUM_TYPE i, 
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       const BGRID_TYPE & boundary_grid,
       VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NUM_TYPE,BBITS_TYPE> & neighbors,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        VTYPE iv = sorted_vertices[i];
        if (boundary_grid.Scalar(iv)) {
          add_boundary_vertex_cube_connect
            (i, sorted_vertices, vertex_loc, neighbors,
             tree, set, compare);
        }
        else {

          set.SetData(i, i);
          tree[iv].SetParent(NULL);
          neighbors.SetVertex(iv);

          for (NUM_TYPE k = 0; k < neighbors.NumNeighbors(); k++) {

            VTYPE iv2 = neighbors.Neighbor(k);

            union_vertices
              (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
          }
        }
      }

      /// Get grid vertices sharing a grid cube with vertex iv
      /// Note: iv is included in the list.
      /// Precondition: Array vlist[] is preallocated to length
      ///   at least num_vertices.
      template<typename DTYPE, typename ATYPE, typename VTYPE, typename VTYPE2,
               typename NUM_TYPE>
      void get_facet_connected_vertices
      (const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE iv, VTYPE2 * vlist, NUM_TYPE & num_vertices)
      {
        //DTYPE coord[dimension];
        //DTYPE coord2[dimension];
		IJK::ARRAY <DTYPE> coord (dimension);
		IJK::ARRAY <DTYPE> coord2 (dimension);

        NUM_TYPE vlist2_length;
        IJK::compute_num_grid_vertices_in_region(dimension, 2, vlist2_length);
        IJK::ARRAY<VTYPE2> vlist2(vlist2_length);
        NUM_TYPE num_cube_connected;

        // Time consuming way to get vertices sharing a facet with iv.

        get_cube_connected_vertices
          (dimension, axis_size, iv, vlist2.Ptr(), num_cube_connected);

        compute_coord(iv, dimension, axis_size, &(coord[0]));

        num_vertices = 0;
        for (NUM_TYPE i = 0; i < num_cube_connected; i++) {
          VTYPE2 iv2 = vlist2[i];

          compute_coord(iv2, dimension, axis_size, &(coord2[0]));
          bool flag_shares_facet = false;
          for (DTYPE d = 0; d < dimension; d++) {
            if (coord[d] == coord2[d]) {
              flag_shares_facet = true;
            }
          }

          if (flag_shares_facet) {
            vlist[num_vertices] = vlist2[i];
            num_vertices++;
          }
        }
				    
      }

      /// Add grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices sharing a grid facet.
      template<typename NUM_TYPE, typename DTYPE, typename ATYPE, typename VTYPE, 
               typename LTYPE, typename ITYPE, typename NODE_TYPE, 
               typename COMPARE_FUNC>
      void add_vertex_facet_connect
      (const NUM_TYPE i, const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       const ITYPE * axis_increment,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        VTYPE iv = sorted_vertices[i];
        set.SetData(i, i);
        tree[iv].SetParent(NULL);

        NUM_TYPE vlist_length;
        IJK::compute_num_grid_vertices_in_region(dimension, 2, vlist_length);

        IJK::ARRAY<VTYPE> vlist(vlist_length);
        NUM_TYPE num_connected_vertices;
        get_facet_connected_vertices
          (dimension, axis_size, iv, vlist.Ptr(), num_connected_vertices);

        for (NUM_TYPE i_list = 0; i_list < num_connected_vertices; i_list++) {
          VTYPE iv2 = vlist[i_list];
          if (iv == iv2) { continue; };

          union_vertices
            (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
        }
      }

      /// \brief Get grid vertices sharing a simplex with vertex iv.
      /// Simplices are d! simplices share edge from lowest to highest grid vertex.
      ///
      /// Vertices iv and iv2 lie in a common simplex iff. all the coordinates of iv 
      ///   are less than or equal to all the coordinates of iv2 or vice versa.
      /// Note: iv is included in the list.
      /// Precondition: Array vlist[] is preallocated to length
      ///   at least num_vertices.
      template<typename DTYPE, typename ATYPE, typename VTYPE, typename VTYPE2,
               typename NUM_TYPE>
      void get_simplex_connected_vertices
      (const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE iv, VTYPE2 * vlist, NUM_TYPE & num_vertices)
      {
        //DTYPE coord[dimension];
        //DTYPE coord2[dimension];
		IJK::ARRAY <DTYPE> coord (dimension);
		IJK::ARRAY <DTYPE> coord2 (dimension);

        NUM_TYPE vlist2_length;
        IJK::compute_num_grid_vertices_in_region(dimension, 2, vlist2_length);
        IJK::ARRAY<VTYPE2> vlist2(vlist2_length);
        NUM_TYPE num_cube_connected;

        // Time consuming way to get vertices sharing a simplex with iv

        get_cube_connected_vertices
          (dimension, axis_size, iv, vlist2.Ptr(), num_cube_connected);

        compute_coord(iv, dimension, axis_size, &(coord[0]));

        num_vertices = 0;
        for (NUM_TYPE i = 0; i < num_cube_connected; i++) {
          VTYPE2 iv2 = vlist2[i];

          compute_coord(iv2, dimension, axis_size, &(coord2[0]));
          bool flag_less = false;
          bool flag_greater = false;
          for (DTYPE d = 0; d < dimension; d++) {
            if (coord[d] < coord2[d]) {
              flag_less = true;
            }
            else if (coord[d] > coord2[d]) {
              flag_greater = true;
            }
          }

          if (!flag_less || !flag_greater) {
            vlist[num_vertices] = vlist2[i];
            num_vertices++;
          }
        }
      }

      /// Add grid vertex to a join/split tree.
      /// Grid vertex is connected to other vertices sharing a simplex.
      template<typename NUM_TYPE, typename DTYPE, typename ATYPE, typename VTYPE, 
               typename LTYPE, typename ITYPE, typename NODE_TYPE, 
               typename COMPARE_FUNC>
      void add_vertex_simplex_connect
      (const NUM_TYPE i, const DTYPE dimension, const ATYPE * axis_size,
       const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
       const ITYPE * axis_increment,
       NODE_TYPE * tree, UNION_FIND<NUM_TYPE,NUM_TYPE> & set,
       COMPARE_FUNC compare)
      {
        VTYPE iv = sorted_vertices[i];
        set.SetData(i, i);
        tree[iv].SetParent(NULL);

        NUM_TYPE vlist_length;
        IJK::compute_num_grid_vertices_in_region(dimension, 2, vlist_length);

        IJK::ARRAY<VTYPE> vlist(vlist_length);
        NUM_TYPE num_connected_vertices;
        get_simplex_connected_vertices
          (dimension, axis_size, iv, vlist.Ptr(), num_connected_vertices);

        for (NUM_TYPE i_list = 0; i_list < num_connected_vertices; i_list++) {
          VTYPE iv2 = vlist[i_list];
          if (iv == iv2) { continue; };

          union_vertices
            (iv, iv2, sorted_vertices, vertex_loc, tree, set, compare);
        }
      }

    }


    /// Construct a join tree from a regular grid.
    /// @param dimension = Grid dimension.
    /// @param axis_size[d] = Number of vertices along grid axis d.
    /// @param sorted_vertices[] = Array of grid vertex indices sorted by increasing scalar value.
    /// @param vertex_loc[iv] = Location of vertex \a iv in array \a sorted_vertics[].
    /// @param num_vertices = Number of grid vertices.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename LTYPE, 
              typename NUM_TYPE, typename NODE_TYPE>
    void construct_join_tree
    (const DTYPE dimension, const ATYPE * axis_size,
     const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
     const NUM_TYPE num_vertices, const CONNECTIVITY connectivity,
     NODE_TYPE * tree)
    {
      //VTYPE axis_increment[dimension];
	  IJK::ARRAY <VTYPE> axis_increment ( dimension );
      VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NUM_TYPE,long> 
        neighbors(dimension, axis_size);
      BOOL_GRID<GRID<DTYPE,ATYPE,VTYPE,NUM_TYPE> > 
        boundary_grid(dimension, axis_size);
      IJK::PROCEDURE_ERROR error("construct_join_tree");

      if (!check_num_vertices(dimension, axis_size, num_vertices, error))
        { throw error; };

      UNION_FIND<NUM_TYPE,NUM_TYPE> set(num_vertices);

      IJK::compute_increment(dimension, axis_size, &(axis_increment[0]));

      compute_boundary_grid(boundary_grid);

      switch(connectivity) {

      case(EDGE_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_edge_connect
            (i, sorted_vertices, vertex_loc, boundary_grid, &(axis_increment[0]),
             tree, set, join_compare);
        }
        break;

      case(CUBE_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_cube_connect
            (i, sorted_vertices, vertex_loc, boundary_grid, neighbors,
             tree, set, join_compare);
        }
        break;

      case(FACET_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_facet_connect
            (i, dimension, axis_size, sorted_vertices, vertex_loc,
             &(axis_increment[0]), tree, set, join_compare);
        }
        break;

      case(SIMPLEX_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_simplex_connect
            (i, dimension, axis_size, sorted_vertices, vertex_loc,
             &(axis_increment[0]), tree, set, join_compare);
        }
        break;

      default:
        error.AddMessage("Programming error.  Unknown vertex connectivity.");
        throw error;
        break;
      }
    }

    /// Construct a split tree from a regular grid.
    /// @param dimension = Grid dimension.
    /// @param axis_size[d] = Number of vertices along grid axis d.
    /// @param sorted_vertices[] = Array of grid vertex indices sorted by increasing scalar value.
    /// @param vertex_loc[iv] = Location of vertex \a iv in array \a sorted_vertics[].
    /// @param num_vertices = Number of grid vertices.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename DTYPE, typename ATYPE, typename VTYPE, typename LTYPE, 
              typename NUM_TYPE, typename NODE_TYPE>
    void construct_split_tree
    (const DTYPE dimension, const ATYPE * axis_size,
     const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
     const NUM_TYPE num_vertices, const CONNECTIVITY connectivity,
     NODE_TYPE * tree)
    {
      //VTYPE axis_increment[dimension];
      IJK:ARRAY <VTYPE> axis_increment ( dimension );

	  VERTEX_NEIGHBORS<DTYPE,ATYPE,VTYPE,NUM_TYPE,long> 
        neighbors(dimension, axis_size);
      BOOL_GRID<GRID<DTYPE,ATYPE,VTYPE,NUM_TYPE> > 
        boundary_grid(dimension, axis_size);
      IJK::PROCEDURE_ERROR error("construct_split_tree");

      if (!check_num_vertices(dimension, axis_size, num_vertices, error))
        { throw error; };

      UNION_FIND<NUM_TYPE,NUM_TYPE> set(num_vertices);

      IJK::compute_increment(dimension, axis_size, &(axis_increment[0]));

      compute_boundary_grid(boundary_grid);

      switch(connectivity) {

      case(EDGE_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_edge_connect
            (num_vertices-1-i, sorted_vertices, vertex_loc, boundary_grid,
             &(axis_increment[0]), tree, set, split_compare);
        }
        break;

      case(CUBE_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_cube_connect
            (num_vertices-1-i, sorted_vertices, vertex_loc, 
             boundary_grid, neighbors,
             tree, set, split_compare);
        }
        break;


      case(FACET_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_facet_connect
            (num_vertices-1-i, dimension, axis_size, 
             sorted_vertices, vertex_loc,
             &(axis_increment[0]), tree, set, split_compare);
        }
        break;

      case(SIMPLEX_CONNECT):
        for (NUM_TYPE i = 0; i < num_vertices; i++) {
          add_vertex_simplex_connect
            (num_vertices-1-i, dimension, axis_size, 
             sorted_vertices, vertex_loc,
             &(axis_increment[0]), tree, set, split_compare);
        }
        break;

      default:
        error.AddMessage("Programming error.  Unknown vertex connectivity.");
        throw error;
        break;
      }
    }

    template <typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE,
              typename NUM_TYPE>
    void construct_augmented_contour_tree
    (const JNODE_TYPE * join_tree, const SNODE_TYPE * split_tree,
     const NUM_TYPE num_nodes, CNODE_TYPE * contour_tree)
    {
      IJK::PROCEDURE_ERROR error("construct_augmented_contour_tree");

      IJK::ARRAY< TREE_NODE_EXT<int> > join_tree_x(num_nodes);
      IJK::ARRAY< TREE_NODE_EXT<int> > split_tree_x(num_nodes);

      copy_tree(join_tree, join_tree_x.Ptr(), num_nodes);
      copy_tree(split_tree, split_tree_x.Ptr(), num_nodes);

      compute_num_children(join_tree_x.Ptr(), num_nodes);
      compute_num_children(split_tree_x.Ptr(), num_nodes);

      std::vector<NUM_TYPE> list;

      for (NUM_TYPE i = 0; i < num_nodes; i++) {

        if (!join_tree_x[i].IsDeleted()) {

          if (split_tree_x[i].IsDeleted()) { 
            error.AddMessage
              ("Programming error. Vertex ", i,
               " deleted from split tree but not from join tree.");
            throw error; 
          };

          if (join_tree_x[i].NumChildren() + split_tree_x[i].NumChildren() == 1)
            { list.push_back(i); };
        }
        else {
          if (!split_tree_x[i].IsDeleted()) { 
            error.AddMessage
              ("Programming error. Vertex ", i,
               " deleted from join tree but not from split tree.");
            throw error; 
          };
        }
      }

      while (list.size() > 1) {
        NUM_TYPE i = list.back();
        list.pop_back();

        NUM_TYPE k;
        if (join_tree_x[i].NumChildren() == 0) {
          join_tree_x[i].ParentCompress();
          if (join_tree_x[i].IsRoot()) { 
            error.AddMessage("Programming error. Vertex ", i,
                             " is root of join tree.");
            throw error; 
          };

          k = ParentIndex(join_tree_x, i);
          contour_tree[i].SetParent(contour_tree+k);
          contour_tree[i].SetJoinEdge();
        }
        else {
          split_tree_x[i].ParentCompress();
          if (split_tree_x[i].IsRoot()) { 
            error.AddMessage("Programming error. Vertex ", i,
                             " is root of split tree.");
            throw error; 
          };

          k = ParentIndex(split_tree_x, i);
          contour_tree[i].SetParent(contour_tree+k);
          contour_tree[i].SetSplitEdge();
        }

        join_tree_x[i].Delete();
        split_tree_x[i].Delete();

        if (join_tree_x[k].NumChildren() + split_tree_x[k].NumChildren() == 1)
          { list.push_back(k); };
      }

    };

    template <typename DTYPE, typename ATYPE, typename NUM_TYPE, typename VTYPE, typename LTYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_augmented_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const NUM_TYPE num_vertices,
     const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
     JNODE_TYPE * join_tree, const CONNECTIVITY join_connectivity,
     SNODE_TYPE * split_tree, const CONNECTIVITY split_connectivity,
     CNODE_TYPE * contour_tree)
    {
      IJK::PROCEDURE_ERROR error("construct_augmented_contour_tree");

      if (!check_num_vertices(dimension, axis_size, num_vertices, error))
        { throw error; };


      construct_join_tree(dimension, axis_size, sorted_vertices, vertex_loc,
                          num_vertices, join_connectivity, join_tree);

      construct_split_tree(dimension, axis_size, sorted_vertices, vertex_loc,
                           num_vertices, split_connectivity, split_tree);

      construct_augmented_contour_tree
        (join_tree, split_tree, num_vertices, contour_tree);
    }

    template <typename DTYPE, typename ATYPE, typename NUM_TYPE, typename VTYPE, typename LTYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_augmented_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const NUM_TYPE num_vertices,
     const VTYPE * sorted_vertices, const LTYPE * vertex_loc,
     std::vector<JNODE_TYPE> & join_tree, 
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree)
    {
      join_tree.resize(num_vertices);
      split_tree.resize(num_vertices);
      contour_tree.resize(num_vertices);

      construct_augmented_contour_tree
        (dimension, axis_size, num_vertices, sorted_vertices, vertex_loc,
         &(join_tree[0]), join_connectivity,
         &(split_tree[0]), split_connectivity,
         &(contour_tree[0]));
    }

    template <typename DTYPE, typename ATYPE, typename STYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_augmented_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     JNODE_TYPE * join_tree, const CONNECTIVITY join_connectivity,
     SNODE_TYPE * split_tree, const CONNECTIVITY split_connectivity,
     CNODE_TYPE * contour_tree)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices = 
        IJK::compute_num_grid_vertices(dimension, axis_size);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      construct_augmented_contour_tree
        (dimension, axis_size, num_vertices,
         sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
         join_tree, join_connectivity, split_tree, split_connectivity,
         contour_tree);
    }

    template <typename DTYPE, typename ATYPE, typename STYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_augmented_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     IJK::ARRAY<JNODE_TYPE> & join_tree, 
     const CONNECTIVITY join_connectivity,
     IJK::ARRAY<SNODE_TYPE> & split_tree, 
     const CONNECTIVITY split_connectivity,
     IJK::ARRAY<CNODE_TYPE> & contour_tree)
    {
      construct_augmented_contour_tree
        (dimension, axis_size, scalar, 
         join_tree.Ptr(), join_connectivity,
         split_tree.Ptr(), split_connectivity, contour_tree.Ptr());
    }

    template <typename DTYPE, typename ATYPE, typename STYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_augmented_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     std::vector<JNODE_TYPE> & join_tree,
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices = 
        IJK::compute_num_grid_vertices(dimension, axis_size);

      join_tree.resize(num_vertices);
      split_tree.resize(num_vertices);
      contour_tree.resize(num_vertices);

      construct_augmented_contour_tree
        (dimension, axis_size, scalar, 
         &(join_tree[0]), join_connectivity,
         &(split_tree[0]), split_connectivity, &(contour_tree[0]));
    }


    /// Delete root from tree making some child of the root, the new root
    /// Precondition: Tree has only one root
    template <typename NODE_TYPE, typename NUM_TYPE>
    void delete_root(NODE_TYPE * tree, const NUM_TYPE num_nodes)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree[i].IsDeleted() && !tree[i].IsRoot()) {
          if (tree[i].Parent()->IsRoot()) {
            NODE_TYPE * old_root = tree[i].Parent();
            tree[i].SetParent(NULL);
            old_root->SetParent(tree+i);
            old_root->Delete();
            return;
          }
        }
      }

      // tree has only root.  Delete it.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree[i].IsDeleted() && tree[i].IsRoot()) {
          tree[i].Delete();
          return;
        }
      }
    }

    /// Delete nodes with one lower and one higher neighbor from contour tree
    template <typename CNODE_TYPE, typename NUM_TYPE>
    void delete_1_1_nodes_from_contour_tree
    (CNODE_TYPE * contour_tree, const NUM_TYPE num_nodes)
    {
      IJK::ARRAY<NUM_TYPE> num_join_children(num_nodes);
      IJK::ARRAY<NUM_TYPE> num_split_children(num_nodes);

      compute_num_children(contour_tree, num_nodes, JOIN_EDGE, 
                           num_join_children.Ptr());
      compute_num_children(contour_tree, num_nodes, SPLIT_EDGE, 
                           num_split_children.Ptr());

      bool flag_delete_root = false;
      NUM_TYPE index_root = 0;
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!contour_tree[i].IsDeleted()) {
          if (num_join_children[i] == 1 && num_split_children[i] == 1)
            if (!contour_tree[i].IsRoot())
              { contour_tree[i].Delete(); }
            else {
              // node i is root of contour tree
              flag_delete_root = true;
            }
        }
      }

      if (flag_delete_root) { delete_root(contour_tree, num_nodes); }

      compress_tree(contour_tree, num_nodes);
    }

   
    template <typename STYPE, typename NODE_TYPE, typename NUM_TYPE>
    void merge_nodes_with_identical_values
    (const STYPE * scalar, NODE_TYPE * tree, const NUM_TYPE num_nodes)
    {
      typedef NODE_TYPE * NODE_PTR;

      for (NUM_TYPE i = 0; i < num_nodes; i++) {

        while (!tree[i].IsRoot() && !tree[i].IsDeleted()) {

          NUM_TYPE j = ParentCompressIndex(tree, i);

          if (scalar[i] == scalar[j]) {

            j = ParentCompressIndex(tree, i);

            if (i >= j) {
              tree[i].Delete();
            }
            else {
              // Swap locations of tree[i] and tree[j] in tree
              tree[i].SetParent(NODE_PTR(tree[j].ParentCompress()));
              tree[j].SetParent(tree+i);

              tree[j].Delete();
            }
          }
          else {
            break;
          }
        }
      }

      compress_tree(tree, num_nodes);
    }

    template <typename STYPE, typename NODE_TYPE, typename NUM_TYPE>
    void merge_weighted_nodes_with_identical_values
    (const STYPE * scalar, NODE_TYPE * tree, const NUM_TYPE num_nodes)
    {
      typedef NODE_TYPE * NODE_PTR;

      for (NUM_TYPE i = 0; i < num_nodes; i++) {

        while (!tree[i].IsRoot() && !tree[i].IsDeleted()) {

          NUM_TYPE j = ParentCompressIndex(tree, i);

          if (scalar[i] == scalar[j]) {

            j = ParentCompressIndex(tree, i);

            typename NODE_TYPE::WEIGHT_TYPE w;
            w = std::max(tree[i].Weight(), tree[j].Weight());

            if (i >= j) {
              tree[i].Delete();
              tree[j].SetWeight(w);
            }
            else {
              // Swap locations of tree[i] and tree[j] in tree
              tree[i].SetParent(NODE_PTR(tree[j].ParentCompress()));
              tree[j].SetParent(tree+i);

              tree[j].Delete();
              tree[i].SetWeight(w);
            }
          }
          else {
            break;
          }
        }
      }

      compress_tree(tree, num_nodes);
    }

    /// Delete nodes with one parent and one child from tree
    template <typename TYPEA, typename NUM_TYPE>
    void delete_1_1_nodes(TYPEA * treeA, const NUM_TYPE num_nodes)
    {
      IJK::ARRAY<NUM_TYPE> num_children(num_nodes);
      IJK::ARRAY<bool> is_root(num_nodes);

      compute_num_children(treeA, num_nodes, num_children.Ptr());
      compute_is_root(treeA, num_nodes, is_root.Ptr());
    
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!is_root[i] && num_children[i] == 1)
          { treeA[i].Delete(); }
      }

      for (NUM_TYPE i = 0; i < num_nodes; i++)
        if (!(treeA[i].IsDeleted())) 
          { treeA[i].ParentCompress(); };
    }

    // local namespace
    namespace {

      /// Copy edge type to nodeB is of type CONTOUR_TREE_NODE<TYPEB>
      template <typename TYPEA, typename TYPEB>
      void copy_edge_type(const TYPEA & nodeA, 
                          CONTOUR_TREE_NODE<TYPEB> & nodeB)
      { nodeB.SetEdgeType(nodeA.EdgeType()); };

      /// Do nothing if nodeB is not of type CONTOUR_TREE_NODE<TYPEB>
      template <typename TYPEA, typename TYPEB>
      void copy_edge_type(const TYPEA & nodeA, TYPEB & nodeB)
      {};
    }

    template <typename TYPEA, typename NUM_TYPE, typename TYPEB>
    void copy_tree_skip_deleted_nodes
    (const TYPEA * treeA, const NUM_TYPE num_nodes, std::vector<TYPEB> & treeB)
    {
      treeB.clear();

      IJK::ARRAY<NUM_TYPE> index_treeB(num_nodes);

      // copy nodes from treeA to treeB.  Set node identifiers.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!(treeA[i].IsDeleted())) {
          TYPEB nodeB;
          NUM_TYPE k = treeB.size();
          treeB.push_back(nodeB);
          treeB.back().SetIdent(i);
          copy_edge_type(treeA[i], treeB.back());
          index_treeB[i] = k;
        }
        else { index_treeB[i] = 0; };
      }

      // Set parents.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!treeA[i].IsDeleted()) {
          if (!treeA[i].IsRoot()) {
            NUM_TYPE k = index_treeB[i];
            NUM_TYPE index_parent = index_treeB[ParentIndex(treeA, i)];
            TYPEB * nodeB_parent = &(treeB[index_parent]);
            treeB[k].SetParent(nodeB_parent);
          }
        }
      }

    }

    template <typename TYPEA, typename NUM_TYPE, typename TYPEB>
    void copy_weighted_tree_skip_deleted_nodes
    (const TYPEA * treeA, const NUM_TYPE num_nodes, std::vector<TYPEB> & treeB)
    {
      treeB.clear();

      IJK::ARRAY<NUM_TYPE> index_treeB(num_nodes);

      // copy nodes from treeA to treeB.  Set node identifiers.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!(treeA[i].IsDeleted())) {
          TYPEB nodeB;
          NUM_TYPE k = treeB.size();
          treeB.push_back(nodeB);
          treeB.back().SetIdent(i);
          treeB.back().SetWeight(treeA[i].Weight());
          copy_edge_type(treeA[i], treeB.back());
          index_treeB[i] = k;
        }
        else { index_treeB[i] = 0; };
      }

      // Set parents.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!treeA[i].IsDeleted()) {
          if (!treeA[i].IsRoot()) {
            NUM_TYPE k = index_treeB[i];
            NUM_TYPE index_parent = index_treeB[ParentIndex(treeA, i)];
            TYPEB * nodeB_parent = &(treeB[index_parent]);
            treeB[k].SetParent(nodeB_parent);
          }
        }
      }

    }

    /// Construct a join tree from a regular grid.
    /// @param dimension = Grid dimension.
    /// @param axis_size[d] = Number of vertices along grid axis d.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename DTYPE, typename ATYPE, typename STYPE, typename NODE_TYPE>
    void construct_join_tree
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const CONNECTIVITY connectivity, std::vector<NODE_TYPE> & join_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      if (augment_flag == AUGMENT_ALL && 
          merge_ident_flag == MERGE_IDENT_NONE) {

        join_tree.resize(num_vertices);

        construct_join_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, connectivity, &(join_tree[0]));

        set_node_identifiers(join_tree);

        return;
      }
      else {

        IJK::ARRAY<JOIN_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_jtree(num_vertices);

        construct_join_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, connectivity, augmented_jtree.Ptr());

        if (merge_ident_flag == MERGE_IDENT_ALL) {
          merge_nodes_with_identical_values
            (scalar, augmented_jtree.Ptr(), num_vertices);
        }

        if (augment_flag == AUGMENT_NONE) {
          delete_1_1_nodes(augmented_jtree.Ptr(), num_vertices);
        }

        copy_tree_skip_deleted_nodes
          (augmented_jtree.Ptr(), num_vertices, join_tree);
      }
    }

    /// Construct a join tree from a regular grid.
    template <typename GTYPE, typename NODE_TYPE>
    void construct_join_tree
    (const GTYPE & grid,  std::vector<NODE_TYPE> & join_tree,
     const CONNECTIVITY connectivity, 
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_join_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         connectivity, join_tree, augment_flag, merge_ident_flag);
    }

    /// Construct a split tree from a regular grid.
    /// @param dimension = Grid dimension.
    /// @param axis_size[d] = Number of vertices along grid axis d.
    /// @param sorted_vertices[] = Array of grid vertex indices sorted by increasing scalar value.
    /// @param vertex_loc[iv] = Location of vertex \a iv in array \a sorted_vertics[].
    /// @param num_vertices = Number of grid vertices.
    /// @pre \li Array tree[] is preallocated to length at least \a num_vertices.
    template <typename DTYPE, typename ATYPE, typename STYPE, typename NODE_TYPE>
    void construct_split_tree
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const CONNECTIVITY connectivity, std::vector<NODE_TYPE> & split_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      if (augment_flag == AUGMENT_ALL && 
          merge_ident_flag == MERGE_IDENT_NONE) {

        split_tree.resize(num_vertices);

        construct_split_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, connectivity, &(split_tree[0]));

        set_node_identifiers(split_tree);

        return;
      }
      else {

        IJK::ARRAY<SPLIT_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_stree(num_vertices);

        construct_split_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, connectivity, augmented_stree.Ptr());

        if (merge_ident_flag == MERGE_IDENT_ALL) {
          merge_nodes_with_identical_values
            (scalar, augmented_stree.Ptr(), num_vertices);
        }

        if (augment_flag == AUGMENT_NONE) {
          delete_1_1_nodes(augmented_stree.Ptr(), num_vertices);
        }

        copy_tree_skip_deleted_nodes
          (augmented_stree.Ptr(), num_vertices, split_tree);
      }
    }

    /// Construct a split tree from a regular grid.
    template <typename GTYPE, typename NODE_TYPE>
    void construct_split_tree
    (const GTYPE & grid,  std::vector<NODE_TYPE> & split_tree,
     const CONNECTIVITY connectivity, 
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_split_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         connectivity, split_tree, augment_flag, merge_ident_flag);
    }

    /// Construct join and split tree.
    /// @param JNODE_TYPE Join tree node type.  Must include member function Ident() and SetIdent().
    /// @param SNODE_TYPE Split tree node type.  Must include member function Ident() and SetIdent().
    template <typename DTYPE, typename ATYPE, typename STYPE,
              typename JNODE_TYPE, typename SNODE_TYPE>
    void construct_join_split_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     std::vector<JNODE_TYPE> & join_tree,
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      if (augment_flag == AUGMENT_ALL && 
          merge_ident_flag == MERGE_IDENT_NONE) {

        join_tree.resize(num_vertices);
        split_tree.resize(num_vertices);

        construct_join_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, join_connectivity, &(join_tree[0]));

        construct_split_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, split_connectivity, &(split_tree[0]));

        set_node_identifiers(join_tree);
        set_node_identifiers(split_tree);

        return;
      }
      else {

        IJK::ARRAY<JOIN_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_jtree(num_vertices);
        IJK::ARRAY<SPLIT_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_stree(num_vertices);

        construct_join_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, join_connectivity, augmented_jtree.Ptr());

        construct_split_tree
          (dimension, axis_size, 
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           num_vertices, split_connectivity, augmented_stree.Ptr());


        if (merge_ident_flag == MERGE_IDENT_ALL) {
          merge_nodes_with_identical_values
            (scalar, augmented_jtree.Ptr(), num_vertices);
          merge_nodes_with_identical_values
            (scalar, augmented_stree.Ptr(), num_vertices);
        }

        if (augment_flag == AUGMENT_NONE) {
          delete_1_1_nodes(augmented_jtree.Ptr(), num_vertices);
          delete_1_1_nodes(augmented_stree.Ptr(), num_vertices);
        }

        copy_tree_skip_deleted_nodes
          (augmented_jtree.Ptr(), num_vertices, join_tree);
        copy_tree_skip_deleted_nodes
          (augmented_stree.Ptr(), num_vertices, split_tree);
      }
    }

    /// Construct join and split tree
    /// @param JNODE_TYPE Join tree node type.  Must include member function Ident() and SetIdent().
    /// @param SNODE_TYPE Split tree node type.  Must include member function Ident() and SetIdent().
    template <typename GTYPE, typename JNODE_TYPE, typename SNODE_TYPE>
    void construct_join_split_tree
    (const GTYPE & grid,
     std::vector<JNODE_TYPE> & join_tree,
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_join_split_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         join_tree, join_connectivity, split_tree, split_connectivity,
         augment_flag, merge_ident_flag);
    }

    /// Construct a contour tree.
    /// @param JNODE_TYPE Join tree node type.  Must include member function Ident() and SetIdent().
    /// @param SNODE_TYPE Split tree node type.  Must include member function Ident() and SetIdent().
    /// @param CNODE_TYPE Contour tree node type.  Must include member function Ident() and SetIdent().
    template <typename DTYPE, typename ATYPE, typename STYPE,
              typename JNODE_TYPE, typename SNODE_TYPE, typename CNODE_TYPE>
    void construct_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     std::vector<JNODE_TYPE> & join_tree,
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      if (augment_flag == AUGMENT_ALL && merge_ident_flag == MERGE_IDENT_NONE) {

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           join_tree, join_connectivity, 
           split_tree, split_connectivity, contour_tree);

        set_node_identifiers(join_tree);
        set_node_identifiers(split_tree);
        set_node_identifiers(contour_tree);

        return;
      }
      else {
        IJK::ARRAY<JOIN_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_jtree(num_vertices);
        IJK::ARRAY<SPLIT_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_stree(num_vertices);
        IJK::ARRAY<CONTOUR_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_ctree(num_vertices);

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           augmented_jtree.Ptr(), join_connectivity,
           augmented_stree.Ptr(), split_connectivity,
           augmented_ctree.Ptr());

        if (merge_ident_flag != MERGE_IDENT_NONE) {
          merge_nodes_with_identical_values
            (scalar, augmented_ctree.Ptr(), num_vertices);
        }

        if (augment_flag == AUGMENT_JS) {
          // Join and split trees will have exact same nodes as contour tree
          //   (modulo application of merging nodes with identical values)
          delete_1_1_nodes_from_contour_tree
            (augmented_ctree.Ptr(), num_vertices);

          // recompute augmented_jtree and augmented_ctree
          construct_join_tree
            (augmented_ctree.Ptr(), 
             sorted_vertices.PtrConst(), vertex_loc.PtrConst(), num_vertices,
             augmented_jtree.Ptr());

          construct_split_tree
            (augmented_ctree.Ptr(), 
             sorted_vertices.PtrConst(), vertex_loc.PtrConst(), num_vertices,
             augmented_stree.Ptr());
        };


        if (merge_ident_flag == MERGE_IDENT_ALL) {
          merge_nodes_with_identical_values
            (scalar, augmented_jtree.Ptr(), num_vertices);
          merge_nodes_with_identical_values
            (scalar, augmented_stree.Ptr(), num_vertices);
        }

        if (augment_flag == AUGMENT_NONE) {
          delete_1_1_nodes_from_contour_tree
            (augmented_ctree.Ptr(), num_vertices);
          delete_1_1_nodes(augmented_jtree.Ptr(), num_vertices);
          delete_1_1_nodes(augmented_stree.Ptr(), num_vertices);
        }

        copy_tree_skip_deleted_nodes
          (augmented_jtree.Ptr(), num_vertices, join_tree);
        copy_tree_skip_deleted_nodes
          (augmented_stree.Ptr(), num_vertices, split_tree);
        copy_tree_skip_deleted_nodes
          (augmented_ctree.Ptr(), num_vertices, contour_tree);
      }
    }

    /// Construct a contour tree.
    /// @param JNODE_TYPE Join tree node type.  Must include member function Ident() and SetIdent().
    /// @param SNODE_TYPE Split tree node type.  Must include member function Ident() and SetIdent().
    /// @param CNODE_TYPE Contour tree node type.  Must include member function Ident() and SetIdent().
    template <typename GTYPE, typename JNODE_TYPE, typename SNODE_TYPE, 
              typename CNODE_TYPE>
    void construct_contour_tree
    (const GTYPE & grid,
     std::vector<JNODE_TYPE> & join_tree,
     const CONNECTIVITY join_connectivity,
     std::vector<SNODE_TYPE> & split_tree,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_contour_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         join_tree, join_connectivity,
         split_tree, split_connectivity,
         contour_tree, augment_flag, merge_ident_flag);
    }

    /// Construct a contour tree.
    /// @param CNODE_TYPE Contour tree node type.  Must include member function Ident() and SetIdent().
    template <typename DTYPE, typename ATYPE, typename STYPE, typename CNODE_TYPE>
    void construct_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     const CONNECTIVITY join_connectivity,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      IJK::ARRAY<JOIN_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
        augmented_jtree(num_vertices);
      IJK::ARRAY<SPLIT_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
        augmented_stree(num_vertices);

      if (augment_flag == AUGMENT_ALL && 
          merge_ident_flag == MERGE_IDENT_NONE) {

        contour_tree.resize(num_vertices);

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           augmented_jtree.Ptr(), join_connectivity, 
           augmented_stree.Ptr(), split_connectivity, 
           &(contour_tree[0]));

        set_node_identifiers(contour_tree);

        return;
      }
      else {
        IJK::ARRAY<CONTOUR_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
          augmented_ctree(num_vertices);

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           augmented_jtree.Ptr(), join_connectivity,
           augmented_stree.Ptr(), split_connectivity,
           augmented_ctree.Ptr());

        if (merge_ident_flag != MERGE_IDENT_NONE) {
          merge_nodes_with_identical_values
            (scalar, augmented_ctree.Ptr(), num_vertices);
        }

        if (augment_flag != AUGMENT_ALL) {
          delete_1_1_nodes_from_contour_tree
            (augmented_ctree.Ptr(), num_vertices);
        }

        copy_tree_skip_deleted_nodes
          (augmented_ctree.Ptr(), num_vertices, contour_tree);
      }
    }

    /// Construct a contour tree.
    /// @param CNODE_TYPE Contour tree node type.  
    /// Must include member function Ident() and SetIdent().
    template <typename GTYPE, typename CNODE_TYPE>
    void construct_contour_tree
    (const GTYPE & grid,
     const CONNECTIVITY join_connectivity,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_contour_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         join_connectivity, split_connectivity,
         contour_tree, augment_flag, merge_ident_flag);
    }

    // Forward declarations.
    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_subtree_sizes
    (const NODE_TYPE * tree, const NUM_TYPE num_nodes, NUM_TYPE * subtree_size);
    template <typename WTYPE, typename NUM_TYPE, typename NODE_TYPE>
    void compute_ctree_weights
    (const WTYPE * tree1_weight, const WTYPE * tree2_weight,
     const NUM_TYPE num_nodes, NODE_TYPE * ctree);


    /// Construct a contour tree with weights at nodes.
    /// @param CNODE_TYPE Contour tree node type.  
    ///   Must include member function Ident(), SetIdent(),
    ///     Weight() and SetWeight().
    template <typename DTYPE, typename ATYPE, typename STYPE, 
              typename CNODE_TYPE>
    void construct_weighted_contour_tree
    (const DTYPE dimension, const ATYPE axis_size, const STYPE * scalar,
     const CONNECTIVITY join_connectivity,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      typedef size_t NUM_TYPE;
      typedef typename CNODE_TYPE::WEIGHT_TYPE WTYPE;

      NUM_TYPE num_vertices;
      IJK::compute_num_grid_vertices(dimension, axis_size, num_vertices);

      IJK::ARRAY<NUM_TYPE> sorted_vertices(num_vertices);
      IJK::ARRAY<NUM_TYPE> vertex_loc(num_vertices);
      IJK::ARRAY<NUM_TYPE> jtree_subtree_size(num_vertices);
      IJK::ARRAY<NUM_TYPE> stree_subtree_size(num_vertices);

      sort_grid_vertices(scalar, num_vertices, sorted_vertices.Ptr());
      list_locate(sorted_vertices.PtrConst(), num_vertices, vertex_loc.Ptr());

      IJK::ARRAY<JOIN_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
        augmented_jtree(num_vertices);
      IJK::ARRAY<SPLIT_TREE_NODE<TREE_NODE_DEL<TREE_NODE> > > 
        augmented_stree(num_vertices);

      if (augment_flag == AUGMENT_ALL && 
          merge_ident_flag == MERGE_IDENT_NONE) {

        contour_tree.resize(num_vertices);

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           augmented_jtree.Ptr(), join_connectivity, 
           augmented_stree.Ptr(), split_connectivity, 
           &(contour_tree[0]));

        set_node_identifiers(contour_tree);

        return;
      }
      else {
        IJK::ARRAY<CONTOUR_TREE_NODE<WEIGHTED_TREE_NODE_DEL<TREE_NODE,WTYPE> > > 
          augmented_ctree(num_vertices);

        construct_augmented_contour_tree
          (dimension, axis_size, num_vertices,
           sorted_vertices.PtrConst(), vertex_loc.PtrConst(),
           augmented_jtree.Ptr(), join_connectivity,
           augmented_stree.Ptr(), split_connectivity,
           augmented_ctree.Ptr());

        compute_subtree_sizes
          (augmented_jtree.PtrConst(), num_vertices, jtree_subtree_size.Ptr());
        compute_subtree_sizes
          (augmented_stree.PtrConst(), num_vertices, stree_subtree_size.Ptr());
        compute_ctree_weights
          (jtree_subtree_size.PtrConst(), stree_subtree_size.PtrConst(),
           num_vertices, augmented_ctree.Ptr());

        if (merge_ident_flag != MERGE_IDENT_NONE) {
          merge_weighted_nodes_with_identical_values
            (scalar, augmented_ctree.Ptr(), num_vertices);
        }

        if (augment_flag != AUGMENT_ALL) {
          delete_1_1_nodes_from_contour_tree
            (augmented_ctree.Ptr(), num_vertices);
        }

        copy_weighted_tree_skip_deleted_nodes
          (augmented_ctree.Ptr(), num_vertices, contour_tree);
      }
    }

    /// Construct a weighted contour tree.
    /// @param CNODE_TYPE Contour tree node type.  
    ///   Must include member function Ident(), SetIdent(),
    ///     Weight() and SetWeight().
    template <typename GTYPE, typename CNODE_TYPE>
    void construct_weighted_contour_tree
    (const GTYPE & grid,
     const CONNECTIVITY join_connectivity,
     const CONNECTIVITY split_connectivity,
     std::vector<CNODE_TYPE> & contour_tree,
     const AUGMENT_FLAG augment_flag = AUGMENT_NONE,
     const MERGE_IDENT_FLAG merge_ident_flag = MERGE_IDENT_CTREE)
    {
      construct_weighted_contour_tree
        (grid.Dimension(), grid.AxisSize(), grid.ScalarPtrConst(),
         join_connectivity, split_connectivity,
         contour_tree, augment_flag, merge_ident_flag);
    }

    // **************************************************
    // Contour tree utility functions
    // **************************************************

    template <typename NODE_TYPE, typename NUM_TYPE>
    void compress_tree(NODE_TYPE * tree, const NUM_TYPE num_nodes)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { tree[i].ParentCompress(); }
    }

    template <typename ATYPE, typename BTYPE>
    void copy_parent(ATYPE * treeA, int i, BTYPE * treeB)
    {
      if (treeA[i].IsRoot())
        { treeB[i].SetParent(NULL); }
      else {
        int j = ParentIndex(treeA, i);
        treeB[i].SetParent(treeB+j);
      }
    }

    template <typename ATYPE, typename BTYPE, typename NTYPE>
    void copy_tree
    (const ATYPE * treeA, BTYPE * treeB, const NTYPE num_nodes)
    {
      for (NTYPE i = 0; i < num_nodes; i++) 
        { copy_parent(treeA, i, treeB); }

      for (NTYPE i = 0; i < num_nodes; i++) 
        if (treeA[i].IsDeleted()) 
          { treeB[i].Delete(); }
    }

    template <typename ATYPE, typename BTYPE, typename NTYPE>
    void copy_contour_tree
    (const ATYPE * treeA, BTYPE * treeB, const NTYPE num_nodes)
    {
      for (NTYPE i = 0; i < num_nodes; i++) {
        copy_parent(treeA, i, treeB);
        treeB[i].SetEdgeType(treeA[i].EdgeType());
      }
    }

    template <typename ATYPE, typename BTYPE, typename NTYPE>
    void copy_weighted_contour_tree
    (const ATYPE * treeA, BTYPE * treeB, const NTYPE num_nodes)
    {
      for (NTYPE i = 0; i < num_nodes; i++) {
        copy_parent(treeA, i, treeB);
        treeB[i].SetEdgeType(treeA[i].EdgeType());
        treeB[i].SetWeight(treeA[i].Weight());
      }
    }

    template <typename NODE_TYPE, typename NUM_TYPE>
    void set_node_identifiers(NODE_TYPE * tree, const NUM_TYPE num_nodes)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { tree[i].SetIdent(i); }
    }

    template <typename NODE_TYPE>
    void set_node_identifiers(std::vector<NODE_TYPE> & tree)
    {
      set_node_identifiers(&(tree[0]), tree.size());
    }

    template <typename NODE_TYPE, typename STYPE>
    inline bool compare_scalar
    (const SCALAR_TREE_NODE<NODE_TYPE,STYPE> & a, 
     const SCALAR_TREE_NODE<NODE_TYPE,STYPE> & b)
    {
      return (a.Scalar() < b.Scalar());
    }

    // **************************************************
    // Compute number of children, root nodes
    // **************************************************

    template <typename NODE_TYPE, typename NUM_TYPE, typename NUM2_TYPE>
    void compute_num_children 
    (const NODE_TYPE * tree_node, const NUM_TYPE num_nodes,
     NUM2_TYPE * num_children)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { num_children[i] = 0; };

      for (NUM_TYPE i = 0; i < num_nodes; i++)
        if (!tree_node[i].IsDeleted() &&
            !tree_node[i].IsRoot()) {
          NUM_TYPE j = ParentIndex(tree_node, i);
          num_children[j]++;
        }
    }

    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_num_children 
    (const std::vector<NODE_TYPE> & tree, NUM_TYPE * num_children)
    {
      compute_num_children(&(tree[0]), tree.size(), num_children);
    }

    template <typename NODE_TYPE, typename NUM_TYPE, typename NUM2_TYPE>
    void compute_num_children 
    (const NODE_TYPE * tree_node, const NUM_TYPE num_nodes, 
     const EDGE_TYPE edge_type, NUM2_TYPE * num_children)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { num_children[i] = 0; };

      for (NUM_TYPE i = 0; i < num_nodes; i++)
        if (!tree_node[i].IsDeleted() &&
            !tree_node[i].IsRoot()) {
          NUM_TYPE j = ParentIndex(tree_node, i);

          if (tree_node[i].EdgeType() == edge_type)
            { num_children[j]++; }
          else
            { num_children[i]++; }
        }
    }

    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_num_children 
    (const std::vector<NODE_TYPE> & tree, 
     const EDGE_TYPE edge_type, NUM_TYPE * num_children)
    {
      compute_num_children(&(tree[0]), tree.size(), edge_type, num_children);
    }

    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_num_children 
    (NODE_TYPE * tree_node, const NUM_TYPE num_nodes)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { tree_node[i].SetNumChildren(0); };

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree_node[i].IsDeleted() &&
            !tree_node[i].IsRoot())
          { tree_node[i].Parent()->IncrementNumChildren(); }
      }
    }

    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_is_root
    (const NODE_TYPE * tree_node, const NUM_TYPE num_nodes,
     bool * is_root)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) 
        { is_root[i] = tree_node[i].IsRoot(); }
    }

    // **************************************************
    // Compute subtree sizes and node weights
    // **************************************************

    /// Compute subtree sizes.
    template <typename NODE_TYPE, typename NUM_TYPE>
    void compute_subtree_sizes
    (const NODE_TYPE * tree, const NUM_TYPE num_nodes, NUM_TYPE * subtree_size)
    {
      std::vector<NUM_TYPE> node_list;
      IJK::ARRAY<NUM_TYPE> num_children(num_nodes);

      // Initialize all sizes to 1.
      for (NUM_TYPE i = 0; i < num_nodes; i++)
        { subtree_size[i] = 1; }

      compute_num_children(tree, num_nodes, num_children.Ptr());

      // Add leaves to node_list.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (num_children[i] == 0) 
          { node_list.push_back(i); }
      }

      while (node_list.size() != 0) {
        NUM_TYPE i = node_list.back();
        node_list.pop_back();
        if (!tree[i].IsRoot()) {
          NODE_TYPE * parent = tree[i].Parent();
          NUM_TYPE k = NUM_TYPE(parent - tree);
          subtree_size[k] += subtree_size[i];
          num_children[k]--;
          if (num_children[k] == 0) 
            { node_list.push_back(k); }
        }
      }
    }

    /// Compute contour tree node weights.
    /// @pre tree1, tree2 and ctree all have the same number of nodes.
    ///    i'th node of ctree corresponds to i'th node of tree1 and tree2.
    template <typename WTYPE, typename NUM_TYPE, typename NODE_TYPE>
    void compute_ctree_weights
    (const WTYPE * tree1_weight, const WTYPE * tree2_weight,
     const NUM_TYPE num_nodes, NODE_TYPE * ctree)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        WTYPE w = std::min(tree1_weight[i], tree2_weight[i]);
        ctree[i].SetWeight(w);
      }
    }

    // **************************************************
    // Simplify contour tree
    // **************************************************

    /// Simplify contour tree by removing leaves
    ///   with scalar difference with parent less than min_scalar_diff.
    /// @pre scalar[contour_tree[i].Ident()] is the scalar value
    ///        of contour tree node i.
    template <typename CNODE_TYPE, typename NUM_TYPE,
              typename STYPE, typename STYPE2,
              typename CNODE_TYPE2>
    void simplify_contour_tree
    (const CNODE_TYPE * contour_tree, const NUM_TYPE num_nodes,
     const STYPE * scalar, const STYPE2 min_scalar_diff,
     std::vector<CNODE_TYPE2> & simplified_contour_tree)
    {      
      IJK::ARRAY<CONTOUR_TREE_NODE<WEIGHTED_TREE_NODE_DEL<TREE_NODE,STYPE2> > >
        ctree2(num_nodes);

      copy_contour_tree(contour_tree, ctree2.Ptr(), num_nodes);

      // Set weight of each node to scalar difference with parent.
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (ctree2[i].IsRoot()) {
          ctree2[i].SetWeight(1);
        }
        else {
          CNODE_TYPE2 * p = contour_tree[i].Parent();

          NUM_TYPE iv0 = contour_tree[i].Ident();
          NUM_TYPE iv1 = p->Ident();
          STYPE s0 = scalar[iv0];
          STYPE s1 = scalar[iv1];
          if (s0 < s1) {
            ctree2[i].SetWeight(s1-s0);
          }
          else {
            ctree2[i].SetWeight(s0-s1);
          }

        }
      }

      delete_light_nodes(ctree2.Ptr(), num_nodes, min_scalar_diff);

      compress_tree(ctree2.Ptr(), num_nodes);

      delete_1_1_nodes_from_contour_tree(ctree2.Ptr(), num_nodes);

      copy_weighted_tree_skip_deleted_nodes
        (ctree2.PtrConst(), num_nodes, simplified_contour_tree);

      // Copy identifiers and weights from contour_tree 
      //   to simplified_contour_tree.
      for (NUM_TYPE i = 0; i < simplified_contour_tree.size(); i++) {
        NUM_TYPE j = simplified_contour_tree[i].Ident();
        simplified_contour_tree[i].SetIdent(contour_tree[j].Ident());
        simplified_contour_tree[i].SetWeight(contour_tree[j].Weight());
      }
    }

    /// Simplify contour tree by removing leaves
    ///   with scalar difference with parent less than min_scalar_diff.
    template <typename CNODE_TYPE, typename STYPE, typename STYPE2,
              typename CNODE_TYPE2>
    void simplify_contour_tree
    (const std::vector<CNODE_TYPE> & contour_tree,
     const STYPE * scalar, const STYPE2 min_scalar_diff,
     std::vector<CNODE_TYPE2> & simplified_contour_tree)
    {
      simplify_contour_tree
        (&(contour_tree[0]), contour_tree.size(), scalar, min_scalar_diff,
         simplified_contour_tree);
    }

    /// Simplify weighted contour tree by removing leaves
    ///   with weights less than min_weight.
    template <typename CNODE_TYPE, typename NUM_TYPE,
              typename WTYPE, typename CNODE_TYPE2>
    void simplify_weighted_contour_tree
    (const CNODE_TYPE * contour_tree, const NUM_TYPE num_nodes,
     const WTYPE min_weight,
     std::vector<CNODE_TYPE2> & simplified_contour_tree)
    {
      IJK::ARRAY<CONTOUR_TREE_NODE<WEIGHTED_TREE_NODE_DEL<TREE_NODE,WTYPE> > >
        ctree2(num_nodes);

      copy_weighted_contour_tree(contour_tree, ctree2.Ptr(), num_nodes);

      delete_light_nodes(ctree2.Ptr(), num_nodes, min_weight);

      compress_tree(ctree2.Ptr(), num_nodes);

      delete_1_1_nodes_from_contour_tree(ctree2.Ptr(), num_nodes);

      copy_weighted_tree_skip_deleted_nodes
        (ctree2.PtrConst(), num_nodes, simplified_contour_tree);

      // Copy identifiers from contour_tree to simplified_contour_tree.
      for (NUM_TYPE i = 0; i < simplified_contour_tree.size(); i++) {
        NUM_TYPE j = simplified_contour_tree[i].Ident();
        simplified_contour_tree[i].SetIdent(contour_tree[j].Ident());
      }
    }

    /// Simplify weighted contour tree by removing leaves
    ///   with weights less than min_weight.
    template <typename CNODE_TYPE, typename WTYPE, typename CNODE_TYPE2>
    void simplify_weighted_contour_tree
    (const std::vector<CNODE_TYPE> & contour_tree, const WTYPE min_weight,
     std::vector<CNODE_TYPE2> & simplified_contour_tree)
    {
      simplify_weighted_contour_tree
        (&(contour_tree[0]), contour_tree.size(), min_weight,
         simplified_contour_tree);
    }

    template <typename ITYPE, typename WTYPE>
    class WEIGHT_PAIR {
    public:
      WEIGHT_PAIR() {};

      ITYPE ident;
      WTYPE weight;
    };

    template <typename ITYPE, typename WTYPE>
    inline bool operator < 
      (const WEIGHT_PAIR<ITYPE,WTYPE> & a, 
       const WEIGHT_PAIR<ITYPE,WTYPE> & b)
    {
      return (a.weight < b.weight);
    }

    /// Delete leaves with weights less than min weight.
    template <typename CNODE_TYPE, typename NUM_TYPE, 
              typename WTYPE>
    void delete_light_nodes
    (CNODE_TYPE * contour_tree, const NUM_TYPE num_nodes,
     const WTYPE min_weight)
    {
      IJK::ARRAY<NUM_TYPE> num_join_children(num_nodes);
      IJK::ARRAY<NUM_TYPE> num_split_children(num_nodes);
      IJK::ARRAY<WEIGHT_PAIR<NUM_TYPE,WTYPE> > weight(num_nodes);

      compute_num_children(contour_tree, num_nodes, JOIN_EDGE,
                           num_join_children.Ptr());
      compute_num_children(contour_tree, num_nodes, SPLIT_EDGE,
                           num_split_children.Ptr());

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        weight[i].ident = i;
        weight[i].weight = contour_tree[i].Weight();
      }

      std::sort(weight.Ptr(), weight.Ptr()+num_nodes);

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (weight[i].weight < min_weight) {
          NUM_TYPE j = weight[i].ident;
          if (!contour_tree[j].IsRoot() && 
              !contour_tree[j].IsDeleted() &&
              (num_join_children[j] + num_split_children[j]) == 1) {
            NUM_TYPE k = ParentIndex(contour_tree, j);
            if (num_join_children[j] == 0) {
              if (num_join_children[k] > 1) {
                contour_tree[j].Delete();
                num_join_children[k]--;
              }
            }
            else if (num_split_children[j] == 0) {
              if (num_split_children[k] > 1) {
                contour_tree[j].Delete();
                num_split_children[k]--;
              }
            }
          }
        }
      }
    }


    // **************************************************
    // Class TREE_NODE_DEL member functions
    // **************************************************

    /// Get parent, skipping deleted nodes.
    template <typename TREE_NODE_TYPE>
    TREE_NODE_DEL<TREE_NODE_TYPE> * 
    TREE_NODE_DEL<TREE_NODE_TYPE>::Parent() const
    {
      TREE_NODE_DEL<TREE_NODE_TYPE> * p = 
        static_cast<TREE_NODE_DEL<TREE_NODE_TYPE> *>(this->parent);

      while (p != NULL && (p->IsDeleted())){

        // Check for error.  Tree should not have a loop.
        if (p == p->parent) {
          IJK::PROCEDURE_ERROR error("TREE_NODE_DEL::Parent()");
          throw error;
        }

        p = static_cast<TREE_NODE_DEL<TREE_NODE_TYPE> *>(p->parent);
      }

      return(p);
    }

    /// Get parent, and compress path to parent
    template <typename TREE_NODE_TYPE>
    TREE_NODE_DEL<TREE_NODE_TYPE> *
    TREE_NODE_DEL<TREE_NODE_TYPE>::ParentCompress()
    {
      TREE_NODE_DEL<TREE_NODE_TYPE> * p = Parent();


      // Compress path to parent
      TREE_NODE_DEL<TREE_NODE_TYPE> * q = this;
      while (q != p) {
        TREE_NODE_DEL * next_q = 
          static_cast<TREE_NODE_DEL *> (q->parent);

        q->SetParent(p);
        q = next_q;
      }

      return(p);
    }

    // **************************************************
    // Class WEIGHTED_TREE_NODE_DEL member functions
    // **************************************************

    /// Get parent, and compress path to parent
    template <typename TREE_NODE_TYPE,typename WTYPE>
    WEIGHTED_TREE_NODE_DEL<TREE_NODE_TYPE,WTYPE> *
    WEIGHTED_TREE_NODE_DEL<TREE_NODE_TYPE,WTYPE>::ParentCompress()
    {
      WEIGHTED_TREE_NODE_DEL<TREE_NODE_TYPE,WTYPE> * q;

      WEIGHTED_TREE_NODE_DEL<TREE_NODE_TYPE,WTYPE> * p = Parent();

      // Get max weight from this to p.
      q = this;
      WTYPE w = q->Weight();
      while (q != p) {
        w = std::max(w, q->Weight());
        q = static_cast<WEIGHTED_TREE_NODE_DEL *> (q->parent);
      }

      // Compress path to parent
      q = this;
      while (q != p) {
        WEIGHTED_TREE_NODE_DEL * next_q = 
          static_cast<WEIGHTED_TREE_NODE_DEL *> (q->parent);

        q->SetParent(p);
        q->SetWeight(w);
        q = next_q;
      }

      return(p);
    }

    // **************************************************
    // Class TREE_NODE_EXT member functions
    // **************************************************

    template <typename NTYPE>
    void TREE_NODE_EXT<NTYPE>::Init()
    {
      num_children = 0;
    }

    /// Delete node.
    /// Change number of children of parent, if necessary.
    /// Precondition: Node has at most one child.
    template <typename NTYPE>
    void TREE_NODE_EXT<NTYPE>::Delete()
    {
      is_deleted = true;

      if (NumChildren() == 0) 
        { Parent()->DecrementNumChildren(); }
    }

    // **************************************************
    // Class UNION_FIND_BASE member functions
    // **************************************************

    template <typename NTYPE>
    void UNION_FIND_BASE<NTYPE>::Init(const NTYPE num_elements)
    {
      parent = new NTYPE[num_elements];
      this->num_elements = num_elements;

      for (NTYPE e = 0; e < num_elements; e++)
        { parent[e] = e; };
    }

    template <typename NTYPE>
    void UNION_FIND_BASE<NTYPE>::FreeAll()
    {
      delete [] parent;
      num_elements = 0;
    }

    /// Find root.  No path compression
    template <typename NTYPE>
    NTYPE UNION_FIND_BASE<NTYPE>::FindRoot(const NTYPE e)
    {
      // find root
      NTYPE e2 = e;
      while (parent[e2] != e2) 
        { e2 = parent[e2]; }
      return(e2);
    }

    template <typename NTYPE>
    NTYPE UNION_FIND_BASE<NTYPE>::Find(const NTYPE e)
    {
      // find root
      NTYPE e_root = FindRoot(e);

      // path compression
      NTYPE e_current = e;
      while (parent[e_current] != e_root) {
        NTYPE e_prev = e_current;
        e_current = parent[e_current];
        parent[e_prev] = e_root;
      }
    
      return(e_root);
    }

    template <typename NTYPE>
    void UNION_FIND_BASE<NTYPE>::Union(const NTYPE e1, const NTYPE e2)
    {
      NTYPE e1_root = Find(e1);
      NTYPE e2_root = Find(e2);
      parent[e1_root] = e2_root;
    }

    // **************************************************
    // Class UNION_FIND member functions
    // **************************************************

    template <typename NTYPE, typename DTYPE>
    void UNION_FIND<NTYPE,DTYPE>::Init(const NTYPE num_elements)
    {
      data = new DTYPE[num_elements];
    }

    template <typename NTYPE, typename DTYPE>
    void UNION_FIND<NTYPE,DTYPE>::FreeAll()
    {
      delete [] data;
      data = NULL;
    }

    template <typename NTYPE, typename DTYPE>
    DTYPE UNION_FIND<NTYPE,DTYPE>::Data(const NTYPE e)
    {
      NTYPE e_root = Find(e);
      return(data[e_root]);
    }

    /// Redeclare Union function to set data.
    template <typename NTYPE, typename DTYPE>
    void UNION_FIND<NTYPE,DTYPE>::Union
    (const NTYPE e1, const NTYPE e2, const DTYPE d)
    {
      UNION_FIND_BASE<NTYPE>::Union(e1, e2);
      NTYPE new_root = FindRoot(e1);
      data[new_root] = d;
    }

  }
}

#endif

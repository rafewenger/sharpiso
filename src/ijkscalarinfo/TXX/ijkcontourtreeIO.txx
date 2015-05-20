/// \file ijkcontourteeIO.txx
/// ijk templates defining contour tree IO routines
/// IO routines write output in graphviz .dot format
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

#ifndef _IJKCONTOURTREEIO_
#define _IJKCONTOURTREEIO_

#include <iostream>
#include <sstream>
#include <vector>

#include "ijkcontourtree.txx"

/// Contour tree IO functions
namespace IJK {
  namespace CONTOUR_TREE {

    // **************************************************
    // Output tree
    // **************************************************

    class DOT_PARAM {

    public:
      typedef enum {DEFAULT_RANK_DIR, TOP_TO_BOTTOM, BOTTOM_TO_TOP} RANK_DIR;

      RANK_DIR rank_dir;
      int first_label;
      bool is_directed;
      bool flag_output_scalar;
      bool flag_output_vertex_index;
      bool flag_output_weight;

      /// If true, try to force vertices to be ranked by scalar values
      bool flag_rank_scalar;

    public:
      DOT_PARAM(): 
        rank_dir(DEFAULT_RANK_DIR), is_directed(false), first_label(0)
      {};

      // get functions
      bool IsDirected() const { return(is_directed); };
    };

    /// Output dot header
    void output_dot_header
    (std::ostream & out, const std::string & graph_name,
     const DOT_PARAM & dot_param)
    {
      using std::endl;

      if (dot_param.is_directed) { out << "digraph "; }
      else {out << "graph "; };
      out << graph_name << " {" << endl;
    };

    /// Output dot graph attibutes
    void output_dot_graph_attributes
    (std::ostream & out, const DOT_PARAM & dot_param)
    {
      using std::endl;

      switch(dot_param.rank_dir) {

      case DOT_PARAM::BOTTOM_TO_TOP: 
        out << "  rankdir=BT;" << endl;
        break;

      case DOT_PARAM::TOP_TO_BOTTOM: 
        out << "  rankdir=TB;" << endl;
        break;
      }

    };

    /// Output dot vertex attribute
    template <typename NODE_TYPE>
    void output_dot_vertex_attribute
    (std::ostream & out, const NODE_TYPE & node, const DOT_PARAM & dot_param)
    {
      using std::endl;
      std::ostringstream label;

      label.str("");
      if (dot_param.flag_output_vertex_index) 
        { label << node.Ident() + dot_param.first_label; }

      out << "  " << node.Ident() 
          << " [label=\"" << label.str() << "\"];" << endl;
    }

    /// Output dot vertex attribute
    template <typename NODE_TYPE, typename SCALAR_TYPE>
    void output_dot_vertex_attribute
    (std::ostream & out, const NODE_TYPE & node, 
     const SCALAR_TYPE * scalar, const DOT_PARAM & dot_param)
    {
      using std::endl;
      std::ostringstream label;

      typename NODE_TYPE::IDENT_TYPE iv = node.Ident();

      label.str("");
      if (dot_param.flag_output_scalar) 
        { label << scalar[iv]; }

      if (dot_param.flag_output_vertex_index) {
        std::string str = label.str();
        label.str("");
        label << iv;
        if (str.length() > 0) {
          label << " (" << str << ")";
        }
      }

      out << "  " << iv
          << " [label=\"" << label.str() << "\"];" << endl;
    }

    /// Output dot weighted vertex attribute
    template <typename NODE_TYPE, typename SCALAR_TYPE>
    void output_dot_weighted_vertex_attribute
    (std::ostream & out, const NODE_TYPE & node, 
     const SCALAR_TYPE * scalar, const DOT_PARAM & dot_param)
    {
      using std::endl;
      std::ostringstream label;

      typename NODE_TYPE::IDENT_TYPE iv = node.Ident();

      label.str("");
      if (dot_param.flag_output_scalar) 
        { label << scalar[iv]; }

      if (dot_param.flag_output_weight) {
        if (label.str().length() > 0)
          { label << ","; }
        label << node.Weight();
      }

      if (dot_param.flag_output_vertex_index) {
        std::string str = label.str();
        label.str("");
        label << iv;
        if (str.length() > 0) {
          label << " (" << str << ")";
        }
      }

      out << "  " << iv
          << " [label=\"" << label.str() << "\"];" << endl;
    }

    /// Output dot vertex attributes
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_dot_vertices
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const DOT_PARAM & dot_param)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) 
        { output_dot_vertex_attribute(out, tree[i], dot_param); }
      out << std::endl;
    };

    /// Output dot vertex attibutes
    template <typename NODE_TYPE, typename NUM_TYPE, typename STYPE>
    void output_dot_vertices
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const STYPE * scalar, const DOT_PARAM & dot_param)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        output_dot_vertex_attribute(out, tree[i], scalar, dot_param);
      }
      out << std::endl;
    };

    /// Output dot weighted vertex attibutes
    template <typename NODE_TYPE, typename NUM_TYPE, typename STYPE>
    void output_dot_weighted_vertices
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const STYPE * scalar, const DOT_PARAM & dot_param)
    {
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        output_dot_weighted_vertex_attribute
          (out, tree[i], scalar, dot_param);
      }
      out << std::endl;
    };

    /// Output dot vertices, ranked by scalar value
    template <typename NODE_TYPE, typename NUM_TYPE, typename STYPE>
    void output_dot_ranked_vertices
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const STYPE * scalar, const DOT_PARAM & dot_param,
     void out_func(std::ostream &, const NODE_TYPE &, 
                   const STYPE *, const DOT_PARAM &),
     std::vector<std::string> & dummy_node_label)
    {
      IJK::ARRAY<SCALAR_TREE_NODE<NODE_TYPE,STYPE> > tree2(num_nodes);
      std::ostringstream label;
      using std::endl;

      if (num_nodes < 1) { return; }

      // Copy tree to tree2
      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        NUM_TYPE iv = tree[i].Ident();
        tree2[i].Copy(tree[i]);
        tree2[i].SetScalar(scalar[iv]);
      }

      // Sort tree2.  
      // Note: Sorting tree2 destroys correctness of tree2[i].Parent().
      std::sort(tree2.Ptr(), tree2.Ptr()+num_nodes, 
                compare_scalar<NODE_TYPE,STYPE>);

      NUM_TYPE i = 0;
      while (i < num_nodes) {
        STYPE s = tree2[i].Scalar();
        label.str("");
        label << "dummy" << s;
        dummy_node_label.push_back(label.str());
        out << "  {" << endl;
        out << "    rank=same;" << endl;
        out << "  ";
        out_func(out, tree2[i], scalar, dot_param);
        i++;
        while (i < num_nodes && s == tree2[i].Scalar()) {
          out << "  ";
          out_func(out, tree2[i], scalar, dot_param);
          i++;
        }
        out << "    " << label.str() << "[style=invis];" << endl;
        out << "  }" << endl;
      }
      out << std::endl;
    };

    /// Output dot edges
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_dot_contour_edges
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes,
     const DOT_PARAM & dot_param)
    {
      using std:: endl;

      std::string edge_string = " -- ";
      if (dot_param.is_directed)
        { edge_string = " -> "; }

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree[i].IsDeleted() && !tree[i].IsRoot()) {
          out << "  ";
          if (tree[i].EdgeType() == JOIN_EDGE)
            { out << i << edge_string << ParentIndex(tree, i); }
          else
            { out << ParentIndex(tree, i) << edge_string << i; }
          out << ";" << endl;
        }
      }
    };

    /// Output dot edges
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_id_dot_contour_edges
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes,
     const DOT_PARAM & dot_param)
    {
      using std:: endl;

      std::string edge_string = " -- ";
      if (dot_param.is_directed)
        { edge_string = " -> "; }

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree[i].IsDeleted() && !tree[i].IsRoot()) {
          out << "  ";
          if (tree[i].EdgeType() == JOIN_EDGE)
            { out << tree[i].Ident() << edge_string << tree[i].Parent()->Ident(); }
          else
            { out << tree[i].Parent()->Ident() << edge_string << tree[i].Ident(); }
          out << ";" << endl;
        }
      }
      out << endl;
    };

    /// Output dot edges
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_id_dot_edges
    (std::ostream & out, const NODE_TYPE * tree, const NUM_TYPE num_nodes,
     const DOT_PARAM & dot_param)
    {
      using std:: endl;

      for (NUM_TYPE i = 0; i < num_nodes; i++) {
        if (!tree[i].IsDeleted() && !tree[i].IsRoot()) {
          out << "  ";
          if (dot_param.is_directed) {
            out << tree[i].Ident() << " -> "
                << tree[i].Parent()->Ident();
          }
          else
            { out << tree[i].Ident() << " -- " 
                  << tree[i].Parent()->Ident(); }
          out << ";" << endl;
        }
      }
    };

    /// Output dot dummy edges
    inline void output_dot_dummy_edges
    (std::ostream & out, std::vector<std::string> & dummy_node_label)
    {
      using std::endl;

      typedef std::vector<std::string>::size_type NUM_TYPE;
      if (dummy_node_label.size() < 2) { return; }

      for (NUM_TYPE i = 0; i+1 < dummy_node_label.size(); i++) {
        out << "  " << dummy_node_label[i] << " -- "
            << dummy_node_label[i+1] << " [style=invis];" << endl;
      }
      out << endl;
    }

    /// Output tree in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_dot(std::ostream & out, const std::string & graph_name,
                    const NODE_TYPE * tree, const NUM_TYPE num_nodes,
                    const DOT_PARAM & dot_param)
    {
      using std:: endl;

      // opening {
      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      output_dot_edges(out, tree, num_nodes, dot_param);

      // closing }
      out << "}" << endl;
    }

    /// Output tree in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_dot_contour
    (std::ostream & out, const std::string & graph_name,
     const NODE_TYPE * tree, const NUM_TYPE num_nodes,
     const DOT_PARAM & dot_param)
    {
      using std:: endl;

      // opening {
      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      output_dot_contour_edges(out, tree, num_nodes, dot_param);

      // closing }
      out << "}" << endl;
    }

    /// Output tree identifiers in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_id_dot(std::ostream & out, const std::string & graph_name,
                       const NODE_TYPE * tree, const NUM_TYPE num_nodes,
                       const DOT_PARAM & dot_param)
    {
      using std:: endl;

      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      output_id_dot_edges(out, tree, num_nodes, dot_param);

      // closing }
      out << "}" << endl;
    }

    /// Output tree identifiers in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE>
    void output_id_dot_contour
    (std::ostream & out, const std::string & graph_name,
     const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const DOT_PARAM & dot_param)
    {
      using std:: endl;

      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      if (!dot_param.flag_output_vertex_index) {
        output_dot_vertices(out, tree, num_nodes, dot_param);
      }
      output_id_dot_contour_edges(out, tree, num_nodes, dot_param);

      // closing }
      out << "}" << endl;
    }


    /// Output tree identifiers in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE, typename STYPE>
    void output_id_dot_contour
    (std::ostream & out, const std::string & graph_name,
     const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const STYPE * scalar,
     const DOT_PARAM & dot_param)
    {
      std::vector<std::string> dummy_node_label;

      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      if (dot_param.flag_rank_scalar) {
        output_dot_ranked_vertices
          (out, tree, num_nodes, scalar, dot_param, 
           output_dot_vertex_attribute<NODE_TYPE,STYPE>, 
           dummy_node_label);
      }
      else {
        output_dot_vertices(out, tree, num_nodes, scalar, dot_param);
      }
      output_id_dot_contour_edges(out, tree, num_nodes, dot_param);
      if (dot_param.flag_rank_scalar) {
        output_dot_dummy_edges(out, dummy_node_label);
      }

      // closing }
      out << "}" << std::endl;
    }

    /// Output tree identifiers in graphviz dot format
    template <typename NODE_TYPE, typename NUM_TYPE, typename STYPE>
    void output_weighted_id_dot_contour
    (std::ostream & out, const std::string & graph_name,
     const NODE_TYPE * tree, const NUM_TYPE num_nodes, 
     const STYPE * scalar,
     const DOT_PARAM & dot_param)
    {
      std::vector<std::string> dummy_node_label;
      output_dot_header(out, graph_name, dot_param);
      output_dot_graph_attributes(out, dot_param);
      if (dot_param.flag_rank_scalar) {
        output_dot_ranked_vertices
          (out, tree, num_nodes, scalar, dot_param, 
           output_dot_weighted_vertex_attribute<NODE_TYPE,STYPE>, 
           dummy_node_label);
      }
      else {
        output_dot_weighted_vertices
          (out, tree, num_nodes, scalar, dot_param);
      }
      output_id_dot_contour_edges(out, tree, num_nodes, dot_param);
      if (dot_param.flag_rank_scalar) {
        output_dot_dummy_edges(out, dummy_node_label);
      }

      // closing }
      out << "}" << std::endl;
    }

  }
}

#endif

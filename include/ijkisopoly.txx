/// \file ijkisopoly.txx
/// ijk templates for extracting an isosurface patch from a polyhedron
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2013 Rephael Wenger

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

#ifndef _IJKISOPOLY_
#define _IJKISOPOLY_

#include <numeric>

#include "ijk.txx"
#include "ijkbits.txx"
#include "ijkcube.txx"

namespace IJK {

  // **************************************************
  // SUBROUTINES FOR EXTRACTING ISOSURFACE PATCH
  // **************************************************

  /// Return true if isosurface intersects polyhedron
  template <typename STYPE, typename STYPE2, typename VTYPE, typename ITYPE, 
            typename NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const VTYPE iv0, const ITYPE * increment, const NTYPE num_poly_vertices)
  // return true if isosurface intersects polyhedron
  // Precondition: num_poly_vertices > 0
  {
    VTYPE iv1 = iv0 + increment[0];
    if (scalar[iv1] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        if (scalar[iv1] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[iv1] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        if (scalar[iv1] < isovalue) { return(true); };
      }
    }

    return(false);
  }

  /// Return true if isosurface intersects polyhedron
  template <typename STYPE, typename STYPE2, typename NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const NTYPE num_poly_vertices)
  // return true if isosurface intersects polyhedron
  // Precondition: num_poly_vertices > 0
  {
    if (scalar[0] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        if (scalar[j] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[0] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        if (scalar[j] < isovalue) { return(true); };
      }
    }

    return(false);
  }

  /// Compute isosurface table index of polyhedron with primary vertex iv0
  template <typename STYPE, typename STYPE2, typename VTYPE, typename INC_TYPE, 
            typename NTYPE, typename ITYPE>
  inline void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const VTYPE iv0, const INC_TYPE * increment,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
      VTYPE iv1 = iv0 + increment[j];
      if (scalar[iv1] >= isovalue) {
        it = (it | (1L << j));
      };
    };
  }

  /// Compute isosurface table index of polyhedron
  template <typename STYPE, typename STYPE2, typename NTYPE, typename ITYPE>
  void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
      if (scalar[j] >= isovalue) {
        it = (it | (1L << j));
      };
    };
  }

  /// Compute isosurface table index for list of cubes.
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename SCALAR_TYPE, typename CTYPE, 
            typename TABLE_INDEX_TYPE>
  void compute_cube_isotable_index
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue, const std::vector<CTYPE> & cube_list,
   std::vector<TABLE_INDEX_TYPE> & table_index)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;
 
    const NUMBER_TYPE num_cubes = cube_list.size();
    const NUMBER_TYPE num_cube_vertices = scalar_grid.NumCubeVertices();

    table_index.resize(num_cubes);

    for (NUMBER_TYPE i = 0; i < cube_list.size(); i++) {
      TABLE_INDEX_TYPE it;
      compute_isotable_index
        (scalar_grid.ScalarPtrConst(), isovalue, cube_list[i],
         scalar_grid.CubeVertexIncrement(), num_cube_vertices, it);

      table_index[i] = it;
    }
  }

  /// Add isosurface simplex vertices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename STYPE, 
            typename VTYPE2, typename VTYPE3>
  inline void add_iso_simplex_vertices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv_primary, const STYPE is, 
   const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint)
  {
    for (VTYPE j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
      VTYPE jv = isotable.SimplexVertex(it, is, j);
      VTYPE je = isotable.IsosurfaceVertex(jv).Face();
      VTYPE poly_v0 = isotable.Polyhedron().EdgeEndpoint(je, 0);
      VTYPE poly_v1 = isotable.Polyhedron().EdgeEndpoint(je, 1);

      VTYPE iv0 = iv_primary + increment[poly_v0];
      VTYPE iv1 = iv_primary + increment[poly_v1];
      if (iv0 > iv1) { std::swap(iv0, iv1); };

      iso_simplices.push_back(endpoint.size()/2);
      endpoint.push_back(iv0);
      endpoint.push_back(iv1);
    };

  }

  /// Add isosurface simplices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
                               increment, iso_simplices, endpoint);
    };
  }

  /// Add isosurface simplices and reverse orientation.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_reverse_orient_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
                               increment, iso_simplices, endpoint);

      // reverse orientation by swapping last two vertices
      NTYPE ilast = iso_simplices.size()-1;
      std::swap(iso_simplices[ilast], iso_simplices[ilast-1]);
    };
  }

  /// Add isosurface simplices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   const bool orientation,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      add_iso_simplices(isotable, it, iv0, increment, iso_simplices,
                        endpoint, num_simplices);
    }
    else {
      add_reverse_orient_iso_simplices
        (isotable, it, iv0, increment, iso_simplices,
         endpoint, num_simplices);
    }
  }

  // **************************************************
  // EXTRACT ISOSURFACE PATCH FROM A POLYHEDRON
  // **************************************************

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param mesh_scalar = Array of scalar values of all mesh vertices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iv0 = Index of first vertex in polyhedron.
  /// @param increment = k'th polyhedron vertex has index iv0+increment[k].
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Numver of simplices added to isosurface.
  template <typename STYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void extract_isopatch_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
                        num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
                             num_poly_vertices, it);

      add_iso_simplices(isotable, it, iv0, increment, 
                        iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract reverse orientation isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param mesh_scalar = Array of scalar values of all mesh vertices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iv0 = Index of first vertex in polyhedron.
  /// @param increment = k'th polyhedron vertex has index iv0+increment[k].
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <typename STYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void extract_isopatch_reverse_orient_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
                        num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
                             num_poly_vertices, it);

      add_reverse_orient_iso_simplices
        (isotable, it, iv0, increment, 
         iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_iso_simplices(isotable, it, 0, poly_vertex,
                        iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract reverse orientation isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isovalue = Isovalue.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_reverse_orient_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_reverse_orient_iso_simplices
        (isotable, it, 0, poly_vertex, iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param orientation = Orientation of isosurface simplices.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   const bool orientation, std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      extract_isopatch_from_poly
        (scalar, poly_vertex, isotable, isovalue, 
         iso_simplices, endpoint, num_simplices);
    }
    else {
      extract_isopatch_reverse_orient_from_poly
        (scalar, poly_vertex, isotable, isovalue, 
         iso_simplices, endpoint, num_simplices);
    }
  }

  // **************************************************
  // EXTRACT DUAL ISOSURFACE POLYTOPE AROUND EDGE
  // **************************************************

  /// Extract dual isosurface polytope around an edge.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ISOV_TYPE>
  inline void extract_dual_isopoly_around_edge
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, std::vector<ISOV_TYPE> & iso_poly)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    
    const NUM_TYPE num_facet_vertices = grid.NumFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (NUM_TYPE k = 0; k < num_facet_vertices; k++) 
      { iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); }
  }

  /// Extract dual isosurface polytope around an edge.
  /// Return location of isosurface vertex on facet.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ISOV_TYPE, typename FACETV_TYPE>
  inline void extract_dual_isopoly_around_edge
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, 
   std::vector<ISOV_TYPE> & iso_poly,   std::vector<FACETV_TYPE> & facet_vertex)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (NUM_TYPE k = 0; k < num_facet_vertices; k++) {
      iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); 
      facet_vertex.push_back(edge_dir*num_facet_vertices+k);
    }
  }

  /// Extract dual isosurface polytope around an edge, reverse orientation.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ISOV_TYPE>
  inline void extract_dual_isopoly_around_edge_reverse_orient
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, std::vector<ISOV_TYPE> & iso_poly)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumFacetVertices();
    const NUM_TYPE half_num_facet_vertices = num_facet_vertices/2;

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    // Add last half_num_facet_vertices vertices first.
    for (NUM_TYPE k = half_num_facet_vertices; k < num_facet_vertices; k++) 
      { iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); }

    for (NUM_TYPE k = 0; k < half_num_facet_vertices; k++) 
      { iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); }
  }

  /// Extract dual isosurface polytope around an edge, reverse orientation.
  /// Return location of isosurface vertex on facet.
  /// @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ISOV_TYPE, typename FACETV_TYPE>
  inline void extract_dual_isopoly_around_edge_reverse_orient
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir,
   std::vector<ISOV_TYPE> & iso_poly, std::vector<FACETV_TYPE> & facet_vertex)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumFacetVertices();
    const NUM_TYPE half_num_facet_vertices = num_facet_vertices/2;

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    // Add last half_num_facet_vertices vertices first.
    for (NUM_TYPE k = half_num_facet_vertices; 
         k < num_facet_vertices; k++) {
      iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); 
      facet_vertex.push_back(edge_dir*num_facet_vertices+k);
    }

    for (NUM_TYPE k = 0; k < half_num_facet_vertices; k++) {
      iso_poly.push_back(grid.FacetVertex(iv0, edge_dir, k)); 
      facet_vertex.push_back(edge_dir*num_facet_vertices+k);
    }
  }

  /// Extract dual isosurface polytope around bipolar edge.
  /// Checks that edge is bipolar.
  /// Returns list of isosurface polytope vertices.
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename DTYPE, typename ISOV_TYPE>
  void inline extract_dual_isopoly_around_bipolar_edge
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE iend0, const DTYPE edge_dir,
   std::vector<ISOV_TYPE> & iso_poly)
  {
    VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopoly_around_edge
        (scalar_grid, iend0, iend1, edge_dir, iso_poly);
    }
    else if (is_end0_positive && !is_end1_positive) {
      extract_dual_isopoly_around_edge_reverse_orient
        (scalar_grid, iend0, iend1, edge_dir, iso_poly);
    }
  }

  /// Extract dual isosurface polytope around bipolar edge.
  /// Checks that edge is bipolar.
  /// Returns list of isosurface polytope vertices.
  /// Return locations of isosurface vertices on each facet.
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param facet_vertex[i] = Edge of cube containing iso_poly[i].
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename DTYPE,
            typename ISOV_TYPE, typename FACETV_TYPE>
  inline void extract_dual_isopoly_around_bipolar_edge
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE iend0, const DTYPE edge_dir,
   std::vector<ISOV_TYPE> & iso_poly,
   std::vector<FACETV_TYPE> & facet_vertex)
  {
    VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopoly_around_edge
        (scalar_grid, iend0, iend1, edge_dir, iso_poly, facet_vertex);
    }
    else if(is_end0_positive && !is_end1_positive) {
      extract_dual_isopoly_around_edge_reverse_orient
        (scalar_grid, iend0, iend1, edge_dir, iso_poly, facet_vertex);
    }
  }


  // **************************************************
  // COUNT SUBROUTINES
  // **************************************************

  /// Count number of vertices with scalar value above/below isovalue
  /// which are in given facet plane and adjacent to the facet vertices.
  /// @param iv0 Index of leftmost/lowest vertex in facet.
  /// @param orth_dir Direction orthogonal to facet.
  /// @param[out] num_neg Number of vertices with scalar value below isovalue.
  /// @param[out] num_pos Number of vertices with scalar value equal to
  ///   or above isovalue.
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename ORTH_DIR_TYPE, typename NUM_TYPE>
  void count_adjacent_vertices_in_facet_plane
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE & isovalue,
   const VTYPE iv0, const ORTH_DIR_TYPE orth_dir,
   NUM_TYPE & num_neg, NUM_TYPE & num_pos)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    NUM_TYPE num_vertices_in_cube_ridge;
    unsigned long bit0, bit1;      // boundary bits
    unsigned long mask;

    const DTYPE dimension = scalar_grid.Dimension();
    VTYPE iv1;

    // Initialize.
    num_neg = 0;
    num_pos = 0;

    // Compute rightmost/highest vertex in facet.
    iv1 = iv0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir)
        { iv1 = scalar_grid.NextVertex(iv0, d); }
    }

    num_vertices_in_cube_ridge = scalar_grid.NumCubeRidgeVertices();

    scalar_grid.ComputeBoundaryBits(iv0, bit0);
    scalar_grid.ComputeBoundaryBits(iv1, bit1);

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {

        mask = (1L << d);
        if ((bit0 & mask) == 0) {
          for (NUM_TYPE k = 0; k < num_vertices_in_cube_ridge; k++) {
            VTYPE iv = scalar_grid.RidgeVertex(iv0, orth_dir, d, k);
            iv = scalar_grid.PrevVertex(iv, d);

            if (scalar_grid.Scalar(iv) < isovalue)
              { num_neg++; }
            else
              { num_pos++; }
          }
        }

        if ((bit1 & mask) == 0) {
          for (NUM_TYPE k = 0; k < num_vertices_in_cube_ridge; k++) {
            VTYPE iv = scalar_grid.RidgeVertex(iv0, orth_dir, d, k);
            iv += (2*scalar_grid.AxisIncrement(d));

            if (scalar_grid.Scalar(iv) < isovalue)
              { num_neg++; }
            else
              { num_pos++; }
          }
        }
      }
    }

  }

  // **************************************************
  // SPLIT SUBROUTINES
  // **************************************************

  /// Construct list of dual isosurface vertices from cube_list and cube_data.
  /// @isodual_table Dual isosurface lookup table.
  /// @param cube_list[] List of cubes.
  /// @param table_index table_index[i] is index in isosurface lookup table
  ///    of configuration for cube cube_list[i].
  template <typename ISODUAL_TABLE, typename CTYPE, typename TABLE_INDEX_TYPE, 
            typename ISOV0_TYPE, typename ISOV1_TYPE>
  void construct_dual_isovert_list
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<CTYPE> & cube_list,
   const std::vector<TABLE_INDEX_TYPE> & table_index,
   std::vector<ISOV0_TYPE> & first_isov,
   std::vector<ISOV1_TYPE> & iso_vlist)
  {
    typedef typename std::vector<CTYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const SIZE_TYPE num_cubes = cube_list.size();
    first_isov.resize(num_cubes);

    SIZE_TYPE total_num_isov = 0;
    for (SIZE_TYPE i = 0; i < num_cubes; i++) {
      TABLE_INDEX_TYPE it = table_index[i];
      total_num_isov += isodual_table.NumIsoVertices(it);
    }

    iso_vlist.resize(total_num_isov);

    SIZE_TYPE k = 0;
    for (SIZE_TYPE i = 0; i < num_cubes; i++) {
      first_isov[i] = k;
      TABLE_INDEX_TYPE it = table_index[i];
      SIZE_TYPE num_isov = isodual_table.NumIsoVertices(table_index[i]);

      for (SIZE_TYPE j = 0; j < num_isov; j++) {
        iso_vlist[k+j].cube_index = cube_list[i];
        iso_vlist[k+j].patch_index = j;
        iso_vlist[k+j].table_index = it;
      }
      k = k+num_isov;
    }
  }

  /// Set isosurface polytope vertices.
  template <typename ISODUAL_TABLE, typename CINDEX_TYPE, typename FACETV_TYPE,
            typename TABLE_INDEX_TYPE, typename ISOV0_INDEX_TYPE,
            typename ISOV1_INDEX_TYPE>
  void set_dual_isopoly_vertices
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<TABLE_INDEX_TYPE> & table_index,
   const std::vector<CINDEX_TYPE> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   const std::vector<ISOV0_INDEX_TYPE> & first_isov,
   std::vector<ISOV1_INDEX_TYPE> & isopoly)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const int dimension = isodual_table.Dimension();
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    const int num_facet_vertices = cube.NumFacetVertices();

    isopoly.resize(isopoly_cube.size());

    for (ISOV1_INDEX_TYPE i = 0; i < isopoly_cube.size(); i++) {
      ISOV1_INDEX_TYPE k = isopoly_cube[i];
      TABLE_INDEX it = table_index[k];

      // Compute index of facet vertex opposite to facet_vertex[i]
      int facet_vertex_i = facet_vertex[i];
      int ifacet = cube.FacetIndex(facet_vertex_i);
      int j = facet_vertex_i - ifacet*num_facet_vertices;
      int opposite_vertex = (num_facet_vertices-1) - j;
      opposite_vertex += (ifacet*num_facet_vertices);

      isopoly[i] = first_isov[k] + 
        isodual_table.IncidentIsoVertex(it, opposite_vertex);
    }
  }

  /// Split isosurface vertex pairs which create non-manifold edges.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename CINDEX_TYPE, typename TABLE_INDEX_TYPE,  typename NTYPE>
  void split_non_manifold_isov_pairs
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const AMBIG_INFO & ambig_info,
   const std::vector<CINDEX_TYPE> & cube_list, 
   std::vector<TABLE_INDEX_TYPE> & table_index,
   NTYPE & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename std::vector<CINDEX_TYPE>::size_type SIZE_TYPE;

    const DTYPE dimension = grid.Dimension();
    const NUM_TYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);
    const NUM_TYPE num_vertices = grid.NumVertices();
    IJK::ARRAY<SIZE_TYPE> index_to_cube_list(num_vertices);

    num_split = 0;

    // Set up index_to_cube_list.
    for (SIZE_TYPE i = 0; i < cube_list.size(); i++) {
      CINDEX_TYPE cube_index = cube_list[i];
      index_to_cube_list[cube_index] = i;
    }

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX_TYPE it0 = table_index[i0];
      NUM_TYPE num_isov0 = isodual_table.NumIsoVertices(it0);

      if (num_isov0 == 1) {
        if (ambig_info.NumAmbiguousFacets(it0) == 1) {
          CINDEX_TYPE cube_index0 = cube_list[i0];
          int facet_set = ambig_info.AmbiguousFacetBits(it0);
          NUM_TYPE kf = get_first_one_bit(facet_set, num_cube_facets);

          DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
          NUM_TYPE side = IJK::cube_facet_side(dimension, kf);

          if (!grid.IsCubeFacetOnGridBoundary(cube_index0, orth_dir, side)) {

            CINDEX_TYPE cube_index1 =
              grid.AdjacentVertex(cube_index0, orth_dir, side);
            SIZE_TYPE i1 = index_to_cube_list[cube_index1];
            TABLE_INDEX_TYPE it1 = table_index[i1];
            NUM_TYPE num_isov1 = isodual_table.NumIsoVertices(it1);
            if (num_isov1 == 1) {
              if (ambig_info.NumAmbiguousFacets(it1) == 1) {
                table_index[i0] = isodual_table.Complement(it0);
                table_index[i1] = isodual_table.Complement(it1);
                num_split += 2;
              }
            }
          }
          else {
            // Split isosurface vertices in cube_index0.
            table_index[i0]= isodual_table.Complement(it0);
            num_split++;
          }
        }
      }
    }

  }

  template <typename ISODUAL_TABLE, typename TABLE_INDEX,
            typename ITYPE, typename INDEX0_TYPE, typename INDEX1_TYPE,
            typename NUM_TYPE>
  void complement_table_indices
  (const ISODUAL_TABLE & isodual_table, 
   const ITYPE i0, const ITYPE i1, 
   const INDEX0_TYPE it0, const INDEX1_TYPE it1,
   std::vector<TABLE_INDEX> & table_index, NUM_TYPE & num_changed)
  {
    table_index[i0] = isodual_table.Complement(it0);
    table_index[i1] = isodual_table.Complement(it1);
    num_changed += 2;
  }

  /// Select which cube has configuration of split isosurface vertices
  /// where adjacent cubes share an ambiguous facet and one will have
  /// one isosurface vertex while the other has two isosurface vertices.
  /// Choosing the configuration improves reconstruction of sharp edges.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE, typename TABLE_INDEX_TYPE,  typename NTYPE>
  void select_split_1_2_ambig
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE> & cube_list, 
   std::vector<TABLE_INDEX_TYPE> & table_index,
   NTYPE & num_changed)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<CINDEX_TYPE>::size_type SIZE_TYPE;

    const DTYPE dimension = grid.Dimension();
    const NUM_TYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    const NUM_TYPE num_cube_facets = compute_num_cube_facets(dimension);
    const NUM_TYPE num_vertices = grid.NumVertices();
    IJK::ARRAY<SIZE_TYPE> index_to_cube_list(num_vertices);

    num_changed = 0;

    // Set up index_to_cube_list.
    for (SIZE_TYPE i = 0; i < cube_list.size(); i++) {
      CINDEX_TYPE cube_index = cube_list[i];
      index_to_cube_list[cube_index] = i;
    }

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX_TYPE it0 = table_index[i0];
      NUM_TYPE num_isov0 = isodual_table.NumIsoVertices(it0);

      if (ambig_info.NumAmbiguousFacets(it0) == 1) {

        CINDEX_TYPE cube_index0 = cube_list[i0];
        int facet_set = ambig_info.AmbiguousFacetBits(it0);
        NUM_TYPE kf = get_first_one_bit(facet_set, num_cube_facets);

        DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
        NUM_TYPE side = IJK::cube_facet_side(dimension, kf);

        if (!grid.IsCubeFacetOnGridBoundary(cube_index0, orth_dir, side)) {

          // Only process facets on top/right to avoid processing a facet twice.
          if (side) {
            CINDEX_TYPE cube_index1 =
              grid.AdjacentVertex(cube_index0, orth_dir, side);
            SIZE_TYPE i1 = index_to_cube_list[cube_index1];
            TABLE_INDEX_TYPE it1 = table_index[i1];
            NUM_TYPE num_isov1 = isodual_table.NumIsoVertices(it1);

            // if (num_isov0 == 1 and num_isov1 == 2) or 
            //    (num_isov1 == 1 and num_isov0 == 2) or 
            if (num_isov0*num_isov1 == 2) {

              if (ambig_info.NumAmbiguousFacets(it1) == 1) {

                NUM_TYPE num_neg_in_plane, num_pos_in_plane;
                VTYPE iv0 = grid.FacetVertex(cube_index0, kf, 0);

                count_adjacent_vertices_in_facet_plane
                  (grid, isovalue, iv0, orth_dir, 
                   num_neg_in_plane, num_pos_in_plane);

                NUM_TYPE num_neg_cube0, num_pos_cube0;
                NUM_TYPE num_neg_cube1, num_pos_cube1;
                count_bits
                  (it0, num_cube_vertices, num_neg_cube0, num_pos_cube0);
                count_bits
                  (it1, num_cube_vertices, num_neg_cube1, num_pos_cube1);

                if (num_neg_cube0 > num_neg_cube1) {
                  if (num_pos_in_plane > num_neg_in_plane) {
                    if (num_isov0 == 1) {
                      complement_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         table_index, num_changed);
                    }
                  }
                  else if (num_pos_in_plane < num_neg_in_plane) {
                    if (num_isov1 == 1) {
                      complement_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         table_index, num_changed);
                    }
                  }
                }
                else if (num_neg_cube0 < num_neg_cube1) {
                  if (num_pos_in_plane > num_neg_in_plane) {
                    if (num_isov0 == 1) {
                      complement_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         table_index, num_changed);
                    }
                  }
                  else if (num_pos_in_plane < num_neg_in_plane) {
                    complement_table_indices
                      (isodual_table, i0, i1, it0, it1, 
                       table_index, num_changed);
                  }
                }
              }
            }
          }
        }
        else {
          if (num_isov0 == 2) {
            // When ambiguous facet is on the boundary, always prefer 
            //   one isosurface patch.
            table_index[i0] = isodual_table.Complement(it0);
            num_changed ++;
          }
        }
      }
    }

  }

  /// Compute number of cubes which have more than one isosurface vertex.
  template <typename ISODUAL_TABLE, typename TABLE_INDEX_TYPE, 
            typename NTYPE>
  void compute_num_split(const ISODUAL_TABLE & isodual_table,
                         const std::vector<TABLE_INDEX_TYPE> & table_index,
                         NTYPE & num_split)
  {
    num_split = 0;
    for (NTYPE i = 0; i < table_index.size(); i++) {
      TABLE_INDEX_TYPE it = table_index[i];
      if (isodual_table.NumIsoVertices(it) > 1) 
        { num_split ++; }
    }
  }
                              

  // **************************************************
  // SPLIT DUAL ISOSURFACE VERTICES: DATA STRUCTURES
  // **************************************************

  template <typename CI_TYPE, typename PI_TYPE, typename TI_TYPE>
  class DUAL_ISOVERT {
  public:
    CI_TYPE cube_index;
    PI_TYPE patch_index;
    TI_TYPE table_index;
  };

  // **************************************************
  // SPLIT DUAL ISOSURFACE VERTICES: ROUTINES
  // **************************************************

  /// Split dual isosurface vertices.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param table_index table_index[i] is index in isosurface lookup table
  ///    of configuration for cube cube_list[i].
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename ISODUAL_TABLE, typename TABLE_INDEX_TYPE,
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<TABLE_INDEX_TYPE> & table_index,
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split)
  {
    std::vector<ISOV_INDEX_TYPE> first_isov(cube_list.size());

    construct_dual_isovert_list
      (isodual_table, cube_list, table_index, first_isov, iso_vlist);

    set_dual_isopoly_vertices
      (isodual_table, table_index, isopoly_cube, facet_vertex, 
       first_isov, isopoly);

    compute_num_split(isodual_table, table_index, num_split);
  }

  /// Split dual isosurface vertices.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, typename NTYPE>
  void split_dual_isovert
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    std::vector<TABLE_INDEX> table_index;

    compute_cube_isotable_index
      (scalar_grid, isodual_table, isovalue, cube_list, table_index);

    split_dual_isovert
      (isodual_table, cube_list, table_index, 
       isopoly_cube, facet_vertex, iso_vlist, isopoly, num_split);
  }

  /// Split dual isosurface vertices.
  /// Split non-manifold isosurface vertex pairs.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE,
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename TABLE_INDEX_TYPE, 
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_manifold
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<TABLE_INDEX_TYPE> & table_index,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split)
  {
    table_index.resize(cube_list.size());

    compute_cube_isotable_index
      (scalar_grid, isodual_table, isovalue, cube_list, table_index);

    split_non_manifold_isov_pairs
      (scalar_grid, isodual_table, ambig_info, cube_list, table_index,
       num_non_manifold_split);

    split_dual_isovert(isodual_table, cube_list, table_index,
                       isopoly_cube, facet_vertex, 
                       iso_vlist, isopoly, num_split);
  }

  /// Split dual isosurface vertices.
  /// Split non-manifold isosurface vertex pairs.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, typename NTYPE>
  void split_dual_isovert_manifold
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    std::vector<TABLE_INDEX> table_index;

    split_dual_isovert_manifold
      (scalar_grid, isodual_table, ambig_info, isovalue, cube_list, 
       isopoly_cube, facet_vertex, table_index,
       iso_vlist, isopoly, num_split, num_non_manifold_split);
  }
  /// Split dual isosurface vertices.
  /// Select which cube has configuration of split isosurface vertices
  /// where adjacent cubes share an ambiguous facet and one will have
  /// one isosurface vertex while the other has two isosurface vertices.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE,
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename TABLE_INDEX_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_select_split
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<TABLE_INDEX_TYPE> & table_index,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_1_2_change)
  {
    table_index.resize(cube_list.size());

    compute_cube_isotable_index
      (scalar_grid, isodual_table, isovalue, cube_list, table_index);

    select_split_1_2_ambig
      (scalar_grid, isodual_table, ambig_info, isovalue, cube_list,
       table_index, num_1_2_change);

    split_dual_isovert(isodual_table, cube_list, table_index,
                       isopoly_cube, facet_vertex, 
                       iso_vlist, isopoly, num_split);
  }

  /// Split dual isosurface vertices.
  /// Select which cube has configuration of split isosurface vertices
  /// where adjacent cubes share an ambiguous facet and one will have
  /// one isosurface vertex while the other has two isosurface vertices.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_select_split
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_1_2_change)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    std::vector<TABLE_INDEX> table_index;

    split_dual_isovert_select_split
      (scalar_grid, isodual_table, ambig_info, isovalue, cube_list, 
       isopoly_cube, facet_vertex, table_index, iso_vlist, isopoly,
       num_split, num_1_2_change);
  }


  /// Split dual isosurface vertices with special processing 
  ///   of ambiguous facets.
  /// Calls both split_non_manifold_isov_pairs and select_split_1_2_ambig.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename TABLE_INDEX_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_ambig
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<TABLE_INDEX_TYPE> & table_index,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split,
   NTYPE & num_1_2_change)
  {
    table_index.resize(cube_list.size());

    compute_cube_isotable_index
      (scalar_grid, isodual_table, isovalue, cube_list, table_index);

    split_non_manifold_isov_pairs
      (scalar_grid, isodual_table, ambig_info, cube_list, table_index,
       num_non_manifold_split);

    select_split_1_2_ambig
      (scalar_grid, isodual_table, ambig_info, isovalue, cube_list,
       table_index, num_1_2_change);

    split_dual_isovert
      (isodual_table, cube_list, table_index, 
       isopoly_cube, facet_vertex, iso_vlist, isopoly, num_split);
  }

  /// Split dual isosurface vertices with special processing 
  ///   of ambiguous facets.
  /// Calls both split_non_manifold_isov_pairs and select_split_1_2_ambig.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename AMBIG_INFO,
            typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename FACETV_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_ambig
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table, const AMBIG_INFO & ambig_info,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split,
   NTYPE & num_1_2_change)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    std::vector<TABLE_INDEX> table_index;

    split_dual_isovert_ambig
      (scalar_grid, isodual_table, ambig_info, isovalue, cube_list, 
       isopoly_cube, facet_vertex, table_index, iso_vlist, isopoly,
       num_split, num_non_manifold_split, num_1_2_change);
  }

}



#endif

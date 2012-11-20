/// \file ijkgraph.txx
/// ijk templates for graph routines
/// Version 0.1.0

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

#ifndef _IJKGRAPH_
#define _IJKGRAPH_

namespace IJK {

  /// Compute the number of connected components.
  /// Grow components by searching elist[]. (Slow.)
  /// @param vlist[] Array of vertex indices.
  /// @param numv Number of vertices in vlist[].
  /// @param max_vertex_index Max index of any vertex in vlist[].
  ///        All vertices in vlist[] have indices from 0 to max_vertex_index.
  /// @param elist[] Array of edges.
  ///    j'th edges is (elist[2*j],elist[2*j+1])
  /// @param nume Number of edges.
  /// @param[out] num_components Number of connected components.
  /// @param visited[]  Temporary boolean array indicating visited vertices.
  /// @pre Array visited[] must be pre-allocated to size 
  ///        at least max_vertex_index+1.
  template <typename VTYPE1, typename VTYPE2, typename ETYPE, 
            typename NTYPE1, typename NTYPE2, typename NTYPE3>
  void compute_num_connected_components_using_elist
  (const VTYPE1 vlist[], const NTYPE1 numv, const VTYPE2 max_vertex_index,
   const ETYPE elist[], const NTYPE2 nume,
   NTYPE3 & num_components, bool visited[])
  {
    bool foundNext;

    num_components = 0;

    for (NTYPE1 i = 0; i < numv; i++) 
      { visited[vlist[i]] = false; }

    for (NTYPE1 i = 0; i < numv; i++) {
      VTYPE1 iv = vlist[i];
      if (visited[iv] == false) {
        visited[iv] = true;
        num_components++;

        do {
          foundNext = false;
          for (NTYPE2 j = 0; j < nume; j++) {
            VTYPE1 jv0 = elist[2*j];
            VTYPE1 jv1 = elist[2*j+1];
            if (visited[jv0] == true) {
              if (visited[jv1] == false) {
                visited[jv1] = true;
                foundNext = true;
              }
            }
            else { // visited[jv0] = false.
              if (visited[jv1] == true) {
                visited[jv0] = true;
                foundNext = true;
              }
            }
          }
        }
        while (foundNext);
      }
    }
  }

  /// Find identifier of component containing i.
  /// Compress component[] to point to tree root.
  template <typename ITYPE, typename CTYPE>
  ITYPE find_compress(const ITYPE i, CTYPE component[])
  {
    ITYPE root = component[i];
    while (root != component[root])
      { root = component[root]; }

    ITYPE i2 = i;
    while (component[i2] != root) {
      ITYPE i3 = component[i2];
      component[i2] = root;
      i2 = i3;
    }

    return (root);
  }
 
  /// Union components containing i1 and i2.
  template <typename ITYPE1, typename ITYPE2, typename CTYPE>
  void union_components(const ITYPE1 i0, const ITYPE2 i1,
                        CTYPE component[])
  {
    ITYPE1 j0 = find_compress(i0, component);
    ITYPE1 j1 = find_compress(i1, component);
    component[j1] = j0;
  }

  /// Compute connected components using union-find.
  /// Use union-find to compute components.
  /// @param vlist[] Array of vertex indices.
  /// @param numv Number of vertices in vlist[].
  /// @param max_vertex_index Max index of any vertex in vlist[].
  ///        All vertices in vlist[] have indices from 0 to max_vertex_index.
  /// @param elist[] Array of edges.
  ///    j'th edges is (elist[2*j],elist[2*j+1])
  /// @param nume Number of edges.
  /// @param component[out]  Array indicating component.
  ///           component[iv] is component number of vertex iv.
  /// @pre Array component[] must be pre-allocated to size 
  ///        at least max_vertex_index+1.
  template <typename VTYPE1, typename ETYPE, 
            typename NTYPE1, typename NTYPE2, typename CTYPE>
  void compute_connected_components
  (const VTYPE1 vlist[], const NTYPE1 numv, 
   const ETYPE elist[], const NTYPE2 nume,
   CTYPE component[])
  {

    // Initialize each vertex in vlist[] to be its own component.
    for (NTYPE1 i = 0; i < numv; i++) {
      VTYPE1 iv = vlist[i];
      component[iv] = iv;
    }

    for (NTYPE2 j = 0; j < nume; j++) {
      VTYPE1 jv0 = elist[2*j];
      VTYPE1 jv1 = elist[2*j+1];

      union_components(jv0,jv1,component);
    }
  }

  /// Count number of components containing elements in list[].
  template <typename LTYPE, typename NTYPE, typename CTYPE,
            typename NUMC_TYPE> 
  void count_num_components
  (const LTYPE list[], const NTYPE numv, CTYPE component[],
   NUMC_TYPE & num_component)
  {
    num_component = 0;
    for (NTYPE i = 0; i < numv; i++) {
      LTYPE list_i = list[i];
      LTYPE icomp = find_compress(list_i, component);
      if (icomp == list_i)
        { num_component++; }
    }
  }

}

#endif

/// \file ijkmesh_c++11.txx
/// ijk templates for handling polyhedral meshes using C++11 stl functions.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2013 Rephael Wenger

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

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "ijkmesh.txx"

#ifndef _IJKMESH_CPP11_
#define _IJKMESH_CPP11_

namespace IJK {

  template <typename T0, typename TRIPLET_TYPE>
  struct HASH_TRIPLET {

    std::hash<T0> hash_func;

    HASH_TRIPLET(){};

    size_t operator() (const TRIPLET_TYPE & key) const
    {
      return(hash_func(std::get<0>(key)+std::get<1>(key)+std::get<2>(key)));
    }
  };

  /// Sort two objects.
  template <typename T>
  inline void sort_two(T & a0, T & a1)
  {
    if (a0 > a1) { std::swap(a0,a1); };
  }

  /// Sort three objects.
  template <typename T>
  inline void sort_three(T & a0, T & a1, T & a2)
  {
    sort_two(a0, a1);
    sort_two(a0, a2);
    sort_two(a1, a2);
  }

  /// Sort three objects in an array
  template <typename T>
  inline void sort_three(T a[])
  {
    sort_three(a[0], a[1], a[2]);
  }

  /// Insert in hash table, using ordered triplet as key.
  template <typename T, typename VALUE_TYPE, typename HASH_TABLE_TYPE>
  void insert_using_ordered_triplet
  (const T triplet[], const VALUE_TYPE value, HASH_TABLE_TYPE & hash_table)
  {
    T ordered_triplet[3];

    ordered_triplet[0] = triplet[0];
    ordered_triplet[1] = triplet[1];
    ordered_triplet[2] = triplet[2];

    sort_three(ordered_triplet);

    std::tuple<T, T, T> key =
      std::make_tuple(ordered_triplet[0], ordered_triplet[1], 
                      ordered_triplet[2]);
    hash_table.insert(typename HASH_TABLE_TYPE::value_type(key, value));
  }
  
  /// Insert in hash table, using ordered triplet as key.
  /// @param istart Index of starting location of triplet.
  template <typename T, typename NUM_TYPE, 
            typename VALUE_TYPE, typename HASH_TABLE_TYPE>
  void insert_using_ordered_triplet
  (const std::vector<T> & list, const NUM_TYPE istart,
   const VALUE_TYPE value, HASH_TABLE_TYPE & hash_table)
  {
    insert_using_ordered_triplet(&(list[istart]), value, hash_table);
  }

  /// Triangulate quadrilateral containing any triangle from triangle list.
  /// @param triangle_hash Hash table of triangles.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename HASH_TABLE_TYPE>
  void triangulate_quad_containing_triangle
  (const VTYPE0 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2,
   HASH_TABLE_TYPE & triangle_hash)
  {
    const NTYPE NUM_VERT_PER_TRI = 3;
    const NTYPE NUM_VERT_PER_QUAD = 4;
    VTYPE1 vert[NUM_VERT_PER_TRI];

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;

      bool flag_split = false;
      for (NTYPE i = 0; i < NUM_VERT_PER_QUAD; i++) {
        NTYPE j0 = 0;
        for (NTYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
          if (j != i) {
            vert[j0] = quad_vert[k+j];
            j0++;
          }
        }

        sort_three(vert);
        typename HASH_TABLE_TYPE::key_type
          key = std::make_tuple(vert[0], vert[1], vert[2]);

        typename HASH_TABLE_TYPE::const_iterator 
          triangle_iter = triangle_hash.find(key);

        if (triangle_iter != triangle_hash.end()) {
          // split quad at j'th vertex.
          NTYPE k2 = tri_vert.size();
          triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert+k, i, tri_vert);
          flag_split = true;

          // add new triangles to triangle_hash.
          insert_using_ordered_triplet(tri_vert, k2, 1, triangle_hash);
          insert_using_ordered_triplet
            (tri_vert, k2+NUM_VERT_PER_TRI, 1, triangle_hash);
          break;
        }
      }

      if (!flag_split) {
        quad_vert2.push_back(quad_vert[k]);
        quad_vert2.push_back(quad_vert[k+1]);
        quad_vert2.push_back(quad_vert[k+2]);
        quad_vert2.push_back(quad_vert[k+3]);
      }
    }
  }

  /// Triangulate quadrilateral containing any triangle from triangle list.
  /// Triangulate quadrilateral so that triangle is not duplicated.
  /// Quadrilateral vertices are listed in clockwise or counter-clockwise order.
  /// Add new triangles to vector tri_vert.
  /// Put quadrilaterals which are not triangulated in quad_vert2.
  /// Note: Triangles created from splitting quadrilaterals can
  ///       cause other quadrilaterals to split.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_containing_triangle
  (const VTYPE0 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    const NTYPE NUM_VERT_PER_TRI = 3;
    const NTYPE NUM_VERT_PER_QUAD = 4;
    const NTYPE num_tri = tri_vert.size()/NUM_VERT_PER_TRI;
    VTYPE1 vert[NUM_VERT_PER_TRI];

    typedef std::tuple<VTYPE1, VTYPE1, VTYPE1> TRIANGLE_VERT;
    typedef std::unordered_map<TRIANGLE_VERT, NTYPE, 
      HASH_TRIPLET<VTYPE1,TRIANGLE_VERT> > TRIANGLE_HASH_TABLE;

    TRIANGLE_HASH_TABLE triangle_hash;

    for (NTYPE itri = 0; itri < num_tri; itri++) {
      NTYPE k = itri*NUM_VERT_PER_TRI;
      insert_using_ordered_triplet(tri_vert, k, 1, triangle_hash);
    }

    NTYPE old_tri_vert_size = tri_vert.size();
    triangulate_quad_containing_triangle
      (quad_vert, num_quad, tri_vert, quad_vert2, triangle_hash);

    while (old_tri_vert_size != tri_vert.size()) {

      // copy quad_vert2 to quad_vert3.
      std::vector<VTYPE2> quad_vert3;
      quad_vert3.resize(quad_vert2.size());
      std::copy(quad_vert2.begin(), quad_vert2.end(), quad_vert3.begin());
      quad_vert2.clear();

      old_tri_vert_size = tri_vert.size();
      triangulate_quad_containing_triangle
        (&(quad_vert3[0]), quad_vert3.size()/NUM_VERT_PER_QUAD, 
         tri_vert, quad_vert2, triangle_hash);
    }
  }

  /// Triangulate quadrilateral containing any triangle from triangle list.
  /// C++ STL vector format for quad_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_containing_triangle
  (const std::vector<VTYPE0> & quad_vert,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (quad_vert.size() > 0) {
      triangulate_quad_containing_triangle
        (&(quad_vert[0]), num_quad, tri_vert, quad_vert2);
    }
  }

  /// Triangulate quadrilateral i if it shares two (or more) edges 
  ///   with another quad.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1,
            typename HASH_TABLE_TYPE>
  void triangulate_quad_i_sharing_multiple_edges_with_quad
  (const VTYPE0 * quad_vert, const NTYPE iquad,
   std::vector<VTYPE1> & tri_vert, HASH_TABLE_TYPE & triplet_hash,
   std::vector<bool> & flag_split)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;
    const NTYPE NUM_VERT_PER_TRI = 3;
    VTYPE1 vert[NUM_VERT_PER_TRI];

    NTYPE k = iquad*NUM_VERT_PER_QUAD;

    if (!flag_split[iquad]) {

      for (NTYPE i = 0; i < NUM_VERT_PER_QUAD; i++) {
        NTYPE j0 = 0;
        for (NTYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
          if (j != i) {
            vert[j0] = quad_vert[k+j];
            j0++;
          }
        }

        sort_three(vert);
        typename HASH_TABLE_TYPE::key_type
          key = std::make_tuple(vert[0], vert[1], vert[2]);

        typename HASH_TABLE_TYPE::const_iterator 
          triangle_iter = triplet_hash.find(key);

        if (triangle_iter != triplet_hash.end()) {
          // split quad at j'th vertex.
          NTYPE k2 = tri_vert.size();
          triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert+k, i, tri_vert);
          flag_split[iquad] = true;

          break;
        }
      }
    }

    for (NTYPE i = 0; i < NUM_VERT_PER_QUAD; i++) {
      NTYPE j0 = 0;
      for (NTYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
        if (j != i) {
          vert[j0] = quad_vert[k+j];
          j0++;
        }
      }

      // add new triangles to triangle_hash.
      insert_using_ordered_triplet(vert, iquad, triplet_hash);
    }

  }


  /// Triangulate quadrilaterals sharing two (or more) edges 
  ///   with another quad.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_sharing_multiple_edges_with_quad
  (const VTYPE0 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    const NTYPE NUM_VERT_PER_TRI = 3;
    const NTYPE NUM_VERT_PER_QUAD = 4;
    VTYPE1 vert[NUM_VERT_PER_TRI];
    std::vector<bool> flag_split(num_quad, false);

    typedef std::tuple<VTYPE1, VTYPE1, VTYPE1> TRIPLET_TYPE;
    typedef std::unordered_map<TRIPLET_TYPE, NTYPE, 
      HASH_TRIPLET<VTYPE1,TRIPLET_TYPE> > TRIPLET_HASH_TABLE;

    TRIPLET_HASH_TABLE triplet_hash;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      triangulate_quad_i_sharing_multiple_edges_with_quad
        (quad_vert, iquad, tri_vert, triplet_hash, flag_split);
    }

    // Process quadrilaterals in reverse order.
    triplet_hash.clear();
    for (NTYPE j = num_quad; j > 0; j--) {
      NTYPE iquad = j-1;
      triangulate_quad_i_sharing_multiple_edges_with_quad
        (quad_vert, iquad, tri_vert, triplet_hash, flag_split);
    }

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      if (!flag_split[iquad]) {
        NTYPE k = iquad*NUM_VERT_PER_QUAD;
        quad_vert2.push_back(quad_vert[k]);
        quad_vert2.push_back(quad_vert[k+1]);
        quad_vert2.push_back(quad_vert[k+2]);
        quad_vert2.push_back(quad_vert[k+3]);
      }
    }
  }

  /// Triangulate quadrilaterals sharing two (or more) edges 
  ///   with another quad.
  /// C++ STL vector format for quad_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_sharing_multiple_edges_with_quad
  (const std::vector<VTYPE0> & quad_vert,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (quad_vert.size() > 0) {
      triangulate_quad_sharing_multiple_edges_with_quad
        (&(quad_vert[0]), num_quad, tri_vert, quad_vert2);
    }
  }

  /// Triangulate quadrilaterals sharing two (or more) edges 
  ///   with another quad or triangle.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_sharing_multiple_edges
  (const VTYPE0 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    std::vector<VTYPE2> quad_vert3;

    triangulate_quad_sharing_multiple_edges_with_quad
      (quad_vert, num_quad, tri_vert, quad_vert3);
    triangulate_quad_containing_triangle
      (quad_vert3, tri_vert, quad_vert2);
  }

  /// Triangulate quadrilaterals sharing two (or more) edges 
  ///   with another quad or triangle.
  /// C++ STL vector format for quad_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_sharing_multiple_edges
  (const std::vector<VTYPE0> & quad_vert,
   std::vector<VTYPE1> & tri_vert, std::vector<VTYPE2> & quad_vert2)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (quad_vert.size() > 0) {
      triangulate_quad_sharing_multiple_edges
        (&(quad_vert[0]), num_quad, tri_vert, quad_vert2);
    }
  }

};

#endif

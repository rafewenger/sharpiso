/// \file ijklist.txx
/// ijk templates for handling lists
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010,2012 Rephael Wenger

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

#ifndef _IJKLIST_
#define _IJKLIST_

#include "ijk.txx"

#include <algorithm>
#include <utility>
#include <vector>

namespace IJK {

  // **************************************************
  // REMAP LIST
  // **************************************************

  /// Find "first" (identifying) element of the set containing x.
  /// Modifies set with path compression.
  /// @pre T is an integer type.
  template <typename T1, typename T2>
  T2 find_set(const T1 x, T2 * set_ident)
  {
    // Find first element
    T1 y = set_ident[x];
    while (y != set_ident[y])
      { y = set_ident[y]; }

    // apply path compression
    T1 z = set_ident[x];
    while (z != set_ident[z]) {
      T1 z_next = set_ident[z];
      set_ident[z] = y;
      z = z_next;
    }

    return(y);
  }

  /// Find "first" (identifying) element of the set containing x.
  /// Modifies set with path compression.
  /// @pre T is an integer type.
  /// @pre set_ident.size() > 0.
  /// C++ vector version of set_ident.
  template <typename T1, typename T2>
  T2 find_set(const T1 x, std::vector<T2> & set_ident)
  {
    return(find_set(x, &(set_ident.front())));
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  /// @pre max_element is a non-negative integer greater than or equal to
  ///      any number in list[]
  template <typename T1, typename T2, typename T3, typename T4,
            typename NTYPE>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  const T4 max_element, T3 * list, const NTYPE list_length)
  {
    IJK::PROCEDURE_ERROR error("remap_list");

    if (list_length <= 0 || remap_pair.size() == 0) { return; }

    if (max_element < 0) {
      error.AddMessage
        ("Programming error.  max_element must be non-negative.");
      error.AddMessage("  max_element = ", max_element, ".");
      throw error;
    };

    std::vector<T3> remap_set(max_element+1);
    for (T3 i = 0; i <= max_element; i++) 
      { remap_set[i] = i; }

    for (T3 i = 0; i < remap_pair.size(); i++) {
      T1 k1 = remap_pair[i].first;
      T2 k2 = remap_pair[i].second;

      if (k1 <= max_element && k2 <= max_element) {
        remap_set[k1] = find_set(k2, remap_set);
      }
    }

    for (T3 i = 0; i < list_length; i++) 
      { list[i] = find_set(list[i], remap_set); }
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  template <typename T1, typename T2, typename T3, typename NTYPE>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  T3 * list, const NTYPE list_length)
  {
    IJK::PROCEDURE_ERROR error("remap_list");

    if (list_length <= 0 || remap_pair.size() == 0) { return; }

    T3 * max_ptr = std::max_element(list, list+list_length);

    if (*max_ptr < 0) {
      error.AddMessage
        ("Programming error.  Only non-negative numbers allowed in list[].");
      error.AddMessage("  List contains ", *max_ptr, ".");
      throw error;
    };

    remap_list(remap_pair, *max_ptr, list, list_length);
  }

  /// Remap members of list[] based on pairs in remap_pair[]
  /// @pre T1, T2 and T3 must be an integer type, i.e. short, int, long.
  /// @pre list[] is a list of non-negative integers.
  template <typename T1, typename T2, typename T3>
  void remap_list(const std::vector< std::pair< T1, T2 > > & remap_pair,
                  std::vector<T3> & list)
  {
    if (list.size() == 0) { return; }

    remap_list(remap_pair, &(list.front()), list.size());
  }

}

#endif

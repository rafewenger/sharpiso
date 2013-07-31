/// \file ijkprint.txx
/// Templates for printing lists, coordinates, etc.
/// Useful for debugging or for routines printing information
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

#ifndef _IJKPRINT_
#define _IJKPRINT_

#include <cmath>
#include <iostream>
#include <vector>

// *** SHOULD BE <ctime> ***
#include <time.h>

#include "ijk.txx"

namespace IJK {


  /// Print values in list.  No parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list_values
  (std::ostream & out, const ETYPE * list, const NTYPE length,
   const char separator)
  {
    for (NTYPE i = 0; i < length; i++) {
      if (i > 0) { out << separator; };
      out << list[i];
    }
  }

  /// Print list with left and right delimiters.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length,
                  const char separator,
                  const char left_delim, const char right_delim)
  {
    out << left_delim;
    print_list_values(out, list, length, separator);
    out << right_delim;
  }

  /// Print list enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length,
                  const char separator)
  {
    print_list(out, list, length, separator, '(', ')');
  }

  /// Print list separated by commas and enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length)
  {
    print_list(out, list, length, ',');
  }

  /// Print list separated by commas and enclosed in parentheses.
  /// C++ STL vector format for list[].
  template <typename ETYPE>
  void print_list(std::ostream & out, const std::vector<ETYPE> & list)
  {
    print_list(out, IJK::vector2pointer(list), list.size());
  }

  /// Print list of tuples.
  template <typename ETYPE, typename N0TYPE, typename N1TYPE>
  void print_list_of_tuples
  (std::ostream & out, const ETYPE * list, 
   const N0TYPE tuple_size, const N1TYPE num_tuples,
   const char separator0=',', const char separator1=' ',
   const char left_delim='(', const char right_delim=')')
  {
    for (N1TYPE i = 0; i < num_tuples; i++) {
      N1TYPE k = i*tuple_size;
      print_list
        (out, list+k, tuple_size, separator0, left_delim, right_delim);
      if (i+1 < num_tuples)
        { out << separator1; }
    }
  }

  /// Print list of tuples.
  template <typename ETYPE, typename NTYPE>
  void print_list_of_tuples
  (std::ostream & out, const std::vector<ETYPE> & list, 
   const NTYPE tuple_size)
  {
    typedef typename std::vector<ETYPE>::size_type SIZE_TYPE;

    if (list.size() <= 0) { return; }

    SIZE_TYPE num_tuples = list.size()/tuple_size;
    print_list_of_tuples(out, &(list[0]), tuple_size, num_tuples);
  }

  template <typename CTYPE>
  void print_coord3D(std::ostream & out, const CTYPE * coord)
  {
    print_list(out, coord, 3);
  }

  template <typename CTYPE>
  void print_coord3D(std::ostream & out, const CTYPE * coord, const char * s)
  {
    print_coord3D(out, coord);
    out << s;
  }

  template <typename GTYPE, typename VTYPE>
  void print_grid_coord(std::ostream & out, const GTYPE & grid, const VTYPE iv)
  {
    typename GTYPE::DIMENSION_TYPE dimension = grid.Dimension();
    VTYPE coord[dimension];

    grid.ComputeCoord(iv, coord);
    print_list(out, coord, dimension);
  }

  template <typename GTYPE, typename VTYPE>
  void print_grid_coord(std::ostream & out, const GTYPE & grid, const VTYPE iv,
                        char * s)
  {
    typename GTYPE::DIMENSION_TYPE dimension = grid.Dimension();
    VTYPE coord[dimension];

    grid.ComputeCoord(iv, coord);
    print_list(out, coord, dimension);
    out << s;
  }

  template <typename GTYPE>
  void print_grid_size(std::ostream & out, const GTYPE & grid, const char * s)
  {
    if (grid.Dimension() < 1) 
      { out << "Grid dimension: " << grid.Dimension(); }
    else {
      out << "Grid ";
      for (int d = 0; d < grid.Dimension(); d++) {
        out << grid.AxisSize(d);
        if (d+1 < grid.Dimension()) 
          { out << "x"; }
      };
    }
    out << s; 
  }

  inline void print_time(std::ostream & out, const char * s, 
                         const clock_t & t_start, const clock_t & t_end)
  {
    out << s;
    out << float(t_end-t_start)/CLOCKS_PER_SEC;
    out << std::endl;
  }

  /// Convert (numerator/denominator) to percentage.
  /// @pre denominator is not zero.
  template <typename T0, typename T1, typename T2>
  void compute_percent
  (const T0 numerator, const T1 denominator, T2 & percent)
  {
    if (denominator == 0) {
      percent = 100*numerator;
    }
    else {
      percent = 100*(double(numerator)/denominator);
    }
  }

  /// Print percentage(if denominator is not zero.)
  /// Percentage is rounded down to nearest integer.
  template <typename T0, typename T1>
  void print_percent
  (std::ostream & out, const T0 numerator, const T1 denominator)
  {
    if (denominator != 0) {
      float percent;
      compute_percent(numerator, denominator, percent);
      percent = std::floor(percent);
      out << "(" << percent << "%)";
    }
  }

}

#endif

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
#include <time.h>


namespace IJK {


  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length)
  {
    out << "(";
    for (NTYPE i = 0; i < length; i++) {
      if (i > 0) { out << ","; };
      out << list[i];
    }
    out << ")";
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

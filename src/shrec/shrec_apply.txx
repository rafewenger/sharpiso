/// \file shrec_apply.txx
/// Apply templates

/*
Copyright (C) 2015 Arindam Bhattacharya and Rephael Wenger

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _SHREC_APPLY_TXX_
#define _SHREC_APPLY_TXX_

#include "shrec_isovert.h"

#include <vector>

namespace SHREC {

  // **************************************************
  // UTILITY TEMPLATES
  // **************************************************

  /// Return (b & (1 << i))
  template <typename BTYPE, typename NTYPE>
  inline BTYPE and_bit(const BTYPE b, const NTYPE i)
  {
    return(b & (BTYPE(1) << i));
  }

  /// Convert cube corner index to vertex bit mask
  template <typename DTYPE, typename ITYPE, typename BTYPE>
  inline void convert2vertex_bit_mask
  (const DTYPE dimension, const ITYPE corner_index, BTYPE & bit_mask)
  {
    bit_mask = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (and_bit(corner_index, d) == 0)
        { bit_mask = (bit_mask | (BTYPE(1) << (2*d))); }
      else
        { bit_mask = (bit_mask | (BTYPE(1) << (2*d+1))); }
    }
  }

  /// Convert cube corner index to edge bit mask
  template <typename DTYPE, typename DTYPE2, typename ITYPE, typename BTYPE>
  inline void convert2edge_bit_mask
  (const DTYPE dimension, const DTYPE2 edge_dir,
   const ITYPE corner_index, BTYPE & bit_mask)
  {
    bit_mask = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != edge_dir) {
        if (and_bit(corner_index, d) == 0)
          { bit_mask = (bit_mask | (BTYPE(1) << (2*d))); }
        else
          { bit_mask = (bit_mask | (BTYPE(1) << (2*d+1))); }
      }
    }
  }

  /// Compute cube which shares i'th corner (and no other corners)
  ///   with icubeA.
  template <typename GTYPE, typename VTYPE, typename ITYPE>
  inline VTYPE compute_vertex_adjacent_cube
  (const GTYPE & grid, const VTYPE icubeA, const ITYPE i)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    VTYPE icubeB = icubeA;
    for (DTYPE d = 0; d < grid.Dimension(); d++) {

      if (and_bit(i, d) == 0)
        { icubeB = grid.PrevVertex(icubeB, d); }
      else
        { icubeB = grid.NextVertex(icubeB, d); }
    }

    return(icubeB);
  }

  /// Compute cube which shares i'th edge in direction edge_dir 
  ///   (and no other edges) with icubeA.
  template <typename GTYPE, typename VTYPE, typename DTYPE2, typename ITYPE>
  inline VTYPE compute_edge_adjacent_cube
  (const GTYPE & grid, const VTYPE icubeA, const DTYPE2 edge_dir,
   const ITYPE i)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    VTYPE icubeB = icubeA;
    for (DTYPE d = 0; d < grid.Dimension(); d++) {

      if (d != edge_dir) {
        if (and_bit(i, d) == 0)
          { icubeB = grid.PrevVertex(icubeB, d); }
        else
          { icubeB = grid.NextVertex(icubeB, d); }
      }
    }

    return(icubeB);
  }

  // **************************************************
  // FACET ADJACENT APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is facet adjacent to icubeA. (4+4 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4>
  void apply_to_cubes_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
        const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4);
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (and_bit(boundary_bits, 2*d) == 0) {
          const VTYPE icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4);
        }

        if (and_bit(boundary_bits, 2*d+1) == 0) {
          const VTYPE icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4);
        }
      }
    }
  }

  /// Apply to each cube which is facet adjacent to icubeA. (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5>
  void apply_to_cubes_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
        const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (and_bit(boundary_bits, 2*d) == 0) {
          const VTYPE icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
        }

        if (and_bit(boundary_bits, 2*d+1) == 0) {
          const VTYPE icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
        }
      }
    }
  }

  /// Apply to each cube which is facet adjacent to icubeA. (4+8 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
        const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (and_bit(boundary_bits, 2*d) == 0) {
          const VTYPE icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }

        if (and_bit(boundary_bits, 2*d+1) == 0) {
          const VTYPE icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
      }
    }
  }

  /// Apply to each cube which is facet adjacent to icubeA. (4+9 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9>
  void apply_to_cubes_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
        const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (and_bit(boundary_bits, 2*d) == 0) {
          const VTYPE icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9);
        }

        if (and_bit(boundary_bits, 2*d+1) == 0) {
          const VTYPE icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9);
        }
      }
    }
  }

  /// Apply to each cube which is facet adjacent to icubeA. (4+10 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9,
            typename ATYPE10>
  void apply_to_cubes_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9, ATYPE10 & arg10)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsF(); j++) {
        const VTYPE icubeB = grid.CubeNeighborF(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
          arg6, arg7, arg8, arg9, arg10);
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (and_bit(boundary_bits, 2*d) == 0) {
          const VTYPE icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9, arg10);
        }

        if (and_bit(boundary_bits, 2*d+1) == 0) {
          const VTYPE icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9, arg10);
        }
      }
    }
  }


  // **************************************************
  // FACET ADJACENT IN PLANE APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube in the plane which is facet adjacent to icubeA. 
  /// (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE, typename DIR_TYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5>
  void apply_to_cubes_in_plane_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits, const DIR_TYPE orth_dir,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.
      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (d != orth_dir) {
          VTYPE icubeB;

          icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);

          icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
        }
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (d != orth_dir) {

          if (and_bit(boundary_bits, 2*d) == 0) {
            const VTYPE icubeB = grid.PrevVertex(icubeA, d);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
          }

          if (and_bit(boundary_bits, 2*d+1) == 0) {
            const VTYPE icubeB = grid.NextVertex(icubeA, d);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
          }
        }
      }
    }
  }


  /// Apply to each cube in the plane which is facet adjacent to icubeA. 
  /// (4+6 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE, typename DIR_TYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6>
  void apply_to_cubes_in_plane_facet_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits, const DIR_TYPE orth_dir,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, 
   ATYPE5 & arg5, ATYPE6 & arg6)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.
      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (d != orth_dir) {
          VTYPE icubeB;

          icubeB = grid.PrevVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);

          icubeB = grid.NextVertex(icubeA, d);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);
        }
      }
    }
    else {

      for (DTYPE d = 0; d < grid.Dimension(); d++) {

        if (d != orth_dir) {

          if (and_bit(boundary_bits, 2*d) == 0) {
            const VTYPE icubeB = grid.PrevVertex(icubeA, d);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);
          }

          if (and_bit(boundary_bits, 2*d+1) == 0) {
            const VTYPE icubeB = grid.NextVertex(icubeA, d);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);
          }
        }
      }
    }
  }


  // **************************************************
  // EDGE ADJACENT APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is edge adjacent to icubeA. (4+4 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4>
  void apply_to_cubes_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
        const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4);
      }
    }
    else {

      for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

        BTYPE mask_d = (BTYPE(1) << edge_dir);
        for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

          if ((mask_d & j) != 0) { continue; }

          BTYPE mask;
          convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

          if ((boundary_bits & mask) == 0) {
            VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4);
          }
        }
      }
    }
  }

  /// Apply to each cube which is edge adjacent to icubeA. (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5>
  void apply_to_cubes_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
        const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
      }
    }
    else {

      for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

        BTYPE mask_d = (BTYPE(1) << edge_dir);
        for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

          if ((mask_d & j) != 0) { continue; }

          BTYPE mask;
          convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

          if ((boundary_bits & mask) == 0) {
            VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
          }
        }
      }
    }
  }

  /// Apply to each cube which is edge adjacent to icubeA. (4+8 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
        const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }
    }
    else {

      for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

        BTYPE mask_d = (BTYPE(1) << edge_dir);
        for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

          if ((mask_d & j) == 1) { continue; }

          BTYPE mask;
          convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

          if ((boundary_bits & mask) == 0) {
            VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
          }
        }
      }
    }
  }

  /// Apply to each cube which is edge adjacent to icubeA. (4+9 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9>
  void apply_to_cubes_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
        const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
      }
    }
    else {

      for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

        BTYPE mask_d = (BTYPE(1) << edge_dir);
        for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

          if ((mask_d & j) == 1) { continue; }

          BTYPE mask;
          convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

          if ((boundary_bits & mask) == 0) {
            VTYPE icubeB = 
              compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
              arg6, arg7, arg8, arg9);
          }
        }
      }
    }
  }

  /// Apply to each cube which is edge adjacent to icubeA. (4+10 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9,
            typename ATYPE10>
  void apply_to_cubes_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9, ATYPE10 & arg10)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsE(); j++) {
        const VTYPE icubeB = grid.CubeNeighborE(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
          arg6, arg7, arg8, arg9, arg10);
      }
    }
    else {

      for (DTYPE edge_dir = 0; edge_dir < grid.Dimension(); edge_dir++) {

        BTYPE mask_d = (BTYPE(1) << edge_dir);
        for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

          if ((mask_d & j) == 1) { continue; }

          BTYPE mask;
          convert2edge_bit_mask(grid.Dimension(), edge_dir, j, mask);

          if ((boundary_bits & mask) == 0) {
            VTYPE icubeB = 
              compute_edge_adjacent_cube(grid, icubeA, edge_dir, j);
            f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
              arg6, arg7, arg8, arg9, arg10);
          }
        }
      }
    }
  }


  // **************************************************
  // EDGE ADJACENT IN PLANE APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is edge adjacent to icubeA. (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE, typename DIR_TYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5>
  void apply_to_cubes_in_plane_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits, const DIR_TYPE orth_dir,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    const BTYPE mask_d = (BTYPE(1) << orth_dir);

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) != 0) { continue; }

        VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, orth_dir, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) != 0) { continue; }

        BTYPE mask;
        convert2edge_bit_mask(grid.Dimension(), orth_dir, j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, orth_dir, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
        }
      }
    }
  }

  /// Apply to each cube which is edge adjacent to icubeA. (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE, typename DIR_TYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6>
  void apply_to_cubes_in_plane_edge_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits, const DIR_TYPE orth_dir,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, 
   ATYPE5 & arg5, ATYPE6 & arg6)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    const BTYPE mask_d = (BTYPE(1) << orth_dir);

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) != 0) { continue; }

        VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, orth_dir, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        if ((mask_d & j) != 0) { continue; }

        BTYPE mask;
        convert2edge_bit_mask(grid.Dimension(), orth_dir, j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_edge_adjacent_cube(grid, icubeA, orth_dir, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6);
        }
      }
    }
  }


  // **************************************************
  // VERTEX ADJACENT APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is vertex adjacent to icubeA. (4+4 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4>
  void apply_to_cubes_vertex_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
        const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        BTYPE mask;
        convert2vertex_bit_mask(grid.Dimension(), j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4);
        }
      }
    }
  }

  /// Apply to each cube which is vertex adjacent to icubeA. (4+5 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5>
  void apply_to_cubes_vertex_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
        const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        BTYPE mask;
        convert2vertex_bit_mask(grid.Dimension(), j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5);
        }
      }
    }
  }

  /// Apply to each cube which is vertex adjacent to icubeA. (4+8 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_vertex_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
        const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        BTYPE mask;
        convert2vertex_bit_mask(grid.Dimension(), j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
      }
    }
  }

  /// Apply to each cube which is vertex adjacent to icubeA. (4+9 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9>
  void apply_to_cubes_vertex_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
        const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
          arg6, arg7, arg8, arg9);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        BTYPE mask;
        convert2vertex_bit_mask(grid.Dimension(), j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9);
        }
      }
    }
  }


  /// Apply to each cube which is vertex adjacent to icubeA. (4+9 arguments.)
  template <typename FTYPE, typename GRID_TYPE, 
            typename VTYPE, typename BTYPE,
            typename ATYPE1, typename ATYPE2, typename ATYPE3,
            typename ATYPE4, typename ATYPE5, typename ATYPE6, 
            typename ATYPE7, typename ATYPE8, typename ATYPE9,
            typename ATYPE10>
  void apply_to_cubes_vertex_adjacent_to
  (FTYPE f, const GRID_TYPE & grid, const VTYPE icubeA, 
   const BTYPE boundary_bits,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9, ATYPE10 & arg10)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    if (boundary_bits == 0) {
      // Cube icubeA is an interior cube.

      for (NUM_TYPE j = 0; j < grid.NumCubeNeighborsV(); j++) {
        const VTYPE icubeB = grid.CubeNeighborV(icubeA, j);
        f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
          arg6, arg7, arg8, arg9, arg10);
      }
    }
    else {

      for (NUM_TYPE j = 0; j < grid.NumCubeVertices(); j++) {

        BTYPE mask;
        convert2vertex_bit_mask(grid.Dimension(), j, mask);

        if ((boundary_bits & mask) == 0) {
          VTYPE icubeB = compute_vertex_adjacent_cube(grid, icubeA, j);
          f(icubeA, icubeB, arg1, arg2, arg3, arg4, arg5, 
            arg6, arg7, arg8, arg9, arg10);
        }
      }
    }
  }


  // **************************************************
  // FACET ADJACENT TO LIST APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is facet adjacent to a cube in gcube_list.
  /// (3 + 4 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2, 
            typename ATYPE3, typename ATYPE4>
  void apply_to_cubes_facet_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_facet_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits, arg1, arg2, arg3, arg4);
    }
  }

  /// Apply to each cube which is facet adjacent to a cube in gcube_list.
  /// (3 + 8 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_facet_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_facet_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }

  /// Apply to each cube which is facet adjacent to a cube in gcube_list.
  /// (3 + 9 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8,
            typename ATYPE9>
  void apply_to_cubes_facet_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_facet_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }


  /// Apply to each cube which is facet adjacent to a cube in gcube_list.
  /// (3 + 10 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8,
            typename ATYPE9, typename ATYPE10>
  void apply_to_cubes_facet_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9, ATYPE10 & arg10)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_facet_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10);
    }
  }


  // **************************************************
  // EDGE ADJACENT TO LIST APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is vertex adjacent to a in gcube_list.
  /// (3 + 4 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2, 
            typename ATYPE3, typename ATYPE4>
  void apply_to_cubes_edge_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_edge_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits, arg1, arg2, arg3, arg4);
    }
  }

  /// Apply to each cube which is edge adjacent to a in gcube_list.
  /// (3 + 8 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_edge_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_edge_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }

  /// Apply to each cube which is edge adjacent to a in gcube_list.
  /// (3 + 9 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8,
            typename ATYPE9>
  void apply_to_cubes_edge_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_edge_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }


  // **************************************************
  // VERTEX ADJACENT TO LIST APPLY TEMPLATES
  // **************************************************

  /// Apply to each cube which is vertex adjacent to a in gcube_list.
  /// (3 + 4 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2, 
            typename ATYPE3, typename ATYPE4>
  void apply_to_cubes_vertex_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_vertex_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits, arg1, arg2, arg3, arg4);
    }
  }


  /// Apply to each cube which is vertex adjacent to a in gcube_list.
  /// (3 + 8 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8>
  void apply_to_cubes_vertex_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_vertex_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    }
  }

  /// Apply to each cube which is vertex adjacent to a in gcube_list.
  /// (3 + 9 arguments.)
  template <typename FTYPE, typename ATYPE1, typename ATYPE2,
            typename ATYPE3,typename ATYPE4, typename ATYPE5,
            typename ATYPE6, typename ATYPE7, typename ATYPE8,
            typename ATYPE9>
  void apply_to_cubes_vertex_adjacent_to_list
  (FTYPE f, const std::vector<NUM_TYPE> & gcube_list, const ISOVERT & isovert,
   ATYPE1 & arg1, ATYPE2 & arg2, ATYPE3 & arg3, ATYPE4 & arg4, ATYPE5 & arg5,
   ATYPE6 & arg6, ATYPE7 & arg7, ATYPE8 & arg8, ATYPE9 & arg9)
  {
    for (NUM_TYPE i = 0; i < gcube_list.size(); i++) {

      NUM_TYPE to_gcube = gcube_list[i];
      VERTEX_INDEX to_cube = isovert.CubeIndex(to_gcube);

      BOUNDARY_BITS_TYPE boundary_bits = 
        isovert.gcube_list[to_gcube].boundary_bits;

      apply_to_cubes_vertex_adjacent_to
        (f, isovert.grid, to_cube, boundary_bits,
         arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    }
  }

}

#endif /* _SHREC_APPLY_H_ */

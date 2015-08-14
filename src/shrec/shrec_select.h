/// \file shrec_select.h
/// Select cubes containing isosurface vertices on sharp edges and corners.

/*
Copyright (C) 2012-2015 Arindam Bhattacharya and Rephael Wenger

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

#ifndef _SHREC_SELECT_H_
#define _SHREC_SELECT_H_

#include "shrec_isovert.h"

namespace SHREC {

  // **************************************************
  // MISMATCH TABLE
  // **************************************************

  class MISMATCH_ENTRY {

  public:
    GRID_COORD_DIFF_TYPE offset_coord[DIM3];
    INDEX_DIFF_TYPE offset_increment;
    std::vector<INDEX_DIFF_TYPE> blocking_location_increment;

    /// Only create a mismatch if the line containing the sharp edge
    ///   is within max_distance of the target cube.
    COORD_TYPE max_distance;

  public:
    MISMATCH_ENTRY(){}

    NUM_TYPE NumBlockingLocations() const
    { return(blocking_location_increment.size()); }
  };

  class MISMATCH_TABLE {

  public:
    SHARPISO_GRID grid;
    std::vector<MISMATCH_ENTRY> entry;

  public:
    /// Constructor.
    MISMATCH_TABLE(const SHARPISO_GRID & sharpiso_grid);

    /// Return maximum absolute value of any offset coordinate.
    inline GRID_COORD_TYPE MaxOffsetCoord() const
    { return(5); }

    /// Add entry.  Return entry number.
    NUM_TYPE AddEntry
    (const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
     const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance);

    /// Add blocking location to entry ientry.
    void AddBlockingLocation
    (const NUM_TYPE ientry, const GRID_COORD_DIFF_TYPE c0, 
     const GRID_COORD_DIFF_TYPE c1, const GRID_COORD_DIFF_TYPE c2);

    /// Compute blocking locations for entry ientry.
    void ComputeBlockingLocations(const NUM_TYPE ientry);

    /// Compute all blocking locations.
    void ComputeAllBlockingLocations();

    /// Add all permutations of c0, c1 and c2.
    void AddPermutations
    (const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
     const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance);

    /// Add all permutations of +/- c0, +/- c1 and +/- c2.
    void AddSignedPermutations
    (const GRID_COORD_DIFF_TYPE c0, const GRID_COORD_DIFF_TYPE c1, 
     const GRID_COORD_DIFF_TYPE c2, const COORD_TYPE max_distance);

    /// Clear mismatch table entries.
    void ClearEntries()
    { entry.clear(); }

    /// Set mismatch grid from given cube.
    /// @pre mismatch_grid dimensions must match grid dimensions.
    void SetMismatchGrid
    (const GRID_CUBE_DATA & gcubeA,
     const ISOVERT & isovert,
     const SHARPISO_BOOL_GRID & selected_grid,
     SHARPISO_BOOL_GRID & mismatch_grid) const;

    /// Return true if (cubeA_index + entry offset) is blocked 
    ///   by some selected cube.
    /// @pre (cubeA_index + entry offset) is contained in the grid.
    bool IsBlocked
    (const NUM_TYPE ientry, const VERTEX_INDEX cubeA_index,
     const SHARPISO_BOOL_GRID & selected_grid) const;
  };


  // **************************************************
  // SELECTION DATA STRUCTURES
  // **************************************************

  struct MOD6_ENTRY {
    int coord[DIM3];
  };

  class MOD6_LIST {

  public:
    std::vector<MOD6_ENTRY> entry;

  public:
    MOD6_LIST() {};

    NUM_TYPE NumEntries() const
    { return(entry.size()); }

    void AddEntry(const int ka, const int kb, const int kc);

    /// Add all permutations of ka, kb, kc.
    void AddPermutations(const int ka, const int kb, const int kc);

    /// Add all permutations of ka, kb, kc.
    void AddSignedPermutations(const int ka, const int kb, const int kc);
  };

  class SELECTION_DATA {

  protected:
    int bin_width;

  public:
    BIN_GRID<VERTEX_INDEX> bin_grid;
    SHARPISO_BOOL_GRID selected_grid;

    /// covered_grid.Scalar(icube) = true if icube is covered.
    SHARPISO_BOOL_GRID covered_grid;

    /// corner_covered_grid.Scalar(icube) = true,
    ///   if icube is covered by a corner cube.
    SHARPISO_BOOL_GRID corner_covered_grid;

    /// Cubes which should be avoided because of 3x3x3 region mismatches.
    SHARPISO_BOOL_GRID mismatch_grid;

    /// Table of offsets of cubes which do not match a given cube.
    MISMATCH_TABLE mismatch_table;

  public:
    SELECTION_DATA
    (const SHARPISO_GRID & grid, const SHARP_ISOVERT_PARAM & isovert_param);

    int BinWidth() const { return(bin_width); }
  };

  class SELECTION_DATA_MOD6:public SELECTION_DATA {

  protected:
    /// Add a cube to gcube_lists_mod6[][][].
    /// @pre 0 <= coord_mod6[d] <= 5 for d = 0,1,2.
    void AddCube(const GRID_COORD_TYPE coord_mod6[DIM3], 
                 const NUM_TYPE gcube_index);

    /// Set up mod6 lists.
    void SetMod6Lists();

  public:
    std::vector<VERTEX_INDEX> gcube_lists_mod6[6][6][6];

    /// List of permutations of (k,0,0) mod 6, sorted by distance to (0,0,0).
    MOD6_LIST k_0_0_mod6_list;

    /// List of permutations of (ka,kb,0) mod 6, sorted by distance to (0,0,0).
    MOD6_LIST ka_kb_0_mod6_list;

    /// List of permutations of (ka,kb,kc) mod 6, sorted by distance to (0,0,0).
    MOD6_LIST ka_kb_kc_mod6_list;

  public:
    SELECTION_DATA_MOD6
    (const SHARPISO_GRID & grid, const SHARP_ISOVERT_PARAM & isovert_param);

    int Modulus() const { return(6); }

    const std::vector<VERTEX_INDEX> & GCubeListsMod6(const MOD6_ENTRY & entry)
    {
      return(gcube_lists_mod6
             [entry.coord[0]][entry.coord[1]][entry.coord[2]]);
    }

    const std::vector<VERTEX_INDEX> & GCubeLists_k_0_0(const int i)
    { return(GCubeListsMod6(k_0_0_mod6_list.entry[i])); }

    const std::vector<VERTEX_INDEX> & GCubeLists_ka_kb_0(const int i)
    { return(GCubeListsMod6(ka_kb_0_mod6_list.entry[i])); }

    const std::vector<VERTEX_INDEX> & GCubeLists_ka_kb_kc(const int i)
    { return(GCubeListsMod6(ka_kb_kc_mod6_list.entry[i])); }

    void AddEdgeCubesToGCubeLists
      (const std::vector<GRID_CUBE_DATA> & gcube_list);
  };


  // **************************************************
  // SELECT ROUTINES
  // **************************************************

  /// Select sharp isosurface vertices.
  void select_sharp_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT & isovertData);

  /// Select sharp isosurface vertices.
  void select_sharp_isovert
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT & isovertData);

  /// Select sharp isosurface vertices using mod3 algorithm.
  void select_sharp_isovert_mod3
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT & isovertData);

  /// Select sharp isosurface vertices using mod6 algorithm.
  void select_sharp_isovert_mod6
  (const SHARPISO_SCALAR_GRID_BASE & scalar_grid,
   const GRADIENT_GRID_BASE & gradient_grid,
   const SCALAR_TYPE isovalue,
   const SHARP_ISOVERT_PARAM & isovert_param,
   ISOVERT & isovert);

}

#endif

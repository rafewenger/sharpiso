/// \file ijkmerge.txx
/// Merge data structures.

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008 Rephael Wenger

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

#ifndef _IJKMERGE_
#define _IJKMERGE_

#include <algorithm>

#include "ijk.txx"

namespace IJK {

  // **************************************************
  // class INTEGER_LIST
  // **************************************************

  /// Base class for list of non-negative integers with no duplicate entries.
  template < class INTEGER_TYPE, class LOC_TYPE> 
  class INTEGER_LIST {

  protected:
    typedef std::vector<INTEGER_TYPE> INTEGER_VECTOR;
    typedef typename INTEGER_VECTOR::iterator INTEGER_VECTOR_ITERATOR;

    INTEGER_TYPE min_int;     // minimum integer
    INTEGER_TYPE max_num_int; 
    // Maximum number of integers that can be inserted.
    // Inserted integers must be in {x: min_int <= x < min_int+num_int}.
    INTEGER_VECTOR list;
    LOC_TYPE * list_loc;     // list location
    bool * in_list;          // in_list[i] = true if i is in list

    void Init(const INTEGER_TYPE min_int, const INTEGER_TYPE max_num_int);
    void Init(const INTEGER_TYPE max_num_int) { Init(0, max_num_int); };
    void FreeListLoc();

    INTEGER_LIST() { Init(0, 0);}; // for derived classes

  public:
    INTEGER_LIST(const INTEGER_TYPE min_int,
                 const INTEGER_TYPE max_num_int) 
    { Init(min_int, max_num_int); };
    // Note:  INTEGER_LIST allocates and initializes memory of size max_num_int

    INTEGER_LIST(const INTEGER_TYPE max_num_int) { Init(max_num_int);};
    ~INTEGER_LIST() { FreeListLoc(); }; // destructor

    // get functions
    INTEGER_TYPE MinInt() const { return(min_int); };
    INTEGER_TYPE MaxNumInt() const { return(max_num_int); };
    inline INTEGER_TYPE Index(const INTEGER_TYPE i) const
    { return(i - min_int); };
    bool InList(const INTEGER_TYPE i) const 
    { return(in_list[Index(i)]); };
    LOC_TYPE ListLoc(const INTEGER_TYPE i) const 
    { return(list_loc[Index(i)]); };
    LOC_TYPE ListLength() const { return(list.size()); };
    INTEGER_TYPE List(const int i) const 
    { return(list[i]); };

    // set functions
    LOC_TYPE Insert(const INTEGER_TYPE i)
    {
      INTEGER_TYPE j = Index(i);

      if (!in_list[j]) {
        LOC_TYPE loc = this->list.size();
        this->list.push_back(i);
        this->list_loc[j] = loc;
        this->in_list[j] = true;
        return(loc);
      }
      else {
        return(this->list_loc[j]);
      }
    };

    void ResetMin(const INTEGER_TYPE min_int);

    void ClearList()
    {
      // set in_list[i] to false for all i
      for (INTEGER_VECTOR_ITERATOR iter = this->list.begin(); 
           iter != this->list.end(); iter++) {
        INTEGER_TYPE j = this->Index(*iter);
        in_list[j] = false;
      }

      this->list.clear();
    };
  };

  /// Initialize INTEGER_LIST.
  /// @param min_int = Minimum value of any integer in list. Precondition: min_in t is non-negative.
  /// @param max_num_int Maximum number of integers that can be inserted.  Precondition: max_num_int is non-negative.
  template < class INTEGER_TYPE, class LOC_TYPE> 
  void INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::
  Init(const INTEGER_TYPE min_int, const INTEGER_TYPE max_num_int)
  {
    IJK::PROCEDURE_ERROR error("INTEGER_LIST::Init");
    if (min_int < 0) { 
      error.AddMessage("Programming error. Parameter min_int must be non-negative.");
      throw error;
    }

    if (max_num_int < 0) {
      error.AddMessage("Programming error. Parameter max_num_int must be non-negative.");
      throw error;
    }

    this->min_int = min_int;
    this->max_num_int = max_num_int;
    list_loc = NULL;
    in_list = NULL;

    if (max_num_int > 0) { 
      list_loc = new LOC_TYPE[max_num_int]; 
      in_list = new bool[max_num_int];
      for (INTEGER_TYPE i = 0; i < max_num_int; i++) 
        { in_list[i] = false; };
    };

  }

  /// Free list loc.
  template <class INTEGER_TYPE, class LOC_TYPE> 
  void INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::FreeListLoc()
  {
    delete [] list_loc;
    list_loc = NULL;
    delete [] in_list;
    in_list = NULL;
    min_int = 0;
    max_num_int = 0;
  }

  /// Reset min integer.
  /// @param min_int = Minimum value of any integer in list. Precondition: min_int is non-negative.
  template < class INTEGER_TYPE, class LOC_TYPE> 
  void INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::
  ResetMin(const INTEGER_TYPE new_min_int)
  {
    ClearList();
    this->min_int = new_min_int;
  }

  // **************************************************
  // class INTEGER_COUNT
  // **************************************************

  /// Count non-negative integers.
  template < class INTEGER_TYPE, class LOC_TYPE, class NTYPE> 
  class INTEGER_COUNT:public INTEGER_LIST<INTEGER_TYPE,LOC_TYPE> {

  protected:
    typedef std::vector<NTYPE> COUNT_VECTOR;
    typedef typename COUNT_VECTOR::iterator COUNT_VECTOR_ITERATOR;

    COUNT_VECTOR count;      // count[i] = number of insertions of list[i]

    void Init(const INTEGER_TYPE min_int, const INTEGER_TYPE max_num_int);
    void Init(const INTEGER_TYPE max_num_int) { Init(0, max_num_int); };

    INTEGER_COUNT(){ Init(0, 0); };  // for derived classes

  public:
    INTEGER_COUNT(const INTEGER_TYPE min_int,
                  const INTEGER_TYPE max_num_int)
    { Init(min_int, max_num_int); };
    // Note:  INTEGER_COUNT allocates and initializes memory 
    //   of size (max_int-min_int+1);
    INTEGER_COUNT(const INTEGER_TYPE max_num_int) { Init(0, max_num_int); };
    ~INTEGER_COUNT() {}; // destructor

    // get functions
    NTYPE Count(const INTEGER_TYPE i) const { return(count[i]); };

    // set functions
    void Insert(const INTEGER_TYPE i)
    {
      LOC_TYPE loc = INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::Insert(i);
      
      if (loc < LOC_TYPE(count.size())) {	count[loc]++; }
      else { count.push_back(1); }
    };

    void ResetMin(const INTEGER_TYPE min_int);

    void ClearList()
    {
      count.clear();
      INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::ClearList();
    };
  };

  /// Initialize INTEGER_COUNT.
  /// @param min_int = Minimum value of any integer in list. Precondition: min_in t is non-negative.
  /// @param max_num_int Maximum number of integers that can be inserted.  Precondition: max_num_int is non-negative.
  template < class INTEGER_TYPE, class LOC_TYPE, class NTYPE> 
  void INTEGER_COUNT<INTEGER_TYPE,LOC_TYPE,NTYPE>::
  Init(const INTEGER_TYPE min_int, const INTEGER_TYPE max_num_int)
  {
    INTEGER_LIST<INTEGER_TYPE,LOC_TYPE>::Init(min_int, max_num_int);
  }

  /// Reset min integer.
  /// @param min_int = Minimum value of any integer in list. Precondition: min_int is non-negative.
  template < class INTEGER_TYPE, class LOC_TYPE, class NTYPE> 
  void INTEGER_COUNT<INTEGER_TYPE,LOC_TYPE,NTYPE>::
  ResetMin(const INTEGER_TYPE new_min_int)
  {
    ClearList();
    this->min_int = new_min_int;
  }

  // **************************************************
  // TEMPLATE merge_identical
  // **************************************************

  /// Merge identical values in an integer list
  /// @param list0 List of integers.
  /// @param list0_length Length of list0.
  /// @param[out] list1_nodup List without any duplicate values.
  /// @param[out] list0_map Mapping from list0 to locations in list1_nodup.
  /// @pre Array list0_map is preallocated to length at least list0_length.
  /// @param int_list Data structure of type INTEGER_LIST.
  ///                 int_list is modified by this procedure.
  /// @pre int_list.MinInt() <= list0[i] < int_list.MinInt()+int_list.MaxNumInt()
  ///      for all elements list0[i] of list0.
  template <typename ITYPE, typename NTYPE, typename MTYPE,
            typename INTEGER_LIST_TYPE>
  void merge_identical
  (const ITYPE * list0, const NTYPE list0_length,
   std::vector<ITYPE> & list1_nodup,
   MTYPE * list0_map, INTEGER_LIST_TYPE & int_list)
  {
    // initialize data structures
    int_list.ClearList();
    list1_nodup.clear();

    NTYPE loc;
    for (NTYPE i = 0; i < list0_length; i++) {
      loc = int_list.Insert(list0[i]);
      list0_map[i] = loc;
    }

    list1_nodup.resize(int_list.ListLength());
    for (NTYPE i = 0; i < int_list.ListLength(); i++)
      { list1_nodup[i] = int_list.List(i); }
  }

  /// Merge identical values in an integer list
  /// Version using std::vector.
  /// @param list0 List of integers.
  /// @param[out] list1_nodup List without any duplicate values.
  /// @param[out] list0_map Mapping from list0 to locations in list1_nodup.
  /// @param int_list Data structure of type INTEGER_LIST.
  ///                 int_list is modified by this procedure.
  /// @pre int_list.MinInt() <= list0[i] < int_list.MinInt()+int_list.MaxNumInt()
  ///      for all elements list0[i] of list0.
  template <typename ITYPE, typename MTYPE, typename INTEGER_LIST_TYPE>
  void merge_identical
  (const std::vector<ITYPE> list0, std::vector<ITYPE> & list1_nodup,
   std::vector<MTYPE> & list0_map, INTEGER_LIST_TYPE & int_list)
  {
    list0_map.resize(list0.size());

    if (list0.size() == 0) { return; };

    merge_identical(vector2pointer(list0), list0.size(),
                    list1_nodup, &(list0_map.front()),int_list);
  }

  // **************************************************
  // TEMPLATE ARRAY_LESS_THAN
  // **************************************************

  /// Comparison function template.
  /// Compare array[i0] with array[i1]
  template <class T> class ARRAY_LESS_THAN {
  protected:
    const T * array;

  public:
    ARRAY_LESS_THAN(const T * a):array(a) {};

    bool operator ()(const int i0, const int i1) const
    { 
      return(array[i0] < array[i1]);
    };
  };

  // **************************************************
  // TEMPLATE separate_pairs
  // **************************************************

  /// Separate array of pairs into array of first and array of second.
  /// @param first_list = array of first entries of pair.
  ///   Precondition: first_list is preallocated to length at least num_pair.
  /// @param second_list = array of second entries of pair.
  ///   Precondition: second_list is preallocated to length at least num_pair.
  template <class ITYPE, class SIZE_TYPE>
  void separate_pairs(const ITYPE * pair_list, const SIZE_TYPE num_pair,
                      ITYPE * first_list, ITYPE * second_list)
  {
    for (SIZE_TYPE i = 0; i < num_pair; i++) {
      first_list[i] = pair_list[2*i];
      second_list[i] = pair_list[2*i+1];
    }
  }


  // **************************************************
  // TEMPLATE compute_bin
  // **************************************************

  /// Compute bin.
  template <class NTYPE, class LTYPE>
  inline NTYPE compute_bin(const NTYPE x, const LTYPE log2_interval)
  { return(x >> log2_interval); }

  /// Compute number of elements in each bin.
  template <class ITYPE, class SIZE_TYPE, 
            class NTYPE, class NTYPE2, class NTYPE3>
  void compute_num_in_bins
  (const ITYPE * list, const SIZE_TYPE list_length, const NTYPE num_bin,
   const NTYPE2 log2_interval, NTYPE3 * count)
  {
    for (NTYPE ibin = 0; ibin < num_bin; ibin++) { count[ibin] = 0; }

    for (SIZE_TYPE i = 0; i < list_length; i++) {
      ITYPE ibin = compute_bin(list[i], log2_interval);
      count[ibin]++;
    }
  }


  /// Compute maximum of first elements in pair.
  /// Return 0 if number of pairs is 0.
  template <class ITYPE, class SIZE_TYPE>
  inline ITYPE compute_max_first
  (const ITYPE * pair_list, const SIZE_TYPE num_pair)
  {
    if (num_pair < 1) { return(0); }

    ITYPE max_first = pair_list[0];
    for (SIZE_TYPE i = 1; i < num_pair; i++) {
      if (max_first < pair_list[2*i]) { max_first = pair_list[2*i]; }
    }

    return(max_first);
  }

  // **************************************************
  // TEMPLATE separate_into_bins
  // **************************************************

  template <class ITYPE, class SIZE_TYPE, class NTYPE, 
            class NTYPE2, class LOC_TYPE>
  void separate_into_bins
  (const ITYPE * pair_list, const SIZE_TYPE num_pair, 
   const NTYPE num_bin, const NTYPE2 log2_interval,
   LOC_TYPE * list, LOC_TYPE * first, LOC_TYPE * count)
  {
    for (NTYPE ibin = 0; ibin < num_bin; ibin++)
      { count[ibin] = 0; };

    for (LOC_TYPE i = 0; i < num_pair; i++) {
      ITYPE ibin = compute_bin(pair_list[2*i], log2_interval);
      count[ibin]++;
    }

    IJK::ARRAY<LOC_TYPE> current(num_bin);

    first[0] = 0;
    for (NTYPE ibin = 1; ibin < num_bin; ibin++) {
      first[ibin] = first[ibin-1] + count[ibin-1];
    }

    std::copy(first, first+num_bin, current.Ptr());

    for (LOC_TYPE i = 0; i < num_pair; i++) {
      ITYPE ibin = compute_bin(pair_list[2*i], log2_interval);

      list[current[ibin]] = i;
      current[ibin]++;
    }

  }

  template <class ITYPE, class SIZE_TYPE, class NTYPE, 
            class NTYPE2, class LOC_TYPE>
  void separate_into_bins2
  (const ITYPE * list, const SIZE_TYPE list_length,
   const NTYPE num_bin, const NTYPE2 log2_interval,
   LOC_TYPE * bin_list, LOC_TYPE * first, LOC_TYPE * count)
  {
    compute_num_in_bins(list, list_length, num_bin, log2_interval, count);

    IJK::ARRAY<LOC_TYPE> current(num_bin);

    first[0] = 0;
    for (NTYPE ibin = 1; ibin < num_bin; ibin++) {
      first[ibin] = first[ibin-1] + count[ibin-1];
    }

    std::copy(first, first+num_bin, current.Ptr());

    for (SIZE_TYPE i = 0; i < list_length; i++) {
      ITYPE ibin = compute_bin(list[i], log2_interval);

      bin_list[current[ibin]] = i;
      current[ibin]++;
    }

  }

  /// Reorder list so that identical data entries are clustered together
  /// @param icount = integer count list data structure
  ///   Precondition: icount.MinInt() <= min data element in list.
  ///   Precondition: icount.MinInt()+icount.NumInt() > max data element in list.
  template <class ITYPE, class LTYPE, class SIZE_TYPE, 
            class NTYPE1, class NTYPE2>
  void cluster_identical
  (const ITYPE * data, const LTYPE * list, const SIZE_TYPE list_length,
   LTYPE * list2, INTEGER_COUNT<ITYPE, NTYPE1, NTYPE2> & icount)
  {
    if (list_length <= 0) { return; };

    // insert data in icount
    for (SIZE_TYPE i = 0; i < list_length; i++) {
      icount.Insert(data[list[i]]);
    }

    SIZE_TYPE nodup_list_length = icount.ListLength();
    // Note: nodup_list_length > 0  since list_length > 0.

    IJK::ARRAY<LTYPE> current(nodup_list_length);
    current[0] = 0;
    for (LTYPE j = 1; j < nodup_list_length; j++) {
      current[j] = current[j-1] + icount.Count(j-1);
    }

    for (SIZE_TYPE i = 0; i < list_length; i++) {
      ITYPE x = data[list[i]];
      LTYPE jloc = icount.ListLoc(x);
      list2[current[jloc]] = list[i];
      current[jloc]++;
    }
  }

  // **************************************************
  // TEMPLATE merge_bins
  // **************************************************

  template <class ITYPE, class SIZE_TYPE, class NTYPE, 
            class NTYPE2, class LOC_TYPE, class CTYPE>
  void merge_bins
  (const ITYPE * first_list, const ITYPE * second_list,
   const SIZE_TYPE num_pair, 
   const NTYPE num_bin, const NTYPE2 log2_interval, 
   const LOC_TYPE * bin_list, const LOC_TYPE * first, const CTYPE * bin_count,
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc)
  {
    ARRAY_LESS_THAN<ITYPE> second_list_less_than(second_list);

    NTYPE2 interval = (1 << log2_interval);

    INTEGER_COUNT<ITYPE, LOC_TYPE, LOC_TYPE> icount(interval);
    IJK::ARRAY<LOC_TYPE> cluster_list(num_pair);

    for (NTYPE ibin = 0; ibin < num_bin; ibin++) {

      if (bin_count[ibin] > 0) {

        ITYPE imin = ibin*interval;
        icount.ResetMin(imin);

        /// Cluster identical first entries.
        cluster_identical(first_list, bin_list+first[ibin], bin_count[ibin],
                          cluster_list.Ptr(), icount);

        /// Sort by second entry.
        LOC_TYPE * cluster_list_ptr = cluster_list.Ptr();
        for (LOC_TYPE j = 0; j < icount.ListLength(); j++) {
          std::sort(cluster_list_ptr, cluster_list_ptr+icount.Count(j),
                    second_list_less_than);
          cluster_list_ptr += icount.Count(j);
        }

        // insert in pair_list_nodup
        LOC_TYPE k0 = cluster_list[0];
        ITYPE x0 = first_list[k0];
        ITYPE y0 = second_list[k0];
        LOC_TYPE iloc = pair_list_nodup.size()/2;
        pair_list_loc[k0] = iloc;
        pair_list_nodup.push_back(x0);
        pair_list_nodup.push_back(y0);

        for (LOC_TYPE i = 1; i < bin_count[ibin]; i++) {
          LOC_TYPE k = cluster_list[i];
          ITYPE x = first_list[k];
          ITYPE y = second_list[k];
          if (x == x0 && y == y0) {
            pair_list_loc[k] = iloc;
          }
          else {
            iloc = pair_list_nodup.size()/2;
            pair_list_loc[k] = iloc;
            pair_list_nodup.push_back(x);
            pair_list_nodup.push_back(y);
            x0 = x;
            y0 = y;
          }
        }
      }
    }

  }

  // **************************************************
  // TEMPLATE merge_pairs_using_bins
  // **************************************************

  /// Merge pairs using bins.
  /// log2_interval = log base 2 of interval covered by each bin.
  ///   Precondition: log2_interval >= 0.
  template <class ITYPE, class STYPE, class STYPE2, class LOC_TYPE>
  void merge_pairs_using_bins
  (const ITYPE * pair_list, const STYPE num_pair, const ITYPE max_first,
   const STYPE2 log2_interval,
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc)
  {
    IJK::PROCEDURE_ERROR error("merge_pairs_using_bins");

    if (num_pair < 1) { return; };

    if (log2_interval < 0) {
      error.AddMessage("Programming error. log2_interval must be non-negative.");
      throw error;
    }

    STYPE num_bin = compute_bin(max_first, log2_interval)+1;

    IJK::ARRAY<ITYPE> first_list(num_pair);
    IJK::ARRAY<ITYPE> second_list(num_pair);
    separate_pairs(pair_list, num_pair, first_list.Ptr(), second_list.Ptr());

    IJK::ARRAY<LOC_TYPE> index_sorted(num_pair);
    for (STYPE i = 0; i < num_pair; i++) { index_sorted[i] = i; }

    IJK::ARRAY<LOC_TYPE> first(num_bin);
    IJK::ARRAY<LOC_TYPE> count(num_bin);
    separate_into_bins2
      (first_list.PtrConst(), num_pair, num_bin, log2_interval,
       index_sorted.Ptr(), first.Ptr(), count.Ptr());

    merge_bins(first_list.Ptr(), second_list.Ptr(), num_pair, 
               num_bin, log2_interval,
               index_sorted.Ptr(), first.Ptr(), count.Ptr(),
               pair_list_nodup, pair_list_loc);
  }

  /// Merge pairs using bins.
  /// log2_interval = log base 2 of interval covered by each bin.
  ///   Precondition: log2_interval >= 0.
  template <class ITYPE, class STYPE, class STYPE2, class LOC_TYPE>
  void merge_pairs_using_bins
  (const ITYPE * pair_list, const STYPE num_pair,
   const STYPE2 log2_interval,
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc)
  {
    if (num_pair < 1) { return; };

    ITYPE max_first = compute_max_first(pair_list, num_pair);

    merge_pairs_using_bins
      (pair_list, num_pair, max_first, log2_interval,
       pair_list_nodup, pair_list_loc);
  }

  // *******************************************************
  // CLASS MERGE_PAIRS_PARAMETERS
  // *******************************************************

  /// Parameters for MERGE_PAIRS algorithms.
  template <class NTYPE>
  class MERGE_PAIRS_PARAMETERS {

  protected:
    NTYPE log2_interval;
    static const NTYPE DEFAULT_LOG2_INTERVAL = 10;

  public:
    MERGE_PAIRS_PARAMETERS(const NTYPE log2_interval)
    { this->log2_interval = log2_interval; };
    MERGE_PAIRS_PARAMETERS()
    { this->log2_interval = DefaultLog2Interval(); };

    // Get functions.
    NTYPE Log2Interval() const { return(log2_interval); };
    static NTYPE DefaultLog2Interval() { return(DEFAULT_LOG2_INTERVAL); };

    // Set functions.
    void SetLog2Interval(const NTYPE log2_interval);
  };

  template <class NTYPE>
  void MERGE_PAIRS_PARAMETERS<NTYPE>::
  SetLog2Interval(const NTYPE log2_interval)
  { 
    IJK::PROCEDURE_ERROR error("MERGE_PAIRS_PARAMETERS::SetLog2Interval");

    if (log2_interval < 0) {
      error.AddMessage("Programming error. log2_interval must be non-negative.");
      throw error;
    }

    this->log2_interval = log2_interval; 
  }

  // **************************************************
  // TEMPLATE merge_pairs
  // **************************************************

  template <class ITYPE, class NTYPE, class LOC_TYPE, class PTYPE>
  void merge_pairs
  (const ITYPE * pair_list, const NTYPE num_pair,   
   std::vector<ITYPE> & pair_list_nodup,
   LOC_TYPE * pair_list_loc, 
   const MERGE_PAIRS_PARAMETERS<PTYPE> & parameters)
  {
    merge_pairs_using_bins
      (pair_list, num_pair, parameters.Log2Interval(),
       pair_list_nodup, pair_list_loc);
  }

  template <class ITYPE, class NTYPE, class LOC_TYPE>
  void merge_pairs
  (const ITYPE * pair_list, const NTYPE num_pair, 
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc)
  {
    MERGE_PAIRS_PARAMETERS<NTYPE> default_parameters;

    merge_pairs(pair_list, num_pair,
                pair_list_nodup, pair_list_loc, default_parameters);

  }

  template <class ITYPE, class LOC_TYPE, class PTYPE>
  void merge_pairs
  (const std::vector<ITYPE> & pair_list,
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc,
   const MERGE_PAIRS_PARAMETERS<PTYPE> & parameters)
  {
    merge_pairs(&(pair_list[0]), pair_list.size()/2,
                pair_list_nodup, pair_list_loc, parameters);
  }

  template <class ITYPE, class LOC_TYPE>
  void merge_pairs
  (const std::vector<ITYPE> & pair_list,
   std::vector<ITYPE> & pair_list_nodup, LOC_TYPE * pair_list_loc)
  {
    MERGE_PAIRS_PARAMETERS<int> default_parameters;

    merge_pairs(pair_list, pair_list_nodup, pair_list_loc, 
                default_parameters);
  }

  // **************************************************
  // check function
  // **************************************************

  /// Check for error in pair_list_nodup and pair_list_loc.
  /// Return false and set error message if error found.
  template <class ITYPE, class SIZE_TYPE, class LOC_TYPE>
  bool check_pair_list
  (const ITYPE * pair_list, const SIZE_TYPE num_pair,
   const ITYPE * pair_list_nodup, const SIZE_TYPE num_noup_pair,
   LOC_TYPE * pair_list_loc, IJK::ERROR & error)
  {
    for (SIZE_TYPE i = 0; i < num_pair; i++) {
      LOC_TYPE iloc = pair_list_loc[i];
      if (pair_list[2*i] != pair_list_nodup[2*iloc] ||
          pair_list[2*i+1] != pair_list_nodup[2*iloc+1]) {
        error.AddMessage("Error in pair_list_loc.");
        error.AddMessage("pair_list_loc[", i, "] = ", pair_list_loc[i], ".");
        error.AddMessage("pair ", i, " = (", pair_list[2*i],
                         ",", pair_list[2*i+1], ").");
        error.AddMessage("nodup pair ", iloc, " = (", 
                         pair_list_nodup[2*iloc],
                         ",", pair_list_nodup[2*iloc+1], ").");

        return(false);
      }
    }

    return(true);
  }

  /// Check for error in pair_list_nodup and pair_list_loc.
  /// Return false and set error message if error found.
  template <class ITYPE, class LOC_TYPE>
  bool check_pair_list
  (const std::vector<ITYPE> & pair_list,
   const std::vector<ITYPE> & pair_list_nodup,
   LOC_TYPE * pair_list_loc, IJK::ERROR & error)
  {
    return(check_pair_list(&(pair_list[0]), pair_list.size()/2,
                           &(pair_list_nodup[0]), pair_list_nodup.size()/2,
                           pair_list_loc, error));
  }

}

#endif

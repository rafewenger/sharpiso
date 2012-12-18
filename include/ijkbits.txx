/// \file ijkbits.txx
/// ijk templates for bit operations
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

#ifndef _IJK_BITS_
#define _IJK_BITS_

namespace IJK {

  // **************************************************
  // BIT OPERATIONS
  // **************************************************

  template <typename ITYPE, typename BTYPE, typename DTYPE,
            typename NTYPE, typename ETYPE>
  void convert2base(const ITYPE ival, const BTYPE base, 
                    DTYPE * digit, const NTYPE max_num_digits,
                    ETYPE & error)
  {
    ITYPE jval = ival;
    for (NTYPE i = 0; i < max_num_digits; i++) {
      digit[i] = jval % base;
      jval = jval/base;
    }

    if (jval != 0) {
      error.AddMessage("Error converting ", ival, " to base ", base, ".");
      error.AddMessage("Output has more than ", max_num_digits, " digits.");

      throw error;
    };
  }

  template <typename ITYPE, typename NTYPE, typename NTYPE2, typename NTYPE3>
  void count_bits(const ITYPE ival, const NTYPE num_bits,
                  NTYPE2 & num_zeros, NTYPE3 & num_ones)
  {
    const ITYPE base = 2;

    num_zeros = 0;
    num_ones = 0;

    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      int bit = jval % base;
      if (bit == 0) { num_zeros++; }
      else { num_ones++; }
      jval = jval/base;
    }
  }

  /// Reverse order of bits in ival.
  template <typename ITYPE, typename NTYPE>
  ITYPE reverse_bits(const ITYPE ival, const NTYPE num_bits)
  {
    const ITYPE base = 2;

    ITYPE reverse_val = 0;
    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      ITYPE bit = jval % base;
      reverse_val = base*reverse_val;
      reverse_val = (reverse_val | bit);
      jval = jval/base;
    }

    return(reverse_val);
  }

  /// Return true if (ival == reverse_bits(ival,num_bits))
  template <typename ITYPE, typename NTYPE>
  bool equals_reverse_bits(const ITYPE ival, const NTYPE num_bits)
  {
    ITYPE reverse_val = reverse_bits(ival, num_bits);
    return((ival == reverse_val));
  }

  /// Return index of first one bit.
  /// Return num_bits if all bits are zero.
  template <typename ITYPE, typename NTYPE>
  NTYPE get_first_one_bit(const ITYPE val, const NTYPE num_bits)
  {
    ITYPE mask = ITYPE(1);
    for (NTYPE i = 0; i < num_bits; i++) {
      if ((val & mask) != 0) { return(i); }
      mask = (mask << ITYPE(1));
    }
    return(num_bits);
  }
    
}

#endif

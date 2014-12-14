/// \file sharpiso_array.txx
/// Array templates for sharp isosurface processing.
/// Version v0.1.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2014 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _SHARPISO_ARRAY_
#define _SHARPISO_ARRAY_

namespace SHARPISO {

  // *******************************************************
  // TEMPLATE FIXED_ARRAY
  // *******************************************************

  template <const int MAX_NUM_ELEMENTS, typename ETYPE, typename NTYPE>
  class FIXED_ARRAY {

  protected:
    ETYPE element[MAX_NUM_ELEMENTS];
    NTYPE num_elements;
  

  public:
    FIXED_ARRAY() 
    { num_elements = 0; };

    int MaxNumElements() const
    { return(MAX_NUM_ELEMENTS); }

    NTYPE NumElements() const
    { return(num_elements); }

    bool Contains(const ETYPE x) const
    {
      for (NTYPE i = 0; i < NumElements(); i++)
        { if (x == element[i]) { return(true); } }
      return(false);
    }
        
    template <typename ITYPE>
    ETYPE & operator [] (const ITYPE i) { return(*(element+i)); }

    template <typename ITYPE>
    ETYPE operator [] (const ITYPE i) const { return(*(element+i)); }

    const ETYPE * PtrConst() const { return(element); }
    ETYPE * Ptr() { return(element); };

    /// Add element and increase num_elements.
    /// @pre num_elements < MAX_NUM_ELEMENTS.  (No check for this condition).
    void PushBack(const ETYPE x) {
      element[num_elements] = x;
      num_elements++;
    }

    void Clear()
    { num_elements = 0; }
  };

}

#endif

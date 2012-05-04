/// \file ijkparse.txx
/// ijk templates for converting strings to values/arrays
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

#ifndef _IJKSTRING_
#define _IJKSTRING_

#include <sstream>
#include <vector>


/// ijk templates for manipulating strings.
namespace IJK {


  // **************************************************
  // Convert strings to numeric values/arrays.
  // **************************************************

  namespace {

    /// Remove trailing blanks from string s
    inline void remove_trailing_blanks(std::string & s)
    {
      size_t pos = 0;
      for (size_t i = 0; i < s.length(); i++) {
        if (!isspace(s[i])) { pos = i+1; }
      }
      if (pos < s.length()) { s.erase(pos); };
        
    }
  };

  /// Convert string to value
  template <typename VTYPE>
  bool string2val(const char * s, VTYPE & val)
  {
    std::istringstream v_string;

    std::string s2 = s;
    remove_trailing_blanks(s2);

    v_string.str(s2);

    v_string >> val;

    if (!v_string.eof()) 
      { return(false); }
    else
      { return(true); }
  }

  /// Convert string to array of elements
  template <typename ETYPE>
  bool string2vector(const char * s, std::vector<ETYPE> & v)
  {
    std::istringstream v_string;

    v.clear();

    std::string s2 = s;
    remove_trailing_blanks(s2);

    v_string.str(s2);
    while (v_string.good()) {
      ETYPE x;
      v_string >> x;
      v.push_back(x);
    }

    if (!v_string.eof()) 
      { return(false); }
    else
      { return(true); }
  }

  // **************************************************
  // Split string into prefix and suffix.
  // **************************************************

  template <typename STRING_TYPE>
  void split_string(const STRING_TYPE & s, const char c,
                    STRING_TYPE & prefix, STRING_TYPE & suffix)
  // split string at last occurrence of character c into prefix and suffix
  {
    typename STRING_TYPE::size_type i = s.rfind(c);
    if (i == STRING_TYPE::npos) {
      prefix = s;
      suffix = "";
    }
    else {
      if (i > 0) { prefix = s.substr(0,i); }
      else { prefix = ""; };

      if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
      else { suffix = ""; };
    }
  }

  template <typename STRING_TYPE>
  void split_string(const char * s, const char c,
                    STRING_TYPE & prefix, STRING_TYPE & suffix)
  {
    split_string(STRING_TYPE(s), c, prefix, suffix);
  }

}

#endif

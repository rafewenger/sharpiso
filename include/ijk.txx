/// \file ijk.txx
/// ijk templates defining general ijk objects, i.e. classes BOX and ERROR
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008-2013 Rephael Wenger

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

#ifndef _IJK_
#define _IJK_

#include <limits>
#include <iostream>
#include <sstream>
#include <vector>
#include <utility>

/// General purpose ijk classes and routines.
/// Box, line segments, array and error classes.
namespace IJK {

  // **************************************************
  // TEMPLATE CLASS BOX
  // **************************************************

  /// \brief Axis-parallel box data structure.  
  /// Represents box by minimum and maximum coordinates.
  template <typename COORD_TYPE> class BOX {
  protected:
    int dimension;
    COORD_TYPE * min_coord;
    COORD_TYPE * max_coord;

    void Init();
    void FreeAll();
    
  public:
    BOX() { Init(); };
    BOX(const int dimension);
    ~BOX() { FreeAll(); };
    BOX(const BOX & box);                 // copy constructor
    const BOX & operator = (const BOX &); // copy assignment

    // get functions
    int Dimension() const { return(dimension); };
    COORD_TYPE MinCoord(const int d) const
    { return(min_coord[d]); };
    const COORD_TYPE * MinCoord() const
    { return(min_coord); };
    COORD_TYPE MaxCoord(const int d) const
    { return(max_coord[d]); };
    const COORD_TYPE * MaxCoord() const
    { return(max_coord); };
    COORD_TYPE AxisSize(const int d) const
    { return(max_coord[d]+1-min_coord[d]); }

    template <typename CTYPE2>
    bool Contains(const CTYPE2 * coord) const;

    // set functions
    void SetDimension(const int d);
    void SetMinCoord(const int d, const COORD_TYPE c)
    { min_coord[d] = c; };
    void SetMaxCoord(const int d, const COORD_TYPE c)
    { max_coord[d] = c; };
    template <typename CTYPE2>                // set min_coord[] to coord[]
    void SetMinCoord(const CTYPE2 * coord);
    template <typename CTYPE2>                // set max_coord[] to coord[]
    void SetMaxCoord(const CTYPE2 * coord);
    template <typename CTYPE2, typename CTYPE3>
    void SetCoord(const CTYPE2 * minc, const CTYPE3 * maxc);
    void SetAllMinCoord(const COORD_TYPE c);  // set all min coord to c
    void SetAllMaxCoord(const COORD_TYPE c);  // set all max coord to c
  };

  // **************************************************
  // CLASS LINE_SEGMENT & ASSOCIATED FUNCTIONS
  // **************************************************

  /// \brief Line segment between two grid vertices
  template <typename VTYPE>
  class LINE_SEGMENT:public std::pair<VTYPE,VTYPE> {  

  public:
    LINE_SEGMENT() {};  // constructor
    LINE_SEGMENT(const VTYPE iv0, const VTYPE iv1)
    { SetEnd(iv0, iv1); };

    // set functions

    /// Set line segment endpoints. 
    void SetEnd(const VTYPE iv0, const VTYPE iv1) 
    { this->first = iv0; this->second = iv1; Order(); };

    /// Order endpoints so that  V0() <= V1().
    void Order()                      
    { if (V0() > V1()) { std::swap(this->first, this->second); }; };


    // get functions
    VTYPE V0() const { return(this->first); };  ///< Return endpoint 0.
    VTYPE V1() const { return(this->second); };  ///< Return endpoint 1.
    template <typename ITYPE> 
    VTYPE V(const ITYPE i) const       ///< Return endpoint i.
    {
      if (i == 0) { return(V0()); }
      else { return(V1()); }
    }

  };

  /// Return true if (V0 < V1) for every pair (V0,V1)
  ///    of line segment endpoints
  template <typename VTYPE>
  bool is_ordered(const std::vector< LINE_SEGMENT<VTYPE> > & list)
  {
    typename std::vector< LINE_SEGMENT<VTYPE> > ::const_iterator pos;
    for (pos = list.begin(); pos != list.end(); ++pos) {
      if (pos->V0() > pos->V1())
        return(false);
    }
    return(true);
  }


  // **************************************************
  // TEMPLATE CLASS ARRAY AND ARRAY_L
  // **************************************************

  /// \brief Simple array class for creating static arrays.
  ///
  /// Making this a class guarantees that when a program leaves
  /// the function where this class is declared, the array memory is freed.
  template <typename ETYPE> class ARRAY {
  protected:
    ETYPE * element;

    template <typename LTYPE>
    void Init(const LTYPE array_length)
    { element = new ETYPE[array_length]; }

    template <typename LTYPE>
    void Init(const LTYPE array_length, const ETYPE init_value)
    { 
      Init(array_length);
      for (LTYPE i = 0; i < array_length; i++) { element[i] = init_value; };
    }

  public:
    template <typename LTYPE>
    ARRAY(const LTYPE array_length) // constructor
    { Init(array_length); }
    template <typename LTYPE>
    ARRAY(const LTYPE array_length, const ETYPE init_value) // constructor
    { Init(array_length, init_value); }
    ~ARRAY();

    // get functions
    ETYPE * Ptr() { return(element); };
    const ETYPE * PtrConst() const { return(element); }
    template <typename ITYPE>
    ETYPE & operator [] (const ITYPE i) { return(*(element+i)); }
    template <typename ITYPE>
    ETYPE operator [] (const ITYPE i) const { return(*(element+i)); }

    // Free function
    void Free();
  };

  /// Class array with length stored
  template <typename ETYPE, typename LTYPE>
  class ARRAY_L:public ARRAY<ETYPE> {
  protected:
    LTYPE length;

    void Init(const LTYPE length) 
    { this->length = length; }

  public:
    ARRAY_L(const LTYPE length):ARRAY<ETYPE>(length) 
      // constructor
    { Init(length); };

    ARRAY_L(const LTYPE length, const ETYPE init_value):
      ARRAY<ETYPE>(length, init_value)
      // constructor
    { Init(length); };

    // get functions
    LTYPE Length() const { return(length); };
    ETYPE * End() { return(this->element+length); };
  };

  // **************************************************
  // TEMPLATE CLASS CONSTANT
  // **************************************************

  /// Class CONSTANT always returns the same value
  template <typename ITYPE, typename CTYPE>
  class CONSTANT {
  protected:
    CTYPE c;           ///< The constant value.

  public:
    CONSTANT(const CTYPE c) { this->c = c; };

    // Get functions.
    CTYPE operator [] (const ITYPE i) const { return(c); };
    CTYPE operator () (const ITYPE i) const { return(c); };
  };

  // **************************************************
  // CONVERT A C++ VECTOR TO A POINTER
  // **************************************************

  /// \brief Returns C++ pointer to C array storing in vector data.
  /// If vector is empty, returns NULL pointser.
  /// If vector is not empty, returns pointer to a C array.
  template <typename T>
  const T * vector2pointer(const std::vector<T> & v)
  {
    if (v.empty()) { return(NULL); }
    else { return(&(v.front())); }
  }

  // **************************************************
  // SET ARRAY
  // **************************************************

  /// Set all elements of array a to x.
  template <typename NTYPE, typename ETYPE0, typename ETYPE1>
  void set_c_array(const NTYPE alength, const ETYPE0 x, ETYPE1 a[])
  {
    for (NTYPE i = 0; i < alength; i++)
      { a[i] = x; }
  }

  /// Set a[i] to x if flag[i] is true.
  template <typename NTYPE, typename ETYPE0, typename ETYPE1,
            typename FTYPE>
  void set_c_array(const NTYPE alength, const ETYPE0 x, 
                   const FTYPE & flag, ETYPE1 a[])
  {
    for (NTYPE i = 0; i < alength; i++)
      if (flag[i]) 
        { a[i] = x; }
  }

  /// Set a[i] to x if flag[i] is true.
  template <typename ETYPE0, typename ETYPE1, typename FTYPE>
  void set_array(const ETYPE0 x, const FTYPE & flag, std::vector<ETYPE1> & a)
  {
    if (a.empty()) { return; };
    set_c_array(a.size(), x, flag, &a.front());
  }

  // **************************************************
  // INTEGER POWER FUNCTION
  // **************************************************

  /// *** DEPRECATED ***
  /// \brief Integer power function.
  /// Return (base)^p.
  template <typename ITYPE0, typename ITYPE1>
  ITYPE0 int_power(const ITYPE0 base, const ITYPE1 p)
  {
    ITYPE0 result = 1;
    for (ITYPE1 k = 0; k < p; k++) {
      result *= base;
    }

    return(result);
  }

  /// \brief Integer power function.
  /// @param[out] result Equals (base)^p.
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void int_power(const ITYPE0 base, const ITYPE1 p, 
                 ITYPE2 & result)
  {
    result = 1;
    for (ITYPE1 k = 0; k < p; k++) {
      result *= base;
    }
  }

  /// \brief Integer power function with error checking.
  /// @param[out] result Equals (base)^p.
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
            typename ETYPE>
  void int_power(const ITYPE0 base, const ITYPE1 p, 
                 ITYPE2 & result, ETYPE & error)
  {
    const ITYPE2 result_bound = (std::numeric_limits<ITYPE2>::max()/base);
    
    result = 1;
    for (ITYPE1 k = 0; k < p; k++) {

      if (result > result_bound) {
        error.AddMessage
          ("Result out of bounds. ", base, "^", p,
           " is larger than ", std::numeric_limits<ITYPE2>::max(), ".");
        throw error;
      }
      result *= base;
    }
  }

  // **************************************************
  // COUNT TEMPLATE FUNCTION
  // **************************************************

  /// Count number of elements ge X in array a[].
  /// @tparam ATYPE Array type. Could be C array or stl vector.
  template <int X, typename ATYPE, typename NTYPE>
  NTYPE count_ge(const ATYPE a, const NTYPE length)
  {
    NTYPE num_ge = 0;
    for (NTYPE i = 0; i < length; i++) {
      if (a[i] >= X) { num_ge++; }
    }
    return(num_ge);
  }

  /// Count number of elements ge X in stl vector v.
  template <int X, typename T>
  typename std::vector<T>::size_type count_ge(const std::vector<T> & v)
  {
    return(count_ge<X>(v, v.size()));
  }

  // **************************************************
  // PUSH MULTIPLE ELEMENTS ON A C++ VECTOR
  // **************************************************

  /// Push two elements on a C++ vector.
  template <typename T>
  inline void push_back(const T & a0, const T & a1, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
  }

  /// Push three elements on a C++ vector.
  template <typename T>
  inline void push_back
  (const T & a0, const T & a1, const T & a2, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
    v.push_back(a2);
  }

  /// Push four elements on a C++ vector.
  template <typename T>
  inline void push_back
  (const T & a0, const T & a1, const T & a2, const T & a3, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
    v.push_back(a2);
    v.push_back(a3);
  }

  // **************************************************
  // ERROR CLASSES
  // **************************************************

  /// \brief Compose string from two elements.
  /// Format string using output string stream.
  template < typename T1, typename T2> 
  std::string compose_string(T1 a1, T2 a2)
  {
    std::ostringstream os0;

    os0 << a1 << a2;
    return(os0.str());
  }

  /// \brief Compose string from three elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3;
    return(os0.str());
  }

  /// \brief Compose string from four elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4;
    return(os0.str());
  }

  /// \brief Compose string from five elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4, T5 a5)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5;
    return(os0.str());
  }

  /// \brief Compose string from seven elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6;
    return(os0.str());
  }

  /// \brief Compose string from seven elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7;
    return(os0.str());
  }

  /// \brief Compose string from eight elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8>
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8;
    return(os0.str());
  }

  /// \brief Compose string from nine elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9;
    return(os0.str());
  }

  /// \brief Compose string from ten elements.
  /// Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9, typename T10 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9,
   T10 a10)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 << a10;
    return(os0.str());
  }

  /// Error class.
  class ERROR {

  protected:
    std::vector<std::string> msg;   // error messages

  public:
    ERROR(){};
    ERROR(const char * error_msg) { AddMessage(error_msg); };
    ERROR(const std::string & error_msg) { AddMessage(error_msg); };

    ERROR(const ERROR & error)                      // copy constructor
    {
      for (int i = 0; i < error.NumMessages(); i++)
        { msg.push_back(error.Message(i)); };
    };

    const ERROR & operator = (const ERROR & right)  // copy assignment
    {
      if (&right != this) {
        for (int i = 0; i < right.NumMessages(); i++)
          { msg.push_back(right.Message(i)); };
      }

      return(*this);
    };

    // get functions
    std::string Message(const int i) const
    {
      if (i >= 0 && i < NumMessages()) { return(msg[i]); }
      else { return(""); }
    }
    int NumMessages() const { return(msg.size()); };

    // get procedure error messages
    std::string ProcMessage(const std::string & procedure_name)
    {
      std::string error_msg = 
        "Error detected in procedure: " +  procedure_name + ".";
      return(error_msg);
    }
    std::string ProcMessage(const char * procedure_name)
    { return(ProcMessage(std::string(procedure_name))); }

    // set functions
    void AddMessage(const std::string & error_msg) ///< Add error message.
    { msg.push_back(error_msg); };
    void AddMessage(const char * error_msg)     
    { AddMessage(std::string(error_msg)); };

    // Add message with two elements.
    template <typename T1, typename T2>
    void AddMessage(T1 m1, T2 m2)
    { AddMessage(compose_string(m1, m2)); };

    // Add message with three elements
    template <typename T1, typename T2, typename T3>
    void AddMessage(T1 m1, T2 m2, T3 m3)
    { AddMessage(compose_string(m1, m2, m3)); };

    // Add message with four elements
    template <typename T1, typename T2, typename T3, typename T4>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4)
    { AddMessage(compose_string(m1, m2, m3, m4)); };

    // Add message with five elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4, T5 m5)
    { AddMessage(compose_string(m1, m2, m3, m4, m5)); };

    // Add message with six elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6)
    { AddMessage(compose_string(m1, m2, m3, m4, m5, m6)); };

    // Add message with seven elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7)
    { AddMessage(compose_string(m1, m2, m3, m4, m5, m6, m7)); };

    // Add message with eight elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8)); };

    // Add message with nine elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8, m9)); };

    // Add message with ten elements
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9, typename T10>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9,
     T10 m10)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)); };

    void AddProcMessage(const char * procname)
    { AddMessage(ProcMessage(procname)); };
    void AddProcMessage(const std::string & procname)
    { AddMessage(ProcMessage(procname)); };
    void SetMessage(const int i, const std::string & error_msg)
    {
      if (i >= 0 && i < NumMessages())
        msg[i] = error_msg;
    };
    void SetMessage(const int i, const char * error_msg)
    { SetMessage(i, std::string(error_msg)); };
    void SetProcMessage(const int i, const char * procname)
    { SetMessage(i, ProcMessage(procname)); };
    void SetProcMessage(const int i, const std::string & procname)
    { SetMessage(i, ProcMessage(procname)); };
    void ClearAll() { msg.clear(); };

    ERROR & operator ()(const std::string & msg1)
    {
      AddMessage(msg1);
      return(*this);
    };

    ERROR & operator ()(const std::string & msg1, const std::string & msg2)
    {
      AddMessage(msg1);
      AddMessage(msg2);
      return(*this);
    }

    // print error messages to ostream out
    void Print(std::ostream & out) const
    {
      for (int i = 0; i < NumMessages(); i++) {
        out << Message(i) << std::endl;
      }
    }

  };

  /// \brief Procedure error class.
  /// Includes procedure name as first error message.
  class PROCEDURE_ERROR:public ERROR {

  protected:
    std::string procedure_name;

  public:
    PROCEDURE_ERROR(const char * procedure_name)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
    };
    PROCEDURE_ERROR(const char * procedure_name, const char * error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };
    PROCEDURE_ERROR(const std::string & procedure_name, 
                    const std::string & error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };

  };

  // **************************************************
  // CHECK
  // **************************************************

  /// Return true if array memory is allocated.
  /// Return false and set error message if array is NULL.
  template <typename T>
  bool check_array_allocated
  (const T * array, const char * array_name, ERROR & error)
  {
    if (array == NULL) {
      error.AddMessage("Programming error. Memory for array ",
                       array_name, "[] not allocated.");
      return(false);
    }

    return(true);
  }

  /// Return true if ptr is NULL
  template <typename T>
  bool check_is_NULL
  (const T * ptr, const char * variable_name, ERROR & error)
  {
    if (ptr != NULL) {
      error.AddMessage
        ("Programming error.  Previously allocated memory for variable ",
         variable_name, " not released.");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS BOX MEMBER FUNCTIONS
  // **************************************************

  template <typename T> void BOX<T>::Init()
  {
    dimension = 0;
    min_coord = 0;
    max_coord = 0;
  }

  template <typename T> void BOX<T>::FreeAll()
  {
    dimension = 0;
    if (min_coord != NULL) { delete [] min_coord; }
    if (max_coord != NULL) { delete [] max_coord; }
    min_coord = NULL;
    max_coord = NULL;
  }

  template <typename T> BOX<T>::BOX(const int dimension)
  { Init(); SetDimension(dimension); }

  template <typename COORD_TYPE>
  template <typename CTYPE2>
  bool BOX<COORD_TYPE>::Contains(const CTYPE2 * coord) const
  {
    for (int d = 0; d < dimension; d++) {
      if (coord[d] < min_coord[d]) { return(false); }
      if (coord[d] > max_coord[d]) { return(false); }
    }
    return(true);
  }

  template <typename COORD_TYPE> 
  void BOX<COORD_TYPE>::SetDimension(const int d)
  {
    FreeAll();
    Init();
    if (d <= 0) return;
    min_coord = new COORD_TYPE[d];
    max_coord = new COORD_TYPE[d];
    dimension = d;
  }

  template <typename COORD_TYPE>
  template <typename CTYPE2>
  void BOX<COORD_TYPE>::
  SetMinCoord(const CTYPE2 * coord)
  {
    for (int d = 0; d < dimension; d++)
      SetMinCoord(d, coord[d]);
  }

  template <typename COORD_TYPE>
  template <typename CTYPE2>
  void BOX<COORD_TYPE>::
  SetMaxCoord(const CTYPE2 * coord)
  {
    for (int d = 0; d < dimension; d++)
      SetMaxCoord(d, coord[d]);
  }

  template <typename COORD_TYPE>
  template <typename CTYPE2, typename CTYPE3>
  void BOX<COORD_TYPE>::SetCoord
  (const CTYPE2 * minc, const CTYPE3 * maxc)
  {
    SetMinCoord(minc);
    SetMaxCoord(maxc);
  }

  template <typename COORD_TYPE>
  void BOX<COORD_TYPE>::SetAllMinCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMinCoord(d, c);
  }

  template <typename COORD_TYPE>
  void BOX<COORD_TYPE>::SetAllMaxCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMaxCoord(d, c);
  }

  // copy constructor template
  template <typename COORD_TYPE> 
  BOX<COORD_TYPE>::BOX(const BOX<COORD_TYPE> & box)
  {
    Init();
    SetDimension(box.Dimension());
    for (int d = 0; d < box.Dimension(); d++) {
      SetMinCoord(d, box.MinCoord(d));
      SetMaxCoord(d, box.MaxCoord(d));
    }
  }

  // copy assigment template
  template <typename COORD_TYPE> 
  const BOX<COORD_TYPE> & BOX<COORD_TYPE>::operator = 
  (const BOX<COORD_TYPE> & right)
  {
    if (&right != this) {         // avoid self-assignment
      FreeAll();
      SetDimension(right.Dimension());
      for (int d = 0; d < right.Dimension(); d++) {
        SetMinCoord(d, right.MinCoord(d));
        SetMaxCoord(d, right.MaxCoord(d));
      }
    }
    return *this;
  }

  // **************************************************
  // TEMPLATE CLASS ARRAY MEMBER FUNCTIONS
  // **************************************************

  template <typename ETYPE> ARRAY<ETYPE>::~ARRAY()
  {
    Free();
  }

  template <typename ETYPE> void ARRAY<ETYPE>::Free()
  {
    delete [] element;
    element = NULL;
  }

}

#endif

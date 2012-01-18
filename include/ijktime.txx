/// \file ijktime.txx
/// Time related templates
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


#ifndef _IJKTIME_
#define _IJKTIME_

#include <ctime>

namespace IJK {

  // **************************************************
  // CONVERT CPU TIME TO SECONDS
  // **************************************************

  template <typename S_TYPE>
  inline void clock2seconds(const clock_t t, S_TYPE & seconds)
  {
    seconds = S_TYPE(t)/CLOCKS_PER_SEC;
  }

}

#endif

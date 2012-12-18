/// \file ijkgrid_macros.txx
/// ijk macros for grid classes.
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

#ifndef _IJKGRID_MACROS_
#define _IJKGRID_MACROS_

// Use "for loop" to create local declaration
#define IJK_LOCAL(_init_statement,_counter0,_counter1)          \
  for (int _counter0=0, _counter1=0; _counter0<1; _counter0++)  \
    for (_init_statement; _counter1<1; _counter1++)

#define IJK_FOR_EACH_GRID_CUBE_LOCAL(_iv,_grid,_VTYPE,_cube_list,_i,_endv) \
  if (_grid.Dimension() > 0 && _grid.AxisSize(0) > 1)                \
    IJK_LOCAL(IJK::FACET0_CUBE_LIST<_VTYPE> _cube_list(_grid),       \
              _k0 ## __LINE__, _k1 ## __LINE__)                      \
      for (_VTYPE _i = 0; _i < _cube_list.NumCubes(); _i++)          \
        for (_VTYPE _iv = _cube_list.CubeIndex(_i),                  \
               _endv = _iv + _grid.AxisSize(0)-1;                    \
             _iv < _endv; _iv++)


#define IJK_FOR_EACH_GRID_CUBE(_iv,_grid,_VTYPE)                     \
  IJK_FOR_EACH_GRID_CUBE_LOCAL(_iv,_grid,_VTYPE,                     \
                               _cube_list ## __LINE__,               \
                               _i ## __LINE__,                       \
                               _endv ## __LINE__)

#define IJK_FOR_EACH_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE,_vlist,_i,_axis_inc,_axis_size,_endv) \
  if (_grid.Dimension() > 0)                                        \
    IJK_LOCAL(IJK::FACET_VERTEX_LIST<_VTYPE>               \
              _vlist(_grid, 0, true),                               \
              _k0 ## __LINE__, _k1 ## __LINE__)                     \
      for (_VTYPE _edge_dir = 0; _edge_dir < _grid.Dimension();     \
           _edge_dir++)                                             \
        if (_grid.AxisSize(_edge_dir) > 0)                          \
          IJK_LOCAL(_vlist.GetVertices(_grid, _edge_dir),           \
                   _k2 ## __LINE__, _k3 ## __LINE__)                \
            for (_VTYPE _i = 0,                                     \
                   _axis_inc = _grid.AxisIncrement(_edge_dir),      \
                   _axis_size = _grid.AxisSize(_edge_dir);          \
                 _i < _vlist.NumVertices(); _i++)                   \
              for (_VTYPE _iend0 = _vlist.VertexIndex(_i),          \
                     _endv = _iend0 + (_axis_size-1)*_axis_inc;     \
                   _iend0 < _endv;                                  \
                   _iend0 +=  _axis_inc)

#define IJK_FOR_EACH_GRID_EDGE(_iend0,_edge_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE, \
                                        _vlist ## __LINE__,         \
                                        _i ## __LINE__,             \
                                        _axis_inc ## __LINE__,      \
                                        _axis_size ## __LINE__,     \
                                        _endv ## __LINE__)


#define IJK_FOR_EACH_INTERIOR_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE,_vlist,_i,_axis_inc,_axis_size,_endv) \
  if (_grid.Dimension() > 0)                                        \
    IJK_LOCAL(IJK::FACET_INTERIOR_VERTEX_LIST<_VTYPE>               \
              _vlist(_grid, 0, true),                               \
              _k0 ## __LINE__, _k1 ## __LINE__)                     \
      for (_VTYPE _edge_dir = 0; _edge_dir < _grid.Dimension();     \
           _edge_dir++)                                             \
        if (_grid.AxisSize(_edge_dir) > 0)                          \
          IJK_LOCAL(_vlist.GetVertices(_grid, _edge_dir),           \
                   _k2 ## __LINE__, _k3 ## __LINE__)                \
            for (_VTYPE _i = 0,                                     \
                   _axis_inc = _grid.AxisIncrement(_edge_dir),      \
                   _axis_size = _grid.AxisSize(_edge_dir);          \
                 _i < _vlist.NumVertices(); _i++)                   \
              for (_VTYPE _iend0 = _vlist.VertexIndex(_i),          \
                     _endv = _iend0 + (_axis_size-1)*_axis_inc;     \
                   _iend0 < _endv;                                  \
                   _iend0 +=  _axis_inc)

#define IJK_FOR_EACH_INTERIOR_GRID_EDGE(_iend0,_edge_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_INTERIOR_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE, \
                                        _vlist ## __LINE__,         \
                                        _i ## __LINE__,             \
                                        _axis_inc ## __LINE__,      \
                                        _axis_size ## __LINE__,     \
                                        _endv ## __LINE__)

#define IJK_FOR_EACH_BOUNDARY_GRID_VERTEX_LOCAL(_iv,_grid,_VTYPE,_vertex_list,_i) \
  IJK_LOCAL(IJK::GRID_BOUNDARY_VERTEX_LIST<_VTYPE> _vertex_list(_grid), \
            _k0 ## __LINE__, _k1 ## __LINE__)                           \
    for (_VTYPE _i = 0, _iv = _vertex_list.VertexIndex(_i);         \
         _i < _vertex_list.NumVertices();                           \
         _i++, _iv = _vertex_list.VertexIndex(_i))

#define IJK_FOR_EACH_BOUNDARY_GRID_VERTEX(_iv,_grid,_VTYPE)          \
  IJK_FOR_EACH_BOUNDARY_GRID_VERTEX_LOCAL(_iv,_grid,_VTYPE,          \
                                          _vertex_list ## __LINE__,  \
                                          _i ## __LINE__)

#endif


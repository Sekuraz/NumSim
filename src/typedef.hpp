/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdint.h>
#include <cstdlib>
#include <iostream>

//------------------------------------------------------------------------------

#ifndef __TYPEDEF_HPP
#define __TYPEDEF_HPP

//------------------------------------------------------------------------------

#define DIM 2
#define REAL_TYPE DOUBLE // FLOAT, DOUBLE or LONG_DOUBLE
#define INDEX_TYPE std::size_t

//------------------------------------------------------------------------------

/// Typedef for reals
#if REAL_TYPE==DOUBLE
typedef double real_t;
#elif REAL_TYPE==LONG_DOUBLE
typedef long double real_t;
#else
typedef float real_t;
#endif

/// Typedef for integers
typedef INDEX_TYPE index_t;

//------------------------------------------------------------------------------

/// Template for array/vector types
template <typename _type, index_t _dim> struct array_t {
  /// Constructor initializes the array to 0.
  array_t() {
    for (index_t i = 0; i < _dim; ++i)
      x[i] = 0;
  }

  /// Constructor initializes the array to v1
  array_t(const _type &v1) {
    for (index_t i = 0; i < _dim; ++i)
      x[i] = v1;
  }

  /// Constructor initializes first value to v1, all others to v2
  array_t(const _type &v1, const _type &v2) {
    x[0] = v1;
    for (index_t i = 1; i < _dim; ++i)
      x[i] = v2;
  }

  /// Copy-constructor from field of values
  array_t(const _type(&cp)[_dim]) {
    for (index_t i = 0; i < _dim; ++i)
      x[i] = cp[i];
  }
  /// Copy-constructor from other array_t
  array_t(const array_t<_type, _dim> &cp) {
    for (index_t i = 0; i < _dim; ++i)
      x[i] = cp.x[i];
  }

  /// Writing access operator
  _type &operator[](index_t i) { return x[i]; }
  /// Reading access operator
  const _type &operator[](index_t i) const { return x[i]; }

  _type x[_dim]; ///< values of the array
  static const index_t dim = _dim; ///< size of the array
};

/// output of the array_t in readable form
template <typename _type, index_t _dim> std::ostream& operator<<(std::ostream& os, const array_t<_type, _dim> value) {
  os << "[" << value.x[0];
  for (index_t i = 1; i < value.dim; i++) {
      os << ", " << value.x[i];
  }
  os << "]";
  return os;
}
/// Typedef for d-dimensional array of reals
typedef array_t<real_t, DIM> multi_real_t;
/// Typedef for d-dimensional array of integer
typedef array_t<index_t, DIM> multi_index_t;

//------------------------------------------------------------------------------

// Forward declaration of classes used
class Communicator;
class Compute;
class Geometry;
class Grid;
class Iterator;
class Parameter;
class Solver;

#endif // __TYPEDEF_HPP

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

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //      u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //      u=0, v=0
  Geometry();

  /// Loads a geometry from a file
  void Load(const char file[]);

  /// Returns the number of cells in each dimension
  inline const multi_index_t &Size() const {
    return this->_size;
  }

  /// Returns the number of u-variables in each dimension
  inline const multi_index_t &SizeU() const {
    return this->_sizeU;
  }
  /// Returns the number of v-variables  in each dimension
  inline const multi_index_t &SizeV() const {
    return this->_sizeV;
  }
  /// Returns the number of p-variables  in each dimension
  inline const multi_index_t &SizeP() const {
    return this->_sizeP;
  }

  /// Returns the length of the domain
  inline const multi_real_t &Length() const {
    return this->_length;
  }
  /// Returns the velocity at the boundary
  inline const multi_real_t &Velocity() const {
    return this->_velocity;
  }
  /// Returns the meshwidth
  inline const multi_real_t &Mesh() const {
    return this->_h;
  }

  /// Updates the velocity field u
  void Update_U(Grid &u) const;
  /// Updates the velocity field v
  void Update_V(Grid &v) const;
  /// Updates the pressure field p
  void Update_P(Grid &p) const;

private:
  multi_index_t _size;
  multi_index_t _sizeU;
  multi_index_t _sizeV;
  multi_index_t _sizeP;
  multi_real_t _length;
  multi_real_t _h;

  multi_real_t _velocity;
  real_t _pressure;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP

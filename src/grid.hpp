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

//------------------------------------------------------------------------------
#ifndef __GRID_HPP
#define __GRID_HPP
//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

/// Class describing a grid of a variable for the computation
class Grid {
public:
  /// Constructs a grid based on a geometry
  Grid(const Geometry &geom);

  /// Constructs a grid based on a geometry with an offset
  /// \param geom   Geometry information
  /// \param offset distance of staggered grid point to cell's anchor point;
  ///               (anchor point = lower left corner)
  Grid(const Geometry &geom, const multi_real_t &offset);

  /// Deletes the grid
  ~Grid();

  /// Initializes the grid with a value
  void Initialize(const real_t &value);

  /// Write access to the grid cell at position [it]
  inline real_t& Cell(const Iterator &it) { return this->_data[it]; }
  /// Read access to the grid cell at position [it]
  inline const real_t& Cell(const Iterator &it) const { return this->_data[it]; }
  /// Interpolate the value at a arbitrary (real) position
  real_t Interpolate(const multi_real_t &pos) const;

  /// Computes the left-sided difference quotient in x-dim at [it]
  real_t dx_l(const Iterator &it) const;
  /// Computes the right-sided difference quotient in x-dim at [it]
  real_t dx_r(const Iterator &it) const;
  /// Computes the left-sided difference quotient in y-dim at [it]
  real_t dy_l(const Iterator &it) const;
  /// Computes the right-sided difference quotient in y-dim at [it]
  real_t dy_r(const Iterator &it) const;
  /// Computes the central difference quotient of 2nd order in x-dim at [it]
  real_t dxx(const Iterator &it) const;
  /// Computes the central difference quotient of 2nd order in y-dim at [it]
  real_t dyy(const Iterator &it) const;

  /// Computes u*du/dx with the donor cell method
  real_t DC_udu_x(const Iterator &it, const real_t &alpha) const;
  /// Computes v*du/dy with the donor cell method
  real_t DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const;
  /// Computes u*dv/dx with the donor cell method
  real_t DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const;
  /// Computes v*dv/dy with the donor cell method
  real_t DC_vdv_y(const Iterator &it, const real_t &alpha) const;

  /// Returns the maximal value of the grid
  real_t Max() const;
  /// Returns the minimal value of the grid
  real_t Min() const;
  /// Returns the absolute maximal value of the grid
  real_t AbsMax() const;
  /// Sets min, max to the minimal and maximal value of the grid
  void MinMax(real_t &min, real_t &max) const;

  /// Returns the Geometry of the grid
  inline const Geometry &Geom() const { return this->_geom; }
  /// Returns the size of the Grid in each dimension
  inline const multi_index_t &Size() const { return this->_geom.SizeP(); }
  /// Returns the total data size of the Grid
  inline const index_t &dataSize() const { return this->_sizeData; }

  /// Returns a pointer to the raw data
  inline real_t *Data() { return _data; };
  /// Returns a pointer to the raw data
  inline const real_t *cData() const { return _data; };

  /// Prints the data of the grid
  void print() const;

  /// Computes the inner L2-product of this and the other grid (over all interior cells)
  real_t InnerProduct(const Grid& other) const;

private:
  real_t *_data;
  const multi_real_t _offset;
  const Geometry &_geom;
  const multi_real_t &_hInv;
  const multi_real_t _hInv2;
  index_t _sizeData;
};
//------------------------------------------------------------------------------
#endif // __GRID_HPP

/*
 * Copyright (C) 2016   Stephan Lunowa, Markus Baur, Jonas Harsch
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
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "typedef.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"

// Constructs a grid based on a geometry
Grid::Grid(const Geometry &geom) : Grid(geom, multi_real_t(0)) {}

// Constructs a grid based on a geometry with an offset
// \param geom   Geometry information
// \param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry &geom, const multi_real_t &offset)
  : _data(nullptr), _offset(offset), _geom(geom), _hInv(geom.invMesh()),
    _hInv2(_hInv[0]*_hInv[0], _hInv[1]*_hInv[1]
#if (DIM == 3)
    , _hInv[2]*_hInv[2]
#endif
    ) {
  const multi_index_t &size = this->Size();
  this->_sizeData = 1;
  for(index_t i = 0; i < DIM; i++) {
    this->_sizeData *= size[i];
  }
  this->_data = new real_t[_sizeData];
}

// Deletes the grid
Grid::~Grid() {
  delete[] _data;
}

// Initializes the grid with a value
void Grid::Initialize(const real_t &value) {
  for(index_t i = 0; i < this->_sizeData; i++) {
    this->_data[i] = value;
  }
}

// Write access to the grid cell at position [it]
real_t &Grid::Cell(const Iterator &it) {
  return this->_data[it];
}
// Read access to the grid cell at position [it]
const real_t &Grid::Cell(const Iterator &it) const {
  return this->_data[it];
}

// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t &pos) const {
  const multi_index_t &size = this->Size();

  // compute index from position
  index_t multSize = 1;
  index_t i = 0;
  multi_real_t delta;
  for(index_t dim = 0; dim < DIM; dim++) {
    delta[dim] = (pos[dim] + this->_offset[dim]) * this->_hInv[dim];
    index_t iDim = (index_t)( delta[dim] );
    delta[dim] -= iDim;
    i += multSize * iDim;
    multSize *= size[dim];
  }

  Iterator it(*this, i);
  real_t value = 0;
  if(it.Valid()) {
    // constant interpolation
    //value = this->_data[it];

    // bilinear interpolation
    // TODO rewrite for n dim
    value = this->_data[it] * (1.0 - delta[0])*(1.0 - delta[1])
          + this->_data[it.Right()] * delta[0]*(1.0 - delta[1])
          + this->_data[it.Top()] * (1.0 - delta[0])*delta[1]
          + this->_data[it.Top().Right()] * delta[0]*delta[1];
  } else { // Error in index computation
    std::cerr << "Error: Grid: Pos out of area. Pos = " << pos[0] << ", " << pos[1]
              << ", it = " << it << std::endl;
  }
  return value;
}

// Computes the left-sided difference quotient in x-dim at [it]
real_t Grid::dx_l(const Iterator &it) const {
  return (this->_data[it]  - this->_data[it.Left()]) * this->_hInv[0];
}
// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const {
  return (this->_data[it.Right()] - this->_data[it]) * this->_hInv[0];
}
// Computes the left-sided difference quotient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const {
  return (this->_data[it] - this->_data[it.Down()]) * this->_hInv[1];
}
// Computes the right-sided difference quotient in y-dim at [it]
real_t Grid::dy_r(const Iterator &it) const {
  return (this->_data[it.Top()] - this->_data[it]) * this->_hInv[1];
}
// Computes the central difference quotient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
  return (this->_data[it.Left()] - 2*this->_data[it] + this->_data[it.Right()]) * this->_hInv2[0];
}
// Computes the central difference quotient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
  return (this->_data[it.Down()] - 2*this->_data[it] + this->_data[it.Top()]) * this->_hInv2[1];
}

// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const {
  real_t uMeanRight = (this->_data[it] + this->_data[it.Right()])/2;
  real_t uMeanLeft = (this->_data[it.Left()] + this->_data[it])/2;
  return ( (uMeanRight*uMeanRight - uMeanLeft*uMeanLeft) * this->_hInv[0]
         + alpha * (-std::fabs(uMeanRight)*this->dx_r(it)/2
                    +std::fabs(uMeanLeft) *this->dx_l(it)/2) );
}
// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const {
  real_t vMeanRight = (v->_data[it] + v->_data[it.Right()])/2;
  real_t vMeanDown = (v->_data[it.Down()] + v->_data[it.Down().Right()])/2;

  return ( ( vMeanRight *(this->_data[it] + this->_data[it.Top()])/2
            -vMeanDown*(this->_data[it.Down()] + this->_data[it])/2 ) * this->_hInv[1]
          + alpha * (-std::fabs(vMeanRight) *this->dy_r(it)/2
                     +std::fabs(vMeanDown)*this->dy_l(it)/2) );
}
// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const {
  real_t uMeanTop = (u->_data[it] + u->_data[it.Top()])/2;
  real_t uMeanLeft = (u->_data[it.Left()] + u->_data[it.Left().Top()])/2;

  return ( ( uMeanTop*(this->_data[it] + this->_data[it.Right()])/2
            -uMeanLeft *(this->_data[it.Left()] + this->_data[it])/2 ) * this->_hInv[0]
          + alpha * (-std::fabs(uMeanTop)*this->dx_r(it)/2
                     +std::fabs(uMeanLeft) *this->dx_l(it)/2) );
}
// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const {
  real_t vMeanTop = (this->_data[it] + this->_data[it.Top()])/2;
  real_t vMeanDown = (this->_data[it.Down()] + this->_data[it])/2;

  return ( (vMeanTop*vMeanTop - vMeanDown*vMeanDown) * this->_hInv[0]
         + alpha * (-std::fabs(vMeanTop) *this->dy_r(it)/2
                    +std::fabs(vMeanDown)*this->dy_l(it)/2) );
}

// Returns the maximal value of the grid
real_t Grid::Max() const {
  // TODO rewrite efficient
  real_t max = this->_data[0];
  for(index_t i = 1; i < this->_sizeData; i++)
    max = std::max(max, this->_data[i]);
  return max;
}
// Returns the minimal value of the grid
real_t Grid::Min() const {
  // TODO rewrite efficient
  real_t min = this->_data[0];
  for(index_t i = 1; i < this->_sizeData; i++)
    min = std::min(min, this->_data[i]);
  return min;
}
// Returns the absolute maximal value
real_t Grid::AbsMax() const {
  // TODO rewrite efficient
  real_t max = std::fabs(this->_data[0]);
  for(index_t i = 1; i < this->_sizeData; i++)
    max = std::max(max, std::fabs(this->_data[i]));
  return max;
}
// Sets min, max to the minimal and maximal value of the grid
void Grid::MinMax(real_t &min, real_t &max) const {
  min = max = this->_data[0];
  for(index_t i = 1; i < this->_sizeData; i++) {
    if(min > this->_data[i]) {
      min = this->_data[i];
    } else if(max < this->_data[i]) {
      max = this->_data[i];
    }
  }
}

void Grid::print() const {
  for(Iterator it(*this); it.Valid(); it.Next()) {
    if(it.Pos()[0] == 0) {
      std::cout << std::endl;
    }
    printf("%+.2f, ", this->Cell(it));
  }
  std::cout << std::endl;
}


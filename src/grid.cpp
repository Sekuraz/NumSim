#include <iostream>
#include <cmath>
#include <algorithm>
#include "typedef.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"

// Constructs a grid of type t based on a geometry
Grid::Grid(const Geometry &geom, const Grid::type &t) : _data(nullptr), _offset(), _geom(geom), _type(t) {
  const multi_index_t &size = this->Size();
  this->_sizeData = 1;
  for(index_t i = 0; i < DIM; i++) {
    this->_sizeData *= size[i];
  }
  this->_data = new real_t[_sizeData];
}

// Constructs a grid based on a geometry with an offset
// \param geom   Geometry information
// \param t      type of the grid (u, v, p)
// \param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry &geom, const Grid::type &t, const multi_real_t &offset)
    : _data(nullptr), _offset(offset), _geom(geom), _type(t) {
  const multi_index_t &size = this->Size();
  this->_sizeData = 1;
  for(index_t i = 0; i < DIM; i++) {
    this->_sizeData *= size[i];
  }
  this->_data = new real_t[_sizeData];
  // TODO: testing
  std::cout << "offset " << offset[0] << " " << offset[1] << std::endl;
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
  const multi_real_t &length = this->_geom.Length();
  const multi_index_t &size = this->Size();
  const multi_real_t &h = this->_geom.Mesh();

  // compute index from position
  index_t multSize = 1;
  index_t i = 0;
  multi_real_t delta;
  for(index_t dim = 0; dim < DIM; dim++) {
    delta[dim] = (pos[dim] + this->_offset[dim])/ h[dim];
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
  return (this->_data[it]  - this->_data[it.Left()]) / this->_geom.Mesh()[0];
}
// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const {
  return (this->_data[it.Right()] - this->_data[it]) / this->_geom.Mesh()[0];
}
// Computes the left-sided difference quotient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const {
  return (this->_data[it] - this->_data[it.Down()]) / this->_geom.Mesh()[1];
}
// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dy_r(const Iterator &it) const {
  return (this->_data[it.Top()] - this->_data[it]) / this->_geom.Mesh()[1];
}
// Computes the central difference quotient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
  const real_t &h = this->_geom.Mesh()[0];
  return (this->_data[it.Left()] - 2*this->_data[it] + this->_data[it.Right()]) /(h*h);
}
// Computes the central difference quotient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
  const real_t &h = this->_geom.Mesh()[1];
  return (this->_data[it.Down()] - 2*this->_data[it] + this->_data[it.Top()]) /(h*h);
}

// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const {
  // TODO implement
  return 0;
}
// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const {
  // TODO implement
  return 0;
}
// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const {
  // TODO implement
  return 0;
}
// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const {
  // TODO implement
  return 0;
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


#include <cstddef>
#include "typedef.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"

// Constructs a grid based on a geometry
Grid::Grid(const Geometry *geom) : _data(nullptr), _offset(), _geom(geom) {
  const multi_index_t &gSize = _geom->Size();
  index_t size = 1;
  for(index_t i = 0; i < DIM; i++) {
    size *= gSize[i];
  }
  this->_data = new real_t[size];
}

// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry *geom, const multi_real_t &offset)
    : _data(nullptr), _offset(offset), _geom(geom) {
  const multi_index_t &gSize = _geom->Size();
  index_t size = 1;
  for(index_t i = 0; i < DIM; i++) {
    size *= gSize[i];
  }
  this->_data = new real_t[size];
}

// Deletes the grid
Grid::~Grid() {
  delete[] _data;
}

// Initializes the grid with a value
void Grid::Initialize(const real_t &value) {
  // TODO
}

// Write access to the grid cell at position [it]
real_t &Grid::Cell(const Iterator &it) {
  return _data[it]; // TODO: + offset ?
}
// Read access to the grid cell at position [it]
const real_t &Grid::Cell(const Iterator &it) const {
  return _data[it]; // TODO: + offset ?
}

// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t &pos) const {
  // TODO
  return 0;
}

// Computes the left-sided difference quotient in x-dim at [it]
real_t Grid::dx_l(const Iterator &it) const {
  return (_data[it]  - _data[it.Left()]) / _geom->Mesh()[0];
}
// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const {
  return (_data[it.Right()] - _data[it]) / _geom->Mesh()[0];
}
// Computes the left-sided difference quotient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const {
  return (_data[it] - _data[it.Down()]) / _geom->Mesh()[1];
}
// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dy_r(const Iterator &it) const {
  return (_data[it.Top()] - _data[it]) / _geom->Mesh()[1];
}
// Computes the central difference quotient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
  const real_t &h =  _geom->Mesh()[0];
  return (_data[it.Left()] - 2*_data[it] + _data[it.Right()]) /(h*h);
}
// Computes the central difference quotient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
  const real_t &h =  _geom->Mesh()[1];
  return (_data[it.Down()] - 2*_data[it] + _data[it.Top()]) /(h*h);
}

// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const {
  // TODO
  return 0;
}
// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const {
  // TODO
  return 0;
}
// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const {
  // TODO
  return 0;
}
// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const {
  // TODO
  return 0;
}

// Returns the maximal value of the grid
real_t Grid::Max() const {
  // TODO
  return 0;
}
// Returns the minimal value of the grid
real_t Grid::Min() const {
  // TODO
  return 0;
}
// Returns the absolute maximal value
real_t Grid::AbsMax() const {
  // TODO
  return 0;
}


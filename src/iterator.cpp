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
#include <omp.h>
#include "iterator.hpp"
#include "geometry.hpp"

Iterator::Iterator(const Geometry &geom) : _geom(geom), _valid(true) {
  const int n = omp_get_num_threads();
  const int i = omp_get_thread_num();
  const index_t &size = geom.DataSize();
  index_t part = size / n;
  this->_begin = i * part;
  this->_value = this->_begin;
  this->_end = (i+1 == n)? size : (i+1) * part;
//  printf("It: n = %d, i = %d, starting at %lu end at %lu, total %lu\n", n, i, _value, _end, size);
}
Iterator::Iterator(const Geometry &geom, const index_t &value) : Iterator(geom) {
  this->_value = std::max(this->_begin, value);
  this->_valid = value < this->_end;
}
Iterator::Iterator(const Geometry &geom, const multi_index_t &pos) : Iterator(geom) {
  this->_value = pos[DIM-1];
  const multi_index_t &size = this->_geom.SizeP();
  for(index_t dim = DIM-1; dim-- > 0; ) {
    this->_value *= size[dim];
    this->_value += pos[dim];
  }
  this->_valid = this->_value < this->_end;
}

// Returns the position coordinates
multi_index_t Iterator::Pos() const {
    // TODO: rewrite for n dimensions
    multi_index_t pos((this->_value % this->_geom.SizeP()[0]),
                      (this->_value / this->_geom.SizeP()[0]));
    return pos;
}

// Sets the iterator to the first element
void Iterator::First() {
    this->_value = this->_begin;
    this->_valid = true;
}

// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
    if (++this->_value >= this->_end) {
        this->_valid = false;
    }
}

// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const {
    if (this->_value % this->_geom.SizeP()[0] == 0) {
        return *this;
    }
    return Iterator(this->_geom, this->_value - 1);
}

// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
    if ((this->_value + 1) % this->_geom.SizeP()[0] == 0) {
        return *this;
    }
    return Iterator(this->_geom, this->_value + 1);
}

// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
    if (this->_value >= this->_geom.SizeP()[0]*(this->_geom.SizeP()[1]-1)) {
        return *this;
    }
    return Iterator(this->_geom, this->_value + this->_geom.SizeP()[0]);
}

// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
    if (this->_value < this->_geom.SizeP()[0]) {
        return *this;
    }
    return Iterator(this->_geom, this->_value - this->_geom.SizeP()[0]);
}

// Returns the value of an Iterator that is located left from this one.
// If we are at the left boundary, the cell sees itself.
index_t Iterator::Left(int __attribute__((unused))) const {
    if (this->_value % this->_geom.SizeP()[0] == 0) {
        return this->_value;
    }
    return this->_value - 1;
}

// Returns the value of an Iterator that is located right from this one.
// If we are at the right boundary, the cell sees itself.
index_t Iterator::Right(int __attribute__((unused))) const {
    if ((this->_value + 1) % this->_geom.SizeP()[0] == 0) {
        return this->_value;
    }
    return this->_value + 1;
}

// Returns the value of an Iterator that is located above this one.
// If we are at the upper domain boundary, the cell sees itself.
index_t Iterator::Top(int __attribute__((unused))) const {
    if (this->_value >= this->_geom.SizeP()[0]*(this->_geom.SizeP()[1]-1)) {
        return this->_value;
    }
    return this->_value + this->_geom.SizeP()[0];
}

// Returns the value of an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
index_t Iterator::Down(int __attribute__((unused))) const {
    if (this->_value < this->_geom.SizeP()[0]) {
        return this->_value;
    }
    return this->_value - this->_geom.SizeP()[0];
}

//------------------------------------------------------------------------------
// Iterator for interior cells

// Sets the iterator to the first element
void InteriorIterator::First() {
  Iterator::First();
  if( !this->_geom.isFluid(*this) ) {
    this->Next();
  }
}

// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
  Iterator::Next();
  while(this->_valid && !this->_geom.isFluid(*this) ) {
    Iterator::Next();
  }
}

//------------------------------------------------------------------------------
// Iterator for domain boundary cells.

// Sets the iterator to the first element
void BoundaryIterator::First() {
  switch(this->_boundary) {
    case BoundaryIterator::boundary::top:
      this->_value = this->_geom.SizeP()[0]*(this->_geom.SizeP()[1]-1);
      break;
    case BoundaryIterator::boundary::right:
      this->_value = this->_geom.SizeP()[0]-1;
      break;
    default: // left, down
      this->_value = 0;
      break;
  }
  this->_valid = true;
}
// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
  switch(this->_boundary) {
    case BoundaryIterator::boundary::top:
      Iterator::Next();
      break;
    case BoundaryIterator::boundary::left:
      this->_value += this->_geom.SizeP()[0];
      if(this->_value >= this->_geom.DataSize()) {
        this->_valid = false;
      }
      break;
    case BoundaryIterator::boundary::down:
      Iterator::Next();
      if(this->_value >= this->_geom.SizeP()[0]) {
        this->_valid = false;
      }
      break;
    case BoundaryIterator::boundary::right:
      this->_value += this->_geom.SizeP()[0];
      if(this->_value >= this->_geom.DataSize()) {
        this->_valid = false;
      }
      break;
    default:
      std::cerr << "Error: BoundaryIterator: boundary " << this->_boundary
                << " does not exist!" << std::endl;
      this->_valid = false;
      break;
  }
}


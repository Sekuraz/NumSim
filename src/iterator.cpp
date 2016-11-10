/*
 * Copyright (C) 2016   Stephan Lunowa, Markus Baur
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
#include "iterator.hpp"
#include "grid.hpp"

Iterator::Iterator(const Grid &grid) : _grid(grid), _value(0), _valid(true) {}
Iterator::Iterator(const Grid &grid, const index_t &value)
    : _grid(grid), _value(value), _valid(value < grid.dataSize()) {}
Iterator::Iterator(const Grid &grid, const multi_index_t &pos) : _grid(grid) {
  this->_value = pos[DIM-1];
  const multi_index_t &size = this->_grid.Size();
  for(index_t dim = DIM-1; dim-- > 0; ) {
    this->_value *= size[dim];
    this->_value += pos[dim];
  }
  this->_valid = this->_value < this->_grid.dataSize();
}

// Returns the current position value
const index_t& Iterator::Value() const {
    return this->_value;
}

// Cast operator to convert Iterators to integers
Iterator::operator const index_t& () const {
    return this->_value;
}

// Returns the position coordinates
multi_index_t Iterator::Pos() const {
    // TODO: rewrite for n dimensions
    multi_index_t pos((this->_value % this->_grid.Size()[0]),
                      (this->_value / this->_grid.Size()[0]));
    return pos;
}

// Sets the iterator to the first element
void Iterator::First() {
    this->_value = 0;
    this->_valid = true;
}

// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
    if (++this->_value >= this->_grid.dataSize()) {
        this->_valid = false;
    }
}

// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
    return this->_valid;
}

// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const {
    if (this->_value % this->_grid.Size()[0] == 0) {
        return *this;
    }
    return Iterator(this->_grid, this->_value - 1);
}

// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
    if ((this->_value + 1) % this->_grid.Size()[0] == 0) {
        return *this;
    }
    return Iterator(this->_grid, this->_value + 1);
}

// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
    if (this->_value >= this->_grid.Size()[0]*(this->_grid.Size()[1]-1)) {
        return *this;
    }
    return Iterator(this->_grid, this->_value + this->_grid.Size()[0]);
}

// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
    if (this->_value < this->_grid.Size()[0]) {
        return *this;
    }
    return Iterator(this->_grid, this->_value - this->_grid.Size()[0]);
}

//------------------------------------------------------------------------------
// Iterator for interior cells

// Sets the iterator to the first element
void InteriorIterator::First() {
  this->_value = this->_grid.Size()[0] + 1;
  this->_valid = true;
}
// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
  Iterator::Next();
  if (this->_value >= this->_grid.Size()[0]*(this->_grid.Size()[1]-1)) {
    this->_valid = false;
  } else if (this->_value % this->_grid.Size()[0] == 0) {
    this->Next();
  } else if ((this->_value + 1) % this->_grid.Size()[0] == 0) {
    this->Next();
  }
}

//------------------------------------------------------------------------------
// Iterator for domain boundary cells.

// Sets the iterator to the first element
void BoundaryIterator::First() {
  switch(this->_boundary) {
    case 1: // top
      this->_value = this->_grid.Size()[0]*(this->_grid.Size()[1]-1);
      break;
    case 4: // right
      this->_value = this->_grid.Size()[0]-1;
      break;
    default: // left, down, all
      this->_value = 0;
      break;
  }
  this->_valid = true;
}
// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
  switch(this->_boundary) {
    case 0: // all
      // TODO implement all boundaries
      std::cerr << "Error: BoundaryIterator: all boundaries not implemented.";
      this->_valid = false;
      break;
    case 1: // top
      Iterator::Next();
      break;
    case 2: // left
      this->_value += this->_grid.Size()[0];
      if(this->_value >= this->_grid.dataSize()) {
        this->_valid = false;
      }
      break;
    case 3: // down
      Iterator::Next();
      if(this->_value >= this->_grid.Size()[0]) {
        this->_valid = false;
      }
      break;
    case 4: // right
      this->_value += this->_grid.Size()[0];
      if(this->_value >= this->_grid.dataSize()) {
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


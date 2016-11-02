#include <iostream>
#include "iterator.hpp"
#include "geometry.hpp"

Iterator::Iterator(const Geometry &geom) : _geom(geom), _value(0), _valid(true) {}
Iterator::Iterator(const Geometry &geom, const index_t &value) : _geom(geom), _value(value), _valid(true) {
  if (this->_value >= this->_geom.Size()[0] * this->_geom.Size()[1]) {
    this->_valid = false;
  }
}


//     Returns the current position value
const index_t& Iterator::Value() const {
    return this->_value;
}

// Cast operator to convert Iterators to integers
Iterator::operator const index_t& (void) const {
    return this->_value;
}

// Returns the position coordinates
multi_index_t Iterator::Pos() const {
    return multi_index_t(this->_value % this->_geom.Size()[0], this->_value / this->_geom.Size()[0]);
}

// Sets the iterator to the first element
void Iterator::First() {
    this->_value = 0;
}

// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
    if (++this->_value >= this->_geom.Size()[0] * this->_geom.Size()[1]) {
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
    if (this->_value % this->_geom.Size()[0] == 0) {
        return *this;
    }
    return Iterator(_geom, this->_value - 1);
}

// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
    if ((this->_value + 1) % this->_geom.Size()[0] == 0) {
        return *this;
    }
    return Iterator(_geom, this->_value + 1);
}

// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
    if (this->_value >= this->_geom.Size()[0]*(this->_geom.Size()[1]-1)) {
        return *this;
    }
    return Iterator(_geom, this->_value + this->_geom.Size()[0]);
}

// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
    if (this->_value < this->_geom.Size()[0]) {
        return *this;
    }
    return Iterator(_geom, this->_value - this->_geom.Size()[0]);
}

//------------------------------------------------------------------------------
// Iterator for interior cells

// Construct a new InteriorIterator
InteriorIterator::InteriorIterator(const Geometry &geom) : Iterator(geom, geom.Size()[0] + 1) {}

// Sets the iterator to the first element
void InteriorIterator::First() {
  this->_value = this->_geom.Size()[0] + 1;
}
// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
  Iterator::Next();
  if (this->_value >= this->_geom.Size()[0]*(this->_geom.Size()[1]-1)) {
    this->_valid = false;
  } else if (this->_value % this->_geom.Size()[0] == 0) {
    this->Next();
  } else if ((this->_value + 1) % this->_geom.Size()[0] == 0) {
    this->Next();
  }
}

//------------------------------------------------------------------------------
// Iterator for domain boundary cells.

// Sets the iterator to the first element
void BoundaryIterator::First() {
  switch(this->_boundary) {
    case 1:
      this->_value = this->_geom.Size()[0]*(this->_geom.Size()[1]-1);
      break;
    case 4:
      this->_value = this->_geom.Size()[0]-1;
      break;
    default:
      this->_value = 0;
      break;
  }
}
// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
  switch(this->_boundary) {
    case 1:
      Iterator::Next();
      break;
    case 2:
      this->_value += this->_geom.Size()[0];
      if(this->_value >= this->_geom.Size()[0]*this->_geom.Size()[1]) {
        this->_valid = false;
      }
      break;
    case 3:
      Iterator::Next();
      if(this->_value >= this->_geom.Size()[0]) {
        this->_valid = false;
      }
      break;
    case 4:
      this->_value += this->_geom.Size()[0];
      if(this->_value >= this->_geom.Size()[0]*this->_geom.Size()[1]) {
        this->_valid = false;
      }
      break;
    default:
      // TODO
      std::cerr << "Error: BoundaryIterator: all boundaries not implemented.";
      this->_value = 0;
      break;
  }
}


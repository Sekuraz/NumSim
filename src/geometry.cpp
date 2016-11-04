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
#include <fstream>
#include <string>
#include "typedef.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

using namespace std;

Geometry::Geometry() : _size(128,128), _sizeU(_size), _sizeV(_size),
    _sizeP(_size), _length(1,1) {
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_sizeU[dim] += (dim == 0)? 1 : 2;
    this->_sizeV[dim] += (dim == 1)? 1 : 2;
    this->_sizeP[dim] += 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
  }
}

/// Loads a geometry from a file
void Geometry::Load(const char file[]) {
  string param;
  real_t value;
  ifstream in(file);
  while(in.good()) {
    in >> ws >> param >> ws;
    in.ignore(); // removes '=' between parameter and value
    in >> ws >> value >> ws;
    if(!param.compare("Size") || !param.compare("size")) {
      this->_size[0] = (index_t)value;
      cout << "Geometry: Load size = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_size[dim] = (index_t)value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Length") || !param.compare("length")) {
      this->_length[0] = value;
      cout << "Geometry: Load length = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_length[dim] = value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Velocity") || !param.compare("velocity")) {
      this->_velocity[0] = value;
      cout << "Geometry: Load velocity = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_velocity[dim] = value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Pressure") || !param.compare("pressure")) {
      cout << "Geometry: Load pressure = " << value << endl;
      this->_pressure = value;
    } else {
      cerr << "Geometry: unknown identifier " << param;
    }
  }
  in.close();

  for(index_t dim = 0; dim < DIM; dim++) {
    this->_sizeU[dim] = this->_size[dim] + ((dim == 0)? 1 : 2);
    this->_sizeV[dim] = this->_size[dim] + ((dim == 1)? 1 : 2);
    this->_sizeP[dim] = this->_size[dim] + 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
  }
}

// Updates the velocity field u
void Geometry::Update_U(Grid &u) const {
    //TODO noop
}

// Updates the velocity field v
void Geometry::Update_V(Grid &v) const {
    //TODO noop
}

// Updates the pressure field p
void Geometry::Update_P(Grid &p) const {
  // homogenous Neumann condition
  BoundaryIterator it(p);
  it.SetBoundary(1);
  for(it.First(); it.Valid(); it.Next()) {
    p.Cell(it) = p.Cell(it.Down());
  }
  it.SetBoundary(2);
  for(it.First(); it.Valid(); it.Next()) {
    p.Cell(it) = p.Cell(it.Right());
  }
  it.SetBoundary(3);
  for(it.First(); it.Valid(); it.Next()) {
    p.Cell(it) = p.Cell(it.Top());
  }
  it.SetBoundary(4);
  for(it.First(); it.Valid(); it.Next()) {
    p.Cell(it) = p.Cell(it.Left());
  }
}

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
#include <fstream>
#include <string>
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

using namespace std;

Geometry::Geometry(const Communicator &comm) : Geometry(comm, multi_index_t {128, 128}) {}

Geometry::Geometry(const Communicator &comm, const multi_index_t& size)
    : _comm(comm), _totalSize(size), _size(), _sizeU(), _sizeV(), _sizeP(), _length(1,1) {
  this->computeSizes();
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
      this->_totalSize[0] = (index_t)value;
      cout << "Geometry: Load size = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_totalSize[dim] = (index_t)value;
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

  this->computeSizes();
}

// Updates the velocity field u
void Geometry::Update_U(Grid &u) const {
  // driven cavity u given at top
  this->_comm.copyBoundary(u);

  BoundaryIterator it(u);
  if(this->_comm.isLeft()) {
    // homogenous Dirichlet condition left
    it.SetBoundary(2);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = 0;
    }
  }
  if(this->_comm.isBottom()) {
    // homogenous Dirichlet condition (mean) down
    it.SetBoundary(3);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = -u.Cell(it.Top());
    }
  }
  if(this->_comm.isRight()) {
    // homogenous Dirichlet condition right
    it.SetBoundary(4);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = 0;
    }
  }
  if(this->_comm.isTop()) {
    // inhomogenous Dirichlet condition (mean) top
    it.SetBoundary(1);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = 2*this->_velocity[0]-u.Cell(it.Down());
    }
  }
}

// Updates the velocity field v
void Geometry::Update_V(Grid &v) const {
  // driven cavity v given at top
  this->_comm.copyBoundary(v);

  BoundaryIterator it(v);
  if(this->_comm.isLeft()) {
    // homogenous Dirichlet condition (mean) left
    it.SetBoundary(2);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = -v.Cell(it.Right());
    }
  }
  if(this->_comm.isBottom()) {
    // homogenous Dirichlet condition down
    it.SetBoundary(3);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = 0;
    }
  }
  if(this->_comm.isRight()) {
    // homogenous Dirichlet condition (mean) right
    it.SetBoundary(4);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = -v.Cell(it.Left());
    }
  }
  if(this->_comm.isTop()) {
    // inhomogenous Dirichlet condition top
    it.SetBoundary(1);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = this->_velocity[1];
    }
  }
}

// Updates the pressure field p
void Geometry::Update_P(Grid &p) const {
  // homogenous Neumann condition at all sides
  this->_comm.copyBoundary(p);

  BoundaryIterator it(p);
  if(this->_comm.isLeft()) {
    it.SetBoundary(2);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Right());
    }
  }
  if(this->_comm.isBottom()) {
    it.SetBoundary(3);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Top());
    }
  }
  if(this->_comm.isRight()) {
    it.SetBoundary(4);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Left());
    }
  }
  if(this->_comm.isTop()) {
    it.SetBoundary(1);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Down());
    }
  }
}

// compute grid sizes
void Geometry::computeSizes() {
  const multi_index_t& numProc = this->_comm.ThreadDim();
  //const multi_index_t& idx = this->_comm.ThreadIdx();
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_size[dim] = this->_totalSize[dim] / numProc[dim];
    if(this->_size[dim] % 2 == 1) {
      this->_size[dim]++;
    }
    // TODO: fix last grid size

    this->_sizeU[dim] = this->_size[dim] + ((dim == 0)? 1 : 2);
    this->_sizeV[dim] = this->_size[dim] + ((dim == 1)? 1 : 2);
    this->_sizeP[dim] = this->_size[dim] + 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
  }
}

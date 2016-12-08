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

Geometry::Geometry(const Communicator &comm) : Geometry(comm, multi_index_t {128, 128}) {}

Geometry::Geometry(const Communicator &comm, const multi_index_t& size)
    : _comm(comm), _totalSize(size), _totalLength(1,1) {
  this->computeSizes();
}

/// Loads a geometry from a file
void Geometry::Load(const char file[]) {
  std::string param;
  real_t value;
  std::ifstream in(file);
  while(in.good()) {
    in >> param >> std::ws;
    if(!std::isdigit(in.peek())) {
      in.ignore(); // removes separator between parameter and value
    }
    in >> value;
    if(!param.compare("Size") || !param.compare("size")) {
      this->_totalSize[0] = (index_t)value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value;
        this->_totalSize[dim] = (index_t)value;
      }
      if(this->_comm.ThreadNum() == 0) {
        std::cout << "Geometry: Load size = " << this->_totalSize << std::endl;
      }
    } else if(!param.compare("Length") || !param.compare("length")) {
      this->_totalLength[0] = value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value;
        this->_totalLength[dim] = value;
      }
      if(this->_comm.ThreadNum() == 0) {
        std::cout << "Geometry: Load length = " << this->_totalLength << std::endl;
      }
    } else if(!param.compare("Velocity") || !param.compare("velocity")) {
      this->_velocity[0] = value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value;
        this->_velocity[dim] = value;
      }
      if(this->_comm.ThreadNum() == 0) {
        std::cout << "Geometry: Load velocity = " << this->_velocity << std::endl;
      }
    } else if(!param.compare("Pressure") || !param.compare("pressure")) {
      this->_pressure = value;
      if(this->_comm.ThreadNum() == 0) {
        std::cout << "Geometry: Load pressure = " << this->_pressure << std::endl;
      }
    } else {
      std::cerr << "Geometry: unknown identifier " << param << std::endl;
    }
    in >> std::ws;
  }
  in.close();

  this->computeSizes();
}

// Updates the velocity field u
void Geometry::Update_U(Grid &u) const {
  // driven cavity u given at top
  this->_comm.copyBoundary(u);

  BoundaryIterator it(u, BoundaryIterator::boundary::left);
  if(this->_comm.isLeft()) {
    // homogenous Dirichlet condition left
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = 0;
    }
  }
  if(this->_comm.isBottom()) {
    // homogenous Dirichlet condition (mean) down
    it.SetBoundary(BoundaryIterator::boundary::down);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = -u.Cell(it.Top());
    }
  }
  if(this->_comm.isRight()) {
    // homogenous Dirichlet condition right
    it.SetBoundary(BoundaryIterator::boundary::right);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it.Left()) = 0;
      u.Cell(it) = 0;
    }
  }
  if(this->_comm.isTop()) {
    // inhomogenous Dirichlet condition (mean) top
    it.SetBoundary(BoundaryIterator::boundary::top);
    for(it.First(); it.Valid(); it.Next()) {
      u.Cell(it) = 2*this->_velocity[0]-u.Cell(it.Down());
    }
  }
}

// Updates the velocity field v
void Geometry::Update_V(Grid &v) const {
  // driven cavity v given at top
  this->_comm.copyBoundary(v);

  BoundaryIterator it(v, BoundaryIterator::boundary::left);
  if(this->_comm.isLeft()) {
    // homogenous Dirichlet condition (mean) left
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = -v.Cell(it.Right());
    }
  }
  if(this->_comm.isBottom()) {
    // homogenous Dirichlet condition down
    it.SetBoundary(BoundaryIterator::boundary::down);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = 0;
    }
  }
  if(this->_comm.isRight()) {
    // homogenous Dirichlet condition (mean) right
    it.SetBoundary(BoundaryIterator::boundary::right);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = -v.Cell(it.Left());
    }
  }
  if(this->_comm.isTop()) {
    // inhomogenous Dirichlet condition top
    it.SetBoundary(BoundaryIterator::boundary::top);
    for(it.First(); it.Valid(); it.Next()) {
      v.Cell(it) = this->_velocity[1];
      v.Cell(it.Down()) = this->_velocity[1];
    }
  }
}

// Updates the pressure field p
void Geometry::Update_P(Grid &p) const {
  // homogenous Neumann condition at all sides
  this->_comm.copyBoundary(p);

  BoundaryIterator it(p, BoundaryIterator::boundary::left);
  if(this->_comm.isLeft()) {
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Right());
    }
  }
  if(this->_comm.isBottom()) {
    it.SetBoundary(BoundaryIterator::boundary::down);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Top());
    }
  }
  if(this->_comm.isRight()) {
    it.SetBoundary(BoundaryIterator::boundary::right);
    for(it.First(); it.Valid(); it.Next()) {
      p.Cell(it) = p.Cell(it.Left());
    }
  }
  if(this->_comm.isTop()) {
    it.SetBoundary(BoundaryIterator::boundary::top);
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
    this->_length[dim] = this->_totalLength[dim] / numProc[dim];

    // TODO: fix last grid size
    // until then fix total size
    this->_totalSize[dim] = this->_size[dim] * numProc[dim];

    this->_sizeP[dim] = this->_size[dim] + 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
    this->_hInv[dim] = 1.0 / this->_h[dim];
  }
}

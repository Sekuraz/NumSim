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
    : _comm(comm), _free(false), _flags(nullptr), _totalSize(size), _totalLength(1,1) {
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
    } else if(!param.compare("Geometry") || !param.compare("geometry")) {
      this->_free = true;
      in.ignore(10000, '\n');  // removes remaining line
      this->_flags = new char[this->_totalSize[0] * this->_totalSize[1]];
      for(index_t j = 0; j < this->_totalSize[1]; j++) {
        in.read(&this->_flags[this->_totalSize[0] * j], this->_totalSize[0]);
        in >> std::ws; // removes newline between lines
      }
      if(this->_comm.ThreadNum() == 0) {
        std::cout << "Geometry: Load free geometry." << std::endl;
        // TODO: output for testing
        for(index_t j = 0; j < this->_totalSize[1]; j++) {
          std::cout.write(&this->_flags[this->_totalSize[0] * j], this->_totalSize[0]);
          std::cout << std::endl;
        }
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
  this->_comm.copyBoundary(u);
  // driven cavity u given at top
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

void Geometry::Update(Grid &u, Grid &v, Grid &p) const {
  for(Iterator it(u); it.Valid(); it.Next()) {
    switch(this->flag(it.Pos())) {
      case '#': // '#' Wall/Obstacle/NoSlip boundary (u = v = 0, dp/dn = 0)
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it.Left()) = 0;
          if(this->isFluid(it.Top().Pos())) {
            v.Cell(it) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            v.Cell(it) = -v.Cell(it.Left());
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top().Pos())) {
          v.Cell(it) = 0;
          if(this->isFluid(it.Right().Pos())) {
            u.Cell(it) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            u.Cell(it) = -u.Cell(it.Top());
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it) = 0;
          if(this->isFluid(it.Down().Pos())) {
            v.Cell(it.Down()) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            v.Cell(it) = -v.Cell(it.Right());
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down().Pos())) {
          v.Cell(it.Down()) = 0;
          if(this->isFluid(it.Left().Pos())) {
            u.Cell(it.Left()) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            u.Cell(it) = -u.Cell(it.Down());
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'I': // General Inflow boundary (u = u_0, v = v_0, dp/dn = 0)
        // TODO: other than const velocities
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it.Left()) = this->_velocity[0];
          if(this->isFluid(it.Top().Pos())) {
            v.Cell(it) = this->_velocity[1];
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Left());
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top().Pos())) {
          v.Cell(it) = this->_velocity[1];
          if(this->isFluid(it.Right().Pos())) {
            u.Cell(it) = this->_velocity[0];
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            u.Cell(it) = 2*this->_velocity[0]-u.Cell(it.Top());
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it) = this->_velocity[0];
          if(this->isFluid(it.Down().Pos())) {
            v.Cell(it.Down()) = this->_velocity[1];
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Right());
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down().Pos())) {
          v.Cell(it.Down()) = this->_velocity[1];
          if(this->isFluid(it.Left().Pos())) {
            u.Cell(it.Left()) = this->_velocity[0];
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Down());
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'H': // Horizontal Inflow boundary (u = u_0, v = ?)
        // TODO: other than const velocities
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it) = this->_velocity[0];
          if(this->isFluid(it.Top().Pos())) {
            v.Cell(it) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            v.Cell(it) = -v.Cell(it.Left());
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top().Pos())) {
          v.Cell(it) = 0;
          if(this->isFluid(it.Right().Pos())) {
            u.Cell(it) = this->_velocity[0];
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            u.Cell(it) = 2*this->_velocity[0]-u.Cell(it.Top());
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it) = this->_velocity[0];
          if(this->isFluid(it.Down().Pos())) {
            v.Cell(it.Down()) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            v.Cell(it) = -v.Cell(it.Right());
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down().Pos())) {
          v.Cell(it.Down()) = 0;
          if(this->isFluid(it.Left().Pos())) {
            u.Cell(it.Left()) = this->_velocity[0];
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Down());
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'V': // Vertical Inflow boundary (u = ?, v = v_0)
        // TODO: other than const velocities
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it) = 0;
          if(this->isFluid(it.Top().Pos())) {
            v.Cell(it) = this->_velocity[1];
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            v.Cell(it) = 2*this->_velocity[1]-v.Cell(it.Left());
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top().Pos())) {
          v.Cell(it) = this->_velocity[1];
          if(this->isFluid(it.Right().Pos())) {
            u.Cell(it) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            u.Cell(it) = -u.Cell(it.Top());
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it) = 0;
          if(this->isFluid(it.Down().Pos())) {
            v.Cell(it.Down()) = this->_velocity[1];
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Right());
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down().Pos())) {
          v.Cell(it.Down()) = this->_velocity[1];
          if(this->isFluid(it.Left().Pos())) {
            u.Cell(it.Left()) = 0;
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            u.Cell(it) = -u.Cell(it.Down());
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'O': // Outflow boundary (du/dx = 0 resp. dv/dy = 0)
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it) = u.Cell(it.Left());
          if(this->isFluid(it.Top().Pos())) {
            v.Cell(it) = v.Cell(it.Top());
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            v.Cell(it) = v.Cell(it.Top());
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top().Pos())) {
          v.Cell(it) = v.Cell(it.Top());
          if(this->isFluid(it.Right().Pos())) {
            u.Cell(it.Left()) = u.Cell(it);
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            u.Cell(it.Left()) = u.Cell(it);
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it.Left()) = u.Cell(it);
          if(this->isFluid(it.Down().Pos())) {
            v.Cell(it) = v.Cell(it.Down());
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            v.Cell(it) = v.Cell(it.Down());
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down().Pos())) {
          v.Cell(it) = v.Cell(it.Down());
          if(this->isFluid(it.Left().Pos())) {
            u.Cell(it) = u.Cell(it.Left());
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            u.Cell(it) = u.Cell(it.Left());
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case '|': // Vertical Slip-boundary (u = 0, dv/dx = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Left().Pos())) {
          u.Cell(it) = 0;
          v.Cell(it) = v.Cell(it.Left());
          p.Cell(it) = this->_pressure;
        } else if(this->isFluid(it.Right().Pos())) {
          u.Cell(it.Left()) = 0;
          v.Cell(it.Left()) = v.Cell(it);
          p.Cell(it) = this->_pressure;
        }
        break;
      case '-': // Horizontal Slip-boundary (du/dy = 0, v = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Top().Pos())) {
          u.Cell(it) = u.Cell(it.Down());
          v.Cell(it) = 0;
          p.Cell(it) = this->_pressure;
        } else if(this->isFluid(it.Down().Pos())) {
          u.Cell(it) = 0;
          v.Cell(it.Down()) = v.Cell(it);
          p.Cell(it) = this->_pressure;
        }
        break;
    }
  }
}

// compute grid sizes
void Geometry::computeSizes() {
  const multi_index_t& numProc = this->_comm.ThreadDim();
  //const multi_index_t& idx = this->_comm.ThreadIdx();
  for(index_t dim = 0; dim < DIM; dim++) {
    if(this->_free) {
      this->_totalSize[dim] -= 2;
    }
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

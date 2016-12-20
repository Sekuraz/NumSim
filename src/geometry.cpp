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


Geometry::Geometry(const Communicator &comm, const multi_index_t& size,
                   const char file[], const bool printinfo)
    : _comm(comm), _free(false), _flags(nullptr), _totalSize(size), _totalLength(1,1),
    _velocity(1,0), _pressure(0.1) {
  if(file != nullptr) {
    this->Load(file, printinfo);
  } else {
    this->computeSizes();
  }
}

/// Loads a geometry from a file
void Geometry::Load(const char file[], const bool printinfo) {
  std::string param;
  real_t value;
  std::ifstream in(file);
  if(!in) {
    std::cerr << "Geometry: Could not open file \"" << file << "\"!" << std::endl;
  } else if(printinfo) {
    std::cout << "Geometry: Loading from file \"" << file << "\"." << std::endl;
  }
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
      if(printinfo) {
        std::cout << "Geometry: Load size = " << this->_totalSize << std::endl;
      }
    } else if(!param.compare("Length") || !param.compare("length")) {
      this->_totalLength[0] = value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value;
        this->_totalLength[dim] = value;
      }
      if(printinfo) {
        std::cout << "Geometry: Load length = " << this->_totalLength << std::endl;
      }
    } else if(!param.compare("Velocity") || !param.compare("velocity")) {
      this->_velocity[0] = value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value;
        this->_velocity[dim] = value;
      }
      if(printinfo) {
        std::cout << "Geometry: Load velocity = " << this->_velocity << std::endl;
      }
    } else if(!param.compare("Pressure") || !param.compare("pressure")) {
      this->_pressure = value;
      if(printinfo) {
        std::cout << "Geometry: Load pressure = " << this->_pressure << std::endl;
      }
    } else if(!param.compare("Geometry") || !param.compare("geometry")) {
      this->_free = true;
      in.ignore(10000, '\n');  // removes remaining line
      if(this->_flags != nullptr) {
        delete[] this->_flags;
      }
      this->_flags = new char[this->_totalSize[0] * this->_totalSize[1]];
      for(index_t j = this->_totalSize[1]; j-- > 0;) {
        in.read(&this->_flags[this->_totalSize[0] * j], this->_totalSize[0]);
        in >> std::ws; // removes newline between lines
      }
      if(printinfo) {
        std::cout << "Geometry: Load free geometry." << std::endl;
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
/*
void Geometry::Update_U(Grid &u) const {
  this->_comm.copyBoundary(u);
  // driven cavity u given at top
  BoundaryIterator it(*this, BoundaryIterator::boundary::left);
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
*/
// Updates the velocity field v
/*void Geometry::Update_V(Grid &v) const {
  // driven cavity v given at top
  this->_comm.copyBoundary(v);

  BoundaryIterator it(*this, BoundaryIterator::boundary::left);
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
*/
// Updates the pressure field p
/*
void Geometry::Update_P(Grid &p) const {
  // homogenous Neumann condition at all sides
  this->_comm.copyBoundary(p);

  BoundaryIterator it(*this, BoundaryIterator::boundary::left);
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
*/
void Geometry::Update_P(Grid &p) const {
  this->_comm.copyBoundary(p);

  for(Iterator it(*this); it.Valid(); it.Next()) {
    switch(this->flag(it)) {
      case '#': // '#' Wall/Obstacle/NoSlip boundary (u = v = 0, dp/dn = 0)
        if(this->isFluid(it.Left())) {
          if(this->isFluid(it.Top())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top())) {
          if(this->isFluid(it.Right())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right())) {
          if(this->isFluid(it.Down())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down())) {
          if(this->isFluid(it.Left())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'I': // General Inflow boundary (u = u_0, v = v_0, dp/dn = 0)
        if(this->isFluid(it.Left())) {
          if(this->isFluid(it.Top())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Left()) + p.Cell(it.Top()));
          } else {
            p.Cell(it) = p.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top())) {
          if(this->isFluid(it.Right())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Top()) + p.Cell(it.Right()));
          } else {
            p.Cell(it) = p.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right())) {
          if(this->isFluid(it.Down())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Right()) + p.Cell(it.Down()));
          } else {
            p.Cell(it) = p.Cell(it.Down());
          }
        } else if(this->isFluid(it.Down())) {
          if(this->isFluid(it.Left())) {
            p.Cell(it) = 0.5 * (p.Cell(it.Down()) + p.Cell(it.Left()));
          } else {
            p.Cell(it) = p.Cell(it.Down());
          }
        }
        break;
      case 'V': // Vertical Inflow boundary (u = u_0, v = v_0, but only fluid right or left)
        if(this->isFluid(it.Left())) {
          p.Cell(it) = p.Cell(it.Left());
        } else if(this->isFluid(it.Right())) {
          p.Cell(it) = p.Cell(it.Right());
        }
        break;
      case 'H': // Horizontal Inflow boundary (u = u_0, v = v_0, but only fluid top or down)
        if(this->isFluid(it.Top())) {
          p.Cell(it) = p.Cell(it.Top());
        } else if(this->isFluid(it.Down())) {
          p.Cell(it) = p.Cell(it.Down());
        }
        break;
      case 'O': // Outflow boundary (d/dn (u,v) = 0)
        p.Cell(it) = 0;
        break;
      case '|': // Vertical Slip-boundary (du/dx = 0, v = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Left())) {
          p.Cell(it) = 2*this->_pressure - p.Cell(it.Left());
        } else if(this->isFluid(it.Right())) {
          p.Cell(it) = 2*this->_pressure - p.Cell(it.Right());
        }
        break;
      case '-': // Horizontal Slip-boundary (u = 0, dv/dy = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Top())) {
          p.Cell(it) = 2*this->_pressure - p.Cell(it.Top());
        } else if(this->isFluid(it.Down())) {
          p.Cell(it) = 2*this->_pressure - p.Cell(it.Down());
        }
        break;
    }
  }
}

void Geometry::Update(Grid &u, Grid &v) const {
  this->_comm.copyBoundary(u);
  this->_comm.copyBoundary(v);

  for(Iterator it(*this); it.Valid(); it.Next()) {
    switch(this->flag(it)) {
      case '#': // '#' Wall/Obstacle/NoSlip boundary (u = v = 0, dp/dn = 0)
        if(this->isFluid(it.Left())) {
          u.Cell(it.Left()) = 0;
          if(this->isFluid(it.Top())) {
            u.Cell(it) = -u.Cell(it.Top());
            v.Cell(it) = 0;
          } else {
            v.Cell(it) = -v.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top())) {
          v.Cell(it) = 0;
          if(this->isFluid(it.Right())) {
            u.Cell(it) = 0;
          } else {
            u.Cell(it) = -u.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right())) {
          u.Cell(it) = 0;
          v.Cell(it) = -v.Cell(it.Right());
          if(this->isFluid(it.Down())) {
            v.Cell(it.Down()) = 0;
          }
        } else if(this->isFluid(it.Down())) {
          v.Cell(it.Down()) = 0;
          u.Cell(it) = -u.Cell(it.Down());
          if(this->isFluid(it.Left())) {
            u.Cell(it.Left()) = 0;
            v.Cell(it) = - v.Cell(it.Left());
          }
        }
        break;
      case 'I': // General Inflow boundary (u = u_0, v = v_0, dp/dn = 0)
        if(this->isFluid(it.Left())) {
          u.Cell(it) = this->_velocity[0];
          if(this->isFluid(it.Top())) {
            v.Cell(it) = this->_velocity[1];
          } else {
            v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Left());
          }
        } else if(this->isFluid(it.Top())) {
          v.Cell(it) = this->_velocity[1];
          if(this->isFluid(it.Right())) {
            u.Cell(it) = this->_velocity[0];
          } else {
            u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Top());
          }
        } else if(this->isFluid(it.Right())) {
          u.Cell(it) = this->_velocity[0];
          if(this->isFluid(it.Down())) {
            v.Cell(it) = this->_velocity[1];
          } else {
            v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Right());
          }
        } else if(this->isFluid(it.Down())) {
          v.Cell(it) = this->_velocity[1];
          if(this->isFluid(it.Left())) {
            u.Cell(it.Left()) = this->_velocity[0];
          } else {
            u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Down());
          }
        }
        break;
      case 'V': // Vertical Inflow boundary (u = u_0, v = v_0, but only fluid right or left)
        // TODO: other than const velocities
        if(this->isFluid(it.Left())) {
          u.Cell(it) = this->_velocity[0];
          v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Left());
        } else if(this->isFluid(it.Right())) {
          u.Cell(it) = parabolic(it);//this->_velocity[0];
          v.Cell(it) = 2*this->_velocity[1] - v.Cell(it.Right());
        }
        break;
      case 'H': // Horizontal Inflow boundary (u = u_0, v = v_0, but only fluid top or down)
        // TODO: other than const velocities
        if(this->isFluid(it.Top())) {
          v.Cell(it) = this->_velocity[1];
          u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Top());
        } else if(this->isFluid(it.Down())) {
          v.Cell(it) = this->_velocity[1];
          u.Cell(it) = 2*this->_velocity[0] - u.Cell(it.Down());
        }
        break;
      case 'O': // Outflow boundary (d/dn (u,v) = 0)
        if(this->isFluid(it.Left())) {
          u.Cell(it) = u.Cell(it.Left());
          v.Cell(it) = v.Cell(it.Left());
        } else if(this->isFluid(it.Top())) {
          v.Cell(it) = v.Cell(it.Top());
          u.Cell(it) = u.Cell(it.Top());
        } else if(this->isFluid(it.Right())) {
          u.Cell(it) = u.Cell(it.Right());
          v.Cell(it) = v.Cell(it.Right());
        } else if(this->isFluid(it.Down())) {
          v.Cell(it) = v.Cell(it.Down());
          u.Cell(it) = u.Cell(it.Down());
        }
        break;
      case '|': // Vertical Slip-boundary (du/dx = 0, v = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Left())) {
          u.Cell(it) = u.Cell(it.Left());
          v.Cell(it) = 0;
        } else if(this->isFluid(it.Right())) {
          u.Cell(it) = u.Cell(it.Right());
          v.Cell(it) = 0;
        }
        break;
      case '-': // Horizontal Slip-boundary (u = 0, dv/dy = 0, dp determined by parameter pressure)
        if(this->isFluid(it.Top())) {
          u.Cell(it) = 0;
          v.Cell(it) = v.Cell(it.Top());
        } else if(this->isFluid(it.Down())) {
          u.Cell(it) = 0;
          v.Cell(it) = v.Cell(it.Down());
        }
        break;
    }
  }
}

// compute grid sizes
void Geometry::computeSizes() {
  const multi_index_t& numProc = this->_comm.ThreadDim();
  const multi_index_t& tIdx = this->_comm.ThreadIdx();

  if(! this->_free) { // initialize Driven Cavity
    this->_flags = new char[this->_totalSize[0] * this->_totalSize[1]];
    for(index_t j = 0; j < this->_totalSize[1]; j++) {
      this->_flags[this->_totalSize[0] * j] = (j == this->_totalSize[1]-1? 'H' : '#');
      for(index_t i = 1; i < this->_totalSize[0]-1; i++) {
        this->_flags[i + this->_totalSize[0] * j] = (j == 0? '#' : (j == this->_totalSize[1]-1? 'H' : ' ') );
      }
      this->_flags[this->_totalSize[0] * (j+1)-1] = (j == this->_totalSize[1]-1? 'H' : '#');
    }
    if(this->_comm.ThreadNum() == 0) {
      std::cout << "Geometry: Load Driven Cavity." << std::endl;
    }
  }

  this->_sizeData = 1;
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_totalSize[dim] -= 2;
    this->_size[dim] = this->_totalSize[dim] / numProc[dim];
    this->_offset[dim] = tIdx[dim] * this->_size[dim];

    // fix last grid size
    if(this->_totalSize[dim] != this->_size[dim] * numProc[dim]) {
      if(tIdx[dim] == numProc[dim]-1) {
        this->_size[dim] += this->_totalSize[dim] % numProc[dim];
        std::cerr << "Grid not divisible without remainder for " << this->_comm.ThreadCnt() << " Threads!" << std::endl;
      }
    }

    this->_length[dim] = this->_totalLength[dim] * this->_size[dim] / (real_t)this->_totalSize[dim];
    this->_sizeP[dim] = this->_size[dim]+2;
    this->_sizeData *= this->_sizeP[dim];
    this->_h[dim] = this->_totalLength[dim] / this->_totalSize[dim];
    this->_hInv[dim] = 1.0 / this->_h[dim];
  }

  // only use local flag field
  char *flags = new char[this->_sizeP[0] * this->_sizeP[1]];
  for(index_t j = 0; j < this->_sizeP[1]; ++j) {
    strncpy(&flags[j * _sizeP[0]],
            &this->_flags[this->_offset[0] + (this->_offset[1]+j)*(_totalSize[0]+2)], _sizeP[0]);
  }
  delete[] this->_flags;
  this->_flags = flags;

  // update eXchange boundaries in local flag field
  BoundaryIterator bit(*this, BoundaryIterator::boundary::left);
  for(bit.First(); bit.Valid(); bit.Next()) {
    if(this->isFluid(bit)) {
      this->_flags[bit] = 'X';
    }
  }
  bit.SetBoundary(BoundaryIterator::boundary::down);
  for(bit.First(); bit.Valid(); bit.Next()) {
    if(this->isFluid(bit)) {
      this->_flags[bit] = 'X';
    }
  }
  bit.SetBoundary(BoundaryIterator::boundary::right);
  for(bit.First(); bit.Valid(); bit.Next()) {
    if(this->isFluid(bit)) {
      this->_flags[bit] = 'X';
    }
  }
  bit.SetBoundary(BoundaryIterator::boundary::top);
  for(bit.First(); bit.Valid(); bit.Next()) {
    if(this->isFluid(bit)) {
      this->_flags[bit] = 'X';
    }
  }

  this->_N = 0;
  // compute local number of fluid cells
  for(InteriorIterator it(*this); it.Valid(); it.Next()) {
    this->_N++;
  }
}

const char& Geometry::flag(const Iterator &it) const {
  return this->_flags[it];
}

// Returns the parabolic velocity profile (u_0*y*(y-L_y)).
// TODO: rewrite for parallelization
real_t Geometry::parabolic(const Iterator &it) const {
  const real_t y = (it.Pos()[1]-0.5) * this->_h[1];
  return this->_velocity[0] * y * (this->_h[1]*this->_size[1] - y);
}

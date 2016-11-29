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
#include <cmath>
#include "typedef.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "solver.hpp"

// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator &it, const Grid &grid, const Grid &rhs) const {
  return (rhs.Cell(it) - grid.dxx(it) - grid.dyy(it));
}

//------------------------------------------------------------------------------
// concrete SOR solver

// Constructs an actual SOR solver
SOR::SOR(const Geometry &geom, const real_t &omega) : Solver(geom), _omega(omega) {
  if(this->_omega <= 0.0 || this->_omega > 2.0) {
    real_t h = std::fmax(this->_geom.Mesh()[0], this->_geom.Mesh()[1]);
    this->_omega = 2.0 / (1.0 + std::sin(M_PI * h));
    std::cout << "SOR: New omega = " << this->_omega << std::endl;
  }
}

// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &h = this->_geom.Mesh();
  const real_t prefactor = this->_omega * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]);

  real_t residual = 0;
  for(InteriorIterator it(grid); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - prefactor * localRes;
  }
  return residual / (h[0] * h[1]);
}

//------------------------------------------------------------------------------
// concrete Red or Black SOR solver

real_t RedOrBlackSOR::RedCycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &h = this->_geom.Mesh();
  const real_t prefactor = this->_omega * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]);

  real_t residual = 0;

  for(InteriorIterator it(grid); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - prefactor * localRes;
    it.Next();
  }

  return residual / (h[0] * h[1]);
}

real_t RedOrBlackSOR::BlackCycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &h = this->_geom.Mesh();
  const real_t prefactor = this->_omega * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]);

  real_t residual = 0;

  InteriorIterator it(grid);
  for(it.Next(); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - prefactor * localRes;
    it.Next();
  }

  return residual / (h[0] * h[1]);
}

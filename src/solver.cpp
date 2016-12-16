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
SOR::SOR(const Geometry &geom, const real_t &omega) : Solver(geom),
    _correction(omega * 0.5 * (geom.Mesh()[0]*geom.Mesh()[0]*geom.Mesh()[1]*geom.Mesh()[1])
                / (geom.Mesh()[0]*geom.Mesh()[0]+geom.Mesh()[1]*geom.Mesh()[1])) {}

// Returns the total residual and executes a solver cycle
real_t SOR::Cycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &hInv = this->_geom.invMesh();
  real_t residual = 0;

  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - this->_correction * localRes;
  }

  return residual * hInv[0] * hInv[1];
}

//------------------------------------------------------------------------------
// concrete Red or Black SOR solver

real_t RedOrBlackSOR::RedCycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &hInv = this->_geom.invMesh();
  real_t residual = 0;

  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - this->_correction * localRes;
    it.Next();
  }

  return residual * hInv[0] * hInv[1];
}

real_t RedOrBlackSOR::BlackCycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &hInv = this->_geom.invMesh();
  real_t residual = 0;

  InteriorIterator it(this->_geom);
  for(it.Next(); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - this->_correction * localRes;
    it.Next();
  }

  return residual * hInv[0] * hInv[1];
}

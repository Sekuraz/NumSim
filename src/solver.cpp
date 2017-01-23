/* vim: set tabstop=2:softtabstop=2:shiftwidth=2:expandtab */
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
#include "comm.hpp"
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
SOR::SOR(const Geometry &geom, const Communicator &comm, const real_t &omega)
  : Solver(geom, comm),
    _correction(omega * 0.5 * (geom.Mesh()[0]*geom.Mesh()[0]*geom.Mesh()[1]*geom.Mesh()[1])
                / (geom.Mesh()[0]*geom.Mesh()[0]+geom.Mesh()[1]*geom.Mesh()[1])) {}

// Returns the total residual and executes a solver cycle
real_t SOR::Cycle(Grid &grid, const Grid &rhs) const {
  real_t residual = 0;

  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) -= this->_correction * localRes;
  }

  residual *= this->_invNumFluid;
  // gather global residual
  return this->_comm.gatherSum(residual);
}

//------------------------------------------------------------------------------
// concrete Red or Black SOR solver

real_t RedOrBlackSOR::Cycle(Grid &grid, const Grid &rhs) const {
  real_t res = 0;
  if(this->_firstRed) {
    // first half-step
    res = this->RedCycle(grid, rhs);
    // exchange boundary values for p
    this->_comm.copyBoundary(grid);
    // second half-step
    res += this->BlackCycle(grid, rhs);
  } else {
    // first half-step
    res = this->BlackCycle(grid, rhs);
    // exchange boundary values for p
    this->_comm.copyBoundary(grid);
    // second half-step
    res += this->RedCycle(grid, rhs);
  }
  // gather global residual
  return this->_comm.gatherSum(res);
}

real_t RedOrBlackSOR::RedCycle(Grid &grid, const Grid &rhs) const {
  real_t residual = 0;

  /*for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) -= this->_correction * localRes;
    it.Next();
  }*/
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    const multi_index_t pos = it.Pos();
    if((pos[0]+pos[1]) % 2 == 0) {
      real_t localRes = Solver::localRes(it, grid, rhs);
      residual += localRes * localRes;
      grid.Cell(it) -= this->_correction * localRes;
    }
  }

  return residual * this->_invNumFluid;
}

real_t RedOrBlackSOR::BlackCycle(Grid &grid, const Grid &rhs) const {
  real_t residual = 0;

  /*InteriorIterator it(this->_geom);
  for(it.Next(); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residual += localRes * localRes;
    grid.Cell(it) -= this->_correction * localRes;
    it.Next();
  }*/
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    const multi_index_t pos = it.Pos();
    if((pos[0]+pos[1]) % 2 == 1) {
      real_t localRes = Solver::localRes(it, grid, rhs);
      residual += localRes * localRes;
      grid.Cell(it) -= this->_correction * localRes;
    }
  }

  return residual * this->_invNumFluid;
}

void CG::reset(const Grid &grid, const Grid &rhs) {
  real_t local;
  for (InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    local = Solver::localRes(it, grid, rhs);

    this->_res.Cell(it) = local;
    this->_direction.Cell(it) = local;
  }
  this->old_residual = this->_res.InnerProduct(this->_res);
  // gather global residual
  this->old_residual = this->_comm.gatherSum(this->old_residual);
}

real_t CG::Cycle(Grid &grid, const Grid &rhs __attribute__((unused))) const {
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    this->_Ad.Cell(it) = this->_direction.dxx(it) + this->_direction.dyy(it);
  }

  real_t dTAd = this->_direction.InnerProduct(this->_Ad);
  real_t alpha = this->old_residual / this->_comm.gatherSum(dTAd);

  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    grid.Cell(it) += alpha * this->_direction.Cell(it);
    this->_res.Cell(it) -= alpha * this->_Ad.Cell(it);
  }

  real_t residual = this->_res.InnerProduct(this->_res);
  // gather global residual
  residual = this->_comm.gatherSum(residual);
  real_t beta = residual / this->old_residual;

  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    this->_direction.Cell(it) *= beta;
    this->_direction.Cell(it) += this->_res.Cell(it);
  }
  this->_comm.copyBoundary(this->_direction);

  this->old_residual = residual;
  return residual * this->_invNumFluid;
}

//------------------------------------------------------------------------------
// Constructs an actual Multigrid solver
MG::MG(const Geometry &geom, const Communicator &comm, const index_t level, const index_t& nu)
    : Solver(geom, comm), _level(level), _nu(nu), _smoother(geom, comm, 1.0),
      _coarse(nullptr), _e(nullptr), _res(nullptr) {
  if(this->_level > 0) {
    Geometry* geom_coarse = geom.coarse();
    this->_e = new Grid(*geom_coarse);
    this->_res = new Grid(*geom_coarse);
    this->_coarse = new MG(*geom_coarse, comm, this->_level-1);
  }
}

MG::~MG() {
  if(this->_level > 0) {
    delete &this->_coarse->_geom;
    delete this->_coarse;
    delete this->_e;
    delete this->_res;
  }
}

// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t MG::Cycle(Grid &p, const Grid &rhs) const {
  this->Smooth(p, rhs);
  if(this->_level > 0) {
    this->Restrict(p, rhs);
    this->_coarse->Cycle(*this->_e, *this->_res);
    this->Interpolate(p);
  } // TODO: else solve
  return this->Smooth(p, rhs);
}

// Restricts the residuals of the solution to the next coarser Grid
void MG::Restrict(const Grid &p, const Grid &rhs) const {
  const Geometry &cGeom = this->_coarse->_geom;
  // compute residuals and restrict
  for(Iterator it(cGeom); it.Valid(); it.Next()) {
    this->_e->Cell(it) = 0;
    multi_index_t pos_fine = it.Pos();
    for(index_t dim = 0; dim < DIM; dim++) {
      pos_fine[dim] = (pos_fine[dim] == 0) ? 0 : 2*pos_fine[dim]-1;
    }
    const Iterator itFine(this->_geom, pos_fine);
    switch(cGeom.flag(it)) {
    case ' ':
      this->_res->Cell(it) = 0.25
        * ( this->localRes(itFine, p, rhs) + this->localRes(itFine.Right(), p, rhs)
          + this->localRes(itFine.Top(), p, rhs) + this->localRes(itFine.Top().Right(), p, rhs) );
      break;
    case '#': // '#' Wall/Obstacle/NoSlip boundary (u = v = 0, dp/dn = 0)
    case 'I': // General Inflow boundary (u = u_0, v = v_0, dp/dn = 0)
      if(cGeom.isFluid(it.Left())) {
        if(cGeom.isFluid(it.Top())) {
          this->_res->Cell(it) = 0.5 * (p.dx_l(itFine) + p.dy_r(itFine));
        } else {
          this->_res->Cell(it) = p.dx_l(itFine);
        }
      } else if(cGeom.isFluid(it.Top())) {
        if(cGeom.isFluid(it.Right())) {
          this->_res->Cell(it) = 0.5 * (p.dy_r(itFine) + p.dx_r(itFine));
        } else {
          this->_res->Cell(it) = p.dy_r(itFine);
        }
      } else if(cGeom.isFluid(it.Right())) {
        if(cGeom.isFluid(it.Down())) {
          this->_res->Cell(it) = 0.5 * (p.dx_r(itFine) + p.dy_l(itFine));
        } else {
          this->_res->Cell(it) = p.dx_r(itFine);
        }
      } else if(cGeom.isFluid(it.Down())) {
        if(cGeom.isFluid(it.Left())) {
          this->_res->Cell(it) = 0.5 * (p.dy_l(itFine) + p.dx_l(itFine));
        } else {
          this->_res->Cell(it) = p.dy_l(itFine);
        }
      }
      break;
    case 'V': // Vertical Inflow boundary (u = u_0, v = v_0, but only fluid right or left)
      if(cGeom.isFluid(it.Left())) {
        this->_res->Cell(it) = p.dx_l(itFine);
      } else if(cGeom.isFluid(it.Right())) {
        this->_res->Cell(it) = p.dx_r(itFine);
      }
      break;
    case 'H': // Horizontal Inflow boundary (u = u_0, v = v_0, but only fluid top or down)
      if(cGeom.isFluid(it.Top())) {
        this->_res->Cell(it) = p.dy_r(itFine);
      } else if(cGeom.isFluid(it.Down())) {
        this->_res->Cell(it) = p.dy_l(itFine);
      }
      break;
    case 'O': // Outflow boundary (d/dn (u,v) = 0)
      this->_res->Cell(it) = p.Cell(itFine) - rhs.Cell(itFine);
      break;
    case '|': // Vertical Slip-boundary (du/dx = 0, v = 0, dp determined by parameter pressure)
      if(cGeom.isFluid(it.Left())) {
        // TODO: ???
        //p.Cell(it) = 2*this->_pressure - p.Cell(it.Left()) + rhs.Cell(it) * this->_h[0];
        this->_res->Cell(it) = 0.5*(p.Cell(itFine)+p.Cell(itFine.Left())) - rhs.Cell(it);
      } else if(cGeom.isFluid(it.Right())) {
        // TODO: ???
        //p.Cell(it) = 2*this->_pressure - p.Cell(it.Right()) - rhs.Cell(it) * this->_h[0];
        this->_res->Cell(it) = 0.5*(p.Cell(itFine)+p.Cell(itFine.Right())) - rhs.Cell(it);
      }
      break;
    case '-': // Horizontal Slip-boundary (u = 0, dv/dy = 0, dp determined by parameter pressure)
      if(cGeom.isFluid(it.Top())) {
        // TODO: ???
        //p.Cell(it) = 2*this->_pressure - p.Cell(it.Top()) - rhs.Cell(it) * this->_h[1];
        this->_res->Cell(it) = 0.5*(p.Cell(itFine)+p.Cell(itFine.Top())) - rhs.Cell(it);
      } else if(cGeom.isFluid(it.Down())) {
        // TODO: ???
        //p.Cell(it) = 2*this->_pressure - p.Cell(it.Down()) + rhs.Cell(it) * this->_h[1];
        this->_res->Cell(it) = 0.5*(p.Cell(itFine)+p.Cell(itFine.Down())) - rhs.Cell(it);
      }
      break;
    }
  }
}

// Interpolates and adds from the coarser solution to this one
void MG::Interpolate(Grid &p) const {
  // add e to p (while interpolating)
  const Geometry &cGeom = this->_coarse->_geom;
  for(Iterator it(cGeom); it.Valid(); it.Next()) {
    multi_index_t pos_fine = it.Pos();
    for(index_t dim = 0; dim < DIM; dim++) {
      pos_fine[dim] = (pos_fine[dim] == 0) ? 0 : 2*pos_fine[dim]-1;
    }
    const Iterator itFine(this->_geom, pos_fine);
    switch(cGeom.flag(it)) {
    case ' ':
      p.Cell(itFine) += this->_e->Cell(it);
      p.Cell(itFine.Right()) += this->_e->Cell(it);
      p.Cell(itFine.Top()) += this->_e->Cell(it);
      p.Cell(itFine.Top().Right()) += this->_e->Cell(it);
      break;
    case '#': // '#' Wall/Obstacle/NoSlip boundary (u = v = 0, dp/dn = 0)
    case 'I': // General Inflow boundary (u = u_0, v = v_0, dp/dn = 0)
      if(cGeom.isFluid(it.Left())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      } else if(cGeom.isFluid(it.Top())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: T/RT ?
      } else if(cGeom.isFluid(it.Right())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      } else if(cGeom.isFluid(it.Down())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: T/RT ?
      }
      break;
    case 'V': // Vertical Inflow boundary (u = u_0, v = v_0, but only fluid right or left)
      if(cGeom.isFluid(it.Left())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      } else if(cGeom.isFluid(it.Right())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      }
      break;
    case 'H': // Horizontal Inflow boundary (u = u_0, v = v_0, but only fluid top or down)
      if(cGeom.isFluid(it.Top())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: T/RT ?
      } else if(cGeom.isFluid(it.Down())) {
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: T/RT ?
      }
      break;
    case 'O': // Outflow boundary (d/dn (u,v) = 0)
      p.Cell(itFine) += this->_e->Cell(it);
      // TODO: R/T/RT ?
      break;
    case '|': // Vertical Slip-boundary (du/dx = 0, v = 0, dp determined by parameter pressure)
      if(cGeom.isFluid(it.Left())) {
        // TODO: ???
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      } else if(cGeom.isFluid(it.Right())) {
        // TODO: ???
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Top()) += this->_e->Cell(it);
        // TODO: R/RT ?
      }
      break;
    case '-': // Horizontal Slip-boundary (u = 0, dv/dy = 0, dp determined by parameter pressure)
      if(cGeom.isFluid(it.Top())) {
        // TODO: ???
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: R/RT ?
      } else if(cGeom.isFluid(it.Down())) {
        // TODO: ???
        p.Cell(itFine) += this->_e->Cell(it);
        p.Cell(itFine.Right()) += this->_e->Cell(it);
        // TODO: R/RT ?
      }
      break;
    }
  }
}

real_t MG::Smooth(Grid &p, const Grid &rhs) const {
  real_t res = 1;
  for(index_t n = 0; n < this->_nu; n++) {
    res = this->_smoother.Cycle(p, rhs);
    this->_geom.Update_P(p, rhs);
  }
  return res;
}

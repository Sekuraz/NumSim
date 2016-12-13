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
#include <algorithm>
#include <list>
#include <cmath>
#include "typedef.hpp"
#include "comm.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"

// Creates a compute instance with given geometry, parameter and communicator
Compute::Compute(const Geometry &geom, const Parameter &param, const Communicator &comm)
    : _t(0),
    _u(new Grid(geom, multi_real_t(0, geom.Mesh()[1]/2))),
    _v(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, 0))),
    _p(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, geom.Mesh()[1]/2))),
    _F(new Grid(geom)), _G(new Grid(geom)), _rhs(new Grid(geom)),
    _velocities(new Grid(geom)), _streamline(new Grid(geom)),
    _vorticity(new Grid(geom)), _particle(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, geom.Mesh()[1]/2))),
    _firstRed(comm.EvenOdd() && (geom.Size()[0] % 2 == 0)), // TODO one size even other odd
    _geom(geom), _param(param), _comm(comm),
    _initPosParticle(param.ParticleInitPos()) {

  // initialize the solver
  this->_solver = new RedOrBlackSOR(geom, param.Omega());

  // initialize u,v,p,particle
  this->_u->Initialize(0);
  this->_v->Initialize(0);
  this->_p->Initialize(0);
  this->_particle->Initialize(0);

  // write first particle pos in list for particle trace and streakline
  _particleTracing.push_back(this->_initPosParticle);
  _streakline.push_back(this->_initPosParticle);

}
// Deletes all grids
Compute::~Compute() {
  delete _u; delete _v; delete _p; delete _F; 
  delete _G; delete _rhs; delete _velocities; delete _solver;
  delete _streamline; delete _vorticity; delete _particle;
}

// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
  // set boundary values for u, v, p
  this->_geom.Update(*(this->_u), *(this->_v), *(this->_p));

  // compute local dt
  // Test CFL and Pr condition
  const multi_real_t &h = this->_geom.Mesh();
  this->_dtlimit = this->_param.Tau()
                 * std::min(std::min(h[0] / this->_u->AbsMax(), h[1] / this->_v->AbsMax()),
                            this->_param.Re() * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]));
  real_t dt = std::min(this->_dtlimit, this->_param.Dt());
  // gather global dt
  dt = this->_comm.gatherMin(dt);

  // compute the momentum equations
  this->MomentumEqu(dt);

  // compute right-hand-side of the poisson equation
  this->RHS(dt);

  // solve the poisson equation for the pressure
  real_t res = 2 * this->_param.Eps() * this->_param.Eps();
  index_t i;
  for(i = 0; (i < this->_param.IterMax()) && (res > this->_param.Eps() * this->_param.Eps()) ; i++) {
    // first half-step
    if(this->_firstRed) {
      res = this->_solver->RedCycle(*(this->_p), *(this->_rhs));
    } else {
      res = this->_solver->BlackCycle(*(this->_p), *(this->_rhs));
    }
    // set boundary values for p
    this->_comm.copyBoundary(*(this->_p));
    // second half-step
    if(this->_firstRed) {
      res += this->_solver->BlackCycle(*(this->_p), *(this->_rhs));
    } else {
      res += this->_solver->RedCycle(*(this->_p), *(this->_rhs));
    }
    // set boundary values for p
    // TODO: only update p
    this->_geom.Update(*(this->_u), *(this->_v), *(this->_p));

    // gather global residual
    res = this->_comm.gatherSum(res);
  }
  if((i == this->_param.IterMax()) && printInfo) {
    std::cerr << "Warning: SOR did not converge! res = " << std::sqrt(res) << std::endl;
  }

  // compute new velocites
  this->NewVelocities(dt);

  // calculate streamlines and vorticity
  this->Stream();
  this->Vort();

  // calculate new particle postitions
  // TODO: rewrite for parallelization
  if (_comm.ThreadCnt() == 1) {
    this->Particle();
    this->Streaklines();
  }

  // print information
  if(printInfo) {
    std::cout.precision(4);
    std::cout << "t = " << this->_t << " \tdt = " << dt << "\titer = " << i
              << "\tres = " << std::scientific << std::sqrt(res)
              << std::defaultfloat << std::endl;
  }

  // update the time
  this->_t += dt;
}

// Computes and returns the absolute velocity
const Grid *Compute::GetVelocity() {
  for(Iterator it(*(this->_velocities)); it.Valid(); it.Next()) {
    real_t uMean = (this->_u->Cell(it) + this->_u->Cell(it.Top()))/2;
    real_t vMean = (this->_v->Cell(it) + this->_v->Cell(it.Right()))/2;
    this->_velocities->Cell(it) = std::sqrt(uMean*uMean + vMean*vMean);
  }
  return _velocities;
}

// Computes and returns the vorticity
void Compute::Vort() {
  for(Iterator it(*(this->_vorticity)); it.Valid(); it.Next()) {
    this->_vorticity->Cell(it) = this->_u->dy_r(it) - this->_v->dx_r(it);
  }
}
// Computes the stream line values
void Compute::Stream() {
  this->_streamline->Data()[0] = 0;
  for(Iterator it(*(this->_streamline), 1); it.Valid(); it.Next()) {
    if(it.Pos()[1] == 0) {
      this->_streamline->Cell(it) = this->_streamline->Cell(it.Left()) - _geom.Mesh()[0] * this->_v->Cell(it);
    } else {
      this->_streamline->Cell(it) = this->_streamline->Cell(it.Down()) + _geom.Mesh()[1] * this->_u->Cell(it);
    }
  }
  real_t offset = this->_comm.copyOffset(*this->_streamline);
  for(Iterator it(*(this->_streamline)); it.Valid(); it.Next()) {
    this->_streamline->Cell(it) += offset;
  }
}

void Compute::Particle() {
  // returns last particlePosition of the list
  multi_real_t lastPos = _particleTracing.back();

  // compute new step without live visualisation
  //this->ParticleStep(lastPos);

  // compute new step with live visualisation
  this->ParticleStepVisu(lastPos);

  // write new position in the list
  _particleTracing.push_back(lastPos);

}

void Compute::Streaklines(){  

  // Cycle the list of particles
  for (std::list<multi_real_t>::iterator it = _streakline.begin(); it != _streakline.end(); ++it) {
    this->ParticleStep((*it));
  }

  // Add new particle at the initial position
  std::list<multi_real_t>::iterator it = _streakline.begin();
  _streakline.insert (_streakline.begin(), this->_initPosParticle);

}	

// ParticleStep() without live visualisation
void Compute::ParticleStep(multi_real_t &lastPos) {

  // calculate new position
  lastPos[0] += _dtlimit * this->_u->Interpolate(lastPos);
  lastPos[1] += _dtlimit * this->_v->Interpolate(lastPos);

}

// ParticleStepVisu() with live visualisation (only valid for particle tracing)
void Compute::ParticleStepVisu(multi_real_t &lastPos) {

  // get cell of the particle
  _particleIndx[0] = (index_t)(lastPos[0]/_geom.Mesh()[0] + 1);
  _particleIndx[1] = (index_t)(lastPos[1]/_geom.Mesh()[1] + 1);

  // set old position to 0 for debug visualization
  Iterator it1(*(this->_particle), _particleIndx);
  this->_particle->Cell(it1) = 0;

  // calculate new position
  lastPos[0] += _dtlimit * this->_u->Interpolate(lastPos);
  lastPos[1] += _dtlimit * this->_v->Interpolate(lastPos);
  // write new position in the list

  // get cell of the particle
  _particleIndx[0] = (index_t)(lastPos[0]/_geom.Mesh()[0] + 1);
  _particleIndx[1] = (index_t)(lastPos[1]/_geom.Mesh()[1] + 1);

  // set new position to 1 for debug visualization
  Iterator it2(*(this->_particle), _particleIndx);
  this->_particle->Cell(it2) = 1;

}

// Compute the new velocites u,v
void Compute::NewVelocities(const real_t &dt) {
  // compute u
  for(InteriorIterator it(*(this->_u)); it.Valid(); it.Next()) {
    this->_u->Cell(it) = this->_F->Cell(it) - dt * this->_p->dx_r(it);
  }
  // compute v
  for(InteriorIterator it(*(this->_v)); it.Valid(); it.Next()) {
    this->_v->Cell(it) = this->_G->Cell(it) - dt * this->_p->dy_r(it);
  }
}

// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt) {
  for(InteriorIterator it(*(this->_F)); it.Valid(); it.Next()) {
    this->_F->Cell(it) = this->_u->Cell(it) + dt * (
        ( this->_u->dxx(it) + this->_u->dyy(it) )/this->_param.Re()
         - this->_u->DC_udu_x(it, this->_param.Alpha())
        - this->_u->DC_vdu_y(it, this->_param.Alpha(), this->_v) );
    this->_G->Cell(it) = this->_v->Cell(it) + dt * (
          ( this->_v->dxx(it) + this->_v->dyy(it) )/this->_param.Re()
          - this->_v->DC_vdv_y(it, this->_param.Alpha())
          - this->_v->DC_udv_x(it, this->_param.Alpha(), this->_u) );
  }

  // boundary values of F and G
  this->_geom.Update(*this->_F, *this->_G, *this->_p);
}

// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt) {
  for(InteriorIterator it(*(this->_p)); it.Valid(); it.Next()) {
    this->_rhs->Cell(it) = (this->_F->dx_l(it) + this->_G->dy_l(it)) / dt;
  }
}


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
#include <iterator>
#include <cmath>
#include <cstring>
#include "typedef.hpp"
#include "comm.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"

// Creates a compute instance with given geometry, parameter and communicator
Compute::Compute(const Geometry &geom, const Parameter &param, const Communicator &comm, const char sol[])
    : _t(0),
    _u(new Grid(geom, multi_real_t(0, geom.Mesh()[1]/2))),
    _v(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, 0))),
    _p(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, geom.Mesh()[1]/2))),
    _F(new Grid(geom)), _G(new Grid(geom)), _rhs(new Grid(geom)),
    _velocities(new Grid(geom)), _streamline(new Grid(geom)),
    _vorticity(new Grid(geom)), _particle(new Grid(geom, multi_real_t(geom.Mesh()[0]/2, geom.Mesh()[1]/2))),
    _geom(geom), _param(param), _comm(comm),
    _initPosParticle(param.ParticleInitPos()), _numParticles(param.ParticleInitPos().size()) {

  // initialize the solver
  if(sol == nullptr || strncmp(sol, "CG", 2) == 0) {
    this->_solver = new CG(geom, comm);
  } else if(strncmp(sol, "MG", 2) == 0) {
    const multi_index_t& size = this->_geom.SizeP();
    const index_t level = (index_t)std::fmax(std::fmin( std::log2(size[0]), std::log2(size[1]))-1, 0);
    this->_solver = new MG(geom, comm, level);
  } else if(strncmp(sol, "RB", 2) == 0) {
    this->_solver = new RedOrBlackSOR(geom, param.Omega(), comm);
  } else {
    this->_solver = new SOR(geom, param.Omega(), comm);
  }

  // initialize u,v,p,particle
  this->_u->Initialize(0);
  this->_v->Initialize(0);
  this->_p->Initialize(0);
  this->_particle->Initialize(0);

  // initialize timestep restriction
  const multi_real_t &h = this->_geom.Mesh();
  this->_dtlimit = std::min(this->_param.Dt(),
      this->_param.Tau() * this->_param.Re() * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]));

#ifndef BLATT4
  if (this->_comm.ThreadCnt() == 1) {
    // write first particle pos in list for particle trace and streakline
    this->_particleTracing.insert(this->_particleTracing.begin(), this->_initPosParticle.begin(), this->_initPosParticle.end());
    this->_streakline.insert(this->_streakline.begin(), this->_initPosParticle.begin(), this->_initPosParticle.end());
  }
#else
  // output for exercise sheet 4
  const Iterator it1(this->_geom, multi_index_t(120, 5));
  const Iterator it2(this->_geom, multi_index_t(64, 64));
  const Iterator it3(this->_geom, multi_index_t(5, 120));
  std::cout << this->_param.Re() << "\t" << this->_t << "\t" << 0.0 << "\t"
            << this->_u->Cell(it1) << "\t" << this->_v->Cell(it1) << "\t"
            << this->_u->Cell(it2) << "\t" << this->_v->Cell(it2) << "\t"
            << this->_u->Cell(it3) << "\t" << this->_v->Cell(it3) << "\t"
            << std::endl;
#endif
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
  // set boundary values for u, v
  this->_geom.Update(*(this->_u), *(this->_v));

  // compute local dt
  // Test CFL and Pr condition
  const multi_real_t &h = this->_geom.Mesh();
  real_t dt = std::min(this->_dtlimit,
      this->_param.Tau() * std::min(h[0] / this->_u->AbsMax(), h[1] / this->_v->AbsMax()));
  // gather global dt
  dt = this->_comm.gatherMin(dt);

  // compute the momentum equations
  this->MomentumEqu(dt);

  // compute right-hand-side of the poisson equation
  this->RHS(dt);

  // set boundary values for p
  this->_geom.Update_P(*(this->_p));
  // solve the poisson equation for the pressure
  real_t res = 2 * this->_param.Eps() * this->_param.Eps();
  index_t i;

  this->_solver->reset(*(this->_p), *(this->_rhs));

  for(i = 0; (i < this->_param.IterMax()) && (res > this->_param.Eps() * this->_param.Eps()) ; i++) {
    res = this->_solver->Cycle(*(this->_p), *(this->_rhs));
    // set boundary values for p
    this->_geom.Update_P(*(this->_p));
  }
  if((i == this->_param.IterMax()) && printInfo) {
    std::cerr << "Warning: Compute: Solver did not converge! res = " << std::sqrt(res) << std::endl;
  }

  // compute new velocites
  this->NewVelocities(dt);

#ifndef BLATT4
  // calculate new particle postitions
  // TODO: rewrite for parallelization
  if (_comm.ThreadCnt() == 1) {
    this->Particle(dt);
    this->Streaklines(dt);
  }

  // print information
  if(printInfo) {
    std::cout.precision(4);
    std::cout << "t = " << this->_t << " \tdt = " << dt << "\titer = " << i
              << "\tres = " << std::scientific << std::sqrt(res)
              << std::defaultfloat << std::endl;
  }
#endif

  // update the time
  this->_t += dt;

#ifdef BLATT4
  // output for exercise sheet 4
  const Iterator it1(this->_geom, multi_index_t(120, 5));
  const Iterator it2(this->_geom, multi_index_t(64, 64));
  const Iterator it3(this->_geom, multi_index_t(5, 120));
  std::cout << this->_param.Re() << "\t" << this->_t << "\t" << std::sqrt(res) << "\t"
            << this->_u->Cell(it1) << "\t" << this->_v->Cell(it1) << "\t"
            << this->_u->Cell(it2) << "\t" << this->_v->Cell(it2) << "\t"
            << this->_u->Cell(it3) << "\t" << this->_v->Cell(it3) << "\t"
            << std::endl;
#endif
}

// Computes and returns the absolute velocity
const Grid *Compute::GetVelocity() {
  for(Iterator it(this->_geom); it.Valid(); it.Next()) {
    real_t uMean = (this->_u->Cell(it) + this->_u->Cell(it.Top()))/2;
    real_t vMean = (this->_v->Cell(it) + this->_v->Cell(it.Right()))/2;
    this->_velocities->Cell(it) = std::sqrt(uMean*uMean + vMean*vMean);
  }
  return _velocities;
}

// Computes and returns the vorticity
const Grid* Compute::GetVorticity() {
  for(Iterator it(this->_geom); it.Valid(); it.Next()) {
    this->_vorticity->Cell(it) = this->_u->dy_r(it) - this->_v->dx_r(it);
  }
  return this->_vorticity;
}
// Computes the stream line values
const Grid* Compute::GetStreamline() {
  this->_streamline->Data()[0] = 0;
  for(Iterator it(this->_geom, 1); it.Valid(); it.Next()) {
    if(it.Pos()[1] == 0) {
      this->_streamline->Cell(it) = this->_streamline->Cell(it.Left()) - _geom.Mesh()[0] * this->_v->Cell(it);
    } else if (this->_geom.noslip(it)) {
      // correct for separated domains because integral is invariant even across process boundaries, only true for pure noslip obstacles
      this->_streamline->Cell(it) = this->_streamline->Cell(it.Down());
    } else {
      this->_streamline->Cell(it) = this->_streamline->Cell(it.Down()) + _geom.Mesh()[1] * this->_u->Cell(it);
    }
  }
  real_t offset = this->_comm.copyOffset(*this->_streamline);
  for(Iterator it(this->_geom); it.Valid(); it.Next()) {
    this->_streamline->Cell(it) += offset;
  }
  return this->_streamline;
}

void Compute::Particle(const real_t &dt) {
  std::list<multi_real_t>::iterator itEnd = this->_particleTracing.end();
  std::list<multi_real_t> append;

  for(std::list<multi_real_t>::iterator it = std::prev(itEnd, this->_numParticles); it != itEnd; ++it) {
    multi_real_t temp = *it;
    this->ParticleStepVisu(temp, dt, this->_numParticles - std::distance(it, itEnd));
    append.push_back(temp);
  }

  // write new positions in the list
  this->_particleTracing.insert(this->_particleTracing.end(), append.begin(), append.end());
}

// TODO: exchange between domains
void Compute::Streaklines(const real_t &dt) {
  // Cycle the list of particles
  for (std::list<multi_real_t>::iterator it = this->_streakline.begin(); it != this->_streakline.end(); ++it) {
    this->ParticleStep((*it), dt);
    // if NAN, outside of domain
    if( std::isnan((*it)[0]) ) {
      it = this->_streakline.erase(it);
      --it;
    }
  }

  // Add new particles at the initial positions
  this->_streakline.insert(this->_streakline.begin(), this->_initPosParticle.begin(), this->_initPosParticle.end());
}

// updates particle position
void Compute::ParticleStep(multi_real_t &pos, const real_t &dt) {
  // calculate new position
  pos[0] += dt * this->_u->Interpolate(pos);
  pos[1] += dt * this->_v->Interpolate(pos);

  if (pos[0] >= this->_geom.Length()[0] || pos[1] >= this->_geom.Length()[1]) {
    pos[0] = std::nan("");
  }
}

// updates particle position with live visualisation (only valid for particle tracing)
// TODO: exchange between domains
void Compute::ParticleStepVisu(multi_real_t &lastPos, const real_t &dt, const index_t numPos) {
  // get cell of the particle
  _particleIndx[0] = (index_t)(lastPos[0]/_geom.Mesh()[0] + 1);
  _particleIndx[1] = (index_t)(lastPos[1]/_geom.Mesh()[1] + 1);

  // set old position to 0 for debug visualization
  this->_particle->Cell(Iterator(this->_geom, _particleIndx)) = 0;

  // calculate new position
  this->ParticleStep(lastPos, dt);

  // if NAN, outside domain, then reset
  if( std::isnan(lastPos[0]) ) {
    // deltes old particle path
    // _particleTracing.clear();
    lastPos = *std::next(this->_initPosParticle.begin(), numPos);
  }

  // get cell of the particle
  _particleIndx[0] = (index_t)(lastPos[0]/_geom.Mesh()[0] + 1);
  _particleIndx[1] = (index_t)(lastPos[1]/_geom.Mesh()[1] + 1);
  // set new position to 1 for debug visualization
  this->_particle->Cell(Iterator(this->_geom, _particleIndx)) = 1;
}

// Compute the new velocites u,v
void Compute::NewVelocities(const real_t &dt) {
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    this->_u->Cell(it) = this->_F->Cell(it) - dt * this->_p->dx_r(it);
    this->_v->Cell(it) = this->_G->Cell(it) - dt * this->_p->dy_r(it);
  }
}

// Compute the temporary velocities F,G
void Compute::MomentumEqu(const real_t &dt) {
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
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
  this->_geom.Update(*this->_F, *this->_G);
}

// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt) {
  for(InteriorIterator it(this->_geom); it.Valid(); it.Next()) {
    this->_rhs->Cell(it) = (this->_F->dx_l(it) + this->_G->dy_l(it)) / dt;
  }
}


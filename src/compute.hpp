/*
 * Copyright (C) 2015   Malte Brunn
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

#include <list>
#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
/// Class for all computations done within one timestep.
class Compute {
public:
  /// Creates a compute instance with given geometry, parameter and communicator
  Compute(const Geometry &geom, const Parameter &param, const Communicator &comm,
          const char sol[] = nullptr);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  /// \param printInfo print information about current solver state (residual etc.)
  void TimeStep(bool printInfo);

  /// Returns the simulated time in total
  inline const real_t &GetTime() const { return this->_t; };

  /// Returns the pointer to U
  inline const Grid *GetU() const { return this->_u; };
  /// Returns the pointer to V
  inline const Grid *GetV() const { return this->_v; };
  /// Returns the pointer to P
  inline const Grid *GetP() const { return this->_p; };
  /// Returns the pointer to RHS
  inline const Grid *GetRHS() const { return this->_rhs; };
  /// Returns the pointer to the Streamline
  const Grid *GetStreamline();
  /// Returns the pointer to the Vorticity
  const Grid *GetVorticity();
  /// Returns the pointer to the Particle Tracing grid
  inline const Grid *GetParticleTrace() const { return this->_particle; };
  /// Returns the pointer to the Particle Tracing list
  inline const std::list<multi_real_t> &GetParticles() const { return this->_particleTracing; };
  /// Returns the pointer to the Particle Tracing list
  inline const std::list<multi_real_t> &GetStreakline() const { return this->_streakline; };

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();

private:
  real_t _t; ///< current timestep
  real_t _dtlimit; ///< donor-cell diffusion condition (p. 27)

  Grid *_u; ///< velocities in x-direction
  Grid *_v; ///< velocities in y-direction
  Grid *_p; ///< pressure

  Grid *_F; /// prel. velocities in x-direction
  Grid *_G; /// prel. velocities in y-direction
  Grid *_rhs; ///< right-hand side of the pressure (poisson) equation

  Grid *_velocities; ///< grid for absolute velocities
  Grid *_streamline; ///< grid for streamfunction
  Grid *_vorticity; ///< grid for the vorticity
  Grid *_particle; ///< grid for the particle tracing

  Solver *_solver; ///< solver for the pressure (poisson) equation

  const Geometry &_geom; ///< geometry of the problem
  const Parameter &_param; ///< parameters for the computation
  const Communicator &_comm; ///< communicator for boundary exchange

  const std::list<multi_real_t> &_initPosParticle; ///< initial positions of the traced particles
  multi_index_t _particleIndx; ///< cell of the traced particle for visu
  std::list<multi_real_t> _particleTracing; ///< list for particle tracing
  std::list<multi_real_t> _streakline; ///< list of streakline points
  const index_t _numParticles; ///< number of traced particles/ streaklines

  /// Compute the new velocities u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocities F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);
  /// Computes the new particles
  void Particle(const real_t &dt);
  /// Computes the new streakline
  void Streaklines(const real_t &dt);
  /// Computes the new position of the particle
  void ParticleStep(multi_real_t &particlePos, const real_t &dt);
  /// Computes new position of the particle with live visualisation
  void ParticleStepVisu(multi_real_t &particlePos, const real_t &dt, const index_t numPos);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP

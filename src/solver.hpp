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
//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "comm.hpp"
//------------------------------------------------------------------------------
#ifndef __SOLVER_HPP
#define __SOLVER_HPP
//------------------------------------------------------------------------------

/// abstract base class for an iterative solver
class Solver {
public:
  /// Constructor of the abstract Solver class
  Solver(const Geometry &geom, const Communicator &comm)
      : _geom(geom), _invNumFluid(1.0/geom.NumFluid()), _comm(comm) {};
  /// Destructor of the Solver Class
  virtual ~Solver() {};

  /// This function must be implemented in a child class
  // @param [in][out] grid current values
  // @param [in]      rhs  right hand side values
  // @returns accumulated residual
  virtual real_t Cycle(Grid &grid, const Grid &rhs) const = 0;

  /// Resets the solver and prepares it for a new time step
  virtual void reset(const Grid &grid __attribute__((unused)), const Grid &rhs __attribute__((unused))) {};


protected:
  const Geometry &_geom;
  const real_t _invNumFluid; ///< The invers of the number of fluid cells
  const Communicator& _comm;

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid &grid, const Grid &rhs) const;
};

//------------------------------------------------------------------------------

/// concrete SOR solver
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry &geom, const Communicator &comm, const real_t &omega);

  /// Returns the total residual and executes a solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  virtual real_t Cycle(Grid &grid, const Grid &rhs) const;

protected:
  const real_t _correction;  ///< The correction factor computed from over-relaxation parameter
};
//------------------------------------------------------------------------------

/// concrete RedBlackSOR solver
class RedOrBlackSOR : public SOR {
public:
  /// Constructs an RedBlackSOR solver
  RedOrBlackSOR(const Geometry &geom, const Communicator &comm, const real_t &omega)
      : SOR(geom, comm, omega), _firstRed(comm.EvenOdd() && (geom.Size()[0] % 2 == 0))  // TODO one size even other odd
      {};

  /// Returns the total residual and executes a red solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

  /// Returns the total residual and executes a red solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  real_t RedCycle(Grid &grid, const Grid &rhs) const;
  /// Returns the total residual and executes a black solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  real_t BlackCycle(Grid &grid, const Grid &rhs) const;

protected:
  const bool _firstRed; ///< Whether first Red of Black Cycle is done
};


class CG : public Solver {
public:
  /// Constructs an Conjugated Gradients solver
  CG(const Geometry &geom, const Communicator &comm)
      : Solver(geom, comm), _res(geom), _direction(geom), _Ad(geom) {};
  /// Destructor
  ~CG() {};

  /// Returns the total residual and executes a red solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

  /// Resets the solver and prepares it for a new time step
  virtual void reset(const Grid &grid, const Grid &rhs);

protected:
  mutable Grid _res;
  mutable Grid _direction;
  mutable Grid _Ad;
  mutable real_t old_residual;
};
//------------------------------------------------------------------------------

/// concrete Multigrid solver
class MG : public Solver {
public:
  /// Constructs an Multigrid solver
  MG(const Geometry &geom, const Communicator &comm, const index_t level, const index_t& nu = index_t(4));
  /// Destructor
  ~MG();

  /// Returns the total residual and executes a solver cycle (V-Cycle)
  // @param grid current pressure values
  // @param rhs right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

protected:
  /// Restricts the residuals of the solution to the next coarser Grid
  /// \param p          The current pressure grid
  /// \param rhs        The current right-hand-side
  void Restrict(const Grid &p, const Grid &rhs) const;
  /// Interpolates and adds from the coarser solution to this one
  /// \param p          The current pressure grid
  void Interpolate(Grid &p) const;
  /// Smoothing cycle
  /// \param p     The pressure grid
  /// \param rhs   The right-hand-side
  real_t Smooth(Grid &p, const Grid &rhs) const;

  const index_t _level;
  const index_t _nu;
  const RedOrBlackSOR _smoother;
  const MG *_coarse;
  Grid *_e;
  Grid *_res;
};
//------------------------------------------------------------------------------
#endif // __SOLVER_HPP

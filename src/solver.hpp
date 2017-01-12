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
  Solver(const Geometry &geom) : _geom(geom) {};
  /// Destructor of the Solver Class
  virtual ~Solver() {};

  /// This function must be implemented in a child class
  // @param [in][out] grid current values
  // @param [in]      rhs  right hand side values
  // @returns accumulated residual
  virtual real_t Cycle(Grid &grid, const Grid &rhs) const = 0;

protected:
  const Geometry &_geom;

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid &grid, const Grid &rhs) const;
};

//------------------------------------------------------------------------------

/// concrete SOR solver
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry &geom, const real_t &omega);
  /// Destructor
  ~SOR() {};

  /// Returns the total residual and executes a solver cycle
  /// \param[in][out] grid current pressure values
  /// \param[in] rhs right hand side
  virtual real_t Cycle(Grid &grid, const Grid &rhs) const;

protected:
  const real_t _correction;  ///< The correction factor computed from over-relaxation parameter
  const real_t _invNumFluid; ///< The invers of the number of fluid cells
};
//------------------------------------------------------------------------------

/// concrete RedBlackSOR solver
class RedOrBlackSOR : public SOR {
public:
  /// Constructs an RedBlackSOR solver
  RedOrBlackSOR(const Geometry &geom, const real_t &omega, const Communicator &comm)
      : SOR(geom,omega), _comm(comm),
      _firstRed(comm.EvenOdd() && (geom.Size()[0] % 2 == 0))  // TODO one size even other odd
      {};
  /// Destructor
  ~RedOrBlackSOR() {};

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
  const Communicator& _comm;
  const bool _firstRed; ///< Whether first Red of Black Cycle is done
};
//------------------------------------------------------------------------------
#endif // __SOLVER_HPP

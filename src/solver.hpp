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

  /// \brief Returns the total residual and executes a solver cycle
  /// This function must be implemented in a child class
  /// \param [in,out] grid current values
  /// \param [in]     rhs  right hand side values
  /// \returns accumulated residual
  virtual real_t Cycle(Grid &grid, const Grid &rhs) const = 0;

  /// Resets the solver and prepares it for a new time step
  virtual void reset(const Grid &grid __attribute__((unused)), const Grid &rhs __attribute__((unused))) {};

protected:
  const Geometry &_geom;     ///< The geometry defining the domain
  const real_t _invNumFluid; ///< The invers of the number of fluid cells
  const Communicator& _comm; ///< The communicator for communication between MPI processes

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid &grid, const Grid &rhs) const;
};

//------------------------------------------------------------------------------

/// concrete SOR solver
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry &geom, const Communicator &comm, const real_t &omega);

  /// \brief Returns the total residual and executes a solver cycle
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
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

  /// \brief Returns the total residual and executes a red and black solver cycle
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

  /// \brief Returns the total residual and executes a red solver cycle
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
  real_t RedCycle(Grid &grid, const Grid &rhs) const;
  /// \brief Returns the total residual and executes a black solver cycle
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
  real_t BlackCycle(Grid &grid, const Grid &rhs) const;

protected:
  const bool _firstRed; ///< Whether first Red of Black Cycle is done
};

/// concrete Conjugated-Gradients (CG) solver
class CG : public Solver {
public:
  /// Constructs an Conjugated Gradients solver
  CG(const Geometry &geom, const Communicator &comm)
      : Solver(geom, comm), _res(geom), _direction(geom), _Ad(geom) {};

  /// \brief Returns the total residual and executes a solver cycle
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

  /// Resets the solver and prepares it for a new time step
  virtual void reset(const Grid &grid, const Grid &rhs);

protected:
  mutable Grid _res;           ///< The pointwise residual
  mutable Grid _direction;     ///< The direction of the correction
  mutable Grid _Ad;            ///< A*_direction where A is the discrete Laplassian
  mutable real_t old_residual; ///< The 2-norm of the old residual
};
//------------------------------------------------------------------------------

/// concrete Multigrid solver
class MG : public Solver {
public:
  /// Constructs an Multigrid solver
  /// \param[in] geom  The geometry defining the domain
  /// \param[in] comm  The communicator for communication between MPI processes
  /// \param[in] level The level of the multigrid solver (0 == coarsest grid)
  /// \param[in] gamma The type of MG-Cycle (1 == V, 2 == W)
  /// \param[in] nu    The number of pre- and postsmooting iterations done per Cycle
  MG(const Geometry &geom, const Communicator &comm, const index_t level,
     const index_t& gamma, const index_t& nu);
  /// Destructor
  ~MG();

  /// \brief Returns the total residual and executes a solver cycle (MG-Cycle)
  /// \param[in,out] grid current pressure values
  /// \param[in]     rhs  right hand side
  real_t Cycle(Grid &grid, const Grid &rhs) const;

protected:
  /// Restricts the residuals of the solution to the next coarser Grid
  /// \param[in] p   The current pressure grid
  /// \param[in] rhs The current right-hand-side
  void Restrict(const Grid &p, const Grid &rhs) const;
  /// Interpolates and adds from the coarser solution to this one
  /// \param[out] p The current pressure grid
  void Interpolate(Grid &p) const;
  /// Smoothing iteration
  /// \param[in,out] p   The pressure grid
  /// \param[in]     rhs The right-hand-side
  real_t Smooth(Grid &p, const Grid &rhs) const;

  const index_t _level;  ///< The level of the multigrid solver (0 == coarsest grid)
  const index_t &_gamma; ///< The type of MG-Cycle (1 == V, 2 == W)
  const index_t &_nu;    ///< The number of pre- and postsmooting iterations done per Cycle
  const RedOrBlackSOR _smoother; ///< The smoother
  const MG *_coarse;     ///< The MG solver on the next coarser level
  Grid *_e;              ///< The error approximation on the next coarser level
  Grid *_res;            ///< The restriction of the residual on the next coarser level
};
//------------------------------------------------------------------------------
#endif // __SOLVER_HPP

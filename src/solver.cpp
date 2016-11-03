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
SOR::SOR(const Geometry &geom, const real_t &omega) : Solver(geom), _omega(omega) {}

// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid &grid, const Grid &rhs) const {
  const multi_real_t &h = this->_geom.Mesh();
  const real_t prefactor = this->_omega * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]);

  real_t residuum = 0;
  for(InteriorIterator it(grid); it.Valid(); it.Next()) {
    real_t localRes = Solver::localRes(it, grid, rhs);
    residuum += localRes * localRes;
    grid.Cell(it) = grid.Cell(it) - prefactor * localRes;
  }
  return std::sqrt(residuum);
}


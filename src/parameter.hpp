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

#include "typedef.hpp"
#include <list>
//------------------------------------------------------------------------------
#ifndef __PARAMETER_HPP
#define __PARAMETER_HPP
//------------------------------------------------------------------------------
/// Class for the parameters used in the simulation
class Parameter {
public:
  /// Constructs a new Parameter set with default values
  // Driven Cavity parameters; see exercise sheet 1
  Parameter(const char file[] = nullptr, const bool printinfo = 0)
      : _re(1000), _omega(1.7), _gamma(1), _nu(4), _alpha(0.9), _dt(1e-2), _tend(50),
        _eps(1e-3), _tau(0.5), _itermax(500), _vtkDt(0.5), _visuDt(0) {
      if(file != nullptr) {
        this->Load(file, printinfo);
      }
  }

  /// Loads the parameter values from a file
  void Load(const char file[], const bool printinfo = 0);

  /// Returns the Reynolds-number
  inline const real_t &Re() const { return this->_re; };
  /// Returns the over-relaxation parameter for SOR and RedBlackSOR
  inline const real_t &Omega() const { return this->_omega; };
  /// Returns the over-relaxation parameter for SOR and RedBlackSOR (write access)
  inline real_t &Omega() { return this->_omega; };
  /// Returns the type of MGCycle for the multigrid solver MG
  inline const index_t &Gamma() const { return this->_gamma; };
  /// Returns the number of smoother iterations for the multigrid solver MG
  inline const index_t &Nu() const { return this->_nu; };
  /// Returns the weighting between Donor-Cell and central differences
  inline const real_t &Alpha() const { return this->_alpha; };
  /// Returns the maximal time-step
  inline const real_t &Dt() const { return this->_dt; };
  /// Returns the final time
  inline const real_t &Tend() const { return this->_tend; };
  /// Returns the maximal number of iterations done by the iterative solver
  inline const index_t &IterMax() const { return this->_itermax; };
  /// Returns the maximal error allowed in the iterative solver
  inline const real_t &Eps() const { return this->_eps; };
  /// Returns the savety-factor for the time-step restriction
  inline const real_t &Tau() const { return this->_tau; };
  /// Returns the time-step after which a vtk-output is written
  inline const real_t &VtkDt() const { return this->_vtkDt; };
  /// Returns the time-step after which the visualization is updated
  inline const real_t &VisuDt() const { return this->_visuDt; };
  /// Returns the initial positions of the particles
  inline const std::list<multi_real_t> &ParticleInitPos() const { return this->_particleInitPos; };

private:
  real_t _re;
  real_t _omega;
  index_t _gamma;
  index_t _nu;
  real_t _alpha;
  real_t _dt;
  real_t _tend;
  real_t _eps;
  real_t _tau;
  index_t _itermax;
  real_t _vtkDt;
  real_t _visuDt;
  std::list<multi_real_t> _particleInitPos;
};
//------------------------------------------------------------------------------
#endif // __PARAMETER_HPP

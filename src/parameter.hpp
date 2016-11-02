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
//------------------------------------------------------------------------------
#ifndef __PARAMETER_HPP
#define __PARAMETER_HPP
//------------------------------------------------------------------------------
class Parameter {
public:
  /// Constructs a new Parameter set with default values
  // Driven Cavity parameters; see exercise sheet 1
  Parameter() : _re(1), _omega(1), _alpha(0.5), _dt(1e-2), _tend(1),
                _eps(1e-3), _tau(0.8), _itermax(100) {};

  /// Loads the parameter values from a file
  void Load(const char file[]);

  /// Getter functions for all parameters
  inline const real_t &Re() const { return this->_re; };
  inline const real_t &Omega() const { return this->_omega; };
  inline const real_t &Alpha() const { return this->_alpha; };
  inline const real_t &Dt() const { return this->_dt; };
  inline const real_t &Tend() const { return this->_tend; };
  inline const index_t &IterMax() const { return this->_itermax; };
  inline const real_t &Eps() const { return this->_eps; };
  inline const real_t &Tau() const { return this->_tau; };

private:
  real_t _re;
  real_t _omega;
  real_t _alpha;
  real_t _dt;
  real_t _tend;
  real_t _eps;
  real_t _tau;
  index_t _itermax;
};
//------------------------------------------------------------------------------
#endif // __PARAMETER_HPP

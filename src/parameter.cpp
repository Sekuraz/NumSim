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
#include <fstream>
#include <string>
#include "typedef.hpp"
#include "parameter.hpp"

using namespace std;

// Loads the parameter values from a file
void Parameter::Load(const char file[]) {
  string param;
  real_t value;
  ifstream in(file);
  while(in.good()) {
    in >> ws >> param >> ws;
    in.ignore(); // removes '=' between parameter and value
    in >> ws >> value >> ws;
    if(!param.compare("Re") || !param.compare("re")) {
      this->_re = value;
      cout << "Parameter: Load Re = " << this->_re << endl;
    } else if(!param.compare("Omega") || !param.compare("omega")) {
      this->_omega = value;
      cout << "Parameter: Load omega = " << this->_omega << endl;
    } else if(!param.compare("Alpha") || !param.compare("alpha")) {
      this->_alpha = value;
      cout << "Parameter: Load alpha = " << this->_alpha << endl;
    } else if(!param.compare("Dt") || !param.compare("dt")) {
      this->_dt = value;
      cout << "Parameter: Load dt = " << this->_dt << endl;
    } else if(!param.compare("Tend") || !param.compare("tend") || !param.compare("tEnd")) {
      this->_tend = value;
      cout << "Parameter: Load tEnd = " << this->_tend << endl;
    } else if(!param.compare("Eps") || !param.compare("eps")) {
      this->_eps = value;
      cout << "Parameter: Load eps = " << this->_eps << endl;
    } else if(!param.compare("Tau") || !param.compare("tau")) {
      this->_tau = value;
      cout << "Parameter: Load tau = " << this->_tau << endl;
    } else if(!param.compare("IterMax") || !param.compare("Itermax") || !param.compare("itermax")) {
      this->_itermax = value;
      cout << "Parameter: Load IterMax = " << this->_itermax << endl;
    } else if(!param.compare("vtkdt") || !param.compare("vtkDt") || !param.compare("VTKdt")) {
      this->_vtkDt = value;
      cout << "Parameter: Load vtkDt = " << this->_vtkDt << endl;
    } else {
      cerr << "Parameter: Unknown identifier " << param << endl;
    }
  }
  in.close();
}


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

// Loads the parameter values from a file
void Parameter::Load(const char file[], const bool printinfo) {
  std::string param;
  real_t value;
  std::ifstream in(file);
  if(printinfo) {
    std::cout << "Parameter: Loading from file " << file << std::endl;
  }
  while(in.good()) {
    in >> param >> std::ws;
    if(!std::isdigit(in.peek())) {
      in.ignore(); // removes separator between parameter and value
    }
    in >> value;
    if(!param.compare("Re") || !param.compare("re")) {
      this->_re = value;
      if(printinfo) {
        std::cout << "Parameter: Load Re = " << this->_re << std::endl;
      }
    } else if(!param.compare("Omega") || !param.compare("omega")) {
      this->_omega = value;
      if(printinfo) {
        std::cout << "Parameter: Load omega = " << this->_omega << std::endl;
      }
    } else if(!param.compare("Alpha") || !param.compare("alpha")) {
      this->_alpha = value;
      if(printinfo) {
        std::cout << "Parameter: Load alpha = " << this->_alpha << std::endl;
      }
    } else if(!param.compare("Dt") || !param.compare("dt")) {
      this->_dt = value;
      if(printinfo) {
        std::cout << "Parameter: Load dt = " << this->_dt << std::endl;
      }
    } else if(!param.compare("Tend") || !param.compare("tend") || !param.compare("tEnd")) {
      this->_tend = value;
      if(printinfo) {
        std::cout << "Parameter: Load tEnd = " << this->_tend << std::endl;
      }
    } else if(!param.compare("Eps") || !param.compare("eps")) {
      this->_eps = value;
      if(printinfo) {
        std::cout << "Parameter: Load eps = " << this->_eps << std::endl;
      }
    } else if(!param.compare("Tau") || !param.compare("tau")) {
      this->_tau = value;
      if(printinfo) {
        std::cout << "Parameter: Load tau = " << this->_tau << std::endl;
      }
    } else if(!param.compare("IterMax") || !param.compare("Itermax") || !param.compare("itermax")) {
      this->_itermax = value;
      if(printinfo) {
        std::cout << "Parameter: Load IterMax = " << this->_itermax << std::endl;
      }
    } else if(!param.compare("vtkdt") || !param.compare("vtkDt") || !param.compare("VTKdt")) {
      this->_vtkDt = value;
      if(printinfo) {
        std::cout << "Parameter: Load vtkDt = " << this->_vtkDt << std::endl;
      }
    } else if(!param.compare("visudt") || !param.compare("visuDt") || !param.compare("VISUDt")) {
      this->_visuDt = value;
      if(printinfo) {
        std::cout << "Parameter: Load visuDt = " << this->_visuDt << std::endl;
      }
    } else if(!param.compare("particleposx") || !param.compare("particlePosX") || !param.compare("PARTICLEPOSX")) {
      this->_particleInitPos[0] = value;
      if(printinfo) {
        std::cout << "Parameter: Load particlePosX = " << this->_particleInitPos[0] << std::endl;
      }
    } else if(!param.compare("particleposy") || !param.compare("particlePosY") || !param.compare("PARTICLEPOSY")) {
      this->_particleInitPos[1] = value;
      if(printinfo) {
        std::cout << "Parameter: Load particlePosY = " << this->_particleInitPos[1] << std::endl;
      }
    } else {
      std::cerr << "Parameter: Unknown identifier " << param << std::endl;
    }
    in >> std::ws;
  }
  in.close();
}


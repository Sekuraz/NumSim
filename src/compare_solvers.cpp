/*
 *  Copyright (C) 2016   Stephan Lunowa, Markus Baur, Jonas Harsch
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "solver.hpp"

int main(int argc, char *argv[]) {
  // set epsilon^2
  const real_t epsilon2 = 1e-12;
  // set maximal Grid size (2^k)
  const index_t maxSize = 128;
  // set number of repetitions
  const index_t nRepeat = 1;

  // Create normal distribution random generator
  const index_t seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<real_t> distribution(0,1);

  Communicator comm(&argc, &argv);
  const bool printinfo = comm.ThreadNum() == 0;

  // print header
  if(printinfo) {
    std::cout << "Size\tSolver\tIterations\tResidual" << std::endl;
  }

  for(index_t n = 4; n <= maxSize; n *= 2) {
    Geometry geom(comm, {n+2,n+2});
    // Set optimal omega for SOR and RBSOR
    const real_t omega = 2.0 / (1.0 + std::sin(M_PI * std::fmax(geom.Mesh()[0], geom.Mesh()[1])));
    // Set optimal level for MG
    const index_t level = (index_t)std::fmax(std::log2(n), 0);

    MG mg(geom, comm, level, 2);
    CG cg(geom, comm);
    RedOrBlackSOR rb(geom, comm, omega);
    SOR sor(geom, comm, omega);

    Grid rhs(geom);
    rhs.Initialize(0);
    Grid pMG(geom), pCG(geom), pRB(geom), pSOR(geom);
    
    for(index_t rep = 0; rep < nRepeat; rep++) {
      for(InteriorIterator it(geom); it.Valid(); it.Next()) {
        const real_t val = distribution(generator);
        pMG.Cell(it) = val;
        pCG.Cell(it) = val;
        pRB.Cell(it) = val;
        pSOR.Cell(it) = val;
      }

      // set boundary values for p
      geom.Update_P(pMG);
      // solve the poisson equation for the pressure
      index_t i;
      real_t res = 2 * epsilon2;
      mg.reset(pMG, rhs);
      for(i = 0; res > epsilon2 ; i++) {
        res = mg.Cycle(pMG, rhs);
        // set boundary values for p
        geom.Update_P(pMG);
      }
      // print information
      if(printinfo) {
        std::cout << n << "\tMG\t" << i << "\t" << std::scientific << std::sqrt(res) << std::endl;
      }

      // set boundary values for p
      geom.Update_P(pCG);
      // solve the poisson equation for the pressure
      res = 2 * epsilon2;
      cg.reset(pCG, rhs);
      for(i = 0; res > epsilon2 ; i++) {
        //if(i % 100 == 0) cg.reset(pCG, rhs);
        res = cg.Cycle(pCG, rhs);
        // set boundary values for p
        geom.Update_P(pCG);
      }
      // print information
      if(printinfo) {
        std::cout << n << "\tCG\t" << i << "\t" << std::scientific << std::sqrt(res) << std::endl;
      }

      // set boundary values for p
      geom.Update_P(pRB);
      // solve the poisson equation for the pressure
      res = 2 * epsilon2;
      rb.reset(pRB, rhs);
      for(i = 0; res > epsilon2 ; i++) {
        res = rb.Cycle(pRB, rhs);
        // set boundary values for p
        geom.Update_P(pRB);
      }
      // print information
      if(printinfo) {
        std::cout << n << "\tRB\t" << i << "\t" << std::scientific << std::sqrt(res) << std::endl;
      }

      // set boundary values for p
      geom.Update_P(pSOR);
      // solve the poisson equation for the pressure
      res = 2 * epsilon2;
      sor.reset(pSOR, rhs);
      for(i = 0; res > epsilon2 ; i++) {
        res = sor.Cycle(pSOR, rhs);
        // set boundary values for p
        geom.Update_P(pSOR);
      }
      // print information
      if(printinfo) {
        std::cout << n << "\tSOR\t" << i << "\t" << std::scientific << std::sqrt(res) << std::endl;
      }
    }
  }

  return 0;
}

/*
 *  Copyright (C) 2016   Malte Brunn, Stephan Lunowa, Markus Baur
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
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#ifdef USE_DEBUG_VISU
  #include "visu.hpp"
#endif
#include "grid.hpp"
#include "solver.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Parameter param;
  param.Load("param.txt");
  Geometry geom(comm, {9,9});

  const multi_real_t &h = geom.Mesh();
  Grid g (geom, multi_real_t(h[0]/2, h[1]/2));
  Grid rhs (geom);
  g.Initialize(0);
  rhs.Initialize(0);

  RedOrBlackSOR rbsor(geom, param.Omega());

#ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  const int hRes = 600;
  const int vRes = (int)(hRes / geom.TotalLength()[0]);
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(hRes, vRes, comm.ThreadNum(), comm.ThreadIdx(), comm.ThreadDim());
#endif

  // insert pertubation
  if(comm.ThreadNum() == 0) {
    multi_index_t half = {geom.SizeP()[0]/2, geom.SizeP()[1]/2};
    g.Cell(Iterator(geom, half)) = 1000;
  }

  geom.Update_P(g);
  bool run = true;
  for (int i = 0; run; i++) {
    real_t res = 0;
    if(comm.EvenOdd() && (geom.Size()[0] % 2 == 0)) {
      res = rbsor.RedCycle(g, rhs);
    } else {
      res = rbsor.BlackCycle(g, rhs);
    }
    comm.copyBoundary(g);
    if(comm.EvenOdd() && (geom.Size()[0] % 2 == 0)) {
      res += rbsor.BlackCycle(g, rhs);
    } else {
      res += rbsor.RedCycle(g, rhs);
    }
    geom.Update_P(g);

    // compute global min and max
    real_t min, max;
    g.MinMax(min, max);
    min = comm.gatherMin(min);
    max = comm.gatherMax(max);
    res = comm.gatherSum(res);

#ifdef USE_DEBUG_VISU
    visu.Render(&g, min, max);
#endif
    if(comm.ThreadNum() == 0) {
      std::cout << "step = " << i << "\tmax-min = " << (max-min)
                << "\tres = " << std::sqrt(res) << std::endl;
    }

    /*for (Iterator it(g); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", g.Cell(it) * 1);
    }
    std::cout << std::endl;*/
    if(std::cin.get() == 'q') run = false;
    run = comm.gatherAnd(run);
  }

  return 0;
}

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
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#ifdef DEBUG_VISU
  #include "visu.hpp"
#endif
#include "grid.hpp"
#include "solver.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Parameter param;
  param.Load("param.txt");
  Geometry geom(comm, {16,16});

  Grid g (geom);
  Grid rhs (geom);
  g.Initialize(0);
  rhs.Initialize(0);

  SOR sor(geom, param.Omega());

#ifdef DEBUG_VISU
  // Create and initialize the visualization
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(600,600);//(800, 800);
#endif

  multi_index_t half = {geom.Size()[0]/2, geom.Size()[1]/2};
  g.Cell(Iterator(geom, half)) = 1000;

  bool run = true;
  for (index_t i = 0; run; i++) {
    std::cout << "step = " << i << " " << std::endl;

    geom.Update_P(g);

    sor.Cycle(g, rhs);
#ifdef DEBUG_VISU
    visu.Render(&g);
#endif
    for (Iterator it(geom); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", g.Cell(it) * 1);
    }
    std::cout << std::endl;

    if(std::cin.get() == 'q') run = false;
    run = comm.gatherAnd(run);
  }
  return 0;
}

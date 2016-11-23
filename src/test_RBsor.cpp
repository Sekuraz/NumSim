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
#include "compute.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"
#include "grid.hpp"
#include "solver.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Parameter param;
  param.Load("param.txt");
  Geometry geom(comm, {16,16});

  Grid g (geom,Grid::type::p);
  Grid rhs (geom,Grid::type::p);
  g.Initialize(0);
  rhs.Initialize(0);

  RedOrBlackSOR rbsor(geom, param.Omega());
  // Create and initialize the visualization
  const int hRes = 600;
  const int vRes = (int)(hRes / geom.TotalLength()[0]);
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(hRes / comm.ThreadDim()[0], vRes / comm.ThreadDim()[1], comm.ThreadNum());
  
  // insert pertubation
  multi_index_t half = {geom.Size()[0]/2, geom.Size()[1]/2};
  g.Cell(Iterator(g, half)) = 1000;

  for (int i = 0; true; i++) {
    std::cout << "step = " << i << " " << std::endl;

    geom.Update_P(g);
    rbsor.RedCycle(g, rhs);
    geom.Update_P(g);
    rbsor.BlackCycle(g,rhs);
    visu.Render(&g);

    /*for (Iterator it(g); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", g.Cell(it) * 1);
    }
    std::cout << std::endl;*/
    std::cin.get();

  }
  return 0;
}

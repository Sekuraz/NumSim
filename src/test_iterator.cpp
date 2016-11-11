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
#include "compute.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"
#include "grid.hpp"
#include "solver.hpp"

int main(int argc, char *argv[]) {
  // Create parameter and geometry instances and load values
  Parameter param;
  param.Load("param.txt");
  Geometry geom;

  const multi_real_t &h = geom.Mesh();
  multi_real_t offset(0, h[1]/2);
  Grid g (geom,Grid::type::u, offset);
  g.Initialize(0);

  // Create and initialize the visualization
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(600,600);//(800, 800);

  for (InteriorIterator it (g); it.Valid(); it.Next()) {

    g.Cell(it) = 100;
    visu.Render(&g);

    for (Iterator it(g); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", g.Cell(it) * 1);
    }
    std::cout << std::endl;
    std::cin.get();

    g.Cell(it) = 0;
  }

/*
    BoundaryIterator it (g);
    for (int i = 0; i < 4; i++) {
        it.SetBoundary(i % 4 + 1);
        for(; it.Valid(); it.Next()) {

            g.Cell(it) = 100;
            visu.Render(&g);

            for (Iterator it(g); it.Valid(); it.Next()) {
                if (it.Pos()[0] == 0) {
                    std::cout << std::endl;
                }
                printf("%+.2f, ", g.Cell(it) * 1);
            }
            std::cout << std::endl;
            std::cin.get();

            g.Cell(it) = 0;
        }
    }
*/
  return 0;
}

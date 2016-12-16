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
#include "grid.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Parameter param;
  param.Load("param.txt");
  Geometry geom (comm, {4, 4});

  const multi_real_t &h = geom.Mesh();
  multi_real_t offset(h[0]/2, 0);
  Grid v (geom, offset);
  offset = {0, h[1]/2};
  Grid u (geom, offset);

  for (Iterator it(geom); it.Valid(); it.Next()) {
      v.Cell(it) = it.Value();
  }
  for (Iterator it(geom); it.Valid(); it.Next()) {
      u.Cell(it) = it.Value();
  }

  for (InteriorIterator it (geom); it.Valid(); it.Next()) {

    u.DC_vdu_y(it, 1.0, &v);

    for (Iterator it(geom); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", u.Cell(it) * 1);
    }
    std::cout << std::endl;
    for (Iterator it(geom); it.Valid(); it.Next()) {
        if (it.Pos()[0] == 0) {
            std::cout << std::endl;
        }
        printf("%+.2f, ", v.Cell(it) * 1);
    }
    std::cout << std::endl;
    std::cin.get();

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

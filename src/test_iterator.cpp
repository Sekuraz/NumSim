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
#include <omp.h>
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#ifdef USE_DEBUG_VISU
  #include "visu.hpp"
#endif
#include "grid.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Parameter param;
  Geometry geom (comm, {6, 6});

  const multi_real_t &h = geom.Mesh();
  multi_real_t offset(h[0]/2, h[1]/2);
  Grid g (geom, offset);
  g.Initialize(0);

#ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(600,600);
#endif

  #pragma omp parallel
  { for (InteriorIterator it (geom); it.Valid(); it.Next()) {
    g.Cell(it) = it;
//    geom.Update_P(g);
#ifdef USE_DEBUG_VISU
    visu.Render(&g);
#endif
//    g.print();
//    std::cin.get();

//    g.Cell(it) = 0;
  } }

  if(omp_get_thread_num() == 0) g.print();
/*
    BoundaryIterator it (geom);
    for (int i = 0; i < 4; i++) {
        it.SetBoundary(i % 4 + 1);
        for(; it.Valid(); it.Next()) {

            g.Cell(it) = 100;
            visu.Render(&g);

            g.print();
            std::cin.get();

            g.Cell(it) = 0;
        }
    }
*/
  return 0;
}

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

#include <mpi.h>
#include <iostream>
#include "typedef.hpp"
#include "comm.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"

int main(int argc, char *argv[]) {
  Communicator comm(&argc, &argv);
  // Create parameter and geometry instances and load values
  Geometry geom (comm, {2, 2});

  const multi_real_t &h = geom.Mesh();
  multi_real_t offset(h[0]/2, h[1]/2);
  Grid g (geom, offset);
  g.Initialize(0);
  for (Iterator it(g); it.Valid(); it.Next()) {
    g.Cell(it) = it * 0.1 + comm.ThreadNum();
  }

  comm.copyBoundary(g);

  for (int i = 0; i < comm.ThreadCnt(); i++) {
    if (i == comm.ThreadNum()) {
      for (Iterator it(g); it.Valid(); it.Next()) {
          if (it.Pos()[0] == 0) {
              std::cout << std::endl;
          }
          printf("%+.2f, ", g.Cell(it) * 1);
      }
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return 0;
}

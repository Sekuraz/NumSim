/*
 *  Copyright (C) 2016   Malte Brunn, Stephan Lunowa, Markus Baur, Jonas Harsch
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

int main(int argc, char *argv[]) {
  // Delete old VTK files, prevents bugs in paraview
  system("exec rm -r ./VTK/*");
  // Create communicator, parameter and geometry instances and load values
  Communicator comm(&argc, &argv);
  Parameter param;
  param.Load("param.txt", (comm.ThreadNum() == 0));
  Geometry geom(comm);
  geom.Load("geometry.txt");

  // TODO: for testing
  if(comm.ThreadNum() == 0) {
    std::cout << "Grids: " << geom.Size() << ", u " << geom.SizeU() << ", v " << geom.SizeV()
              << ", p " << geom.SizeP() << std::endl;
  }

  // Create the fluid solver
  Compute comp(geom, param, comm);

#ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  const int hRes = 600;
  const int vRes = (int)(hRes / geom.TotalLength()[0]);
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(hRes / comm.ThreadDim()[0], vRes / comm.ThreadDim()[1], comm.ThreadNum());
#endif // USE_DEBUG_VISU

  // Create a VTK generator
  VTK vtk(geom.Mesh(), geom.Size());

  const Grid *visugrid;
  bool run = true;

  visugrid = comp.GetVelocity();

  real_t nextTimeVTK = 0;
  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {
#ifdef USE_DEBUG_VISU
    // compute global min and max
    real_t min, max;
    visugrid->MinMax(min, max);
    min = comm.gatherMin(min);
    max = comm.gatherMax(max);
    // Render and check if window is closed
    switch (visu.Render(visugrid, min, max)) {
    case -1:
      run = false;
      break;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    default:
      break;
    };
#endif // USE_DEBUG_VISU

    if(comp.GetTime() >= nextTimeVTK ) {
      // Create a VTK File in the folder VTK (must exist)
      vtk.Init("VTK/field");
      // TODO: get all data from the processes
      vtk.AddField("Velocity", comp.GetU(), comp.GetV());
      vtk.AddScalar("Pressure", comp.GetP());
      vtk.AddScalar("x-Velocity", comp.GetU());
      vtk.AddScalar("y-Velocity", comp.GetV());
      vtk.AddScalar("Vorticity", comp.GetVorticity());
      vtk.AddScalar("Streamlines", comp.GetStreamline());
      vtk.Finish();
      nextTimeVTK += param.VtkDt();
    }

#ifdef USE_DEBUG_VISU
    comp.TimeStep(comm.ThreadNum() == 0);
#else
    comp.TimeStep(false);
#endif

#ifdef USE_DEBUG_VISU
    // Gather if one process stopped
    run = comm.gatherAnd(run);
#endif
  }
  return 0;
}

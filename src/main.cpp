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
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "typedef.hpp"
#include "comm.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#ifdef USE_DEBUG_VISU
  #include "visu.hpp"
#endif
#include "vtk.hpp"

/// reads arguments from command line
void parseCommandLine(const int argc, char *argv[], char* &paramPath, char* &geomPath) {
  for(int i = 1; i < argc; i++) {
    if(strncmp(argv[i], "-G",2)==0) { // geometry path
      // test whether path at flag or afterwards
      if (strlen(argv[i])>2) {
        geomPath = &argv[i][2];
      } else {
        geomPath = argv[++i];
      }
    } else if(strncmp(argv[i], "-g",2)==0) { // geometry standard
      // test whether number at flag or afterwards
      int number;
      if (strlen(argv[i])>2) {
        number = atoi(&argv[i][2]);
      } else {
        number = atoi(argv[++i]);
      }
      geomPath = new char[20];
      switch(number) {
        case 1:
          strcpy(geomPath, "data/pressure.geom");
          break;
        case 2:
          strcpy(geomPath, "data/step.geom");
          break;
        case 3:
          strcpy(geomPath, "data/karman.geom");
          break;
        default:
          strcpy(geomPath, "data/channel.geom");
          break;
      }
    } else if(strncmp(argv[i], "-P",2)==0) { // parameter path
      // test whether path at flag or afterwards
      if (strlen(argv[i])>2) {
        paramPath = &argv[i][2];
      } else {
        paramPath = argv[++i];
      }
    }
  }
}

int main(int argc, char *argv[]) {
  char *paramPath = nullptr, *geomPath = nullptr;

  // Create communicator, parameter and geometry instances and load values
  Communicator comm(&argc, &argv);

  // Delete old VTK files, prevents bugs in paraview
  if(comm.ThreadNum() == 0) {
    int status = system("exec rm -r ./VTK/*");
    if (status < 0) {
      std::cerr << "Error while deleting old VTK files!" << std::endl;
    }
  }

  parseCommandLine(argc, argv, paramPath, geomPath);

  Parameter param;
  param.Load((paramPath != nullptr? paramPath : "param.txt"), (comm.ThreadNum() == 0));
  Geometry geom(comm);
  geom.Load((geomPath != nullptr? geomPath : "geometry.txt"), (comm.ThreadNum() == 0));
  if(param.Omega() <= 0.0 || param.Omega() > 2.0) {
    real_t h = std::fmax(geom.Mesh()[0], geom.Mesh()[1]);
    param.Omega() = 2.0 / (1.0 + std::sin(M_PI * h));
    std::cout << "Set new omega = " << param.Omega() << std::endl;
  }

  // Create the fluid solver
  Compute comp(geom, param, comm);

#ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  const int hRes = 600;
  const int vRes = (int)(hRes * geom.TotalLength()[1] / geom.TotalLength()[0]);
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(hRes, vRes, comm.ThreadNum(), comm.ThreadIdx(), comm.ThreadDim());
#endif // USE_DEBUG_VISU

  // Create a VTK generator
  multi_real_t offset;
  for(index_t dim = 0; dim < DIM; dim++) {
    offset[dim] = comm.ThreadIdx()[dim] * geom.Size()[dim] * geom.Mesh()[dim];
  }
  VTK vtk(geom.Mesh(), geom.SizeP(), offset, geom.TotalSize(), comm.ThreadNum(),
          comm.ThreadIdx(), comm.ThreadDim());

  const Grid *visugrid;
  bool run = true;

  visugrid = comp.GetVelocity();

  real_t nextTimeVTK = 0;
  real_t nextTimeVisu = 0;
  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {

#ifdef USE_DEBUG_VISU
    if(comp.GetTime() >= nextTimeVisu) {
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
      case 4:
        visugrid = comp.GetVorticity();
        break;
      case 5:
        visugrid = comp.GetStreamline();
        break;
      case 6:
        visugrid = comp.GetParticleTrace();
        break;
      default:
        break;
      };
    nextTimeVisu += param.VisuDt();
    }
#endif // USE_DEBUG_VISU

    if(comp.GetTime() >= nextTimeVTK ) {
      // Create a VTK File in the folder VTK (must exist)
      vtk.Init("VTK/field");
      vtk.AddField("Velocity", comp.GetU(), comp.GetV());
      vtk.AddScalar("Pressure", comp.GetP());
      vtk.AddScalar("Vorticity", comp.GetVorticity());
      vtk.AddScalar("Streamlines", comp.GetStreamline());
      vtk.Finish();

      // Create VTK File for Particle tracing in the folder VTK (must exist)
      vtk.InitParticles("VTK/particle");
      vtk.AddParticles("Particle", comp.GetParticles());
      vtk.FinishParticles();

      // Create VTK File for Streaklines in the folder VTK (must exist)
      vtk.InitParticles("VTK/streakline");
      vtk.AddParticles("Streakline", comp.GetStreakline());
      vtk.FinishParticles();

      nextTimeVTK += param.VtkDt();
    }

/*
    std::cout << std::endl << "U: " << comp.GetU()->Size();
    comp.GetU()->print();
    std::cout << std::endl << "V: " << comp.GetV()->Size();
    comp.GetV()->print();
    std::cout << std::endl << "P: " << comp.GetP()->Size();
    comp.GetP()->print();
    std::cin.get();
*/
    comp.TimeStep(comm.ThreadNum() == 0);
#ifdef USE_DEBUG_VISU
    // Gather if one process stopped
    run = comm.gatherAnd(run);
#endif
  }

  return 0;
}

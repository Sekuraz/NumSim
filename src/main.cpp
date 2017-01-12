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
#ifdef DEBUG_VISU
  #include "visu.hpp"
#endif
#ifdef WRITE_VTK
  #include "vtk.hpp"
#endif

/// reads arguments from command line
void parseCommandLine(const int argc, char *argv[], char* &paramPath, char* &geomPath) {
  for(int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "-h", 2) == 0 or strncmp(argv[i], "--help", 6) == 0) {
      std::cout << std::endl << "Numsim - Fluid simulation." << std::endl << std::endl
         << "Usage: numsim [arguments]\tSolves the Navier-Stokes equation for the given geometry and parameters." << std::endl << std::endl
         << "General arguments:" << std::endl
         << "Arguments must be seperated by a space from their additional parameters." << std::endl
         << "Starting without any argument solves a 128x128 driven cavity." << std::endl
         << "\t-h, --help\tPrints this message." << std::endl
         << "\t-g <num>\tChoose among some ready-made geometries." << std::endl
         << "\t\t\t\t0: Driven Cavity" << std::endl
         << "\t\t\t\t1: Simple channel" << std::endl
         << "\t\t\t\t2: Pressure driven channel" << std::endl
         << "\t\t\t\t3: Channel with step" << std::endl
         << "\t\t\t\t4: Channel with barrier" << std::endl
         << "\t-G <path>\tUse the specified geometry file" << std::endl
         << "\t-P <path>\tUse the specified parameter file" << std::endl
         << "\t-I <path>\tUse the specified parameter and geometry file" << std::endl
         << "\t\t\tEquivalent to -G <path>.geom and -P <path>.param" << std::endl
         << std::endl;
      exit(0);
    } else if(strncmp(argv[i], "-G",2)==0) { // geometry path
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
      paramPath = new char[20];
      switch(number) {
        case 1:
          strcpy(geomPath, "data/channel.geom");
          strcpy(paramPath, "data/channel.param");
          break;
        case 2:
          strcpy(geomPath, "data/pressure.geom");
          strcpy(paramPath, "data/pressure.param");
          break;
        case 3:
          strcpy(geomPath, "data/step.geom");
          strcpy(paramPath, "data/step.param");
          break;
        case 4:
          strcpy(geomPath, "data/karman.geom");
          strcpy(paramPath, "data/karman.param");
          break;
        default:
          strcpy(geomPath, "data/drivencavity.geom");
          strcpy(paramPath, "data/drivencavity.param");
          break;
      }
    } else if(strncmp(argv[i], "-P",2)==0) { // parameter path
      // test whether path at flag or afterwards
      if (strlen(argv[i])>2) {
        paramPath = &argv[i][2];
      } else {
        paramPath = argv[++i];
      }
    } else if(strncmp(argv[i], "-I",2)==0) { // input paths
      geomPath = new char[120];
      paramPath = new char[120];
      // test whether path at flag or afterwards
      if (strlen(argv[i])>2) {
        sprintf(paramPath, "%s.param", &argv[i][2]);
        sprintf(geomPath, "%s.geom", &argv[i][2]);
      } else {
        sprintf(paramPath, "%s.param", argv[++i]);
        sprintf(geomPath, "%s.geom", argv[i]);
      }
    }
  }
}

#ifdef WRITE_VTK
void writeVTK(VTK &vtk, Compute &comp, bool writeParticles) {
  if (writeParticles) {
    // Create VTK File for Particle tracing in the folder VTK (must exist)
    vtk.InitParticles("VTK/particle");
    vtk.AddParticles("Particle", comp.GetParticles());
    vtk.FinishParticles();
    // Create VTK File for Streaklines in the folder VTK (must exist)
    vtk.InitParticles("VTK/streakline");
    vtk.AddParticles("Streakline", comp.GetStreakline());
    vtk.FinishParticles();
  }
  // Create a VTK File in the folder VTK (must exist)
  vtk.Init("VTK/field");
  vtk.AddField("Velocity", comp.GetU(), comp.GetV());
  vtk.AddScalar("Pressure", comp.GetP());
  vtk.AddScalar("Vorticity", comp.GetVorticity());
  vtk.AddScalar("Streamfunction", comp.GetStreamline());
  vtk.Finish();
}
#endif

int main(int argc, char *argv[]) {
  char *paramPath = nullptr, *geomPath = nullptr;

  parseCommandLine(argc, argv, paramPath, geomPath);

  // Create communicator, parameter and geometry instances and load values
  Communicator comm(&argc, &argv);
  bool printInfo = comm.ThreadNum() == 0;
  #ifdef BLATT4
    printInfo = false;
  #endif

#ifdef WRITE_VTK
  // Delete old VTK files, prevents bugs in paraview
  if(comm.ThreadNum() == 0) {
    int status = system("rm -r ./VTK/* &>/dev/null");
    if (status < 0) {
      std::cerr << "Error while deleting old VTK files!" << std::endl;
    }
  }
#endif // WRITE_VTK

  Parameter param(paramPath, printInfo);
  Geometry geom(comm, geomPath, printInfo);
  if(param.Omega() <= 0.0 || param.Omega() > 2.0) {
    real_t h = std::fmax(geom.Mesh()[0], geom.Mesh()[1]);
    param.Omega() = 2.0 / (1.0 + std::sin(M_PI * h));
    #ifndef BLATT4
      std::cout << "Set new omega = " << param.Omega() << std::endl;
    #endif
  }

  // Create the fluid solver
  Compute comp(geom, param, comm);

#ifdef DEBUG_VISU
  // Create and initialize the visualization
  const int hRes = 800;
  const int vRes = (int)(hRes * geom.TotalLength()[1] / geom.TotalLength()[0]);
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(hRes, vRes, comm.ThreadNum(), comm.ThreadIdx(), comm.ThreadDim());

  const Grid *visugrid;
  visugrid = comp.GetVelocity();
  real_t nextTimeVisu = 0;
#endif // DEBUG_VISU

#ifdef WRITE_VTK
  // Create a VTK generator
  VTK vtk(geom.Mesh(), geom.SizeP(), geom.Offset(), geom.TotalSize(), comm.ThreadNum(),
          comm.ThreadIdx(), comm.ThreadDim());
  real_t nextTimeVTK = 0;
#endif // WRITE_VTK

  bool run = true;
  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {

#ifdef DEBUG_VISU
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
#endif // DEBUG_VISU

#ifdef WRITE_VTK
    if(comp.GetTime() >= nextTimeVTK ) {
      writeVTK(vtk, comp, comm.ThreadCnt() == 1);
      nextTimeVTK += param.VtkDt();
    }
#endif // WRITE_VTK

/*
    std::cout << std::endl << "U: " << comp.GetU()->Size();
    comp.GetU()->print();
    std::cout << std::endl << "V: " << comp.GetV()->Size();
    comp.GetV()->print();
    std::cout << std::endl << "P: " << comp.GetP()->Size();
    comp.GetP()->print();
    std::cin.get();
*/
    comp.TimeStep(printInfo);
#ifdef DEBUG_VISU
    // Gather if one process stopped
    run = comm.gatherAnd(run);
#endif //DEBUG_VISU
  }

#ifdef WRITE_VTK
  if(comp.GetTime() >= nextTimeVTK ) {
    writeVTK(vtk, comp, comm.ThreadCnt() == 1);
  }
#endif // WRITE_VTK

  return 0;
}

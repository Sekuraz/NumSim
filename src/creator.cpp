/*
 * Copyright (C) 2016   Stephan Lunowa, Markus Baur, Jonas Harsch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <iostream>
#include <cstring>
#include "typedef.hpp"
using namespace std;

/// Returns the number of argv which contains the option
int getOpt(const char *opt, int &argc, char **&argv) {
  for(int i = 1; i < argc ; ++i) {
    if(strcmp(argv[i], opt) == 0) return i;
  }
  return 0;
}

struct param_t {
  real_t  re;
  real_t  omega;
  real_t  alpha;
  real_t  dt;
  real_t  tend;
  real_t  eps;
  real_t  tau;
  index_t itermax;
  real_t  vtkDt;
  real_t  visuDt;
  real_t  particlePosX;
  real_t  particlePosY;

  param_t() : re(1000), omega(1.7), alpha(0.9),  dt(0), tend(50), eps(1e-3),
              tau(0.5), itermax(100), vtkDt(0.2), visuDt(0.1), particlePosX(0.25),
              particlePosY(0.25) {}
};

struct geom_t {
  multi_index_t size;
  multi_real_t  length;
  multi_real_t  velocity;
  real_t        pressure;
  int           type;

  geom_t() : size(128), length(1), velocity(1,0), pressure(0), type(0) {}
};

int main (int argc, char **argv) {
  int pos;
  index_t count = 0;
  char filename[100] = "drivencavity";
  param_t param;
  geom_t geom;

  pos = getOpt("-h",argc,argv);
  if (!pos) pos = getOpt("--help",argc,argv);
  if (pos) {
    cout << endl << "Creator - geometry and parameter file creator." << endl << endl
         << "Usage: creator [arguments]\tGenerates a geometry and a parameter file" << endl << endl
         << "General arguments:" << endl
         << "Arguments must be seperated by a space from their additional parameters." << endl
         << "\t-h, --help\t\tPrints this message." << endl
         << "\t-o <name>\t\tSaves the geometry and the parameters in \"<name>.param\"" << endl
         << "\t\t\tand \"<name>.geom\". The default filenames are" << endl
         << "\t\t\t\"drivencavity.param\" and \"drivencavity.geom\"." << endl << endl
         << "Parameters:" << endl
         << "If some of these arguments are given, the \".param\" file is written." << endl
         << "\t-re <float>\t\tReynolds number. Default: 1000" << endl
         << "\t-omg <float>\tSOR parameter. Default: is 1.7" << endl
         << "\t-alpha <float>\tDonor-Cell weighting parameter. Default: 0.9" << endl
         << "\t-dt <float>\t\tTime-step. Default: 0 (no absolute restriction)" << endl
         << "\t-eps <float>\tTolerance of the solver. Default: 0.001" << endl
         << "\t-tau <float>\tSavety factor in time-step restriction. Default: 0.5" << endl
         << "\t-tend <float>\tFinal (end) time. Default: 50.0" << endl
         << "\t-iter <int>\t\tMaximal number of solver iterations per time-step. Default: 100" << endl
         << "\t-vtkDt <float>\tTime-step for writing vtk-files. Default: 0.2" << endl
         << "\t-visuDt <float>\tTime-step for real time visualization. Default: 0.1" << endl
         << "\t-ppos <float>xfloat\tInital position of the particle for particle tracing." << endl
         << "\t\t\tDefault: 0.25x0.25" << endl << endl
         << "Geometry:" << endl
         << "If some of these arguments are given, the \".geom\" file is written." << endl
         << "The default is a driven cavity. Pre-defined geometries must be set before" << endl
         << "changing other parameters." << endl
         << "\t-pre <num>\t\tChoose among some ready-made geometries." << endl
         << "\t\t\t1: Simple channel" << endl
         << "\t\t\t2: Pressure driven channel" << endl
         << "\t\t\t3: Channel with step" << endl
         << "\t\t\t4: Channel with barrier" << endl
         << "\t\t\tThese worlds might change the default values." << endl
         << "\t-length <x>x<y>\tReal real-world geometry size. Default 1.0x1.0" << endl
         << "\t-pressure <float>\tPressure at dirichlet pressure boundary. Default: 0.0" << endl
         << "\t-size <x>x<y>\tSize of the grid. Default: 128x128" << endl
         << "\t-vel <x>x<y>\tVelocity of the fluid at dirichlet boundary. Default: 1.0x0.0" << endl
         << endl;
    return 0;
  }
  pos = getOpt("-o",argc,argv);
  if (pos) sscanf(argv[pos+1],"%s",filename);

  pos = getOpt("-alpha",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.alpha);
  }
  pos = getOpt("-dt",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.dt);
  }
  pos = getOpt("-eps",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.eps);
  }
  pos = getOpt("-iter",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lu",&param.itermax);
  }
  pos = getOpt("-omg",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.omega);
  }
  pos = getOpt("-re",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.re);
  }
  pos = getOpt("-tau",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.tau);
  }
  pos = getOpt("-tend",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.tend);
  }
  pos = getOpt("-visuDt",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.visuDt);
  }
  pos = getOpt("-vtkDt",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&param.vtkDt);
  }
  pos = getOpt("-ppos",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lfx%lf", &param.particlePosX,&param.particlePosY);
    sscanf(argv[pos+1],"%lf",&param.particlePosX);
  }

  if(count > 0) {
    FILE* handle;
    char name[100];
    sprintf(name,"%s.param",filename);
    handle = fopen(name,"w");
    fprintf(handle,"re = %lf\n",param.re);
    fprintf(handle,"omega = %lf\n",param.omega);
    fprintf(handle,"alpha = %lf\n",param.alpha);
    fprintf(handle,"dt = %lf\n",param.dt);
    fprintf(handle,"tend = %lf\n",param.tend);
    fprintf(handle,"eps = %lf\n",param.eps);
    fprintf(handle,"tau = %lf\n",param.tau);
    fprintf(handle,"itermax = %lu\n",param.itermax);
    fprintf(handle,"vtkDt = %lf\n",param.vtkDt);
    fprintf(handle,"visuDt = %lf\n",param.visuDt);
    fprintf(handle,"particlePosX = %lf\n",param.particlePosX);
    fprintf(handle,"particlePosY = %lf\n",param.particlePosY);
    fclose(handle);
    count = 0;
  }

  pos = getOpt("-pre",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%i",&geom.type);
  }
  switch (geom.type) {
  case 1:
    geom.length[0] = 5;
    geom.size[0] = 160;
    geom.size[1] = 32;
    break;
  case 2:
  case 3:
  case 4:
    geom.length[0] = 5;
    geom.size[0] = 160;
    geom.size[1] = 32;
    geom.velocity[0] = 0;
    geom.pressure = 0.1;
    break;
  default:
    break;
  };
  pos = getOpt("-length",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lfx%lf",&geom.length[0],&geom.length[1]);
  }
  pos = getOpt("-size",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lux%lu",&geom.size[0],&geom.size[1]);
  }
  pos = getOpt("-vel",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lfx%lf",&geom.velocity[0],&geom.velocity[1]);
  }
  pos = getOpt("-pressure",argc,argv);
  if (pos) {
    count++;
    sscanf(argv[pos+1],"%lf",&geom.pressure);
  }

  if(count > 0) {
    FILE* handle;
    char name[120];
    sprintf(name,"%s.geom",filename);
    handle = fopen(name,"w");
    fprintf(handle, "size = %lu %lu\n", geom.size[0], geom.size[1]);
    fprintf(handle, "length = %lf %lf\n", geom.length[0], geom.length[1]);
    fprintf(handle, "velocity = %lf %lf\n", geom.velocity[0], geom.velocity[1]);
    fprintf(handle, "pressure = %lf\n", geom.pressure);
    fprintf(handle, "geometry = %ifree\n", geom.type);
    // print geometry shape
    fprintf(handle,"#");
    for (index_t i = 0; i < geom.size[0]-2; ++i) {
      switch (geom.type) {
      case 1:
      case 2:
      case 3:
      case 4:
        fprintf(handle,"#");
        break;
      default:
        fprintf(handle,"I");
        break;
      }
    }
    fprintf(handle,"#\n");
    for (index_t j = 0; j < geom.size[1]-2; ++j) {
      switch (geom.type) {
      case 1:
        fprintf(handle,"V");
        for (index_t i = 0; i < geom.size[0]-2; ++i) fprintf(handle," ");
        fprintf(handle,"O\n");
        break;
      case 2:
        fprintf(handle,"|");
        for (index_t i = 0; i < geom.size[0]-2; ++i) fprintf(handle," ");
        fprintf(handle,"O\n");
        break;
      case 3:
        if (j+1 < geom.size[1]/2)
          fprintf(handle,"|");
        else
          fprintf(handle,"#");
        for (index_t i = 0; i < geom.size[0]-2; ++i) {
          if (j+1 < geom.size[1]/2 || i > geom.size[1]/2)
            fprintf(handle," ");
          else
            fprintf(handle,"#");
        }
        fprintf(handle,"O\n");
        break;
      case 4:
        fprintf(handle,"|");
        for (index_t i = 0; i < geom.size[0]-2; ++i) {
          if (j+1 < geom.size[1]/4 || j+1 > geom.size[1]*3/4)
            fprintf(handle," ");
          else if (i == j+geom.size[1]/4 || i == j+geom.size[1]/4+1)
            fprintf(handle,"#");
          else if (j == geom.size[1]/4 && i+1 == j+geom.size[1]/4)
            fprintf(handle,"#");
          else if (j+2 == geom.size[1]*3/4 && i-2 == j+geom.size[1]/4)
            fprintf(handle,"#");
          else
            fprintf(handle," ");
        }
        fprintf(handle,"O\n");
        break;
      default:
        fprintf(handle,"#");
        for (index_t i = 0; i < geom.size[0]-2; ++i) fprintf(handle," ");
        fprintf(handle,"#\n");
        break;
      }
    }
    for (index_t i = 0; i < geom.size[0]; ++i) fprintf(handle,"#");
    fclose(handle);
  }
  return 0;
}

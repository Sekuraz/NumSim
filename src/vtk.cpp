/*
 * Copyright (C) 2015   Malte Brunn
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

#include "vtk.hpp"
#include <cstring>
#include <cstdio>
//------------------------------------------------------------------------------
index_t VTK::_cnt = 0;
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_index_t &size)
    : VTK(h, size, multi_index_t(0.0)) {}
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_index_t &size, const multi_index_t &offset)
    : VTK(h, size, offset, size, 0, multi_index_t(0), multi_index_t(1)) {}
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_index_t &size,
      const multi_index_t &offset, const multi_index_t &globalSize, const int &rank,
      const multi_index_t &tidx, const multi_index_t &tdim)
    : _h(h), _size(size), _offset(offset), _globalSize(globalSize), _rank(rank),
    _tidx(tidx), _tdim(tdim), _handle(nullptr), _phandle(nullptr) {}
//------------------------------------------------------------------------------
void VTK::Init(const char *path) {
  if(this->_handle) {
    return;
  }
  int flength = strlen(path) + 30;
  char *filename;
  filename = new char[flength];
  sprintf(filename, "%s_r%02lu-%02lu-%02lu_%05lu.vts", (strlen(path) > 0)? path : "field",
          this->_tidx[0], this->_tidx[1], (DIM == 3? this->_tidx[2] : 0), this->_cnt);
  this->_handle = fopen(filename, "w");

  if(this->_rank == 0) {
    sprintf(filename, "%s_%05lu.pvts",  (strlen(path) > 0)? path : "field", this->_cnt);
    this->_phandle = fopen(filename, "w");
  }
  delete[] filename;

  fprintf(this->_handle, "<?xml version=\"1.0\"?>\n");
  fprintf(this->_handle, "<VTKFile type=\"StructuredGrid\">\n");
  fprintf(this->_handle, "<StructuredGrid WholeExtent=\"%lu %lu %lu %lu %lu %lu\" GhostLevel=\"0\">\n",
          this->_offset[0], this->_offset[0] + (this->_size[0] - 2),
          this->_offset[1], this->_offset[1] + (this->_size[1] - 2),
          (DIM == 3 ? this->_offset[2] : 0),
          (DIM == 3 ? (this->_offset[2] + (this->_size[2] - 2)) : 0));
  fprintf(this->_handle, "<Piece Extent=\"%lu %lu %lu %lu %lu %lu \">\n",
          this->_offset[0], this->_offset[0] + (this->_size[0] - 2),
          this->_offset[1], this->_offset[1] + (this->_size[1] - 2),
          (DIM == 3 ? (this->_offset[2]) : 0),
          (DIM == 3 ? (this->_offset[2] + (this->_size[2] - 2)) : 0));
  fprintf(this->_handle, "<Points>\n");
  fprintf(this->_handle, "<DataArray type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n");

  for(index_t z = 0; z < (DIM == 3 ? this->_size[2]-1 : 1); ++z) {
    for (index_t y = 0; y < this->_size[1]-1; ++y) {
      for (index_t x = 0; x < this->_size[0]-1; ++x) {
        fprintf(this->_handle, "%le %le %le\n",
                (double)(x + this->_offset[0]) * this->_h[0],
                (double)(y + this->_offset[1]) * this->_h[1],
                (DIM == 3 ? (double)(z + this->_offset[2]) * this->_h[2] : 0));
      }
    }
  }

  fprintf(this->_handle, "</DataArray>\n");
  fprintf(this->_handle, "</Points>\n");
  fprintf(this->_handle, "<PointData>\n");

  // write header information in master file
  if(this->_rank == 0) {
    fprintf(this->_phandle, "<?xml version=\"1.0\"?>\n");
    fprintf(this->_phandle, "<VTKFile type=\"PStructuredGrid\">\n");

    // define whole domain extent
    fprintf(this->_phandle,
        "<PStructuredGrid WholeExtent=\"0 %lu 0 %lu 0 %lu\" GhostLevel=\"0\">\n",
        this->_globalSize[0], this->_globalSize[1], (DIM == 3 ? this->_globalSize[2] : 0));

    // announce arrays in 3 dimensions
    fprintf(this->_phandle, "<PPoints>\n");
    fprintf(this->_phandle, "<PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
    fprintf(this->_phandle, "</PPoints>\n");

    // announce pieces and their respective extents
    for(index_t z = 0; z < ((DIM == 3)? this->_tdim[2] : 1); ++z) {
      for(index_t y = 0; y < this->_tdim[1]; ++y) {
        for (index_t x = 0; x < this->_tdim[0]; ++x) {
          // create filename string and build file handle from this.
          filename = new char[flength];
          sprintf(filename, "%s_r%02lu-%02lu-%02lu_%05lu.vts", (strlen(path) > 0)? path : "field",
                  x, y, z, this->_cnt);
          fprintf(this->_phandle, "<Piece Extent=\"%lu %lu %lu %lu %lu %lu\" Source=\"%s\"/>\n",
                  x * (this->_size[0] - 2),
                  ((x == this->_tdim[0]-1)? this->_globalSize[0] : (x + 1) * (this->_size[0] - 2)),
                  y * (this->_size[1] - 2),
                  ((y == this->_tdim[1]-1)? this->_globalSize[1] : (y + 1) * (this->_size[1] - 2)),
                  (DIM == 3 ? (z * (this->_size[2] - 2)) : 0),
                  (DIM == 3 ? ((z == this->_tdim[2]-1)? this->_globalSize[2] : (z + 1) * (this->_size[2] - 2)) : 0),
                  &filename[4]);
          delete[] filename;
        }
      }
    }

    // begin announcing payload
    fprintf(this->_phandle, "<PPointData>\n");
  }
}
//------------------------------------------------------------------------------
void VTK::Finish() {
  if(!this->_handle) {
    return;
  }

  fprintf(this->_handle, "</PointData>\n");
  fprintf(this->_handle, "</Piece>\n");
  fprintf(this->_handle, "</StructuredGrid>\n");
  fprintf(this->_handle, "</VTKFile>\n");

  fclose(this->_handle);
  this->_handle = nullptr;
  this->_cnt++;

  if(this->_rank == 0) {
    fprintf(this->_phandle, "</PPointData>\n");
    fprintf(this->_phandle, "</PStructuredGrid>\n");
    fprintf(this->_phandle, "</VTKFile>\n");

    fclose(this->_phandle);
    this->_phandle = nullptr;
  }
}
//------------------------------------------------------------------------------
void VTK::AddScalar(const char *title, const Grid *grid) {
  if (!this->_handle || (this->_rank == 0 && !this->_phandle)) {
    return;
  }

  fprintf(this->_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", title);

  multi_real_t pos;
  for (index_t z = 0; z < (DIM == 3 ? this->_size[2]-1 : 1); ++z) {
#if DIM == 3
    pos[2] = (double)z * _h[2];
#endif // DIM
    for (index_t y = 0; y < this->_size[1]-1; ++y) {
      pos[1] = (double)y * this->_h[1];
      for (index_t x = 0; x < this->_size[0]-1; ++x) {
        pos[0] = (double)x * this->_h[0];
        fprintf(this->_handle, "%le ", (double)grid->Interpolate(pos));
#if DIM == 3
        fprintf(this->_handle, "\n");
#endif // DIM
      }
#if DIM == 2
      fprintf(this->_handle, "\n");
#endif // DIM
    }
  }

  fprintf(this->_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if(this->_rank == 0) {
    fprintf(this->_phandle, "<PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"/>\n", title);
  }
}
//------------------------------------------------------------------------------
void VTK::AddField(const char *title, const Grid *v1, const Grid *v2) {
  if (!this->_handle || (this->_rank == 0 && !this->_phandle)) {
    return;
  }

  fprintf(this->_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                         "NumberOfComponents=\"3\">\n", title);

  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (index_t y = 0; y < this->_size[1]-1; ++y) {
    pos[1] = (double)y * this->_h[1];
    for (index_t x = 0; x < this->_size[0]-1; ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(this->_handle, "%le %le 0\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos));
    }
  }

  fprintf(this->_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if(this->_rank == 0) {
    fprintf(this->_phandle, "<PDataArray type=\"Float64\" Name=\"%s\" "
            "format=\"ascii\" NumberOfComponents=\"3\"/>\n", title);
  }
}
//------------------------------------------------------------------------------
void VTK::AddField(const char *title, const Grid *v1, const Grid *v2,
                   const Grid *v3) {
  if (!this->_handle || (this->_rank == 0 && !this->_phandle)) {
    return;
  }

  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n", title);

  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (index_t y = 0; y < this->_size[1]-1; ++y) {
    pos[1] = (double)y * _h[1];
    for (index_t x = 0; x < this->_size[0]-1; ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(this->_handle, "%le %le %le\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos), (double)v3->Interpolate(pos));
    }
  }

  fprintf(this->_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if(this->_rank == 0) {
    fprintf(this->_phandle, "<PDataArray type=\"Float64\" Name=\"%s\" "
            "format=\"ascii\" NumberOfComponents=\"3\"/>\n", title);
  }
}

//------------------------------------------------------------------------------
void VTK::InitParticles(const char *path){
  if (_handle) {
    return;
  }
  int flength = strlen(path) + 30;
  char *filename;
  filename = new char[flength];
  sprintf(filename, "%s_%05lu.particles", (strlen(path) > 0)? path : "trace", _cnt);
  this->_handle = fopen(filename, "w");
  delete[] filename;

  fprintf(this->_handle, "<?xml version=\"1.0\"?>\n");
  fprintf(this->_handle, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(this->_handle, "<PolyData>\n");
}
//------------------------------------------------------------------------------
void VTK::FinishParticles(){
  if (!_handle) {
    return;
  }

  fprintf(_handle, "</PolyData>\n");
  fprintf(_handle, "</VTKFile>\n");

  fclose(_handle);

  _handle = NULL;
}
//------------------------------------------------------------------------------

void VTK::AddParticles(const char *title, const std::list<multi_real_t> particles){
  if (!_handle)  {
    return;
  }
  fprintf(_handle, "<Piece NumberOfPoints=\"%lu \" NumberOfVerts=\"0\" NumberOfLines= "
                   "\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n", particles.size());
  fprintf(_handle, "<Points>\n");
  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",title);

  std::list<multi_real_t>::const_iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {
    multi_real_t data = *it;
    fprintf( _handle, "%le %le %le\n", data[0], data[1], (DIM == 3 ? data[2] : 0) );
  }

  fprintf(_handle, "</DataArray>\n");
  fprintf(_handle, "</Points>\n");
  fprintf(_handle, "</Piece>\n");
}

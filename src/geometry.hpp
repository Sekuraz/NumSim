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

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
/// The geometry of the underlying physical domain.
class Geometry {
public:
  /// Constructs a default geometry:
  /// driven cavity with 128 x 128 grid, no-slip boundary conditions
  /// as shown below
  /// <pre>
  ///       u=1, v=0
  ///     -------------
  ///     |           |
  /// u=0 |           | u=0
  /// v=0 |           | v=0
  ///     |           |
  ///     -------------
  ///       u=0, v=0
  /// </pre>
  Geometry(const Communicator &comm, const char file[] = nullptr, const bool printinfo = false)
      : Geometry(comm, multi_index_t(128), file, printinfo) {}
  Geometry(const Communicator &comm, const multi_index_t& size,
           const char file[] = nullptr, const bool printinfo = false);
  ~Geometry() { if(this->_flags != nullptr) delete[] this->_flags; }

  /// \brief Loads a geometry from a file
  /// Expects free geometry in the following character coded format
  /// * ' ' Fluid (interior cell)
  /// * '#' Wall/Obstacle/NoSlip boundary (u = v = 0)
  /// * 'I' General Inflow boundary (u = u_0, v = v_0)
  /// * 'H' Horizontal Inflow boundary (u = u_0, v = ?)
  /// * 'V' Vertical Inflow boundary (u = ?, v = v_0)
  /// * 'O' Outflow boundary (du/dx = 0 resp. dv/dy = 0)
  /// * '|' Vertical Slip-boundary (u = 0, dv/dx = 0, dp determined by parameter pressure)
  /// * '-' Horizontal Slip-boundary (du/dy = 0, v = 0, dp determined by parameter pressure)
  void Load(const char file[], const bool printinfo = false);

  /// Returns the global offset in each dimension
  inline const multi_index_t& Offset() const { return this->_offset; }
  /// Returns the total number of cells in each dimension
  inline const multi_index_t& TotalSize() const { return this->_totalSize; }
  /// Returns the number of cells in each dimension
  inline const multi_index_t &Size() const { return this->_size; }
  /// Returns the number of variables (grid-size) in each dimension
  inline const multi_index_t &SizeP() const { return this->_sizeP; }
  /// Returns the total data size of the Grid
  inline const index_t &DataSize() const { return this->_sizeData; }

  /// Returns the local length of the domain
  inline const multi_real_t &Length() const { return this->_length; }
  /// Returns the local length of the domain
  inline const multi_real_t &TotalLength() const { return this->_totalLength; }
  /// Returns the velocity at the boundary
  inline const multi_real_t &Velocity() const { return this->_velocity; }
  /// Returns the meshwidth
  inline const multi_real_t &Mesh() const { return this->_h; }
  /// Returns the invers meshwidth
  inline const multi_real_t &invMesh() const { return this->_hInv; }
  /// Returns whether the geometry is free.
  inline bool isFree() const { return this->_free; }
  /// Returns the local number of fluid cells
  inline const index_t &NumFluid() const { return this->_N; }

  /// Returns whether the cell at the position of the iterator is fluid cell.
  inline bool isFluid(const Iterator &it) const { return ( this->flag(it) == ' ' ); }

  /// Updates the pressure field p for a free geometry
  void Update_P(Grid &p) const;
  /// Updates the velocity fields u and v for a free geometry
  void Update(Grid &u, Grid &v) const;

  bool noslip(const Iterator &it) const;

  /// \brief Creates a coarsed Geometry for multigrid solver.
  /// The returned Geometry must be deleted after use.
  Geometry* coarse(void) const;

private:
  /// Computes the sizes of the grids.
  void computeSizes();
  /// Returns the grid flag of the cell at the position of it.
  const char& flag(const Iterator &it) const;

  /// Returns the parabolic velocity profile (1/2*Re*dp/L_x*y*(y-L_y)).
  real_t parabolic(const Iterator &it) const;

  const Communicator &_comm; ///< Communicator for boundary exchange and local sizes of Grids
  bool _free;                ///< Whether use free geometry given in loaded file or driven cavity
  char *_flags;              ///< flag field indicating the type of the cells
  multi_index_t _offset;     ///< global offset of the field
  multi_index_t _totalSize;  ///< cartesian total number of inner cells
  multi_index_t _size;       ///< cartesian local number of cells
  multi_index_t _sizeP;      ///< cartesian size of the local Grid
  index_t _sizeData;         ///< number of cells in a Grid
  multi_real_t _totalLength; ///< total length of the physical domain in each dimension
  multi_real_t _length;      ///< local length of the physical domain in each dimension
  multi_real_t _h;           ///< mesh width of the Grids
  multi_real_t _hInv;        ///< the invers of the mesh width
  index_t _N;                ///< local number of fluid cells
  multi_real_t _velocity;    ///< constant boundary velocities (u,v) at upper boundary
  real_t _pressure;          ///< constant pressure difference between left and right
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP

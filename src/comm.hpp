/*
 * Copyright (C) 2015   Malte Brunn, Markus Baur, Jonas Harsch, Stephan Lunowa
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

#include <mpi.h>
#include "grid.hpp"
#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __COMM_HPP
#define __COMM_HPP
//------------------------------------------------------------------------------

/// Class for all data exchange between the processes (via MPI).
class Communicator {
public:
  /// \brief Communicator constructor. Initializes MPI Environment.
  /// \param [in] argc Number of arguments program was started with
  /// \param [in] argv Arguments passed to the program on start
  Communicator(int* argc, char*** argv);
  /// Communicator destructor. Finalizes MPI Environment.
  ~Communicator();

  /// Returns the cartesian position of the current process with respect to the lower left corner.
  inline const multi_index_t& ThreadIdx() const { return this->_tidx; };
  /// Returns the cartesian sizes the domain is partitioned among all processes.
  inline const multi_index_t& ThreadDim() const { return this->_tdim; };

  /// Returns the process' number in continous counting.
  inline const int& ThreadNum() const { return this->_rank; };
  /// Returns the total number of processes.
  inline const int& ThreadCnt() const { return this->_size; };

  /// Returns whether this process is a red or a black field.
  inline const bool& EvenOdd() const { return this->_evenodd; };

  /// \brief Returns if the values of all processes are true.
  /// \param [in] val The data over which the logical operation is to be calculated
  /// \return Whether the values of all processes are true.
  bool gatherAnd(const bool& val) const;
  /// \brief Gets the sum of all values and distributes the result among all processes.
  /// \param [in] val The data over which the sum is to be calculated
  /// \return The sum of the values of all processes
  real_t gatherSum(const real_t& val) const;
  /// \brief Gets the minimum of all values and distributes the result among all processes.
  /// \param [in] val The data over which the minimum is to be calculated
  /// \return The minimum of the values of all processes
  real_t gatherMin(const real_t& val) const;
  /// \brief Gets the maximum of all values and distributes the result among all processes.
  /// \param [in] val The data over which the maximum is to be calculated
  /// \return The maximum of the values of all processes
  real_t gatherMax(const real_t& val) const;

  /// \brief Synchronizes the ghost layers of the grid
  /// \param [in] grid  The grid of which the values are synced
  void copyBoundary(Grid& grid) const;

  /// \brief Transfers the offset of the streamfunction and returns the own offset
  /// \param [in] grid  The grid of the streamfunction
  real_t copyOffset(const Grid& grid) const;

  /// Returns if the domain of this process is at the left boundary of the whole domain
  inline bool isLeft() const { return (this->_tidx[0] == 0); }
  /// Returns if the domain of this process is at the right boundary of the whole domain
  inline bool isRight() const { return (this->_tidx[0] == this->_tdim[0]-1); }
  /// Returns if the domain of this process is at the upper boundary of the whole domain
  inline bool isTop() const { return (this->_tidx[1] == this->_tdim[1]-1); }
  /// Returns if the domain of this process is at the lower boundary of the whole domain
  inline bool isBottom() const { return (this->_tidx[1] == 0); }

private:
  multi_index_t _tidx; ///< The cartesian position of the process.
  multi_index_t _tdim; ///< The cartesian sizes in the dimensions.
  MPI_Comm _mpi_cart_comm; ///< MPI_Comm with cartesian coordinates.
  int _rank; ///< The continous rank of this process.
  int _size; ///< The total amount of processes.
  bool _evenodd; ///< whether this process has a black or red domain.

  /// Synchronizes the left ghost layer of the grid
  /// \param [in] grid  The grid of which the values are synced
  bool copyLeftBoundary (Grid& grid) const;
  /// Synchronizes the right ghost layer of the grid
  /// \param [in] grid  The grid of which the values are synced
  bool copyRightBoundary (Grid& grid) const;
  /// Synchronizes the upper ghost layer of the grid
  /// \param [in] grid  The grid of which the values are synced
  bool copyTopBoundary (Grid& grid) const;
  /// Synchronizes the lower ghost layer of the grid
  /// \param [in] grid  The grid of which the values are synced
  bool copyBottomBoundary (Grid& grid) const;
};
//------------------------------------------------------------------------------
#endif // __COMM_HPP
//------------------------------------------------------------------------------

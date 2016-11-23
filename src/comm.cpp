#include <mpi.h>
#include "comm.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "typedef.hpp"

Communicator::Communicator (int* argc, char*** argv) : _tidx(), _tdim(), _mpi_cart_comm() {
  MPI_Init(argc, argv);

  // Get the rank of the current thread
  MPI_Comm_rank(MPI_COMM_WORLD, &this->_rank);
  // Get total number of threads
  MPI_Comm_size(MPI_COMM_WORLD, &this->_size);

  // compute idx and dim from rank and size
  int dims[DIM], periodic[DIM];
  for(index_t dim = 0; dim < DIM; dim++) {
    dims[dim] = 0;
    periodic[dim] = 0;
  }
  MPI_Dims_create(this->_size, (int)DIM, dims);
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_tdim[dim] = dims[dim];
  }

  MPI_Cart_create(MPI_COMM_WORLD, (int)DIM, dims, periodic, 0, &this->_mpi_cart_comm);
  int coords[DIM];
  MPI_Cart_get(this->_mpi_cart_comm, (int)DIM, dims, periodic, coords);
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_tidx[dim] = coords[dim];
  }

  // TODO: for testing
  printf("Dims: %lu %lu Rank: %i Idx: %lu %lu \n", _tdim[0], _tdim[1], _rank, _tidx[0], _tidx[1]);
}

Communicator::~Communicator () {
  MPI_Finalize();
}

real_t Communicator::gatherSum(const real_t& val) const {
  real_t sum;
  MPI_Allreduce(&val, &sum, 1, MPI_REAL_TYPE, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}

real_t Communicator::gatherMin(const real_t& val) const {
  real_t min;
  MPI_Allreduce(&val, &min, 1, MPI_REAL_TYPE, MPI_MIN, MPI_COMM_WORLD);
  return min;
}

real_t Communicator::gatherMax(const real_t& val) const {
  real_t max;
  MPI_Allreduce(&val, &max, 1, MPI_REAL_TYPE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}

void Communicator::copyBoundary(Grid& grid) const {
  // copy boundaries depending on evenodd such that no deadlock occurs
  if(this->EvenOdd()) {
    if(!this->isLeft()) {
      this->copyLeftBoundary(grid);
    }
    if(!this->isRight()) {
      this->copyRightBoundary(grid);
    }
    if(!this->isTop()) {
      this->copyTopBoundary(grid);
    }
    if(!this->isBottom()) {
      this->copyBottomBoundary(grid);
    }
  } else {
    if(!this->isRight()) {
      this->copyRightBoundary(grid);
    }
    if(!this->isLeft()) {
      this->copyLeftBoundary(grid);
    }
    if(!this->isBottom()) {
      this->copyBottomBoundary(grid);
    }
    if(!this->isTop()) {
      this->copyTopBoundary(grid);
    }
  }
}

bool Communicator::copyLeftBoundary(Grid& grid) const {
  index_t size = grid.Size()[1];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid, 2);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Right());
    i++;
  }

  int source, dest;
  MPI_Cart_shift(this->_mpi_cart_comm, 0, -1, &source, &dest);
  MPI_Sendrecv_replace(buffer, i, MPI_REAL_TYPE, dest, 0, dest, 0, this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return false;
}
bool Communicator::copyRightBoundary(Grid& grid) const {
  index_t size = grid.Size()[1];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid, 4);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Left());
    i++;
  }

  int source, dest;
  MPI_Cart_shift(this->_mpi_cart_comm, 0, 1, &source, &dest);
  MPI_Sendrecv_replace(buffer, i, MPI_REAL_TYPE, dest, 0, dest, 0, this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return false;
}
bool Communicator::copyTopBoundary(Grid& grid) const {
  index_t size = grid.Size()[0];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid, 1);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Down());
    i++;
  }

  int source, dest;
  MPI_Cart_shift(this->_mpi_cart_comm, 1, 1, &source, &dest);
  MPI_Sendrecv_replace(buffer, i, MPI_REAL_TYPE, dest, 0, dest, 0, this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return false;
}
bool Communicator::copyBottomBoundary(Grid& grid) const {
  index_t size = grid.Size()[0];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid, 3);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Top());
    i++;
  }

  int source, dest;
  MPI_Cart_shift(this->_mpi_cart_comm, 1, -1, &source, &dest);
  MPI_Sendrecv_replace(buffer, i, MPI_REAL_TYPE, dest, 0, dest, 0, this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return false;
}


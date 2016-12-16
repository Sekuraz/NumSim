#include <mpi.h>
#include "comm.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "typedef.hpp"

Communicator::Communicator (int* argc, char*** argv) : _tidx(), _tdim(), _mpi_cart_comm() {
  MPI_Init(argc, argv);

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

  MPI_Cart_create(MPI_COMM_WORLD, (int)DIM, dims, periodic, 1, &this->_mpi_cart_comm);
  int coords[DIM];
  MPI_Cart_get(this->_mpi_cart_comm, (int)DIM, dims, periodic, coords);
  this->_evenodd = true;
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_tidx[dim] = coords[dim];
    this->_evenodd ^= ((coords[dim] % 2) == 1);
  }

  // Get the rank of the current thread
  MPI_Comm_rank(this->_mpi_cart_comm, &this->_rank);

  // TODO: for testing
  printf("Dims: %lu %lu Rank: %i Idx: %lu %lu odd: %i\n", _tdim[0], _tdim[1], _rank, _tidx[0], _tidx[1], _evenodd);
}

Communicator::~Communicator () {
  MPI_Finalize();
}

bool Communicator::gatherAnd(const bool& val) const {
  bool ret;
  MPI_Allreduce(&val, &ret, 1, MPI_CXX_BOOL, MPI_LAND, this->_mpi_cart_comm);
  return ret;
}

real_t Communicator::gatherSum(const real_t& val) const {
  real_t sum;
  MPI_Allreduce(&val, &sum, 1, MPI_REAL_TYPE, MPI_SUM, this->_mpi_cart_comm);
  return sum;
}

real_t Communicator::gatherMin(const real_t& val) const {
  real_t min;
  MPI_Allreduce(&val, &min, 1, MPI_REAL_TYPE, MPI_MIN, this->_mpi_cart_comm);
  return min;
}

real_t Communicator::gatherMax(const real_t& val) const {
  real_t max;
  MPI_Allreduce(&val, &max, 1, MPI_REAL_TYPE, MPI_MAX, this->_mpi_cart_comm);
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

// Transfers the offset of the streamfunction and returns the own offset
// \param [in] grid  The grid of the streamfunction
real_t Communicator::copyOffset(const Grid& grid) const {
  int source, dest;
  real_t offset = 0.0;
  if(this->isBottom()) {
    if(!this->isLeft()) {
      real_t off;
      MPI_Cart_shift(this->_mpi_cart_comm, 0, -1, &source, &dest);
      MPI_Recv(&off, 1, MPI_REAL_TYPE, dest, 1, this->_mpi_cart_comm, MPI_STATUS_IGNORE);
      offset += off;
    }
    if(!this->isRight()) {
      real_t off = offset + grid.Cell(Iterator(grid.Geom(), grid.Size()[0]-1));
      MPI_Cart_shift(this->_mpi_cart_comm, 0, 1, &source, &dest);
      MPI_Send(&off, 1, MPI_REAL_TYPE, dest, 1, this->_mpi_cart_comm);
    }
  }
  if(!this->isBottom()) {
    real_t off;
    MPI_Cart_shift(this->_mpi_cart_comm, 1, -1, &source, &dest);
    MPI_Recv(&off, 1, MPI_REAL_TYPE, dest, 1, this->_mpi_cart_comm, MPI_STATUS_IGNORE);
    offset += off;
  }
  if(!this->isTop()) {
    real_t off = offset + grid.Cell(Iterator(grid.Geom(), grid.Size()[0]*(grid.Size()[1]-1)));
    MPI_Cart_shift(this->_mpi_cart_comm, 1, 1, &source, &dest);
    MPI_Send(&off, 1, MPI_REAL_TYPE, dest, 1, this->_mpi_cart_comm);
  }

  return offset;
}

bool Communicator::copyLeftBoundary(Grid& grid) const {
  const index_t& size = grid.Size()[1];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid.Geom(), BoundaryIterator::boundary::left);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Right());
    i++;
  }

  int source, dest, result;
  MPI_Cart_shift(this->_mpi_cart_comm, 0, -1, &source, &dest);
  result = MPI_Sendrecv_replace(buffer, size, MPI_REAL_TYPE, dest, 0, dest, 0,
                                this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return (result == MPI_SUCCESS);
}
bool Communicator::copyRightBoundary(Grid& grid) const {
  const index_t& size = grid.Size()[1];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid.Geom(), BoundaryIterator::boundary::right);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Left());
    i++;
  }

  int source, dest, result;
  MPI_Cart_shift(this->_mpi_cart_comm, 0, 1, &source, &dest);
  result = MPI_Sendrecv_replace(buffer, size, MPI_REAL_TYPE, dest, 0, dest, 0,
                                this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return (result == MPI_SUCCESS);
}
bool Communicator::copyTopBoundary(Grid& grid) const {
  const index_t& size = grid.Size()[0];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid.Geom(), BoundaryIterator::boundary::top);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Down());
    i++;
  }

  int source, dest, result;
  MPI_Cart_shift(this->_mpi_cart_comm, 1, 1, &source, &dest);
  result = MPI_Sendrecv_replace(buffer, size, MPI_REAL_TYPE, dest, 0, dest, 0,
                                this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return (result == MPI_SUCCESS);
}
bool Communicator::copyBottomBoundary(Grid& grid) const {
  const index_t& size = grid.Size()[0];
  real_t* buffer = new real_t[size];

  index_t i = 0;
  BoundaryIterator it (grid.Geom(), BoundaryIterator::boundary::down);
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid.Cell(it.Top());
    i++;
  }

  int source, dest, result;
  MPI_Cart_shift(this->_mpi_cart_comm, 1, -1, &source, &dest);
  result = MPI_Sendrecv_replace(buffer, size, MPI_REAL_TYPE, dest, 0, dest, 0,
                                this->_mpi_cart_comm, MPI_STATUS_IGNORE);

  i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    grid.Cell(it) = buffer[i];
    i++;
  }

  delete[] buffer;

  return (result == MPI_SUCCESS);
}


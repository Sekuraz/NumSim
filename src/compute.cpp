#include <cstddef>
#include "typedef.hpp"
#include "compute.hpp"

// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry *geom, const Parameter *param)
    : _geom(geom), _param(param) {
  // TODO: create grids
};
// Deletes all grids
Compute::~Compute() {
  // TODO: delete grids
};

// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
  // TODO: timestep
}

// Computes and returns the absolute velocity
const Grid *Compute::GetVelocity() { return nullptr; };
// Computes and returns the vorticity
const Grid *Compute::GetVorticity() { return nullptr; };
// Computes and returns the stream line values
const Grid *Compute::GetStream() { return nullptr; };

// Compute the new velocites u,v
void Compute::NewVelocities(const real_t &dt) {
  // TODO
}
// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt) {
  // TODO
}
// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt) {
  // TODO
}


#include <iostream>
#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"

// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry &geom, const Parameter &param)
    : _t(0), _dtlimit(param.Dt()), _epslimit(param.Eps()),
    _F(new Grid(geom)), _G(new Grid(geom)), _rhs(new Grid(geom)),
    _tmp(new Grid(geom)), _geom(geom), _param(param) {

  // initialize the solver
  this->_solver = new SOR(geom, param.Omega());

  // initialize u,v,p
  const multi_real_t &h = this->_geom.Mesh();
  // TODO: Compute offsets
  multi_real_t offset(h[0], h[1]/2);
  _u = new Grid(geom, offset);
  offset[0] = h[0]/2;
  offset[1] = h[1];
  _v = new Grid(geom, offset);
  offset[1] = h[1]/2;
  _p = new Grid(geom, offset);
  this->_u->Initialize(0);
  this->_v->Initialize(0);
  this->_p->Initialize(0);
  this->_F->Initialize(0);
  this->_G->Initialize(0);
  this->_rhs->Initialize(0);
  this->_tmp->Initialize(0);
}
// Deletes all grids
Compute::~Compute() {
  delete _u; delete _v; delete _p; delete _F; delete _G;
  delete _rhs; delete _tmp; delete _solver;
}

// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
  // TODO: timestep

  real_t dt = _param.Dt();

  // TODO: test
  Iterator it(_geom);
  while(it.Valid()){
    _u->Cell(it) = it;
    _v->Cell(it) = it;
    _p->Cell(it) = it;
    it.Next();
  }
//  if(it.Valid()) _u->Cell(it) = it;
//  if(_t == 0) _u->Data()[0] = .1;

/*  this->RHS(dt);
  for(index_t i = 0; i < this->_param.IterMax(); i++) {
    real_t residual = this->_solver->Cycle(*(this->_p), *(this->_rhs));
    if(printInfo) {
      std::cout << "Error in iteration " << i << ":\t" << residual << std::endl;
    }
  }
*/
  this->_t += dt;
}

// Computes and returns the absolute velocity
const Grid *Compute::GetVelocity() { return nullptr; }
// Computes and returns the vorticity
const Grid *Compute::GetVorticity() { return nullptr; }
// Computes and returns the stream line values
const Grid *Compute::GetStream() { return nullptr; }

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


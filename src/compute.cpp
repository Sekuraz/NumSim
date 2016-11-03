#include <iostream>
#include <algorithm>
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
    _F(new Grid(geom, Grid::type::u)), _G(new Grid(geom, Grid::type::v)),
    _rhs(new Grid(geom, Grid::type::p)),
    _tmp(new Grid(geom, Grid::type::p)), _geom(geom), _param(param) {

  // initialize the solver
  this->_solver = new SOR(geom, param.Omega());

  // initialize u,v,p
  const multi_real_t &h = this->_geom.Mesh();
  multi_real_t offset(0, h[1]/2);
  _u = new Grid(geom, Grid::type::u, offset);
  offset[0] = h[0]/2;
  offset[1] = 0;
  _v = new Grid(geom, Grid::type::v, offset);
  offset[1] = h[1]/2;
  _p = new Grid(geom, Grid::type::p, offset);
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
  // TODO: test
  Iterator u_it(*_u);
  while(u_it.Valid()){
    _u->Cell(u_it) = u_it;
    u_it.Next();
  }
  Iterator v_it(*_v);
  while(v_it.Valid()){
    _v->Cell(v_it) = v_it;
    v_it.Next();
  }
  if(_t < _param.Dt()) {
    const multi_index_t &size = _p->Size();
    real_t* data = _p->Data();
    data[size[0]*(size[1]/2-1)+size[0]/2-1] = 1;
    data[size[0]*(size[1]/2-1)+size[0]/2] = 1;
    data[size[0]*size[1]/2+size[0]/2-1] = 1;
    data[size[0]*size[1]/2+size[0]/2] = 1;
  }
  // TODO: End test

  // compute dt
  real_t dt = this->_param.Dt();
  // Test CFL condition
  const multi_real_t &h = this->_geom.Mesh();
  real_t dtCond = this->_param.Tau() * std::min(h[0] / this->_u->AbsMax(), h[1] / this->_v->AbsMax());
  if(dtCond < dt) {
    dt = dtCond;
    std::cerr << "Warning: Compute: CFL > tau. New dt = " << dt << "!" << std::endl;
  }
  // Test Pr condition
  dtCond = this->_param.Re() * 0.5 * (h[0]*h[0]*h[1]*h[1])/(h[0]*h[0]+h[1]*h[1]);
  if(dtCond < dt) {
    dt = dtCond;
    std::cerr << "Warning: Compute: Pr > tau. New dt = " << dt << "!" << std::endl;
  }

  // set boundary values for u, v, F, G
  this->_geom.Update_U(*(this->_u));
  this->_geom.Update_V(*(this->_v));
  // TODO: update F,G
//  this->_geom.Update_F(*(this->_F));
//  this->_geom.Update_G(*(this->_G));

  // compute the momentum equations
  this->MomentumEqu(dt);

  // compute right-hand-side of the poisson equation
  this->RHS(dt);

  // solve the poisson equation
  for(index_t i = 0; i < this->_param.IterMax(); i++) {
    // set boundary values for p
    this->_geom.Update_P(*(this->_p));
    real_t res = this->_solver->Cycle(*(this->_p), *(this->_rhs));
    if(printInfo) {
      std::cout << "\tError in SOR-iteration " << i << ":\t" << res << std::endl;
    }
    if(res < this->_epslimit) break;
  }

  // compute new velocites
  this->NewVelocities(dt);

  // update the time
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
  // compute u
  for(InteriorIterator it(*(this->_u)); it.Valid(); it.Next()) {
    // TODO: dp/dx (it)
    //this->_u->Cell(it) = this->_F->Cell(it) - dt * this->_p->dx_r();
  }
  // compute v
  for(InteriorIterator it(*(this->_v)); it.Valid(); it.Next()) {
    // TODO: dp/dy (it)
    //this->_v->Cell(it) = this->_G->Cell(it) - dt * this->_p->dy_r();
  }
}
// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt) {
  // compute F
  for(InteriorIterator it(*(this->_F)); it.Valid(); it.Next()) {
    this->_F->Cell(it) = this->_u->Cell(it) + dt * (
        ( this->_u->dxx(it) + this->_u->dyy(it) )/this->_param.Re()
        - this->_u->DC_udu_x(it, this->_param.Alpha())
        - this->_u->DC_udv_x(it, this->_param.Alpha(), this->_v) );
  }
  // compute G
  for(InteriorIterator it(*(this->_G)); it.Valid(); it.Next()) {
    this->_G->Cell(it) = this->_v->Cell(it) + dt * (
        ( this->_v->dxx(it) + this->_v->dyy(it) )/this->_param.Re()
        - this->_v->DC_vdv_y(it, this->_param.Alpha())
        - this->_v->DC_vdu_y(it, this->_param.Alpha(), this->_u) );
  }
}

// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt) {
  const multi_real_t &h = this->_geom.Mesh();
  for(InteriorIterator it(*(this->_p)); it.Valid(); it.Next()) {
    // TODO F(it), G(it)
    // this->_rhs->Cell(it) = ( (this->_F->Cell() - this->_F->Cell())/h[0]
    //                        + (this->_G->Cell() - this->_G->Cell())/h[1]) / dt;
  }
}


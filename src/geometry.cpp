#include <iostream>
#include <fstream>
#include <string>
#include "geometry.hpp"

using namespace std;

Geometry::Geometry() : _size(128,128), _sizeU(_size), _sizeV(_size),
    _sizeP(_size), _length(1,1) {
  for(index_t dim = 0; dim < DIM; dim++) {
    this->_sizeU[dim] += (dim == 0)? 1 : 2;
    this->_sizeV[dim] += (dim == 1)? 1 : 2;
    this->_sizeP[dim] += 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
  }
}

/// Loads a geometry from a file
void Geometry::Load(const char file[]) {
  string param;
  real_t value;
  ifstream in(file);
  while(in.good()) {
    in >> ws >> param >> ws;
    in.ignore(); // removes '=' between parameter and value
    in >> ws >> value >> ws;
    if(!param.compare("Size") || !param.compare("size")) {
      this->_size[0] = (index_t)value;
      cout << "Geometry: Load size = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_size[dim] = (index_t)value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Length") || !param.compare("length")) {
      this->_length[0] = value;
      cout << "Geometry: Load length = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_length[dim] = value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Velocity") || !param.compare("velocity")) {
      this->_velocity[0] = value;
      cout << "Geometry: Load velocity = " << value;
      for(index_t dim = 1; dim < DIM; dim++) {
        in >> value >> ws;
        this->_velocity[dim] = value;
        cout << " " << value;
      }
      cout << endl;
    } else if(!param.compare("Pressure") || !param.compare("pressure")) {
      cout << "Geometry: Load pressure = " << value << endl;
      this->_pressure = value;
    } else {
      cerr << "Geometry: unknown identifier " << param;
    }
  }
  in.close();

  for(index_t dim = 0; dim < DIM; dim++) {
    this->_sizeU[dim] = this->_size[dim] + (dim == 0)? 1 : 2;
    this->_sizeV[dim] = this->_size[dim] + (dim == 1)? 1 : 2;
    this->_sizeP[dim] = this->_size[dim] + 2;
    this->_h[dim] = this->_length[dim] / this->_size[dim];
  }
}

// Updates the velocity field u
void Geometry::Update_U(Grid &u) const {
    //TODO noop
}

// Updates the velocity field v
void Geometry::Update_V(Grid &v) const {
    //TODO noop
}

// Updates the pressure field p
void Geometry::Update_P(Grid &p) const {
    //TODO noop
}

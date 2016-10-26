#include "geometry.hpp"


Geometry::Geometry() : _size(128,128), _length(1,1), _h(1.0/128, 1.0/128) {}

/// Loads a geometry from a file
void Geometry::Load(const char *file) {
    //TODO noop
}

/// Updates the velocity field u
void Geometry::Update_U(Grid *u) const {
    //TODO noop
}

/// Updates the velocity field v
void Geometry::Update_V(Grid *v) const {
    //TODO noop
}

/// Updates the pressure field p
void Geometry::Update_P(Grid *p) const {
    //TODO noop
}

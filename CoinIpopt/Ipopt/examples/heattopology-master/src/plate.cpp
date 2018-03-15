#include "plate.hpp"
#include "print.hpp"
#include "solver.hpp"
#include <iostream>
#include <math.h>

namespace wit {
Plate::Plate(unsigned int nx, unsigned int ny, double sizex, double sizey) {
  _Temp.resize(nx, ny);
  _Material.resize(nx, ny);
  _sizex = sizex;
  _sizey = sizey;
  _nx = nx;
  _ny = ny;
}

void Plate::solve() {
  initialize(_Material);
  find_minimum(_Material, _Temp);
}

void Plate::print() {
  wit::print_heat(_Temp);
  wit::print_material(_Material);
}

Plate::~Plate() {}
}

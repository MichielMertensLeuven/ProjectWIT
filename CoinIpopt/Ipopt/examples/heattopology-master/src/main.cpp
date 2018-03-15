#include "plate.hpp"
#include <iostream>

#define SIZE_MESH 500
#define DIMENSION_X 1
#define DIMENSION_Y 1

bool tests();

int main(int args, char *argv[]) {
  // tests
  if (tests()) {
    // solve problem
    wit::Plate Plaatje(SIZE_MESH, SIZE_MESH, DIMENSION_X, DIMENSION_Y);
    Plaatje.solve();
    Plaatje.print();
  } else {
    std::cout << "Error" << std::endl;
    return -1;
  }
  return 0;
}

bool tests() {
  // tests
  return true;
}

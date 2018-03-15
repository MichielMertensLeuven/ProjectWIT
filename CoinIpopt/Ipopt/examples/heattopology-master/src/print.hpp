#ifndef wit_print_hpp
#define wit_print_hpp

#include <Eigen/Dense>
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;

namespace wit {
void print_heat(MatrixXd &Temp);
void print_material(MatrixXd &Material);
}

#endif

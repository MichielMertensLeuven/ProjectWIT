//#include "opti.hpp"
#include "solver.hpp"

using Eigen::MatrixXd;

namespace wit {
void initialize(MatrixXd &mat) { mat.setRandom(); }
bool find_minimum(MatrixXd &Material, MatrixXd &Temp) { return false; }
}

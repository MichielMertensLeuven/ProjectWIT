#ifndef tws_solver_hpp
#define tws_solver_hpp

#include <Eigen/Dense>

namespace wit {
void initialize(Eigen::MatrixXd &mat);
bool find_minimum(Eigen::MatrixXd &Material, Eigen::MatrixXd &Temp);
}
#endif

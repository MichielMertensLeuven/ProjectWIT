#ifndef wit_matrix_hpp
#define wit_matrix_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace wit {
void generate_s(int n, double size, Eigen::MatrixXd &k,
                Eigen::SparseMatrix<double> &S);

void generate_q(int n, double T, double Q, double size, Eigen::VectorXd &q);

void generate_k(int n, Eigen::MatrixXd &A, int p, double k_plastic,
                double k_steel, Eigen::MatrixXd &k);

void generate_dsdpi(int n, double size, int i, int j, Eigen::MatrixXd &k,
                    Eigen::MatrixXd &A, int p, double k_plastic, double k_steel,
                    Eigen::SparseMatrix<double> &dSdpi);
void generate_dydp(int n, double size, Eigen::MatrixXd &k, Eigen::MatrixXd &A,
                   int p, double k_plastic, double k_steel, Eigen::VectorXd &u,
                   Eigen::SparseMatrix<double> &S, Eigen::VectorXd &dydp);
void generate_u(int n, double size, int p, double k_plastic, double k_steel,
                double T, double Q, Eigen::MatrixXd &A, Eigen::VectorXd &u);
}

#endif

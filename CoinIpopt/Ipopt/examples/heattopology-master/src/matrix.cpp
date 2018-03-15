#include "matrix.hpp"

namespace wit {

// Generate the left half side of the equation
// n: amount of grid points in each direction
// size: size of the plate in meters
// k: material parameter in Watt per Kelvin per meter per meter
// S: the system matrix which is generated
void generate_s(int n, double size, Eigen::MatrixXd &k,
                Eigen::SparseMatrix<double> &S) {
  S.reserve(Eigen::VectorXi::Constant(n * n, 5));
  double mu = (n - 1) * (n - 1) / (size * size); // 1/dx^2
  int ref = (int)(0.3 * n); // index where Dirichlet boundary starts
  double k_iplus12;
  double k_jplus12;
  double k_imin12;
  double k_jmin12;
  for (int i = 1; i < n - 1; ++i) {
    for (int j = 1; j < n - 1; ++j) { // General case
      k_iplus12 =
          2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
      k_jplus12 =
          2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
      k_imin12 =
          2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
      k_jmin12 =
          2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
      S.insert(i + j * n, i + (j)*n) =
          -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
      S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
      S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
      S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
      S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
    }
  }

  int j = 0; // Lower edge
  for (int i = 1; i < n - 1; ++i) {
    k_iplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
    k_jplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
    k_imin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
    k_jmin12 = 0;                                              // k(i,j-1/2)
    S.insert(i + j * n, i + (j)*n) =
        -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
    S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
    S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
    S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
  }
  j = n - 1; // Upper edge
  for (int i = 1; i < n - 1; ++i) {
    k_iplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
    k_jplus12 = 0;                                             // k(i,j+1/2)
    k_imin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
    k_jmin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
    S.insert(i + j * n, i + (j)*n) =
        -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
    S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
    S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
    S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
  }
  int i = 0; // Left edge
  for (int j = 1; j < n - 1; ++j) {
    k_iplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
    k_jplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
    k_imin12 = 0;                                              // k(i-1/2,j)
    k_jmin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
    S.insert(i + j * n, i + (j)*n) =
        -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
    S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
    S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
    S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
  }
  i = n - 1; // Right edge, lower part
  for (int j = 1; j < ref; ++j) {
    k_iplus12 = 0; // k(i+1/2,j)
    k_jplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
    k_imin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
    k_jmin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
    S.insert(i + j * n, i + (j)*n) =
        -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
    S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
    S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
    S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
  }
  for (int j = ref; j < n - ref; ++j) { // Dirichlet boundary condition
    S.insert(i + j * n, i + (j)*n) = -1;
  }
  for (int j = n - ref; j < n - 1; ++j) { // Right edge, upper part
    k_iplus12 = 0;                        // k(i+1/2,j)
    k_jplus12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
    k_imin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
    k_jmin12 =
        2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
    S.insert(i + j * n, i + (j)*n) =
        -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
    S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
    S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
    S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
  }

  i = 0; // LL corner
  j = 0;
  k_iplus12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
  k_jplus12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
  k_imin12 = 0;                                              // k(i-1/2,j)
  k_jmin12 = 0;                                              // k(i,j-1/2)
  S.insert(i + j * n, i + (j)*n) =
      -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
  S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
  S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
  i = n - 1;     // LR corner
  k_iplus12 = 0; // k(i+1/2,j)
  k_jplus12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j + 1)); // k(i,j+1/2)
  k_imin12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
  k_jmin12 = 0;                                              // k(i,j-1/2)
  S.insert(i + j * n, i + (j)*n) =
      -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
  S.insert(i + j * n, i + (j + 1) * n) = k_jplus12 * mu;
  S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
  j = n - 1;     // UR corner
  k_iplus12 = 0; // k(i+1/2,j)
  k_jplus12 = 0; // k(i,j+1/2)
  k_imin12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i - 1, j)); // k(i-1/2,j)
  k_jmin12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
  S.insert(i + j * n, i + (j)*n) =
      -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
  S.insert(i + j * n, i - 1 + (j)*n) = k_imin12 * mu;
  S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;
  i = 0; // UL corner
  k_iplus12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i + 1, j)); // k(i+1/2,j)
  k_jplus12 = 0;                                             // k(i,j+1/2)
  k_imin12 = 0;                                              // k(i-1/2,j)
  k_jmin12 =
      2 / (1 / k.coeffRef(i, j) + 1 / k.coeffRef(i, j - 1)); // k(i,j-1/2)
  S.insert(i + j * n, i + (j)*n) =
      -(k_iplus12 + k_imin12 + k_jplus12 + k_jmin12) * mu;
  S.insert(i + j * n, i + 1 + (j)*n) = k_iplus12 * mu;
  S.insert(i + j * n, i + (j - 1) * n) = k_jmin12 * mu;

  S = -S;
}

// Generate the right half side of the equation
// n: amount of grid points in each direction
// T: external temperature in Kelvin
// Q: input heat in Watt
// size: size of the plate in meters
// k: heat conductivity matrix in Watt per Kelvin per meter per meter
// q: the system vector which is generated
void generate_q(int n, double T, double Q, double size, Eigen::VectorXd &q) {
  for (int i = 0; i < n * n; ++i) {
    q.coeffRef(i) = Q / size / size;
  }
  int ref = (int)(0.3 * n);
  int i = n - 1; // Dirichlet boundary condition
  for (int j = ref; j < n - ref; ++j) {
    q.coeffRef(i + j * n) = T;
  }
}

// Generate the thermal conductivity matrix in function of the parameter matrix
// A
// n: amount of grid points in each direction
// A: material parameter matrix, 0 means plastic, 1 means steel
// p: penalization factor
// k_plastic: thermal conductivity of plastic in Watt per Kelvin per meter per
// meter
// k_steel  : thermal conductivity of steel in Watt per Kelvin per meter per
// meter
// k: the thermal conductivity matrix which is generated
void generate_k(int n, Eigen::MatrixXd &A, int p, double k_plastic,
                double k_steel, Eigen::MatrixXd &k) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      k.coeffRef(i, j) =
          k_plastic + (k_steel - k_plastic) * pow(A.coeffRef(i, j), p);
    }
  }
}

// Generate the derivative of the system with respect to a certain parameter
// n: amount of grid points in each direction
// size: size of the plate in meters
// i: index in the x direction of the parameter to which the system should be
// derived
// j: index in the y direction of the parameter to which the system should be
// derived
// k: heat conductivity matrix in Watt per Kelvin per meter per meter
// A: material parameter matrix, 0 means plastic, 1 means steel
// p: penalization factor
// k_plastic: thermal conductivity of plastic in Watt per Kelvin per meter per
// meter
// k_steel  : thermal conductivity of steel in Watt per Kelvin per meter per
// meter
// dSdpi: the derivative of the system matrix with respect to a certain
// parameter
void generate_dsdpi(int n, double size, int i, int j, Eigen::MatrixXd &k,
                    Eigen::MatrixXd &A, int p, double k_plastic, double k_steel,
                    Eigen::SparseMatrix<double> &dSdpi) {
  double mu = (n - 1) * (n - 1) / (size * size);
  double k_iplus;
  double k_jplus;
  double k_imin;
  double k_jmin;
  double dkdp = p * (k_steel - k_plastic) * pow(A.coeffRef(i, j), p - 1);
  int ref = (int)(0.3 * n); // index where Dirichlet boundary starts

  if (i == 0 && j == 0) { // LL corner
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
              pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
  } else if (i == n - 1 && j == 0) { // LR corner
    k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
              pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
             pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;
    dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;

    dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
    dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
  } else if (i == 0 && j == n - 1) { // UL corner
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
             pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
  } else if (i == n - 1 && j == n - 1) { // UR corner
    k_iplus = 0;
    k_jplus = 0;
    k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
             pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
             pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

    dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
  } else if (j == 0) { // Lower edge
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
              pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
             pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jmin = 0;
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;
    dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
    dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
  } else if (j == n - 1) { // Upper edge
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jplus = 0;
    k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
             pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
             pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
  } else if (i == 0) { // Left edge
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
              pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    k_imin = 0;
    k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
             pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
    dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
  } else if (i == n - 1) {         // Right edge
    if (j >= ref && j < n - ref) { // Dirichlet
      dSdpi.setZero();
    } else {
      k_iplus = 0;
      k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
                pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
      k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
               pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
      k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
               pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
      dSdpi.insert(i + j * n, i + (j)*n) =
          -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
      dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;
      dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;
      dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

      dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
      dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
      dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
    }
  } else { // General case
    k_iplus = 2 * pow(k.coeffRef(i + 1, j), 2) /
              pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jplus = 2 * pow(k.coeffRef(i, j + 1), 2) /
              pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    k_imin = 2 * pow(k.coeffRef(i + 1, j), 2) /
             pow((k.coeffRef(i + 1, j) + k.coeffRef(i, j)), 2);
    k_jmin = 2 * pow(k.coeffRef(i, j + 1), 2) /
             pow((k.coeffRef(i, j + 1) + k.coeffRef(i, j)), 2);
    dSdpi.insert(i + j * n, i + (j)*n) =
        -(k_iplus + k_imin + k_jplus + k_jmin) * mu * dkdp;
    dSdpi.insert(i + j * n, i + 1 + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j + 1) * n) = k_jplus * mu * dkdp;
    dSdpi.insert(i + j * n, i - 1 + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + j * n, i + (j - 1) * n) = k_jmin * mu * dkdp;

    dSdpi.insert(i + 1 + j * n, i + (j)*n) = k_imin * mu * dkdp;
    dSdpi.insert(i + (j + 1) * n, i + (j)*n) = k_jmin * mu * dkdp;
    dSdpi.insert(i - 1 + j * n, i + (j)*n) = k_iplus * mu * dkdp;
    dSdpi.insert(i + (j - 1) * n, i + (j)*n) = k_jplus * mu * dkdp;
  }
}

// Generate the derivative of the objective function with respect to all
// parameters
// n: amount of grid points in each direction
// size: size of the plate in meters
// k: heat conductivity matrix in Watt per Kelvin per meter per meter
// A: material parameter matrix, 0 means plastic, 1 means steel
// p: penalization factor
// k_plastic: thermal conductivity of plastic in Watt per Kelvin per meter per
// meter
// k_steel  : thermal conductivity of steel in Watt per Kelvin per meter per
// meter
// solver: solver to solve system S*dudpi=-dSdpi*u
// u: vector containing the heat solution
// S: system used to find u
// dydp: the derivative of the objective function with respect to all parameters
void generate_dydp(int n, double size, Eigen::MatrixXd &k, Eigen::MatrixXd &A,
                   int p, double k_plastic, double k_steel, Eigen::VectorXd &u,
                   Eigen::SparseMatrix<double> &S, Eigen::VectorXd &dydp) {
  Eigen::SparseMatrix<double> dSdpi(n * n, n * n);
  Eigen::VectorXd part_obj_func(n * n); // dydpi=part_obj_func.dot(dudpi)
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(S);
  solver.factorize(S);
  part_obj_func = (Eigen::VectorXd::Ones(n * n) + 2 * u -
                   2 * Eigen::VectorXd::Ones(n * n) *
                       Eigen::VectorXd::Ones(n * n).transpose() * u / (n * n)) /
                  (2 * n * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      generate_dsdpi(n, size, i, j, k, A, p, k_plastic, k_steel, dSdpi);
      dydp.coeffRef(i + j * n) = part_obj_func.dot(solver.solve(-dSdpi * u));
    }
  }
}

// Generate the derivative of the objective function with respect to all
// parameters
// n: amount of grid points in each direction
// size: size of the plate in meters
// p: penalization factor
// k_plastic: thermal conductivity of plastic in Watt per Kelvin per meter per
// meter
// k_steel  : thermal conductivity of steel in Watt per Kelvin per meter per
// meter
// T: external temperature in Kelvin
// Q: input heat in Watt
// A: material parameter matrix, 0 means plastic, 1 means steel
// u: the solution to the heat equation
void generate_u(int n, double size, int p, double k_plastic, double k_steel,
                double T, double Q, Eigen::MatrixXd &A, Eigen::VectorXd &u) {
  Eigen::MatrixXd k(n, n);
  Eigen::VectorXd q(n);
  Eigen::SparseMatrix<double> S(n * n, n * n);
  S.reserve(Eigen::VectorXi::Constant(n * n, 5));
  generate_k(n, A, p, k_plastic, k_steel, k);
  generate_q(n, T, Q, size, q);
  generate_s(n, size, k, S);
  S.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(S);
  solver.factorize(S);
  u = solver.solve(q);
}
}

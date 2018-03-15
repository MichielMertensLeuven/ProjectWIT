#ifndef wit_plate_hpp
#define wit_plate_hpp

#include <Eigen/Dense>

using Eigen::MatrixXd;

namespace wit {
class Plate {
private:
  unsigned int _nx, _ny;
  double _sizex, _sizey;
  MatrixXd _Temp, _Material;

public:
  Plate(unsigned int nx, unsigned int ny, double sizex, double sizey);
  void solve();
  void print();
  ~Plate();
};
}

#endif

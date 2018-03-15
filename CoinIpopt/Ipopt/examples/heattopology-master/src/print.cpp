#include "print.hpp"

namespace wit {

/*struct rgb {
private:
  int _nx, _ny;
  uint8_t *pr, *pg, *pb;

public:
  rgb(int nx, int ny) {
    _nx = nx;
    _ny = ny;
    uint8_t r[nx * ny], g[nx * ny], b[nx * ny];
    pr = r;
    pg = g;
    pb = b;
  }
  setr(int x, int y, double r) {
    r[(_nx - 1) * y + x]->(uint8_t)round(r);
  }
}; */

void print_heat(MatrixXd &Temp) {
  const double maxTemp = Temp.maxCoeff();
  const double minTemp = Temp.minCoeff();

  int nx = Temp.rows();
  int ny = Temp.cols();

  /* This is where we'll store the pixels */
  uint8_t r[nx][ny], g[nx][ny], b[nx][ny];

  /* Run the width and height of the image */
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      double factor = (Temp(x, y) - minTemp) / (maxTemp - minTemp);
      if (factor < 0.5) {
        r[y][x] = (uint8_t)0.;
        g[y][x] = (uint8_t)round(2. * factor * 255);
        b[y][x] = (uint8_t)round((1. - 2. * factor) * 255);
      } else {
        r[y][x] = (uint8_t)round(2. * (factor - .5) * 255);
        g[y][x] = (uint8_t)round((1. - 2. * (factor - .5)) * 255);
        b[y][x] = (uint8_t)0.;
      }
    }
  }

  /* Write the ppm-formatted file */
  FILE *out = fopen("plate_heat.ppm", "w");
  /* Header:
   * Magic number, width, height, maxvalue
   */
  fprintf(out, "P6 %d %d %d\n", nx, ny, 255);
  /* Rows. 2-byte values if max > 255, 1-byte otherwise.
   * These loops explicitly show where every single byte
   * goes; in practice, it would be faster and shorter to
   * interleave the arrays and write bigger blocks of
   * contiguous data.
   */
  for (size_t i = 0; i < (size_t)nx; i++)
    for (size_t j = 0; j < (size_t)ny; j++) {
      fwrite(&r[i][j], sizeof(uint8_t), 1, out); /* Red */
      fwrite(&g[i][j], sizeof(uint8_t), 1, out); /* Green */
      fwrite(&b[i][j], sizeof(uint8_t), 1, out); /* Blue */
    }
  fclose(out);
}

void print_material(MatrixXd &Material) {
  const double maxMaterial = Material.maxCoeff();
  const double minMaterial = Material.minCoeff(); // should be zero

  int nx = Material.rows();
  int ny = Material.cols();

  /* This is where we'll store the pixels */
  uint8_t r[nx][ny], g[nx][ny], b[nx][ny];

  /* Run the width and height of the image */
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      double factor =
          (Material(x, y) - minMaterial) / (maxMaterial - minMaterial);
      r[y][x] = (uint8_t)round((1. - factor) * 255);
      g[y][x] = (uint8_t)round((1. - factor) * 255);
      b[y][x] = (uint8_t)round((1. - factor) * 255);
    }
  }

  /* Write the ppm-formatted file */
  FILE *out = fopen("plate_material.ppm", "w");
  /* Header:
   * Magic number, width, height, maxvalue
   */
  fprintf(out, "P6 %d %d %d\n", nx, ny, 255);
  /* Rows. 2-byte values if max > 255, 1-byte otherwise.
   * These loops explicitly show where every single byte
   * goes; in practice, it would be faster and shorter to
   * interleave the arrays and write bigger blocks of
   * contiguous data.
   */
  for (size_t i = 0; i < (size_t)nx; i++)
    for (size_t j = 0; j < (size_t)ny; j++) {
      fwrite(&r[i][j], sizeof(uint8_t), 1, out); /* Red */
      fwrite(&g[i][j], sizeof(uint8_t), 1, out); /* Green */
      fwrite(&b[i][j], sizeof(uint8_t), 1, out); /* Blue */
    }
  fclose(out);
}
}

#include <cmath>
#include <limits>
#include "solute.h"

void Solute :: init(int n) {
  num = n;
  q = new double[num];
  sig = new double[num];
  eps = new double[num];
  r = new double[num * 4];
}

double * Solute :: centering() {
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double * shift;
  shift = new double[3];
  
  xmin = ymin = zmin = std::numeric_limits<double>::max();
  xmax = ymax = zmax = std::numeric_limits<double>::lowest();

  for (int n = 0; n < num; ++n) {
    int i = n * 4;
    if (xmin > r[i]) xmin = r[i];
    if (ymin > r[i + 1]) ymin = r[i + 1];
    if (zmin > r[i + 2]) zmin = r[i + 2];
    if (xmax < r[i]) xmax = r[i];
    if (ymax < r[i + 1]) ymax = r[i + 1];
    if (zmax < r[i + 2]) zmax = r[i + 2];
  }

  shift[0] = round(- (xmax - xmin) / 2 - xmin);
  shift[1] = round(- (ymax - ymin) / 2 - ymin);
  shift[2] = round(- (zmax - zmin) / 2 - zmin);

  for (int n = 0; n < num; ++n) {
    int i = n * 4;
    r[i] += shift[0];
    r[i + 1] += shift[1];
    r[i + 2] += shift[2];
  }

 return shift;
}

void Solute :: setup_cuda() {
  cudaMalloc(&dq, num * sizeof(double));
  cudaMalloc(&dr, num * sizeof(double4));
  cudaMemcpyAsync(dq, q, num * sizeof(double), cudaMemcpyDefault);
  cudaMemcpyAsync(dr, r, num * sizeof(double4), cudaMemcpyDefault);
}

void Solute :: free_cuda() {
  cudaFree(dq);
  cudaFree(dr);
}


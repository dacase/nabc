#include <iostream>
#include "rism3d.h"

void RISM3D :: cal_LJ() {
  __global__ void LJ(double * du, const double * __restrict__ dsig, 
		     const double * __restrict__ deps, 
		     const double4 * __restrict__ dru,
		     double cut2, double ikbt, double bx, double by, double bz,
		     int nx, int ny, int nz, int natu, int iv);

  const double cut = 1.0e-2;
  const double cut2 = cut * cut;

  cout << "tabulating solute Lennard-Jones potential ..." << endl;

  cudaMalloc(&du, ce -> ngrid * sv -> natv * sizeof(double));
  cudaMalloc(&dsig, su -> num * sv -> natv * sizeof(double));
  cudaMalloc(&deps, su -> num * sv -> natv * sizeof(double));
  cudaMemset(du, 0.0, ce -> ngrid * sv -> natv * sizeof(double));

  siguv = new double[su -> num * sv -> natv];
  epsuv = new double[su -> num * sv -> natv];

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int iu = 0; iu < su -> num; ++iu) {
      int ip = iu + su -> num * iv;
      siguv[ip] = (su -> sig[iu] + sv -> sigv[iv]) * 0.5;
      epsuv[ip] = sqrt (su -> eps[iu] * sv -> epsv[iv] * kcal2J);
    }
  }

  cudaMemcpyAsync(dsig, siguv, su -> num * sv -> natv * sizeof(double),
		  cudaMemcpyDefault);
  cudaMemcpyAsync(deps, epsuv, su -> num * sv -> natv * sizeof(double),
		  cudaMemcpyDefault);

  double iKbT = 1.0 / (avogadoro * boltzmann * sv -> temper);
  for (int iv = 0; iv < sv -> natv; ++iv) {
    LJ <<< g, b >>> (du + (iv * ce -> ngrid), dsig, deps, su -> dr, 
		      cut2, iKbT, ce -> dr[0], ce -> dr[1], ce -> dr[2], 
		      ce -> grid[0], ce -> grid[1], ce -> grid[2], 
		      su -> num, iv);
  }
}

__global__ void LJ(double * du, const double * __restrict__ dsig, 
		   const double * __restrict__ deps, 
		   const double4 * __restrict__ dru,
                   double cut2, double ikbt, double bx, double by, double bz,
                   int nx, int ny, int nz, int natu, int iv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rx = ((int)threadIdx.x - nx / 2) * bx;
  double ry = ((int)blockIdx.x - ny / 2) * by;
  double rz = ((int)blockIdx.y - nz / 2) * bz;
  for (int iu = 0; iu < natu; ++iu) {
    int iuv = iu + natu * iv;
    double dx = rx - dru[iu].x;
    double dy = ry - dru[iu].y;
    double dz = rz - dru[iu].z;
    double r2 = dx * dx + dy * dy + dz * dz ;

    if (r2 < cut2) r2 = cut2;

    double irs2 = dsig[iuv] * dsig[iuv] / r2;

    double irs6 = irs2 * irs2 * irs2;
    du[ip] += deps[iuv] * 4.0 * irs6 * (irs6 - 1.0);
  }
  du[ip] *= ikbt;
}

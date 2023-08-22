#include <iostream>
#include "fft3d.h"

void FFT3D :: initialize (Cell * ce) {
  __global__ void set_kf(double2 * dkf, int nx, int ny, int nz);
  __global__ void set_ir(int * dir);

  ngrid = ce -> ngrid;
  volf = ce -> dv;
  volb = 1.0 / ce -> volume;
  b.x = ce -> grid[0];
  g.x = ce -> grid[1];
  g.y = ce -> grid[2];
  cudaMalloc(&dkf, ngrid * sizeof(double2));
  cudaMalloc(&dir, ngrid * sizeof(int));
  cufftPlan3d(&plan, ce -> grid[0], ce -> grid[1], ce -> grid[2], CUFFT_Z2Z);
  set_kf <<< g, b >>> (dkf, ce -> grid[0], ce -> grid[1], ce -> grid[2]);
  set_ir <<< g, b >>> (dir);
}


void FFT3D :: execute (double2 * da, int key) {
  __global__ void timeirvol(double2 *, const int * __restrict__, double);
  __global__ void timekf(double2 *, const double2 * __restrict__, 
			 const int * __restrict__);
  __global__ void timekb(double2 *, const double2 * __restrict__, 
			 const int * __restrict__, double);

  if (key == - 1) {
    timekf <<< g, b >>> (da, dkf, dir);
    cufftExecZ2Z(plan, da, da, CUFFT_FORWARD);
    timeirvol <<< g, b >>> (da, dir, volf);
  } else {
    timeirvol <<< g, b >>> (da, dir, 1.0);
    cufftExecZ2Z(plan, da, da, CUFFT_INVERSE);
    timekb <<< g, b >>> (da, dkf, dir, volb);
  }
}


__global__ void timekf(double2 * da, const double2 * __restrict__ dkf, 
		       const int * __restrict__ dir) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double tmpr = da[ip].x * dir[ip];
  double tmpi = da[ip].y * dir[ip];
  da[ip].x = tmpr * dkf[ip].x - tmpi * dkf[ip].y;
  da[ip].y = tmpi * dkf[ip].x + tmpr * dkf[ip].y;
}


__global__ void timekb(double2 * da, const double2 * __restrict__ dkf,
		       const int * __restrict__ dir, double vol) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double tmpr = da[ip].x * dir[ip] * vol;
  double tmpi = da[ip].y * dir[ip] * vol;
  da[ip].x = tmpr * dkf[ip].x + tmpi * dkf[ip].y;
  da[ip].y = tmpi * dkf[ip].x - tmpr * dkf[ip].y;
}


__global__ void timeirvol(double2 * da, const int * __restrict__ dir, double vol) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  da[ip].x *= dir[ip] * vol;
  da[ip].y *= dir[ip] * vol;
}


__global__ void set_kf(double2 * dkf, int nx, int ny, int nz) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  int x = threadIdx.x - nx / 2.0;
  int y = blockIdx.x - ny / 2.0;
  int z = blockIdx.y - nz / 2.0;
  double dkx = M_PI / nx;
  double dky = M_PI / ny;
  double dkz = M_PI / nz;
  double dkr = dkx * x + dky * y + dkz * z;
  dkf[ip].x = cos(dkr);
  dkf[ip].y = - sin(dkr);
}


__global__ void set_ir(int * dir) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  if ((threadIdx.x + blockIdx.x + blockIdx.y) % 2 == 0) {
    dir[ip] = 1;
  } else {
    dir[ip] = -1;
  }
}

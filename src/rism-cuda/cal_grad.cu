#include <thrust/device_vector.h>
#include "rism3d.h"

__constant__ double3 dv;
__constant__ int3 grid;
__constant__ int nu;

void RISM3D :: cal_grad(double * & du) {
  __global__ void grad(double * ds, double2 * dguv, double * dsig, 
		       double * deps,  double4 * dr, double * qu, double gv,
		       int natu, int iv, int iu);

  cudaMemcpyToSymbol(dv, ce -> dr, sizeof(double3));
  cudaMemcpyToSymbol(grid, ce -> grid, sizeof(int3));
  double * ds3;
  cudaMalloc(&ds3, g.x * g.y * 3 * sizeof(double));
  thrust::device_ptr<double> ds3_ptr(ds3);

#pragma omp parallel for
  for (int iu = 0; iu < su -> num * 3; ++iu) {
    du[iu] = 0.0;
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int iu = 0; iu < su -> num; ++iu) {
      grad <<< g, b, b.x * 3 * sizeof(double) >>> 
	(ds3, dguv + (iv * ce -> ngrid), dsig, deps, su -> dr, su -> dq,
	 sv -> qv[iv], su -> num, iv, iu);
      for (int i = 0; i< 3; ++i) {
	double s = thrust::reduce(ds3_ptr + (g.x * g.y) * i, 
				  ds3_ptr + (g.x * g.y) * (i + 1));
	du[iu * 3 + i] += s * sv -> rhov[iv];
      }
    }
  }
}


__global__ void grad(double * ds, double2 * dguv, double * dsig, 
		     double * deps,  double4 * dr, double * qu,
		     double qv, int natu, int iv, int iu) {
  extern __shared__ double sdata[];
  const double cc = hartree * bohr * avogadoro;

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  int iuv = iu + iv * natu;

  double dx = ((int)threadIdx.x - grid.x / 2) * dv.x - dr[iu].x;
  double dy = ((int)blockIdx.x - grid.y / 2) * dv.y - dr[iu].y;
  double dz = ((int)blockIdx.y - grid.z / 2) * dv.z - dr[iu].z;
  double r2 = dx * dx + dy * dy + dz * dz;
  double r1 = sqrt(r2);

  if (r1 < dsig[iuv] * 0.5) {
    sdata[threadIdx.x] = 0.0;
    sdata[threadIdx.x + blockDim.x] = 0.0;
    sdata[threadIdx.x + blockDim.x * 2] = 0.0;
  } else {
    double rs2i = dsig[iuv] * dsig[iuv] / r2;
    double rs6i = rs2i * rs2i * rs2i;
    double ulj = deps[iuv] * 24.0 * rs6i / r2 * (2.0 * rs6i - 1.0) 
      * dguv[ip].x;
    double uco = qu[iu] * qv / (r2 * r1) * cc * dguv[ip].x;
    sdata[threadIdx.x] = (ulj + uco) * dx;
    sdata[threadIdx.x + blockDim.x] = (ulj + uco) * dy;
    sdata[threadIdx.x + blockDim.x * 2] = (ulj + uco) * dz;
  }
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x + blockDim.x]
        += sdata[threadIdx.x + blockDim.x + s];
      sdata[threadIdx.x + blockDim.x * 2]
        += sdata[threadIdx.x + blockDim.x * 2 + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
    ds[blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y] =
      sdata[blockDim.x];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 2] =
      sdata[blockDim.x * 2];
  }
}

#include <iostream>
#include "rism3d.h"

void RISM3D :: cal_Coulomb () {
  __global__ void coulomb(double * de, double * dfr,
			  double4 * dru, double * dqu,
			  double dx, double dy, double dz,
			  int nx, int ny, int nz, int natu);
  __global__ void fk(double2 *, const double4 * __restrict__ , 
		     const double4 * __restrict__ , const double * __restrict__, 
		     int);
  __global__ void beta(double * de, double * dfr, double2 * dfk, double ubeta);

  cout << "synthesizing solute Coulomb potential ..." << endl;
  
  cudaMalloc(&de, ce -> ngrid * sizeof(double));
  cudaMalloc(&dfr, ce -> ngrid * sizeof(double));
  cudaMalloc(&dfk, ce -> ngrid * sizeof(double2));
  cudaMemset(de, 0.0, ce -> ngrid * sizeof(double));
  cudaMemset(dfr, 0.0, ce -> ngrid * sizeof(double));
  cudaMemset(dfk, 0.0, ce -> ngrid * sizeof(double2));

  coulomb <<< g, b >>> (de, dfr, su -> dr, su -> dq,
			ce -> dr[0], ce -> dr[1], ce -> dr[2], 
			ce -> grid[0], ce -> grid[1], ce -> grid[2], su -> num);

  fk <<< g, b >>> (dfk, dgv, su -> dr, su -> dq, su -> num);

  double ubeta = hartree * bohr / (boltzmann * sv -> temper);
  beta <<< g, b >>> (de, dfr, dfk, ubeta);
} 


__global__ void coulomb(double * de, double * dfr,
                        double4 * dru, double * dqu,
                        double bx, double by, double bz,
                        int nx, int ny, int nz, int natu) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rx = ((int)threadIdx.x - nx / 2) * bx;
  double ry = ((int)blockIdx.x - ny / 2) * by;
  double rz = ((int)blockIdx.y - nz / 2) * bz;
  for (int iu = 0; iu < natu; ++iu) {
    double delx = rx - dru[iu].x;
    double dely = ry - dru[iu].y;
    double delz = rz - dru[iu].z;
    double ra = sqrt(delx * delx + dely * dely + delz * delz) ;
    if (ra >= 1.0e-5) {
      double qr = dqu[iu] / ra ;
      de[ip] += qr ;
      dfr[ip] += qr * (1 - exp(- ra)) ;
    } else {
      dfr[ip] += dqu[iu] ;
    }
  }
}


__global__ void fk(double2 * dfk, const double4 * __restrict__ dgv, 
		   const double4 * __restrict__ dru, 
		   const double * __restrict__ dqu, int natu) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rk2 = dgv[ip].x * dgv[ip].x
    + dgv[ip].y * dgv[ip].y + dgv[ip].z * dgv[ip].z;
  double rk4i = 1.0 / (rk2 * (rk2 + 1.0));
  for (int iu = 0; iu < natu; ++iu) {
    double ruk = dgv[ip].x * dru[iu].x 
      + dgv[ip].y * dru[iu].y + dgv[ip].z * dru[iu].z;
    double tmp = 4.0 * M_PI * dqu[iu] * rk4i;
    dfk[ip].x += tmp * cos(ruk);
    dfk[ip].y -= tmp * sin(ruk);
  }
}


__global__ void beta(double * de, double * dfr, double2 * dfk, double ubeta) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  de[ip] *= ubeta;
  dfr[ip] *= ubeta;
  dfk[ip].x *= ubeta;
  dfk[ip].y *= ubeta;
}

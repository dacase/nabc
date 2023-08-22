#include "rism3d.h"

void RISM3D :: cal_potential() {
  __global__ void set_du(double * du, double * de, double q);

  cal_LJ();
  cal_Coulomb();

  //  for (int iv = 0; iv < sv -> natv; ++iv) {
  //    set_du <<< g, b >>> (du + (iv * ce -> ngrid), de, sv -> qv[iv]);
  //  }
  cudaFree(dgv);
  //  cudaFree(de);
}


__global__ void set_du(double * du, double * de, double q) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  du[ip] += q * de[ip];
}

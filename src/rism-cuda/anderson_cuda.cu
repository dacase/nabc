__global__ void newdt0 (double * dt, const double * __restrict__ dtr,
                        double * dtp, double * drp) {
    unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
      + blockIdx.y * blockDim.x * gridDim.x;
    dtp[ip] = dt[ip];
    drp[ip] = dtr[ip];
    dt[ip] += dtr[ip];
}


__global__ void newdt(double * dt, const double * __restrict__ dtr, 
		      double * dtp, double * drp,
		      double s1, double s2, double m, int niv, int biv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double u = dt[ip] + s1 * (dtp[ip] - dt[ip]) + s2 * (dtp[ip + niv] - dt[ip]);
  double v = dt[ip] + dtr[ip]
    + s1 * (dtp[ip] + drp[ip] - dt[ip] - dtr[ip])
    + s2 * (dtp[ip + niv] + drp[ip + niv] - dt[ip] - dtr[ip]);
  dtp[ip + biv] = dt[ip];
  drp[ip + biv] = dtr[ip];
  dt[ip] = u + m * (v - u);
}


__global__ void theta30(double * ds, const double * __restrict__ dtr, 
			const double * __restrict__ drp, int niv) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;

  double t1 = dtr[ip] - drp[ip];
  double t2 = dtr[ip] - drp[ip + niv];
  sdata[threadIdx.x] = dtr[ip] * t1;
  sdata[threadIdx.x + blockDim.x] = dtr[ip] * t2;
  sdata[threadIdx.x + blockDim.x * 2] = t1 * t1;
  sdata[threadIdx.x + blockDim.x * 3] = t1 * t2;
  sdata[threadIdx.x + blockDim.x * 4] = t2 * t2;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x + blockDim.x] += sdata[threadIdx.x + blockDim.x + s];
      sdata[threadIdx.x + blockDim.x * 2]
        += sdata[threadIdx.x + blockDim.x * 2 + s];
      sdata[threadIdx.x + blockDim.x * 3]
        += sdata[threadIdx.x + blockDim.x * 3 + s];
      sdata[threadIdx.x + blockDim.x * 4]
        += sdata[threadIdx.x + blockDim.x * 4 + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
    ds[blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y] =
      sdata[blockDim.x];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 2] =
      sdata[blockDim.x * 2];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 3] =
      sdata[blockDim.x * 3];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 4] =
      sdata[blockDim.x * 4];
  }
}


__global__ void theta31(double * ds2, double * ds) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x;

  sdata[threadIdx.x] = ds[ip];
  sdata[threadIdx.x + blockDim.x] = ds[ip + blockDim.x * gridDim.x];
  sdata[threadIdx.x + blockDim.x * 2] = ds[ip + blockDim.x * gridDim.x * 2];
  sdata[threadIdx.x + blockDim.x * 3] = ds[ip + blockDim.x * gridDim.x * 3];
  sdata[threadIdx.x + blockDim.x * 4] = ds[ip + blockDim.x * gridDim.x * 4];
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x + blockDim.x] += sdata[threadIdx.x + blockDim.x + s];
      sdata[threadIdx.x + blockDim.x * 2]
        += sdata[threadIdx.x + blockDim.x * 2 + s];
      sdata[threadIdx.x + blockDim.x * 3]
        += sdata[threadIdx.x + blockDim.x * 3 + s];
      sdata[threadIdx.x + blockDim.x * 4]
        += sdata[threadIdx.x + blockDim.x * 4 + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds2[blockIdx.x] = sdata[0];
    ds2[blockIdx.x +  gridDim.x] = sdata[blockDim.x];
    ds2[blockIdx.x +  gridDim.x * 2] = sdata[blockDim.x * 2];
    ds2[blockIdx.x +  gridDim.x * 3] = sdata[blockDim.x * 3];
    ds2[blockIdx.x +  gridDim.x * 4] = sdata[blockDim.x * 4];
  }
}


__global__ void theta32(double * ds, double * ds2) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x;

  sdata[ip] = ds2[ip];
  sdata[ip + blockDim.x] = ds2[ip + blockDim.x];
  sdata[ip + blockDim.x * 2] = ds2[ip + blockDim.x * 2];
  sdata[ip + blockDim.x * 3] = ds2[ip + blockDim.x * 3];
  sdata[ip + blockDim.x * 4] = ds2[ip + blockDim.x * 4];
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (ip < s) {
      sdata[ip] += sdata[ip + s];
      sdata[ip + blockDim.x] += sdata[ip + blockDim.x + s];
      sdata[ip + blockDim.x * 2] += sdata[ip + blockDim.x * 2 + s];
      sdata[ip + blockDim.x * 3] += sdata[ip + blockDim.x * 3 + s];
      sdata[ip + blockDim.x * 4] += sdata[ip + blockDim.x * 4 + s];
    }
    __syncthreads();
  }
  if (ip == 0) {
    ds[0] = sdata[0];
    ds[1] = sdata[blockDim.x];
    ds[2] = sdata[blockDim.x * 2];
    ds[3] = sdata[blockDim.x * 3];
    ds[4] = sdata[blockDim.x * 4];
  }
}

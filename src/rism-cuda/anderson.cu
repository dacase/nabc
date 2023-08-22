#include <iostream>
#include <thrust/device_vector.h>
#include "anderson.h"

AN2 :: ~AN2 () {
  delete[] a, c, x;
  cudaFree(dtp);
  cudaFree(drp);
  cudaFree(ds);
}


void AN2 :: initialize (Cell * ce, Solvent * sv) {
  ngrid = ce -> ngrid;
  niv = sv -> natv;
  binary = 0;
  cudaMalloc(&dtp, ngrid * niv * 2 * sizeof(double));
  cudaMalloc(&drp, ngrid * niv * 2 * sizeof(double));
  cudaMalloc(&ds, ce -> grid[1] * ce -> grid[2] * 5 * sizeof(double));
  b.x = ce -> grid[0];
  g.x = ce -> grid[1];
  g.y = ce -> grid[2];

  s = new double[5];
  a = new double[4];
  c = new double[2];
  x = new double[2];
  irho = new double[niv];
  for (int iv = 0; iv < niv; ++iv) {
    irho[iv] = 1.0 / sv -> rhov[iv];
  }
}


void AN2 :: calculate (double * dt, double * dtr) {
  __global__ void newdt0(double *, const double * __restrict__,
			 double *, double *);
  __global__ void newdt(double *, const double * __restrict__,
		       double *, double *,
		      double s1, double s2, double m, int niv, int biv);
  void leg(double * a, double * x, double * b, int n, int nn);

  if (count > 0) {
    --count;
    for (int iv = 0; iv < niv; ++iv) {
      int biv = iv + binary * niv;
      newdt0 <<< g, b >>> (dt + (iv * ngrid), dtr + (iv * ngrid),
			    dtp + (biv * ngrid), drp + (biv * ngrid));
    }
  } else {
    cal_theta(dt, dtr);
    a[2] = a[1];
    leg(a, x, c, 2, 2);
    s1 = x[0] + mp * ((binary == 0)? 0 : 1);
    s2 = x[1] + mp * ((binary == 0)? 1 : 0);
    for (int iv = 0; iv < niv; ++iv) {
      newdt <<< g, b >>> (dt + (iv * ngrid), dtr + (iv * ngrid),
			   dtp + (iv * ngrid), drp + (iv * ngrid),
			   s1, s2, m, niv * ngrid, binary * niv * ngrid);
    }
  }
  binary = (binary == 0)? 1 : 0;
}


void AN2 :: cal_theta (double * dt, double * dtr) {
  __global__ void theta30(double *, const double * __restrict__,
                          const double * __restrict__, int);
  c[0] = c[1] = a[0] = a[1] = a[3] = 0.0;
  for (int iv = 0; iv < niv; ++iv) {
    theta30 <<< g, b, b.x * 5 * sizeof(double) >>>
      (ds, dtr + (iv * ngrid), drp + (iv * ngrid), niv * ngrid);
    thrust::device_ptr<double> ds_ptr(ds);
    for (int i = 0; i < 5; ++i) {
      s[i] = thrust::reduce(ds_ptr + i * g.x * g.y, 
			    ds_ptr + (i + 1) * g.x * g.y);
    }
    c[0] += s[0] * irho[iv];
    c[1] += s[1] * irho[iv];
    a[0] += s[2] * irho[iv];
    a[1] += s[3] * irho[iv];
    a[3] += s[4] * irho[iv];
  }
}


void leg(double * a, double * x, double * b, int n, int nn) {
  int i, j, k, k1;
  double s, p, q, r;
  double *ai, *ak;
  for (k = 0, ak = a; k < n -1; ++k, ak += nn) {
    k1 = k + 1;
    p = ak[k];
    for (j = k1; j < n; ++j)
      ak[j] /= p;
    r = b[k] /= p;
    for (i = k1, ai = ak + nn; i < n; ++i, ai += nn) {
      q = ai[k];
      for (j = k1; j < n; ++j)
        ai[j] -= q * ak[j];
      b[i] -= q * r;
    }
  }
  x[n - 1] = b[n - 1] / ak[n - 1];
  for (k = n - 2, ak = a + nn * (n - 2); k >= 0; --k, ak -= nn) {
    k1 = k + 1;
    s = b[k];
    for (j = k1; j < n; ++j)
      s -= ak[j] * x[j];
    x[k] = s;
  }
}

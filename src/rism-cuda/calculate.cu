#include <iostream>
#include "rism3d.h"

void RISM3D :: calculate (double cf) {
  __global__ void kh(double *, const double * __restrict__, 
		     const double * __restrict__, const double * __restrict__,
		     double);
  __global__ void hnc(double * dtr, const double * __restrict__,
		      const double * __restrict__, const double * __restrict__,
		      double);
  __global__ void trm1mt(double2 *, const double * __restrict__, 
			 const double * __restrict__, const double * __restrict__, 
			 double);
  __global__ void mqvfk(double2 *, const double2 * __restrict__, double);
  __global__ void oz(double2 *, const double2 * __restrict__, 
		     const double * __restrict__, int);
  __global__ void tr(double2 *, double *, const double2 * __restrict__);

  int ng = ce -> ngrid;

  if (clos == 0) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      kh <<< g, b >>> (dtr + (iv * ng), dt + (iv * ng), 
		       du + (iv * ng), de, sv -> qv[iv] * cf);
    }
  } else if (clos == 1) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      hnc <<< g, b >>> (dtr + (iv * ng), dt + (iv * ng), 
			du + (iv * ng),  de, sv -> qv[iv] * cf);
    }
  } 

  for (int iv = 0; iv < sv -> natv; ++iv) {
    trm1mt <<< g, b >>> (dguv + (iv * ng), dtr + (iv * ng),
			  dt + (iv * ng), dfr, sv -> qv[iv] * cf);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dguv + (iv * ng), - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    mqvfk <<< g, b >>> (dguv + (iv * ng), dfk, sv -> qv[iv] * cf);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    oz <<< g, b >>> (dhuv + (iv * ng), dguv,
		      sv -> dx + (iv * sv -> natv * ng), sv -> natv);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    tr <<< g, b >>> (dguv + (iv * ng), dtr + (iv * ng), dhuv + (iv * ng));
  }
} 

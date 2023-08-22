#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include "rism3d.h"

template <typename T> struct square {
  __host__ __device__ T operator()(const T &x) const { 
    return x * x;
  }
};


double RISM3D :: cal_rms () {
  square<double> uop;
  thrust::plus<double> bop;
  thrust::device_ptr<double> dtr_ptr(dtr);

  double rms = thrust::transform_reduce(dtr_ptr, dtr_ptr 
				+ sv -> natv * ce -> ngrid, uop, 0.0, bop);
  rms = sqrt (rms / (ce -> ngrid * sv -> natv));
  return rms;
}

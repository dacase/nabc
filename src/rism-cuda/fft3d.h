#ifndef FFT3D_H
#define FFT3D_H
#include <cufft.h>
#include "cell.h"
using namespace std;

class FFT3D {
 public:
  FFT3D () {}
  ~FFT3D () {}
  void initialize (Cell *);
  void execute (double2 *, int);
 private:
  cufftHandle plan;
  dim3 g, b;
  double2 * dkf;
  int * dir;
  double volf, volb;
  int ngrid;
};

#endif 

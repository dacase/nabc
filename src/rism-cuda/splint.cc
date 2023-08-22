#include <vector>
using namespace std;

double splint (double * & xp, double * &  yp, 
	       vector <double *> & coe, int np, double x) {
  int klo = 0;
  int khi = np - 1;

  while (khi - klo > 1) {
    int k = (khi + klo) >> 1;
    if (xp[k] > x) 
      khi = k;
    else 
      klo = k;
  }
  int k = klo;

  double t = x - xp[k];
  return (((coe[2][k] * t + coe[1][k]) * t + coe[0][k]) * t + yp[k]);
}

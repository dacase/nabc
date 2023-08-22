#ifndef SOLVENT_H
#define SOLVENT_H
#include <stdlib.h>
#include <vector>
#include <string>
using namespace std;

class Solvent {
 public:
  Solvent () {}
  ~Solvent ();
  void read (string);
  void spline (vector <double> &, int * &, int, int);
  vector <vector <double *> > xvva;
  double * dx;
  double * rhov;
  double * qv;
  double * sigv;
  double * epsv;
  double temper;
  double xikt;
  int natv;
 private:
  vector <vector <double *> > xvv;
  double * ttab;
  int * multv;
  int ntab;
};

#endif

#include <vector>
#include <complex>
using namespace std;

void calloc2D (vector <complex <double> *> & array2, int r1, int r2) {
  for (int r = 0 ; r < r1 ; ++r) {
    complex <double> * temp1 = new complex <double>[r2];
    array2.push_back(temp1);
  }
}

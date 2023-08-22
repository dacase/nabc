#include <vector>
using namespace std;

void alloc3D (vector < vector <double* > > & array3, int r1, int r2, int r3) {
  for (int r = 0 ; r < r1 ; ++r) {
    vector <double *> temp2;
    for (int r = 0 ; r < r2 ; ++r) {
      double * temp1 = new double[r3];
      temp2 . push_back (temp1) ;
    }
    array3 . push_back (temp2) ;
  }
}

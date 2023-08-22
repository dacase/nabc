#include <vector>
using namespace std;

void alloc2D (vector < double * > & array2, int r1, int r2) {
  for (int r = 0 ; r < r1 ; ++r) {
    double * temp1 = new double[r2];
    array2 . push_back(temp1);
  }
}

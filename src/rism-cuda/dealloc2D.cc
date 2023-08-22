#include <vector>
using namespace std;

void dealloc2D (vector <double *> & array2) {
  for (int r = 0; r < array2.size(); ++r) {
    delete[] array2[r];
  }
  array2.clear();
}

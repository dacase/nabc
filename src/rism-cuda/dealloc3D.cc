#include <vector>
using namespace std;

void dealloc3D (vector < vector < double *> > & array3) {
  for (int r1 = 0 ; r1 < array3.size(); ++r1) {
    for (int r2 = 0 ; r2 < array3[r1].size(); ++r2) {
      delete[] array3[r1][r2];
    }
    array3[r1].clear();
  }
  array3.clear();
}

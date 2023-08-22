#include <vector>
using namespace std;

void spline (double * & xp, double * & yp, int np, vector <double *> & coe) {
  double a[np];
  double b[np];
  double c[np];
  double d[np];
  double s[np];

  int n = np - 2;

  for (int k = 1; k <= n; ++k) {
    double hl = xp[k] - xp[k - 1];
    double hr = xp[k + 1] - xp[k];
    a[k] = hl / (hl + hr);
    b[k] = 1 - a[k];
    double c1 = (yp[k + 1] - yp[k]) / hr;
    double c2 = (yp[k] - yp[k-1]) / hl;
    c[k] = 6.0 * ((c1 - c2) / (hl + hr));
    d[k] = 2.0;
  }

  for (int k = 1; k < n; ++k) {
    b[k] = b[k] / d[k];
    c[k] = c[k] / d[k];
    d[k + 1] = d[k + 1] - a[k + 1] * b[k];
    c[k + 1] = c[k + 1] - a[k + 1] * c[k];
  }
  s[n] = c[n] / d[n];

  for (int k = n - 1; k >= 1; --k) {
    s[k] = c[k] - b[k] * s[k + 1];
  }

  s[0] = s[n + 1] = 0.0;
  for (int k = 0; k <= n; ++k) {
    double h = xp[k + 1] - xp[k];
    coe[1][k] = s[k] / 2.0;
    coe[2][k] = (s[k + 1] - s[k]) / (6.0 * h);
    coe[0][k] = (yp[k + 1] - yp[k]) / h - coe[1][k] * h - coe[2][k] * h * h;
  }
}

void index(double * & a, int * & b, int n) {
  int i, j, k;
  double x;
    
  for (i = 0; i < n; ++i) {
    b[i] = i; 
  }

  for (k = n / 2 - 1; k >= 0; k--) {
    i = k;  x = a[b[i]];
    while ((j = 2 * i + 1) < n) {
      if (j < n - 1 && a[b[j]] < a[b[j + 1]]) j++;
      if (x >= a[b[j]]) break;
      b[i] = b[j];  i = j;
    }
    b[i] = k;
  }
  --n;
  while (n > 0) {
    k = b[n] ; x = a[b[n]];  b[n] = b[0];  n--;
    i = 0;
    while ((j = 2 * i + 1) < n) {
      if (j < n - 1 && a[b[j]] < a[b[j + 1]]) j++;
      if (x >= a[b[j]]) break;
      b[i] = b[j];  i = j;
    }
    b[i] = k;
  }
}

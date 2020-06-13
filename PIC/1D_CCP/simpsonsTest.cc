#include <iostream>
#include "solverModules.h"
#include <math.h>

using namespace std;

int main(void) {
  int n = 1000001;
  double u_approx, u_exact;
  double *f = new double[n];
  double a = 0.0; 
  double b = 100.0;
  double dx = (b-a)/(n-1.0);

  double x;
  u_exact = 0.892979511416;

  for (int i = 0; i < n; ++i) {
    x = a+dx*i;
    f[i] = pow(x,1.0/3.0)*exp(-x);
  }

  u_approx = simpsons(f, a, b,  n);
  cout << u_approx << " " << u_exact << endl;
  cout << abs(u_approx-u_exact) << endl;
}

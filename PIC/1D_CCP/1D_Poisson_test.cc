#include <math.h>
#include "1D_Poisson.h"
#include <fstream>


int main(void) {

for (int nn = 16 + 1; nn < 512 + 2; nn = 2*nn - 1) {
  //int nn = 201;
  double dx = 1.0/(nn-1.0);

  // Electrode info
  // Bias electrode left, grounded electrode right
  int *elec_range = new int[4]; 
  elec_range[0] = 0; // Left bias electrode
  elec_range[1] = 0; // right side bias
  elec_range[2] = nn-1; // left side ground
  elec_range[3] = nn-1; // right side ground
  double phi_left = 1.0;
  double phi_right = -1.0;

  double *f = new double[nn];
  double *u_exact = new double[nn];
  double *u_approx = new double[nn];

  for (int i = 0; i < nn; ++i) {
    u_approx[i] = 0.0;
    u_exact[i] = cos(M_PI*i*dx);
    f[i] = -M_PI*M_PI*u_exact[i];
  }

  //jacobi_Update(u_approx, f, elec_range, phi_left, phi_right,
	//	  dx, nn, 5e5, 1e-6);

  triDiagSolver(u_approx, f, elec_range, phi_left, phi_right, dx, nn);
  double residual = 0.0;
  double norm_exact = 0.0;
  for (int i = 0; i < nn; ++i) {
    residual += pow((u_approx[i]-u_exact[i]),2.0);
    norm_exact += pow(u_exact[i],2.0);

    //cout << u_approx[i] << "   " << u_exact[i] << endl;
  }
  residual = sqrt(residual)/sqrt(norm_exact);

  cout << "nn = " << nn << " Residual = " << residual << endl;

  delete(elec_range);
  delete(f);
  delete(u_exact);
  delete(u_approx);
}
  return 0;
}


#include <math.h>
#include "solverModules.h"
#include <fstream>


int main(void) {

for (int nn = 16 + 1; nn < 512 + 2; nn = 2*nn - 1) {
  int n_cell = nn-1;
  double dx = 1.0/(n_cell);

  // Electrode info
  // Bias electrode left, grounded electrode right
  int *elec_range = new int[4]; 
  elec_range[0] = 0; // Left bias electrode
  elec_range[1] = 0; // right side bias
  elec_range[2] = nn-1; // left side ground
  elec_range[3] = nn-1; // right side ground
  double phi_left = 1.0;
  double phi_right = -1.0;

// Electric potential
  int nn_inner = elec_range[2] - elec_range[1] - 1; // Interior nodes
  double *RHS = new double[n_cell];
  double *a = new double[n_cell];
  double *b = new double[n_cell];
  double *c = new double[n_cell];

  double *u_exact = new double[n_cell];
  double *u_approx = new double[n_cell];

  for (int i = 0; i < n_cell; ++i) {
    u_approx[i] = 0.0;
    u_exact[i] = cos(M_PI*(i+0.5)*dx);
    RHS[i] = -M_PI*M_PI*u_exact[i];
  }

  /*
    for (int i = elec_range[0]; i <= elec_range[1]; ++i) {
      u_approx[i] = phi_left;
      if (i < elec_range[1]) {
        /phi_cc[i] = phi_left;
      }
    }
    for (int i = elec_range[2]; i <= elec_range[3]; ++i) {
      u_approx[i] = phi_right;
      if (i < elec_range[3]) {
        /phi_cc[i] = phi_right;
      }
    }
   */


    for (int i = 0; i < n_cell; ++i) {
      a[i] = 1.0;
      b[i] = -2.0;
      c[i] = 1.0;
      RHS[i] *= dx*dx;
    }
    a[0] = 0.0;
    b[0] -= 1.0;
  
    c[n_cell-1] = 0.0;
    b[n_cell-1] -= 1.0;

    RHS[0] -= 2.0*phi_left;
    RHS[n_cell-1] -= 2.0*phi_right;
    
    triDiagSolver(u_approx, a, b, c, 
		  RHS, n_cell);


  double residual = 0.0;
  double norm_exact = 0.0;
  for (int i = 0; i < n_cell; ++i) {
    residual += pow((u_approx[i]-u_exact[i]),2.0);
    norm_exact += pow(u_exact[i],2.0);

    //cout << u_approx[i] << "   " << u_exact[i] << endl;
  }
  residual = sqrt(residual)/sqrt(norm_exact);

  cout << "nn = " << nn << " Residual = " << residual << endl;

  delete(elec_range);
  delete(RHS);
  delete(u_exact);
  delete(u_approx);
}
  return 0;
}


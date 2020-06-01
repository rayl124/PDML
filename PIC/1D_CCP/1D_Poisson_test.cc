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

// Electric potential
  int nn_inner = elec_range[2] - elec_range[1] - 1; // Interior nodes
  double *RHS = new double[nn];
  double *a = new double[nn_inner];
  double *b = new double[nn_inner];
  double *c = new double[nn_inner];

  double *u_exact = new double[nn];
  double *u_approx = new double[nn];

  for (int i = 0; i < nn; ++i) {
    u_approx[i] = 0.0;
    u_exact[i] = cos(M_PI*i*dx);
    RHS[i] = -M_PI*M_PI*u_exact[i];
  }

    for (int i = elec_range[0]; i <= elec_range[1]; ++i) {
      u_approx[i] = phi_left;
      /*if (i < elec_range[1]) {
        /phi_cc[i] = phi_left;
      }*/
    }
    for (int i = elec_range[2]; i <= elec_range[3]; ++i) {
      u_approx[i] = phi_right;
      /*if (i < elec_range[3]) {
        /phi_cc[i] = phi_right;
      }*/
    }
    
    RHS[elec_range[1] + 1] -= phi_left/(dx*dx);
    RHS[elec_range[2] - 1] -= phi_right/(dx*dx);

    for (int i = 0; i < nn_inner; ++i) {
      a[i] = 1.0;
      b[i] = -2.0;
      c[i] = 1.0;
      RHS[elec_range[1]+1+i] *= dx*dx;
    }
    a[0] = 0.0;
    c[nn_inner-1] = 0.0;
    
    triDiagSolver(&u_approx[elec_range[1] + 1], a, b, c, 
		    &RHS[elec_range[1] + 1], nn_inner);


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
  delete(RHS);
  delete(u_exact);
  delete(u_approx);
}
  return 0;
}


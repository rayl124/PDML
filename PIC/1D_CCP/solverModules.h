#ifndef _SOLVERMODULES_H
#define _SOLVERMODULES_H

#include <iostream>
#include <algorithm>

using namespace std;

// Uses Jacobi method to solve it
// solves laplace(phi) = -rho/epsilon0 = f
void jacobi_Update(double *phi,
		double *RHS,
		int *elec_range,
		double phi_left,
		double phi_right,
		double dx,
		int nn,
		const int jacobi_max_iter,
		const double tol) {

  double *phi_new = new double[nn];
  //double *RHS0 = new double[nn];
  
  int jacobi_iter = 0;
  double residual;
  
  // Add dirichlet boundaries
  for (int i = elec_range[0]; i <= elec_range[1]; ++i) {
    phi[i] = phi_left;
  }

  for (int i = elec_range[2]; i <= elec_range[3]; ++i) {
    phi[i] = phi_right;
  }

  copy_n(phi, nn, phi_new);
  //copy_n(RHS, nn, RHS0);

  while (jacobi_iter <= jacobi_max_iter) {
    ++jacobi_iter;
    residual = 0.0;

    if (jacobi_iter%1000 == 0) {
      std::cout << "Jacobi iter = " << jacobi_iter << std::endl;
    }
/*
    // Update RHS term to include electron terms
    for (int i = 0; i < nn; ++i) {
      RHS[i] = RHS0[i] - n0*exp(q*(phi[i]-phi0)/(k*Te));
      RHS[i] = -q*RHS[i]/epsilon0;
    }
    */
    // Run Jacobi update for interior nodes
    for (int i = elec_range[1] + 1; i < elec_range[2]; ++i) {
      phi_new[i] = 0.5*(phi[i-1] + phi[i+1] - pow(dx,2.0)*RHS[i]);
    }

    // Calculate the residual
    for (int i = 0; i < nn; ++i) {
      residual += pow(phi_new[i] - phi[i], 2.0);
    }
    residual = sqrt(residual);

    // Update phi to be phi_new
    copy_n(phi_new, nn, phi);

    if (residual < tol) {
      break;
    }
    if (jacobi_iter == jacobi_max_iter) {
      cout << "Jacobi could not converge" << endl;
    }
  }
  delete(phi_new);
}

// Solves a tridiagonal matrix
// Currently used for laplace(phi) = -rho/epsilon0
void triDiagSolver(double *phi,
		double *a, double *b, double *c,
		double *d, int nn) {
 
  double *c_prime = new double[nn];
  double *d_prime = new double[nn];

  for (int i = 0; i < nn; ++i) {
    if (i == (0)) {
      c_prime[i] = c[i]/b[i];
      d_prime[i] = d[i]/b[i];
    } else {
      c_prime[i] = c[i]/(b[i]-a[i]*c_prime[i-1]);
      d_prime[i] = (d[i] - a[i]*d_prime[i-1])/
	      (b[i] - a[i]*c_prime[i-1]);
    }
  }
  for (int i = nn - 1; i > -1; --i) {
    if (i == (nn-1)) {
      phi[i] = d_prime[i];
    } else {
      phi[i] = d_prime[i] - c_prime[i]*phi[i+1];
    }
  }

  delete(c_prime);
  delete(d_prime);
}

// Solves dn/dt + d(-D dn/dx)/dx = n_dot
// Assumes D is a constant and neumann BC
// Explicit time stepping
void driftDiffusionFVExplicit(double *u, double *RHS,
		      double gamma_L, double gamma_R,
		      double D, double dx, double dt,
		      int n_cell) 
{
  double *u_new = new double[n_cell];
  double A = D*dt/(dx*dx);

  for (int i = 0; i < n_cell; ++i) {
    if (i == 0) {
      u_new[i] = RHS[i]*dt + gamma_L*dt/dx + 
		 A*u[i+1] + (1.0-A)*u[i];
    } else if (i == n_cell-1) {
      u_new[i] = RHS[i]*dt - gamma_R*dt/dx + 
		 (1.0-A)*u[i] + A*u[i-1];
    } else {
      u_new[i] = RHS[i]*dt + A*u[i+1] + 
	         (1.0-2.0*A)*u[i] + A*u[i-1];
    }
  }
  
  copy_n(u_new, n_cell, u);

  delete(u_new);
}

// Solves dn/dt + d(-D dn/dx)/dx = n_dot
// Assumes D is a constant and neumann BC
// Crank Nicholson
void driftDiffusionFV_CN(double *u, double *RHS_j, double *RHS_j1,
		      double gamma_L, double gamma_R,
		      double D, double dx, double dt,
		      int n_cell) 
{
  double *a = new double[n_cell];
  double *b = new double[n_cell];
  double *c = new double[n_cell];
  double *d = new double[n_cell];

  double A = D*dt/(2.0*dx*dx);

  for (int i = 0; i < n_cell; ++i) {
    if (i == 0) {
      a[i] = 0.0;
      b[i] = 1.0+A;
      c[i] = -A;
      d[i] = A*u[i+1] + (1.0-A)*u[i] + gamma_L*dt/dx +
	      dt*0.5*(RHS_j1[i]+RHS_j[i]);
    }
    else if (i == n_cell-1) {
      a[i] = A;
      b[i] = 1.0-A;
      c[i] = 0.0;
      d[i] = A*u[i-1] + (1.0-A)*u[i]  - gamma_R*dt/dx +
	      dt*0.5*(RHS_j1[i]+RHS_j[i]);
    }
    else {
      a[i] = -A;
      b[i] = 1.0+2.0*A;
      c[i] = -A;
      d[i] =  A*u[i+1] + (1-2.0*A)*u[i] + A*u[i-1] +
	      dt*0.5*(RHS_j1[i]+RHS_j[i]);
    }
  }
  triDiagSolver(u, a, b, c, d, n_cell);

}
#endif

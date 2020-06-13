#include <math.h>
#include <iostream>
#include <algorithm>
#include "solverModules.h"
#include <fstream>

using namespace std;

int main(void) {
  int n_cell = 64;
  //int ts = 3000;

  double t_end = 1.0;  
  double x_end = 2.0;


  double residual = 0.0;
  double norm = 0.0;
  double dx = x_end/(n_cell);
  //double dt = t_end/ts;
  double D = 1e-7;
  double x, t, gamma_L, gamma_R;

  double CFL = 0.1;
  double dt =  CFL*dx/abs(D);
  if (dt > t_end) {
    dt = t_end/10000.0;
  }
  int ts = round(t_end/dt);
  ts = ts/10;

  double *u_k = new double[n_cell];
  double *u_exact = new double[n_cell*ts];

  double *RHS = new double[n_cell*ts];
  double *f = new double[(n_cell+1)*ts];
  
  ofstream testFile("Results/experimentFiles/sin3.txt");
  testFile << "x / t / u_approx / u_exact" << endl;
  
  for (int k = 0; k < ts; ++k) {
    t = dt*k;
    for (int i = 0; i < n_cell; ++i) {
      x = dx*(i+0.5);	    
      u_exact[i + k*n_cell] = cos(M_PI*t)*cos(M_PI*x)+100.0;
      RHS[i + k*n_cell] = -M_PI*sin(M_PI*t)*cos(M_PI*x) +
	      		  M_PI*M_PI*(3.0e2*x)*(100.0*cos(M_PI*x) + cos(M_PI*t))*cos(M_PI*t)/
			  pow(cos(M_PI*t)*cos(M_PI*x)+100.0,2.0);
      
      if (k == 0) {
        testFile <<  dx*(i+0.5) << " " << dt*k << " " <<  u_exact[i+k*n_cell] << " " << u_exact[i+k*n_cell] << endl;
      }
    }
  }

  copy_n(u_exact, n_cell, u_k);

  for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < n_cell+1; ++i) {
	x = i*dx;
	f[i + k*n_cell] = 3.0e2*x;
      }
  }
  
  double n_L, n_R;

  for (int k = 1; k < ts; ++k) {
    residual = 0.0;
    norm = 0.0;
    t = dt*k;
    double x_ghost1 = -0.5*dx;
    double x_ghost2 = (x_end+0.5*dx);

    n_L = cos(M_PI*(t-dt))*cos(M_PI*x_ghost1)+100.0;
    n_R = cos(M_PI*(t-dt))*cos(M_PI*x_ghost2)+100.0;

    driftDiffusionFVExplicit2(u_k, &RHS[(k-1)*n_cell], &f[(k-1)*(n_cell+1)],
		    n_L, n_R, dx, dt, n_cell);

    //driftDiffusionFVExplicit(u_k, &RHS[(k-1)*n_cell], gamma_L, gamma_R,
    //		    D, dx, dt, n_cell);
    //
    //gamma_L = 0.5*(-D*t-D*(t-dt));
    //gamma_R = 0.5*(-D*t-D*(t-dt));
    //driftDiffusionFV_CN(u_k, &RHS[(k-1)*n_cell], &RHS[k*n_cell], gamma_L,
	//	    gamma_R, D, dx, dt, n_cell);
    for (int i = 0; i < n_cell; ++i) {
      //cout << "Exact: " << u_exact[i+k*n_cell] << ", Approx: " << u_k[i] << endl;
      residual += pow((u_k[i] - u_exact[i+k*n_cell]),2.0);
      norm += pow(u_exact[i+k*n_cell],2.0);

      testFile <<  dx*(i+0.5) << " " << dt*k << " " <<  u_k[i] << " " << u_exact[i+k*n_cell] << endl;
    }

    cout << "Time: " << dt*k << " Error is: " << sqrt(residual)/sqrt(norm) << endl;
   
  }


  
  delete(u_k);
  delete(u_exact);
  delete(RHS);
  return 1;
}

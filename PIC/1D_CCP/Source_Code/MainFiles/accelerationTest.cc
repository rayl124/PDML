#include <math.h>
#include <iostream>
#include <algorithm>
#include "solverModules.h"
#include <fstream>

using namespace std;

int main(void) {
  int n_cell = 128;
  //int ts = 3000;

  double t_end = 20.0;  
  double x_end = 320.0;


  double residual = 0.0;
  double norm = 0.0;
  double dx = x_end/(n_cell);
  //double dt = t_end/ts;
  double D = 1.0;
  double x, t, gamma_L, gamma_R;

  double CFL = 0.1;
  double dt =  CFL*dx/abs(D);
  //if (dt > t_end) {
  dt = 0.005;
  //}
  int ts = round(t_end/dt);

  double *u_k = new double[n_cell];
  double *u_exact = new double[n_cell*ts];
  double *u_ss = new double[n_cell];

  double *RHS = new double[n_cell*ts];
  double *f = new double[(n_cell+1)*ts];
  
  ofstream testFile("test.txt");
  testFile << "x / t / u_approx / u_exact" << endl;
  
  for (double C = 1.0; C <= 50.0; ++C) {
    for (int k = 0; k < ts; ++k) {
      t = dt*(k+1);
      for (int i = 0; i < n_cell; ++i) {
        x = dx*(i+0.5);	    
        u_exact[i + k*n_cell] = 1.0e6*exp(-t) + pow(x, 2.0) + 10.0;
        RHS[i + k*n_cell] = -1.0e6*exp(-t) - 2.0*D;
        u_ss[i] = pow(x,2.0) + 10.0;
      
        if (k == 0) {
          testFile <<  dx*(i+0.5) << " " << dt*k << " " <<  u_exact[i+k*n_cell] << " " << u_exact[i+k*n_cell] << endl;
        }
      }
    }

    copy_n(u_exact, n_cell, u_k);
  
    double n_L, n_R, flux_L, flux_R;
    int show = 0;
    //int ss_count = 0;
    for (int k = 0; k < ts*(n_cell+1); ++k) {
      f[k] = D;
    }
  
    int ss_count = 0;
    for (int k = 1; k < ts; ++k) {
      residual = 0.0;
      norm = 0.0;
      t = dt*(k+1);
      double x_ghost1 = -0.5*dx;
      double x_ghost2 = (x_end+0.5*dx);
    
      n_L = 1.0e6*exp(-t) + pow(x_ghost1, 2.0) + 10.0;
      n_R = 1.0e6*exp(-t) + pow(x_ghost2, 2.0) + 10.0;

      if (k < 0.05*ts && k > 0.01*ts) {
        //if (show == 0) {cout << "Acceleration" << endl; show = 1;} 
        driftDiffusionFVExplicit2(u_k, &RHS[(k-1)*n_cell], &f[(k-1)*(n_cell+1)], n_L, n_R, dx, C*dt, n_cell, &flux_L, &flux_R);
        k += (C-1);
      } else {
        driftDiffusionFVExplicit2(u_k, &RHS[(k-1)*n_cell], &f[(k-1)*(n_cell+1)], n_L, n_R, dx, 1.0*dt, n_cell, &flux_L, &flux_R);
      }

      for (int i = 0; i < n_cell; ++i) {
        //cout << "Exact: " << u_exact[i+k*n_cell] << ", Approx: " << u_k[i] << endl;
        residual += pow((u_k[i] - u_ss[i]),2.0);
        norm += pow(u_ss[i],2.0);
      
        if (k%100 == 0) {
          testFile <<  dx*(i+0.5) << " " << dt*k << " " <<  u_k[i] << " " << u_exact[i+k*n_cell] << " " << u_ss[i] << endl;
        }
      }

      if (k%100 == 0) {
        //cout << "Time: " << dt*k << " Error is: " << sqrt(residual)/sqrt(norm) << endl;
      }

      if (sqrt(residual)/sqrt(norm) < 0.10) {++ss_count;} 
      else {ss_count = 0;}
      if (ss_count == 100) {cout << "C = " << C << ", Exited at time " << dt*k <<endl; 
        for (int i = 0; i < n_cell; ++i) {
          testFile <<  dx*(i+0.5) << " " << dt*k << " " <<  u_k[i] << " " << u_exact[i+k*n_cell] << " " << u_ss[i] << endl;
        }
        break;
      }
    }
    if (ss_count != 100) {cout << "C = " << C << ", no convergence" << endl;}
  } 

  delete(u_ss);
  delete(u_k);
  delete(u_exact);
  delete(RHS);
  return 1;
}

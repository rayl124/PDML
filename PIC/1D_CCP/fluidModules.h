#include "solverModules.h"
#include <iostream>
#include <math.h>


double getMomentumCS(double beta) {
  double *c_HEKhrapak = new double[4]; // Coefficients from Khrapak(2014)
  double *c_LEKhrapak = new double[4];

  c_HEKhrapak[0] = -0.692;
  c_HEKhrapak[1] = 9.594;
  c_HEKhrapak[2] = -8.284;
  c_HEKhrapak[3] = -2.355;
  c_LEKhrapak[0] = -0.019;
  c_LEKhrapak[1] = 0.038;
  c_LEKhrapak[2] = -0.049;
  c_LEKhrapak[3] = 0.015;

  double f = 1.0;

  double sigma;

  if (beta < 0.506) {
    for (int i = 0; i < 4; ++i) {
      f += c_HEKhrapak[i]*pow(beta,(i+1));
    }
    sigma = 4.507*pow(beta,1.0/6.0)*f;
  } else {
    for (int i = 0; i < 4; ++i) {
      f += c_LEKhrapak[i]*pow(beta,-(i+1));
    }
    sigma = 9.866*pow(beta,1.0/3.0)*f;
  }

  return sigma;
}


/*
// Epsilon in K
double diffusionCoeff(double T, double epsilon, double d,
		         double P, double m) 
{
  int n = 1000001;
  double k_B = 1.381e-23;

  double T_star = T/(epsilon); // Unitless

  double int1, int2, x1, x2;
  double bound = pow(1.1012*T_star,-1.0);

  double *f1 = new double[n];
  double *f2 = new double[n];
  double dx1 = (200.0-bound)/(n-1);
  double dx2 = (bound)/(n-1);

  double beta1, beta2;

  for (int i = 0; i < n; ++i) {
    x1 = bound + dx1*i;
    x2 = dx2*i + 0.00001;

    beta1 = pow(2.0*T_star*x1,-1.0);
    beta2 = pow(2.0*T_star*x2,-1.0);

    f1[i] = x1*x1*exp(-x1)*g:wqetMomentumCS(beta1);

    f2[i] = x2*x2*exp(-x2)*getMomentumCS(beta2);
  }

  int1 = simpsons(f1, bound, 100.0, n);
  int2 = simpsons(f2, 0.0, bound, n);

  double Omega = 0.5*(int1 + int2);
  double D = 3.0/sqrt(M_PI)/(8.0);
  //double D = 3.0/(16.0*sqrt(M_PI));

  D *= 100.0*100.0*pow(T, 1.5)/(P*sqrt(m)*d*d*Omega);

  delete(f1);
  delete(f2);
  return D;
} */

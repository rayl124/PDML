#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(void) {
  double const m = 1.0;
  double const q = 1.0;
  double const omega_p = 1.0/(1.0*2.0*M_PI);

  // Size of time history
  int n = 256;
  double *t = new double[n];
  double *t_half = new double[n];
  double *x_FD = new double[n];
  double *x_LF = new double[n];
  double *x_Exact = new double[n];
  double *v_FD = new double[n];
  double *v_LF = new double[n];
  double *v_Exact = new double[n];

  double E;
  
  double t_end = 25.0;
  double dt = t_end/(n-1);

  // Boundary conditions
  x_Exact[0] = 1.0;
  x_FD[0] = 1.0;
  x_LF[0] = 1.0;
  
  v_Exact[0] = 0.0;
  v_FD[0] = 0.0;

  E = -(m/q)*pow(omega_p, 2.0)*x_LF[0];
  v_LF[0] = 0.0 - 0.5*(q/m)*E*dt;
  
  t[0] = 0.0;
  t_half[0] = -0.5*dt;


  ofstream MyFile("vel_int.txt");
  MyFile << "Time Time_0.5 x_Exact v_Exact x_FD v_FD x_LF v_LF" << std::endl;
  MyFile << t[0] << " " << t_half[0] << " " << x_Exact[0] << " " << v_Exact[0] << " ";
  MyFile << x_FD[0] << " " << v_FD[0] << " " << x_LF[0] << " " << v_LF[0];
  MyFile << std::endl;

  // Time vector and exact solution
  for (int i = 1; i < n; ++i){
    // Time Vectors
    t[i] = i*dt;
    t_half[i] = t[i] - 0.5*dt;

    // Exact solutions
    x_Exact[i] = cos(omega_p * t[i]);
    v_Exact[i] = -omega_p * sin(omega_p*t[i]);

    // Forward difference
    E = -(m/q)*pow(omega_p, 2.0)*x_FD[i-1];
    v_FD[i] = v_FD[i-1] + q/m*E*dt;
    x_FD[i] = x_FD[i-1] + v_FD[i-1]*dt;

    // Leapfrog method
    E = -(m/q)*pow(omega_p, 2.0)*x_LF[i-1];
    v_LF[i] = v_LF[i-1] + (q/m)*E*dt;
    x_LF[i] = x_LF[i-1] + v_LF[i];

    MyFile << t[i] << " " << t_half[i] << " " << x_Exact[i] << " " << v_Exact[i] << " ";
    MyFile << x_FD[i] << " " << v_FD[i] << " " << x_LF[i] << " " << v_LF[i];
    MyFile << std::endl;
  }

  delete t;
  delete t_half;
  delete x_FD;  
  delete x_LF;   
  delete x_Exact;   
  delete v_FD; 
  delete v_LF;
  delete v_Exact;

  return 0;
}

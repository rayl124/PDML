#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(void) {
  double const m = 1.0;
  double const q = 1.0;
  double const omega_p = 2.0*M_PI/5.0;

  // Size of time history
  int n = 512;
  double t;
  double x_FD_n1;
  double x_FD_n;
  double x_LF_n1;
  double x_LF_n;
  double x_Exact;

  double v_FD_n1;
  double v_FD_n;
  double v_LF_n1;
  double v_LF_n;
  double v_Exact;
 
  double t_end = 25.0;
  double dt = t_end/(n-1);

  double resid_FD_x = 0.0;
  double resid_FD_v = 0.0;
  double resid_LF_x = 0.0;
  double resid_LF_v = 0.0;

  // Boundary conditions
  x_Exact = 1.0;
  x_FD_n = 1.0;
  x_LF_n = 1.0;
  
  v_Exact= 0.0;
  v_FD_n = 0.0;
  v_LF_n = 0.0 + 0.5*pow(omega_p, 2.0)*x_LF_n*dt;
  
  t = 0.0;

  ofstream MyFile("vel_int.txt");
  MyFile << "Time x_Exact v_Exact x_FD v_FD x_LF v_LF" << std::endl;
  MyFile << t << " " << x_Exact << " " << v_Exact << " ";
  MyFile << x_FD_n << " " << v_FD_n << " " << x_LF_n << " " << v_LF_n;
  MyFile << std::endl;

  // Time vector and exact solution
  for (int i = 1; i < n; ++i){
    // Time Vectors
    t = t+dt;

    // Exact solutions
    x_Exact = cos(omega_p * t);
    v_Exact = -omega_p * sin(omega_p*t);

    // Forward difference
    v_FD_n1 = v_FD_n - pow(omega_p, 2.0)*x_FD_n*dt;
    x_FD_n1 = x_FD_n + v_FD_n*dt;

    // Leapfrog method
    v_LF_n1 = v_LF_n - pow(omega_p, 2.0)*x_LF_n*dt;
    x_LF_n1 = x_LF_n + v_LF_n1*dt;

    // x_n+1 becomes x_n, v_n+1 becomes v_n
    x_FD_n = x_FD_n1;
    v_FD_n = v_FD_n1;
    x_LF_n = x_LF_n1;
    v_LF_n = v_LF_n1;

    // Error Calculation
    resid_FD_x += (x_FD_n-x_Exact)*(x_FD_n-x_Exact);
    resid_FD_v += (v_FD_n-v_Exact)*(v_FD_n-v_Exact);
    resid_LF_x += (x_LF_n-x_Exact)*(x_LF_n-x_Exact);
    resid_LF_v += (v_LF_n-v_Exact)*(v_LF_n-v_Exact);


    MyFile << t << " " << x_Exact << " " << v_Exact << " ";
    MyFile << x_FD_n << " " << v_FD_n << " " << x_LF_n << " " << v_LF_n;
    MyFile << std::endl;
  }

  resid_FD_x = sqrt(resid_FD_x);
  resid_FD_v = sqrt(resid_FD_v);
  resid_LF_x = sqrt(resid_LF_x);
  resid_LF_v = sqrt(resid_LF_v);
  
  MyFile << resid_FD_x << " " << resid_FD_v << " " << resid_LF_x << " ";
  MyFile << resid_LF_v << std::endl;
  MyFile.close();

  return 0;
}

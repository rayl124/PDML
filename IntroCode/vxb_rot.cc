#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
int main(void) {
  int N = 512;

  double t;
  double *x_Exact = new double[2];

  // Index 0 = x direction, index 1 = y direction
  double *x_FD_n = new double[2]; // Forward difference
  double *x_FD_n1 = new double[2];
  double *x_TI_n = new double[2]; // Tajima Implicit
  double *x_TI_n1 = new double[2];
  double *x_TE_n = new double[2]; // Tajima Explicit
  double *x_TE_n1 = new double[2];
  double *x_B_n = new double[2]; // Boris
  double *x_B_n1 = new double[2]; 

  double *v_FD_n = new double[2];
  double *v_FD_n1 = new double[2];
  double *v_TI_n = new double[2]; // Represents v_n-1/2
  double *v_TI_n1 = new double[2]; // v_n+1/2
  double *v_TE_n = new double[2];
  double *v_TE_n1 = new double[2];
  double *v_B_plus = new double[2]; // No electric field, v- = v_n-1/2
  double *v_B_minus = new double[2]; // No electric field, v+ = v_n+1/2
  double *v_prime = new double[2];

  // Properties of the problem
  double B = 0.01; // Magnetic field T
  double m = 9.10938356e-31; // Electron mass kg
  double q = -1.602176634e-19; // Electron charge C
  double v_tan = 1.0e5; // Initial tangential velocity m/s
  double r_L =  m*v_tan/(abs(q)*B); // Larmor radius m
  double omega_c = q*B/m; // Cyclotron frequency rad/s

  double t_end = 8*2*M_PI/abs(omega_c); // 5 periods
  double dt = t_end/(N-1.0);
  double epsilon = omega_c*dt/(2.0); // Scaling factor
  double epsilon2 = epsilon*epsilon; // Epsilon^2
  double det = 1.0/(1.0 + epsilon2); // det(I-Repsilon)

  // Boundaries
  x_FD_n[0] = 0.0;
  x_FD_n[1] = -r_L;
  x_TI_n[0] = 0.0;
  x_TI_n[1] = -r_L;
  x_TE_n[0] = 0.0;
  x_TE_n[1] = -r_L;
  x_B_n[0] = 0.0;
  x_B_n[1] = -r_L;

  v_FD_n[0] = v_tan;
  v_FD_n[1] = 0.0;

  // Shift by half a time step
  v_TI_n[0] = v_tan; // v_tan + (q/m)*B*0.0*dt;
  v_TI_n[1] = 0.0 - 0.5*(q/m)*B*v_tan*dt;

  v_TE_n[0] = v_tan;
  v_TE_n[1] = 0.0 - 0.5*(q/m)*B*v_tan*dt;
 
  v_B_minus[0] = v_tan;
  v_B_minus[1] = 0.0 - 0.5*(q/m)*B*v_tan*dt;

  x_Exact[0] = 0.0;
  x_Exact[1] = -r_L;
  
  ofstream MyFile("vxb_rot.txt");
  MyFile << "Time x_Exact y_Exact x_FD y_FD x_TI y_TI" << std::endl;

  MyFile << 0.0 << " " << x_Exact[0] << " " <<  x_Exact[1] << " " << x_FD_n[0];
  MyFile << " " << x_FD_n[1] << " " << x_TI_n[0] << " " <<  x_TI_n[1];
  MyFile << " " << x_TE_n[0] << " " << x_TE_n[1] << " " << x_B_n[0] << " " << x_B_n[1];

  MyFile << std::endl;

  for(int i = 1; i < N; ++i) {
    t += dt;

    // Exact solution
    x_Exact[0] = m*v_tan/(q*B)*sin(omega_c*t);
    x_Exact[1] = m*v_tan/(q*B)*cos(omega_c*t);
    

    // Forward difference
    v_FD_n1[0] = v_FD_n[0] + (q/m)*B*v_FD_n[1]*dt;
    v_FD_n1[1] = v_FD_n[1] - (q/m)*B*v_FD_n[0]*dt;

    x_FD_n1[0] = x_FD_n[0] + v_FD_n1[0]*dt;
    x_FD_n1[1] = x_FD_n[1] + v_FD_n1[1]*dt;

    // Tajima Implicit
    v_TI_n1[0] = det*((1.0-epsilon2) * v_TI_n[0] + (epsilon + epsilon) *
		    v_TI_n[1]);
    v_TI_n1[1] = det*((-epsilon - epsilon) * v_TI_n[0] + (1.0 - epsilon2) *
		    v_TI_n[1]);

    x_TI_n1[0] = x_TI_n[0] + v_TI_n1[0]*dt;
    x_TI_n1[1] = x_TI_n[1] + v_TI_n1[1]*dt;

    // Tajima Explicit
    v_TE_n1[0] = v_TE_n[0] + epsilon*v_TE_n[1];
    v_TE_n1[1] = -epsilon*v_TE_n[0] + v_TE_n[1];

    x_TE_n1[0] = x_TE_n[0] + v_TE_n1[0]*dt;
    x_TE_n1[1] = x_TE_n[1] + v_TE_n1[1]*dt;

    // Boris
    v_prime[0] = v_B_minus[0] + epsilon*v_B_minus[1];
    v_prime[1] = v_B_minus[1] - epsilon*v_B_minus[0];
    v_B_plus[0] = v_B_minus[0] + (2*epsilon/(1+epsilon2)) * v_prime[1];
    v_B_plus[1] = v_B_minus[1] - (2*epsilon/(1+epsilon2)) * v_prime[0];

    x_B_n1[0] = x_B_n[0] + v_B_plus[0]*dt;
    x_B_n1[1] = x_B_n[1] + v_B_plus[1]*dt; 

    // Update n+1 to be n
    x_FD_n[0] = x_FD_n1[0];
    x_FD_n[1] = x_FD_n1[1];
    v_FD_n[0] = v_FD_n1[0];
    v_FD_n[1] = v_FD_n1[1];
    

    x_TI_n[0] = x_TI_n1[0];
    x_TI_n[1] = x_TI_n1[1];
    v_TI_n[0] = v_TI_n1[0];
    v_TI_n[1] = v_TI_n1[1];

    x_TE_n[0] = x_TE_n1[0];
    x_TE_n[1] = x_TE_n1[1];
    v_TE_n[0] = v_TE_n1[0];
    v_TE_n[1] = v_TE_n1[1];

    x_B_n[0] = x_B_n1[0];
    x_B_n[1] = x_B_n1[1];
    v_B_minus[0] = v_B_plus[0];
    v_B_minus[1] = v_B_plus[1];


    MyFile << t << " " << x_Exact[0] << " " <<  x_Exact[1] << " " << x_FD_n[0];
    MyFile << " " << x_FD_n[1] << " " << x_TI_n[0] << " " <<  x_TI_n[1];
    MyFile << " " << x_TE_n[0] << " " << x_TE_n[1] << " " << x_B_n[0] << " " << x_B_n[1];
    MyFile << std::endl;
  }

  MyFile.close();

  delete x_Exact;
  delete x_FD_n1;
  delete x_FD_n;
  delete x_TI_n1;
  delete x_TI_n;
  delete x_TE_n1;
  delete x_TE_n;
  delete x_B_n1;
  delete x_B_n;

  delete v_FD_n1;
  delete v_FD_n;
  delete v_TI_n1;
  delete v_TI_n;
  delete v_TE_n1;
  delete v_TE_n;
  delete v_B_plus;
  delete v_prime;
  delete v_B_minus;

  return 0;
}

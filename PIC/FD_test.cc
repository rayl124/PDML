#include <math.h>
#include "NumericalSchemes.h"
#include <fstream>
#include <iostream>

using namespace std;

double test_u(double x, double y) {
  return sin(M_PI * x)*sin(M_PI * y);
}

double f_rhs(double x, double y) {
  return 2.0*M_PI*M_PI*test_u(x, y);
}

int main(void) {
  for (int interval = 8; interval < 128 + 1; interval *= 2) { 
    // Number of intervals in each direction
    int M = interval;
    int N = interval;

    double *u_old = new double[(M+1)*(N+1)];
    double *u_new = new double[(M+1)*(N+1)];
    double *f = new double[(M+1)*(N+1)];
    double u_ij;

    double dx = 1.0/M;
    double dy = 1.0/N;

    double residual_Jacobi = 0.0;
    double residual_Exact = 0.0;
    double exact_norm = 0.0;

    // Initial Guess
    for (int i = 0; i < M+1; ++i) {
      for (int j = 0; j < N + 1; ++j) {
        u_old[ij_to_index(i, j, N)] = 0.0;
	f[ij_to_index(i, j, N)] = f_rhs(i * dx, j * dy);
      }
    }

    // Iteration parameters for Jaocbi method
    int iter = 0;
    int max_iter = 1e6;
    double tol = 1e-6;

    // Use Dirichlet boundary conditions
    // Boundaries = 0
    // Run Jacobi iterations until potential field converges
    while(iter < max_iter) {
      ++iter;
      PoissonInteriorFDJacobi(u_old, u_new, f, dx, dy, M, N);
      
      residual_Jacobi = 0.0;

      for (int i = 1; i < M; ++i) {
        for (int j = 1; j < N; ++j) {
          u_ij = u_new[ij_to_index(i, j, N)];
	  // See how much u changes between iterations
	  residual_Jacobi += pow((u_ij - u_old[ij_to_index(i, j, N)]), 2.0); 
  
	  // Update u_old
	  u_old[ij_to_index(i, j, N)] = u_ij;

        }
      }
      residual_Jacobi = sqrt(residual_Jacobi);
      if (residual_Jacobi < tol) {
        for (int i = 1; i < M; ++i) {
	  for (int j = 1; j < N; ++j) {
            residual_Exact += pow(u_old[ij_to_index(i, j, N)] -
 				  test_u(dx * i, dy * j), 2.0);
	    exact_norm += pow(test_u(dx * i, dy * j), 2.0);
	  }
	}
      break;
     }
    }

    residual_Exact = sqrt(residual_Exact)/sqrt(exact_norm);
    
    std::cout << "N = " << N << " Residual = " << residual_Exact;
    std::cout << " Jacobi Iter = " << iter << std::endl;
    delete(u_new);
    delete(u_old);
    delete(f);
  }
  return 0;

}

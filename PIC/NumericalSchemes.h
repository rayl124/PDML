#ifndef _NUMERICALSCHEMES_H
#define _NUMERICALSCHEMES_H

/* LIST OF HELPER FUNCTIONS THAT RUN NUMERICAL SCHEMES
 * Last modified May 4, 2020
 * NOTE: some of these may have to be modified depending on the problem
 *
 * ij_to_index: converts i,j coordinates to indices in array
 * PoissonInteriorFDJacobi: Solves poisson equation for interior nodes
 * Leapfrog: Runs leapfrog scheme for velocity integration
 */

// ij_to_index
// Convert i,j coordinate to location in memory
// N = number of columns
int ij_to_index(int i, int j, int N) {
  return i*N + j;
}

// PoissonInteriorFDJacobi
// Solves interior nodes of laplace(u) = -f
// Uses Jacobi method to solve finite difference
// u_old, u_new solution arrays, f RHS array, 
// dx, dy discretizations, M, N number of rows/columns
//
void PoissonInteriorFDJacobi(double *u_old,
		double *u_new,
		double *f,
		double dx,
		double dy,
		int M,
		int N) {
  double A;
  double B;

  double scale = 1.0/(2.0/(dx*dx)+2.0/(dy*dy));

  // Run Jacobi update for interior nodes
  for (int i = 1; i < M; ++i) {
    for (int j = 1; j < N; ++j) {
      A = u_old[ij_to_index(i-1, j, N)] + u_old[ij_to_index(i+1, j, N)];
      A *= 1.0/(dx*dx);

      B = u_old[ij_to_index(i, j-1, N)] + u_old[ij_to_index(i, j+1, N)];
      B *= 1.0/(dy*dy);

      u_new[ij_to_index(i, j, N)] = (A + B + f[ij_to_index(i, j, N)])*scale;
    }
  }

  // Update exterior nodes requires boundary conditions
  // Write this in the main code
}
/*
void LeapFrogVel() {
    // Leapfrog method
    v_LF_n1 = v_LF_n - pow(omega_p, 2.0)*x_LF_n*dt;
    x_LF_n1 = x_LF_n + v_LF_n1*dt;

    // x_n+1 becomes x_n, v_n+1 becomes v_n
    x_FD_n = x_FD_n1;
    v_FD_n = v_FD_n1;
    x_LF_n = x_LF_n1;
    v_LF_n = v_LF_n1;
}
*/

#endif

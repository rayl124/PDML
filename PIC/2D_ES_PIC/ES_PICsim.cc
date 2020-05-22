#include <iostream>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <cstdlib>


using namespace std;

// ij_to_index
// Convert i,j coordinate to location in memory
// N = number of columns
int ij_to_index(int i, int j, int N) {
  return i*N + j;
}

// Get the charge density/EF weights  at the four local nodes
//
//        2 -------------- 3
//        |                |
//        |                |
//        |		   |
//        0----------------1
//
// node_index[0] = x1 coord of node 0
// node_index[1] = x2 coord of node 0
void get_Weights(double *part_x,
			double *weights,
			double *node_index,
			double dh) {
  
    node_index[0] = floor(part_x[0]/dh);
    node_index[1] = floor(part_x[1]/dh);
    double hx = (part_x[0] - (node_index[0])*dh)/dh;
    double hy = (part_x[1] - (node_index[1])*dh)/dh;

    weights[0] = (1.0-hx)*(1.0-hy);
    weights[1] = hx*(1.0-hy);
    weights[2] = (1.0-hx)*hy;
    weights[3] = hx*hy;
}



// solves laplace(phi) = -rho/epsilon0 = f
void jacobi_Update(double *phi,
		double *RHS,
		int *box_range,
		double dx1,
		double dx2,
		int nx1,
		int nx2,
		const int jacobi_max_iter,
		const double tol,
		const double n0,
		const double phi0,
		const double q,
		const double epsilon0,
		const double k,
		const double Te) {

  double *phi_new = new double[nx1*nx2];
  double *RHS0 = new double[nx1*nx2];
  copy_n(phi, nx1*nx2, phi_new);
  copy_n(RHS, nx1*nx2, RHS0);

  
  int jacobi_iter = 0;
  int index_ij;
  double A;
  double B;
  double scale;
  double residual;
  
  while (jacobi_iter <= jacobi_max_iter) {
    ++jacobi_iter;
    residual = 0.0;

    if (jacobi_iter%1000 == 0) {
      std::cout << "Jacobi iter = " << jacobi_iter << std::endl;
    }

    // Update RHS term to include electron terms
    for (int i = 0; i < nx1; ++i) {
      for (int j = 0; j < nx2; ++j) {
        RHS[ij_to_index(i, j, nx2)] = RHS0[ij_to_index(i, j, nx2)] - 
			n0*exp((phi[ij_to_index(i, j, nx2)] - phi0)/(Te));
	RHS[ij_to_index(i, j, nx2)] = -q*RHS[ij_to_index(i, j, nx2)]/epsilon0;
      }
    }

    scale = 1.0/(2.0/(dx1*dx1)+2.0/(dx2*dx2));
    // Run Jacobi update for interior nodes
    for (int i = 1; i < nx1 - 2; ++i) {
      for (int j = 2; j < nx2 - 2; ++j) {
        // Box values stay constant
        if (i >= box_range[0] && i <= box_range[1] && 
	    j >= box_range[2] && j <= box_range[3]) {
          phi_new[ij_to_index(i, j, nx2)] = phi[ij_to_index(i, j, nx2)];
        } 
        else {
          A = phi[ij_to_index(i-1, j, nx2)] + phi[ij_to_index(i+1, j, nx2)];
          A *= 1.0/(dx1*dx1);

          B = phi[ij_to_index(i, j-1, nx2)] + phi[ij_to_index(i, j+1, nx2)];
          B *= 1.0/(dx2*dx2);

          phi_new[ij_to_index(i, j, nx2)] = (A + B - 
		      RHS[ij_to_index(i, j, nx2)])*scale;
        }
      }
    }

    // Run Jacobi update for Neumann boundaries
  
    // Top
    scale = 1.0/(2.0/(dx1*dx1) + 1.0/(dx2*dx2));
    for (int i = 1; i < nx1-2; ++i) {
      A = phi[ij_to_index(i-1, nx2-2, nx2)] + phi[ij_to_index(i+1, nx2-2, nx2)];
      A *= 1.0/(dx1*dx1);

      B = phi[ij_to_index(i,nx2-3, nx2)]/(dx2*dx2);

      phi_new[ij_to_index(i, nx2-2, nx2)] = (A + B -
  		      RHS[ij_to_index(i, nx2-2, nx2)])*scale;
      phi_new[ij_to_index(i, nx2-1, nx2)] = phi_new[ij_to_index(i, nx2-2, nx2)];
    }

    // Bottom
    scale = 1.0/(2.0/(dx1*dx1) + 1.0/(dx2*dx2));
    for (int i = 1; i < nx1-2; ++i) {
      if (i >= box_range[0] && i <= box_range[1]) {
          phi_new[ij_to_index(i, 1, nx2)] = phi[ij_to_index(i, 1, nx2)];
          phi_new[ij_to_index(i, 0, nx2)] = phi[ij_to_index(i, 0, nx2)];
      } 
      else {
        A = phi[ij_to_index(i-1, 1, nx2)] + phi[ij_to_index(i+1, 1, nx2)];
        A *= 1.0/(dx1*dx1);

        B = phi[ij_to_index(i, 2, nx2)]/(dx2*dx2);

        phi_new[ij_to_index(i, 1, nx2)] = (A + B -
		     RHS[ij_to_index(i, 1, nx2)])*scale;
        phi_new[ij_to_index(i, 0, nx2)] = phi_new[ij_to_index(i, 1, nx2)];
      }
    }

    // Right
    scale = 1.0/(1.0/(dx1*dx1) + 2.0/(dx2*dx2));
    for (int j = 2; j < nx2-2; ++j) {
      A = phi[ij_to_index(nx1-3, j, nx2)]/(dx1*dx1);

      B = phi[ij_to_index(nx1-2, j-1, nx2)] + phi[ij_to_index(nx1-2, j+1, nx2)];
      B *= 1.0/(dx2*dx2);

      phi_new[ij_to_index(nx1-2, j, nx2)] = (A + B -
		    RHS[ij_to_index(nx1-2, j, nx2)])*scale;
      phi_new[ij_to_index(nx1-1, j, nx2)] = phi_new[ij_to_index(nx1-2, j, nx2)];
    }

    // Corners
    scale = 1.0/(1.0/(dx1*dx1) + 1.0/(dx2*dx2));
    phi_new[ij_to_index(nx1-2, nx2-2, nx2)] = (phi[ij_to_index(nx1-3, nx2-2, nx2)]/(dx1*dx1) +
  		  			  phi[ij_to_index(nx1-2, nx2-3, nx2)]/(dx2*dx2))*scale;
    phi_new[ij_to_index(nx1-2, 1, nx2)] = (phi[ij_to_index(nx1-3, 1, nx2)]/(dx1*dx1) +
		  			  phi[ij_to_index(nx1-2, 2, nx2)]/(dx2*dx2))*scale;

    phi_new[ij_to_index(nx1-1, nx2-2, nx2)] = phi_new[ij_to_index(nx1-2, nx2-2, nx2)];
    phi_new[ij_to_index(nx1-2, nx2-1, nx2)] = phi_new[ij_to_index(nx1-2, nx2-2, nx2)];
    phi_new[ij_to_index(nx1-1, nx2-1, nx2)] = phi_new[ij_to_index(nx1-2, nx2-2, nx2)];
    phi_new[ij_to_index(nx1-1, 1, nx2)] = phi_new[ij_to_index(nx1-2, 1, nx2)];
    phi_new[ij_to_index(nx1-2, 0, nx2)] = phi_new[ij_to_index(nx1-2, 1, nx2)];
    phi_new[ij_to_index(nx1-1, 0, nx2)] = phi_new[ij_to_index(nx1-2, 1, nx2)];


    // Calculate the residual
    for (int i = 0; i < nx1; ++i) {
      for (int j = 0; j < nx2; ++j) {
        index_ij = ij_to_index(i, j, nx2);
        residual += pow(phi_new[index_ij] - phi[index_ij], 2.0);
      }
    }
    
    residual = sqrt(residual);
    // Update phi to be phi_new
    copy_n(phi_new, nx1*nx2, phi);

    if (residual < tol) {
      break;
    }
    if (jacobi_iter == jacobi_max_iter) {
      cout << "Jacobi could not converge" << endl;
    }
  }
  delete(phi_new);
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              Main Loop                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////

int main(void) {

  // Constants
  const double epsilon0 = 8.854e-12; // Permittivity of free space
  const double q = 1.602e-19; // Elementary charge [C]
  const double k = 1.381e-23; // Boltzman constant [J/K]
  const double AMU = 1.661e-27; // Atomic Mass Unit [kg]
  const double M = 32.0*AMU; // Ion mass of O2 [kg]

  // Input settings
  double n0 = 1.0e12; // Electron density [#/m^3]
  double phi0 = 0.0; // Reference potential
  double Te = 1.0; // Electron temp [eV]
  double Ti = 0.1; // Ion temp [eV]
  double v_drift = 7.0e3; // Ion injection velocity [m/s]
  double phi_p = -5.0; // Plate potential [V]

  // Plasma parameters
  // ******************DOUBLE CHECK THESE FORMULAS*********
  double lambdaD = sqrt(epsilon0*k*Te/(q*q*n0));
  double vth = sqrt(2.0 * q * Ti/M);

  // Problem discretization
  int nx1 = 16; // # of x1 nodes
  int nx2 = 10; // # of x2 nodes
  int ts = 200; // # of time steps
  double dh = lambdaD; // length of each cell
  int np_insert = (nx2-1)*15; // Put 15 macroparticles per time step

  // Other helpful values
  int nn = nx1*nx2; // total nodes
  double dt = 0.1*dh/v_drift;
  double Lx1 = (nx1-1.0)*dh; // X Domain length
  double Lx2 = (nx2-1.0)*dh; // Y Domain length

  // Insert plate 1/3 of the way into x domain
  int *box_range = new int[4]; // 0,1 is x1 range, 2,3 is x2 range
  box_range[0] = floor(nx1/3) - 1;
  box_range[1] = floor(nx1/3) + 1;
  box_range[2] = 0;
  box_range[3] = floor(nx2/2) - 1;

  // Particle info
  int np_real = round(n0*v_drift*Lx2*dt); // Number of real particles entering
  double sp_wt = np_real/np_insert; // Real/Macroparticle
  double mp_q = 1; // Macroparticle charge
  int max_part = np_insert*ts; // max_particles if none leave during time history

  // Particle data
  double *x_part = new double[2*max_part]; // Even indices = x1 values
  double *v_part = new double[2*max_part]; // Odd indices = x2 values
  double *rho = new double[nn];

  // Various helper variabes
  // Charge density
  double *weights = new double[4];
  double *node_index = new double[2*max_part];

  // Electric potential
  double *RHS = new double[nn];
  int jacobi_max_iter = 1e5;
  double tol = 1.0e-4;
  int index_ij;

  // Electric field
  // Even indices x1, odd indices x2
  double *E_field = new double[2*nn];

  // Generate particles
  int np = 0;

  // Move Particles
  double *E_part = new double[2];
  srand(time(NULL));
  
  // Set up phi boundaries and initial values
  double *phi = new double[nn];

  for (int i = 0; i < nx1; ++i) {
    for (int j = 0; j < nx2; ++j) {
      index_ij = ij_to_index(i, j, nx2);
      
      if (i >= box_range[0] && i <= box_range[1] && 
          j >= box_range[2] && j <= box_range[3]) {
        phi[index_ij] = phi_p;
      } 
      else {
        phi[index_ij] = 0.0;
      }
    }
  }
 
  // Output files
  ofstream DataFile("Results/ESFieldData.txt");
  ofstream ParticleFile("Results/ParticleInfo.txt");
  ofstream NumFile("Results/NumParticles.txt");

  DataFile << "Iteration / Ion Density / Electric Potential / Electric Field x1 / ";
  DataFile << "Electric Field x2" << endl;

  ParticleFile << "Iteration / x1 / x2 / v1 / v2";
  ParticleFile << endl;

  NumFile << "Iteration / Number of Macroparticles" << endl;
  // Main Loop
  for (int iter = 0; iter < ts; ++iter) {     
    // Reset variables
    for (int i = 0; i<nn; ++i) {
      rho[i] = 0.0;
      E_field[i] = 0.0;
      RHS[i] = 0.0;
    }
    
    cout << "Iter: " << iter << endl;
    cout << "Computing ion charge density..." << endl;
    //Compute ion charge density 
    
    ///////////////////////////
    for (int p = 0; p < np; ++p) {
      get_Weights(&x_part[2*p], weights, &node_index[2*p], dh);
      
      int i = node_index[2*p];
      int j = node_index[2*p+1];

      rho[ij_to_index(i, j, nx2)] += weights[0];
      rho[ij_to_index(i + 1, j, nx2)] += weights[1];
      rho[ij_to_index(i, j + 1, nx2)] += weights[2];
      rho[ij_to_index(i + 1, j + 1, nx2)] += weights[3];
    }

    for (int i = 0; i < nx1; ++i) {
      for (int j = 0; j < nx2; ++j) {
        index_ij = ij_to_index(i, j, nx2);	
	rho[index_ij] *= sp_wt*mp_q/(dh*dh);

	// If on the boundary, half the volume.
	// If on the corner, quarter the volume.
	if (i == 0 || i == nx1-1) {
	  rho[index_ij] *= 2.0;
	}
	if (j == 0 || j == nx2-1) {
	  rho[index_ij] *= 2.0;
	}
	// Charge density for ions only, does not include electron term
	RHS[index_ij] = rho[index_ij];;
      }
    }
    
    cout << "Computing electric potential..." << endl;
    // Compute electric potential
    
    ////////////////////////////////
    jacobi_Update(phi, RHS, box_range, dh, dh, nx1, nx2, jacobi_max_iter, tol,
		    n0, phi0, q, epsilon0, k, Te);

    cout << "Computing electric field..." << endl;
    // Compute electric field
    
    ////////////////////////////////////////
    for (int i = 0; i < nx1; ++i) {
      for (int j = 0; j < nx2; ++j) {
        index_ij = 2*ij_to_index(i, j, nx2);

	// Left
	if (i == 0) {
	  E_field[index_ij] = -(phi[ij_to_index(1, j, nx2)] -
			  phi[ij_to_index(0, j, nx2)])/dh;
	  if (j == 0) {
            E_field[index_ij + 1] = -(phi[ij_to_index(i, j + 1, nx2)] - 
			  phi[ij_to_index(i, j, nx2)])/(dh);

	  } else if (j == nx2 - 1) {
            E_field[index_ij + 1] = -(phi[ij_to_index(i, j, nx2)] - 
			  phi[ij_to_index(i, j - 1, nx2)])/(dh);
 
	  } else {
	    E_field[index_ij + 1] = -(phi[ij_to_index(i, j + 1, nx2)] - 
			  phi[ij_to_index(i, j - 1, nx2)])/(2.0*dh);
	  }
	} // Right
       	else if (i == nx1-1) {
	  E_field[index_ij] = -(phi[ij_to_index(i, j, nx2)] -
			  phi[ij_to_index(i-1, j, nx2)])/dh;
	  if (j == 0) {
            E_field[index_ij + 1] = -(phi[ij_to_index(i, j + 1, nx2)] - 
			  phi[ij_to_index(i, j, nx2)])/(dh);

	  } else if (j == nx2 - 1) {
            E_field[index_ij + 1] = -(phi[ij_to_index(i, j, nx2)] - 
			  phi[ij_to_index(i, j - 1, nx2)])/(dh);
 
	  } else {
	    E_field[index_ij + 1] = -(phi[ij_to_index(i, j + 1, nx2)] - 
			  phi[ij_to_index(i, j - 1, nx2)])/(2.0*dh);
	  }
	} // Bottom
	else if (j == 0) {
	  E_field[index_ij] = -(phi[ij_to_index(i+1, j, nx2)] - 
			  phi[ij_to_index(i-1, j, nx2)])/(2.0*dh);
	  E_field[index_ij + 1] = -(phi[ij_to_index(i, j+1, nx2)] -
			  phi[ij_to_index(i, j,  nx2)])/(dh);
	} // Top
	else if (j == nx2-1) {
	  E_field[index_ij] = -(phi[ij_to_index(i+1, j, nx2)] - 
			  phi[ij_to_index(i-1, j, nx2)])/(2.0*dh);
	  E_field[index_ij + 1] = -(phi[ij_to_index(i, j, nx2)] -
			  phi[ij_to_index(i, j-1,  nx2)])/(dh);
	}
	else {
          E_field[index_ij] = -(phi[ij_to_index(i + 1, j, nx2)] - 
			  phi[ij_to_index(i - 1, j, nx2)])/(2.0*dh);
	  E_field[index_ij + 1] = -(phi[ij_to_index(i, j + 1, nx2)] - 
			  phi[ij_to_index(i, j - 1, nx2)])/(2.0*dh);
	}
      }
    }

    cout << "Moving Particles..." << endl;
    // Move particles
    
    ///////////////////////////////////////////
    // Uses same weights as charge density
    for (int p = 0; p < np; ++p) {      
      int i = node_index[2*p];
      int j = node_index[2*p+1];

      get_Weights(&x_part[2*p], weights, &node_index[2*p], dh);

      E_part[0] = E_field[2*ij_to_index(i, j, nx2)]*weights[0] +
	          E_field[2*ij_to_index(i + 1, j, nx2)]*weights[1] +
		  E_field[2*ij_to_index(i, j + 1, nx2)]*weights[2] +
		  E_field[2*ij_to_index(i + 1, j + 1, nx2)]*weights[3];
      E_part[1] = E_field[2*ij_to_index(i, j, nx2) + 1]*weights[0] +
	          E_field[2*ij_to_index(i + 1, j, nx2) + 1]*weights[1] +
		  E_field[2*ij_to_index(i, j + 1, nx2) + 1]*weights[2] +
		  E_field[2*ij_to_index(i + 1, j + 1, nx2) + 1]*weights[3];

      v_part[2*p] = v_part[2*p] +(q/M)*E_part[0]*dt;
      v_part[2*p+1] = v_part[2*p+1] + (q/M)*E_part[1]*dt;

      x_part[2*p] = x_part[2*p] + v_part[2*p]*dt;
      x_part[2*p+1] = x_part[2*p+1] + v_part[2*p+1]*dt;

      // Check for boundaries
      // Reflective: elastic collision, bounces off the wall
      // 	     vx1 = vx1, vx2 = -vx2
      // 	     
      if (x_part[2*p+1] < 0) {
        x_part[2*p+1] = -1.0*x_part[2*p+1];
	v_part[2*p+1] = -1.0*v_part[2*p+1];
      }
      // Open: Remove particle
      // Move the last particle to replace the removed particle
      if (x_part[2*p] < 0 ||
	  x_part[2*p] > Lx1 ||
	  x_part[2*p+1] > Lx2 ||
	  (x_part[2*p] > box_range[0]*dh && x_part[2*p] < box_range[1]*dh &&
	  x_part[2*p+1] > box_range[2]*dh && x_part[2*p+1] < box_range[3]*dh)) {
      
	x_part[2*p] = x_part[2*(np-1)];
	x_part[2*p+1] = x_part[2*(np-1)+1];
	v_part[2*p] = v_part[2*(np-1)];
	v_part[2*p+1] = v_part[2*(np-1)+1];
;
	np -= 1;
      }
    }

    cout << "Adding new particles..." << endl;
    // Generate particles
    for (int new_p = 0; new_p < np_insert; ++new_p) {
      x_part[2*(np + new_p)] = dh*(double(rand())/RAND_MAX);
      x_part[2*(np + new_p) + 1] = Lx2*(double(rand())/RAND_MAX);

      v_part[2*(np + new_p)] = v_drift + 2.0*(double(rand())/RAND_MAX
		     		 + double(rand())/RAND_MAX  
			         + double(rand())/RAND_MAX 
			       - 1.5)*vth;
      v_part[2*(np + new_p) + 1] = 2.0*(double(rand())/RAND_MAX
		      + double(rand())/RAND_MAX + double(rand())/RAND_MAX
		      - 1.5)*vth;
    }
    np += np_insert;

    cout << "End of iteration, np = " << np << endl << endl;
    
    // Output info
    if ((iter+1)%25 == 0) {
      for (int i = 0; i < nn; ++i) {
        DataFile << iter << " " << rho[i] << " " << phi[i] << " " << E_field[2*i] << " ";
	DataFile << E_field[2*i + 1] << endl;
      }
      for (int i = 0; i < np; ++i) {
        ParticleFile << iter << " " << x_part[2*i] << " " << x_part[2*i+1] << " " ;
	ParticleFile << v_part[2*i] << " " << v_part[2*i+1] << endl;
      }

      NumFile << iter << " " << np << endl;
    }
  }

  delete(box_range);
  delete(x_part);
  delete(v_part);
  delete(rho);
  delete(weights);
  delete(node_index);
  delete(RHS);
  delete(E_field);
  delete(E_part);
  delete(phi);

  return 0;
}

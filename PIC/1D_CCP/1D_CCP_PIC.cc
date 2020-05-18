#include <iostream>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <cstdlib>


using namespace std;

void get_Weights(double part_x,
		 double *weights,
		 double *node_index,
		 double dx) {

*node_index = floor(part_x/dx);
double hx = (part_x -  (*node_index) * hx)/dx;
weights[0] = 1 - hx;
weights[1] = hx;
}

// solves laplace(phi) = -rho/epsilon0 = f
void jacobi_Update(double *phi,
		double *RHS,
		double dx,
		int nn,
		const int jacobi_max_iter,
		const double tol,
		const double n0,
		const double phi0,
		const double q,
		const double epsilon0,
		const double k,
		const double Te) {

  double *phi_new = new double[nn];
  double *RHS0 = new double[nn];
  copy_n(phi, nn, phi_new);
  copy_n(RHS, nn, RHS0);

  
  int jacobi_iter = 0;
  double residual;
  
  while (jacobi_iter <= jacobi_max_iter) {
    ++jacobi_iter;
    residual = 0.0;

    if (jacobi_iter%1000 == 0) {
      std::cout << "Jacobi iter = " << jacobi_iter << std::endl;
    }

    // Update RHS term to include electron terms
    for (int i = 0; i < nn; ++i) {
      RHS[i] = RHS0[i] - n0*exp((phi[i]-phi0)/(Te));
      RHS[i] = -q*RHS[i]/epsilon0;
    }

    // Run Jacobi update for interior nodes
    for (int i = 1; i < nn - 1; ++i) {
      phi_new[i] = 0.5*pow(dx,2)*(phi[i-1] + phi[i+1] - RHS[i]);
    }

    // Run Jacobi update for Neumann boundaries

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

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              Main Loop                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////

int main(void) {

  // Constantsa
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

  // DOUBLE CHECK THESE FORMULAS !!!!!!!///
  // Plasma parameters
  double lambdaD = sqrt(epsilon0*Te/(q*n0));
  double vth = sqrt(2.0 * q * Ti/M);

  // Problem discretization
  int nn = 191; // # of x1 nodes
  int ts = 4e6; // # of time steps
  double dx = 5.0e-4; // length of each cell
  int np_insert = (nx2-1)*15; // Put 15 macroparticles per time step

  // Other helpful values
  double dt = 1.0/(27.12e6*200.0); 
  double L = 0.095; // Domain length m

  // Electrode info
  // Bias electrode left, grounded electrode right
  int *elec_range = new int[4]; 
  elec_range[0] = 0; // Left bias electrode
  elec_range[1] = round(0.01/dx); // right side bias
  elec_range[2] = round(0.085/dx); // left side ground
  elec_range[3] = 190; // right side ground
  double f_hf = 27.12e6; // Hz
  double f_lf = 271.2e3; // Hz
  double V_hf = 110.0; // V
  double V_lf = 281.9; // V

  // Particle info
  int np_real = round(n0*v_drift*Lx2*dt); // Number of real particles entering
  double sp_wt = np_real/np_insert; // Real/Macroparticle
  double mp_q = 1; // Macroparticle charge
  int max_part = 4e6; // max_particles if none leave during time history

  // Particle data
  double *part_x = new double[max_part]; // Even indices = x1 values
  double *part_v = new double[max_part]; // Odd indices = x2 values
  double *rho = new double[nn];

  // Various helper variabes
  // Charge density
  double *weights = new double[2];
  double *node_index = new double[max_part];

  // Electric potential
  double *RHS = new double[nn];
  int jacobi_max_iter = 1e5;
  double tol = 1.0e-6;

  // Electric field
  double *E_field = new double[nn];

  // Generate particles
  int np = 0; 
  srand(time(NULL));


  // Move Particles
  double *part_E;
  
  // Set up phi boundaries and initial values
  double *phi = new double[nn];

  for (int i = 0; i < nn; ++i) {
    phi[i] = 0.0;
  }
 
  // Output files
  ofstream DataFile("Results/ESFieldData.txt");
  ofstream ParticleFile("Results/ParticleInfo.txt");
  ofstream NumFile("Results/NumParticles.txt");

  DataFile << "Iteration / Ion Density / Electric Potential /";
  DataFile << "Electric Field" << endl;

  ParticleFile << "Iteration / x1 / x2 / v1 / v2";
  ParticleFile << endl;

  NumFile << "Iteration / Number of Macroparticles" << endl;
  // Main Loop
  for (int iter = 0; iter < ts; ++iter) {     
    // Reset variables
    for (int i = 0; i< nn; ++i) {
      rho[i] = 0.0;
      E_field[i] = 0.0;
      RHS[i] = 0.0;
    }
    
    cout << "Iter: " << iter << endl;
    cout << "Computing ion charge density..." << endl;
    //Compute ion charge density 
    
    ///////////////////////////
    for (int p = 0; p < np; ++p) {
      get_Weights(x_part[p], weights, &node_index[p], dx);
      
      rho[node_index[p]] += weights[0];
      rho[node_index[p+1]] += weights[1];
    }

    for (int i = 0; i < nn; ++i) {
      rho[i] *= sp_wt*mp_q/(dx*dx);

      // If on the boundary, half the volume.
      if (i == 0 || i == nn-1) {
	rho[i] *= 2.0;
      }
      // Charge density for ions only, does not include electron term
      RHS[i] = rho[i];
    }
    
    cout << "Computing electric potential..." << endl;
    // Compute electric potential
    
    ////////////////////////////////
    jacobi_Update(phi, RHS, box_range, dh, dh, nx1, nx2, jacobi_max_iter, tol,
		    n0, phi0, q, epsilon0, k, Te);

    cout << "Computing electric field..." << endl;
    // Compute electric field
    
    ////////////////////////////////////////
    for (int i = 0; i < nn; ++i) {
      // Left
      if (i == 0) {
	E_field[index_ij] = -(phi[1] - phi[0])/dx;
      } // Right
      else if (i == nx1-1) {
	E_field[index_ij] = -(phi[i] -  phi[i-1])/dx;
      }
      else {
	E_field[index_ij] = -(phi[i+1] - phi[i-1])/(2.0*dh);
      }
    }

    cout << "Moving Particles..." << endl;
    // Move particles
    
    ///////////////////////////////////////////
    // Uses same weights as charge density
    for (int p = 0; p < np; ++p) {      

      part_E = E_field[node_index[p]]*weights[0] +
	          E_field[node_index[p] + 1]*weights[1];

      part_v[p] = part_v[p] +(q/M)*part_E*dt;

      part_x[p] = part_x[p] + part_v[p]*dt;
 
      /*
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
      } */
    }

    cout << "Adding new particles..." << endl;
    // Generate particles
    for (int new_p = 0; new_p < np_insert; ++new_p) {
      x_part[2*(np + new_p)] = dh*(double(rand())/RAND_MAX);
      x_part[2*(np + new_p) + 1] = Lx2*(double(rand())/RAND_MAX);

      v_part[2*(np + new_p)] = v_drift + 0.5*(double(rand())/RAND_MAX
		     		 + double(rand())/RAND_MAX  
			         + double(rand())/RAND_MAX 
			       - 1.5)*vth;
      v_part[2*(np + new_p) + 1] = 0.5*(double(rand())/RAND_MAX
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
  delete(part_x);
  delete(part_v);
  delete(rho);
  delete(weights);
  delete(node_index);
  delete(RHS);
  delete(E_field);
  delete(part_E);
  delete(phi);

  return 0;
}

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
		int *elec_range,
		double phi_left,
		double phi_right,
		double dx,
		int nn,
		const int jacobi_max_iter,
		const double tol,
		const double n0,
		const double phi0,
		const double e,
		const double epsilon0,
		const double k,
		const double Te) {

  double *phi_new = new double[nn];
  double *RHS0 = new double[nn];
  
  int jacobi_iter = 0;
  double residual;
  
  // Add dirichlet boundaries
  for (int i = elec_range[0]; i <= elec_range[1]; ++i) {
    phi[node] = phi_left;
  }

  for (int i = elec_range[2]; i <= elec_range[3]; ++i) {
    phi[node] = phi_right;
  }

  copy_n(phi, nn, phi_new);
  copy_n(RHS, nn, RHS0);

  while (jacobi_iter <= jacobi_max_iter) {
    ++jacobi_iter;
    residual = 0.0;

    if (jacobi_iter%1000 == 0) {
      std::cout << "Jacobi iter = " << jacobi_iter << std::endl;
    }

    // Update RHS term to include electron terms
    for (int i = 0; i < nn; ++i) {
      RHS[i] = RHS0[i] - n0*exp(q*(phi[i]-phi0)/(k*Te));
      RHS[i] = -q*RHS[i]/epsilon0;
    }

    // Run Jacobi update for interior nodes
    for (int i = elec_range[1] + 1; i < elec_range[2]; ++i) {
      phi_new[i] = 0.5*(phi[i-1] + phi[i+1] - pow(dx,2.0)*RHS[i]);
    }

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

  // Constants
  const double epsilon0 = 8.854e-12; // Permittivity of free space
  const double e = 1.602e-19; // Elementary charge [C]
  const double k = 1.381e-23; // Boltzman constant [J/K]
  const double AMU = 1.661e-27; // Atomic Mass Unit [kg]
  const double M = 32.0*AMU; // Ion mass of O2 [kg]

  // Input settings
  double n0 = 1.0e12; // Electron density [#/m^3]
  double phi0 = 0.0; // Reference potential
  double Te = 3000; // Electron temp [K], figure this out later
  double Ti = 300; // Ion temp [K]
  double v_drift = 7.0e3; // Ion injection velocity [m/s]
  double phi_p = -5.0; // Plate potential [V]

  // DOUBLE CHECK THESE FORMULAS !!!!!!!///
  // Plasma parameters
  //double lambdaD = sqrt(epsilon0*Te/(q*n0));
  double vth = sqrt(2.0 * k * Ti/M);

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
  double phi_left;
  double phi_right = 0.0;

  // Particle info
  int np_real = round(n0*v_drift*Lx2*dt); // Number of real particles entering
  int max_part = 4e6; // max_particles if none leave during time history
  
  // Particle data
  double *part_q = new double[max_part]; // Macroparticle charge
  double *part_spwt = np_real/np_insert; // Particle weight
  double *part_x = new double[max_part]; // Position
  double *part_v = new double[max_part]; // Velocity

  // Various helper variabes
  // Charge density
  double *weights = new double[2];
  double *node_index = new double[max_part];
  double *rho = new double[nn]; // Density at each node

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
  
  // Calculate collisions setup
  // Initialize for Ar^+ + Ar -> Ar^+ + Ar
  double *coll_back = new double[2*176]; // Energy for backscattering for first 176,
  					 // then corresponding cross section after
  double *coll_iso = new double[2*176]; // Energy and CS for isotropic
  // Electron-impact collisions
  double *coll_elec = new double[];


  // Set up phi boundaries and initial values
  double *phi = new double[nn];
  double t;

  for (int i = 0; i < nn; ++i) {
     phi[i] = 0.0;
  }
 
  // Output files
  ofstream FieldFile("Results/ESFieldData.txt");
  ofstream ParticleFile("Results/ParticleInfo.txt");
  ofstream NumFile("Results/NumParticles.txt");

  FieldFile << "Iteration / Ion Density / Electric Potential /";
  FieldFile << "Electric Field" << endl;

  ParticleFile << "Iteration / x / v / spwt / q";
  ParticleFile << endl;

  NumFile << "Iteration / Number of Macroparticles" << endl;
  // Main Loop
  for (int iter = 0; iter < ts; ++iter) {
    // Electrode Potential
    t = iter*dt;
    phi_left = V_hf*sin(2*M_PI*f_hf*t) +  
	    V_lf*sin(2*M_PI*f_lf*t);

    // Reset variables
    for (int i = 0; i< nn; ++i) {
      rho[i] = 0.0;
      E_field[i] = 0.0;
      RHS[i] = 0.0;
    }
    
    cout << "Iter: " << iter << endl;
    cout << "Computing ion charge density..." << endl;

    /////////////////////////////////////////////////

    //Compute ion charge density 
    
    /////////////////////////////////////////////////
    
    // Get charge from particles
    for (int p = 0; p < np; ++p) {
      get_Weights(x_part[p], weights, &node_index[p], dx);
      
      rho[node_index[p]] += part_spwt[p]*part_q[p]*weights[0];
      rho[node_index[p]+1] += part_spwt[p]*part_q[p]*weights[1];
    }

    // Accont for volume
    for (int i = 0; i < nn; ++i) {
      rho[i] /= dx*dx;

      // If on the boundary, half the volume.
      if (i == 0 || i == nn-1) {
	rho[i] *= 2.0;
      }
      // Charge density for ions only, does not include electron term
      RHS[i] = rho[i];
    }
    
    cout << "Computing electric potential..." << endl;
    //////////////////////////////////////////////////

    // Compute electric potential
    
    /////////////////////////////////////////////////
    jacobi_Update(phi, RHS, elec_range, phi_left, phi_right,
		    dx, nn, jacobi_max_iter, tol, 
		    n0, phi0, e, epsilon0, k, Te);

    cout << "Computing electric field..." << endl;
    /////////////////////////////////////////////////

    // Compute electric field
    
    //////////////////////////////////////////////////
    for (int i = 0; i < nn; ++i) {
      // Left
      if (i == 0) {
	E_field[index_ij] = -(phi[1] - phi[0])/dx;
      } // Right
      else if (i == nn-1) {
	E_field[index_ij] = -(phi[i] -  phi[i-1])/dx;
      }
      else {
	E_field[index_ij] = -(phi[i+1] - phi[i-1])/(2.0*dh);
      }
    }

    cout << "Moving Particles..." << endl;
    //////////////////////////////////////////////////

    // Move particles
    
    //////////////////////////////////////////////////
    // Uses same weights as charge density
    for (int p = 0; p < np; ++p) {      

      part_E = E_field[node_index[p]]*weights[0] +
	          E_field[node_index[p] + 1]*weights[1];

      part_v[p] = part_v[p] + (part_q[p]/M)*part_E*dt;
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
        FieldFile << iter << " " << rho[i] << " " << phi[i] << " " << E_field[i];
	FieldFile << endl;
      }
      for (int i = 0; i < np; ++i) {
        ParticleFile << iter << " " << part_x[i] << " " << part_v[i] << " ";
	ParticleFile << part_spwt[i] << " " << part_q[i] << endl;
      }

      NumFile << iter << " " << np << endl;
    }
  }

  delete(elec_range);
  delete(part_x);
  delete(part_v);
  delete(part_q);
  delete(part_spwt);
  delete(rho);
  delete(weights);
  delete(node_index);
  delete(RHS);
  delete(E_field);
  delete(part_E);
  delete(phi);

  return 0;
}

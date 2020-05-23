#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <cstdlib>

#include <1D_Poisson.h>
#include <collisionModules.h>

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


////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              Main Loop                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////

int main(void) {

  // Constants
  const double epsilon0 = 8.854e-12; // Permittivity of free space
  const double e = 1.602e-19; // Elementary charge [C]
  const double k_B = 1.381e-23; // Boltzman constant [J/K]
  const double AMU = 1.661e-27; // Atomic Mass Unit [kg]
  const double m_n = 39.948*AMU; // Ion mass of Ar [kg]
  const double m_e = 9.109e-31; // Electron mass [kg]

  // Input settings
  double n0 = 1.0e12; // Electron density [#/m^3]
  double phi0 = 0.0; // Reference potential
  double T_n = 300; // Neutral temp, [K]

  // Problem discretization
  int nn = 191; // # of x1 nodes
  int ts = 4e6; // # of time steps
  double dx = 5.0e-4; // length of each cell
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
  //int np_real = round(n0*v_drift*Lx2*dt); // Number of real particles entering
  int max_part = 4e6; // max_particles if none leave during time history
  
  // Electron particle data
  double *part_e_spwt = new double[max_part]; // Particle weight
  double *part_e_x = new double[max_part]; // Position
  double *part_e_vx = new double[max_part]; // Velocity
  double *part_e_vy = new double[max_part];
  double *part_e_vz = new double[max_part];
  // To convert J to eV divide by e (elementary charge)
  double *part_e_epsilon = new double[max_part]; // Energy [J]
  int np_e;

  // Ion particle data
  double *part_i_spwt = new double[max_part]; // Particle weight
  double *part_i_x = new double[max_part]; // Position
  double *part_i_vx = new double[max_part]; // Velocity
  double *part_i_vy = new double[max_part];
  double *part_i_vz = new double[max_part];
  // To convert J to eV divide by e (elementary charge)
  double *part_i_epsilon = new double[max_part]; // Energy [J]
  int np_i;

  double part_q;
  double max_epsilon = 0.0;

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
  srand(time(NULL));

  // Move Particles
  double *part_E;
  
  // Calculate collisions setup
  // See cstrs.dat.txt for what each column is
  int *N_coll = new int[4]; // Types of collisions for each reaction
  N_coll[0] = 4;
  N_coll[1] = 2;
  N_coll[2] = 2;
  N_coll[3] = 1;

  double *CS_energy = new double[5996]; // energy values
  double *el_n_CS = new double[N_coll[0]*5996]; // electron-neutral
  double *el_ex_CS = new double[N_coll[1]*5996]; // electron-excited
  double *i_n_CS = new double [N_coll[2]*5996]; // ion-neutral

  ifstream coll_data("crtrs.dat.txt");
  // Ignore the first 2 lines
  coll_data.ignore(1e5, '\n');
  coll_data.ignore(1e5, '\n');
  
  for (int i = 0; i < 5996; ++i) {
    coll_data >> CS_energy[i];
    for (int j = 0; j < N_coll[0]; ++j) {
      coll_data >> el_n_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[1]; ++j) {
      coll_data >> el_ex_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[2]; ++j) {
      coll_data >> i_n_CS[j*5996 + i];
    }
  }

  coll_data.close();

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
    // Electrons
    for (int p = 0; p < np_e; ++p) {
      get_Weights(part_e_x[p], weights, &node_index[p], dx);

      rho[node_index[p]] += part_spwt[p]*(-1.0*e)*weights[0];
      rho[node_index[p]+1] += part_spwt[p]*(-1.0*e)*weights[1];
    }
    // Ions
    for (int p = 0; p < np_i; ++p) {
      get_Weights(part_i_x[p], weights, &node_index[p], dx);

      rho[node_index[p]] += part_spwt[p]*(e)*weights[0];
      rho[node_index[p]+1] += part_spwt[p]*(e)*weights[1];
    }

    // Account for volume
    for (int i = 0; i < nn; ++i) {
      rho[i] /= dx*dx;

      // If on the boundary, half the volume.
      if (i == 0 || i == nn-1) {
	rho[i] *= 2.0;
      }
      // Right hand side of Poisson equation
      RHS[i] = -rho[i]/epsilon0;
    }
    
    cout << "Computing electric potential..." << endl;

    //////////////////////////////////////////////////

    // Compute electric potential
    
    /////////////////////////////////////////////////
    
    triDiagSolver(phi, RHS, elec_range, phi_left, phi_right, dx, nn);

    cout << "Computing electric field..." << endl;

    /////////////////////////////////////////////////

    // Compute electric field
    
    //////////////////////////////////////////////////
    
    for (int i = 0; i < nn; ++i) {
      // Left
      if (i == 0) {
	E_field[i] = -(phi[1] - phi[0])/dx;
      } // Right
      else if (i == nn-1) {
	E_field[i] = -(phi[i] -  phi[i-1])/dx;
      }
      else {
	E_field[i] = -(phi[i+1] - phi[i-1])/(2.0*dh);
      }
    }

    cout << "Moving Particles..." << endl;

    //////////////////////////////////////////////////

    // Move particles
    
    //////////////////////////////////////////////////

    // Uses same weights as charge density
    for (int p = 0; p < np_e; ++p) {      
      get_Weights(part_e_x[p], weights, &node_index[p], dx);

      part_E = E_field[node_index[p]]*weights[0] +
	          E_field[node_index[p] + 1]*weights[1];

      part_e_vx[p] = part_e_vx[p] + (-1.0*e)/(m_e)*part_E*dt;
      part_e_x[p] = part_e_x[p] + part_e_vx[p]*dt;

      part_e_epsilon[p] = 0.5*m_e*(part_e_vx[p]*part_e_vx[p]+
		      part_e_vy[p]*part_e_vy[p] +
		      part_e_vz[p]*part_e_vz[p]);
      if (part_e_epsilon[p] > max_e_epsilon) {
        max_e_epsilon = part_e_epsilon[p];
      }
    }
    for (int p = 0; p < np_i; ++p) {      
      get_Weights(part_i_x[p], weights, &node_index[p], dx);

      part_E = E_field[node_index[p]]*weights[0] +
	          E_field[node_index[p] + 1]*weights[1];

      part_i_vx[p] = part_i_vx[p] + (e)/(m_n)*part_E*dt;
      part_i_x[p] = part_i_x[p] + part_i_vx[p]*dt;

      part_i_epsilon[p] = 0.5*m_e*(part_i_vx[p]*part_i_vx[p]+
		      part_i_vy[p]*part_i_vy[p] +
		      part_i_vz[p]*part_i_vz[p]);

     if (part_i_epsilon[p] > max_i_epsilon) {
        max_i_epsilon = part_i_epsilon[p];
      }
    }

 
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

    cout << "Calculating collisions..." << endl;

    //////////////////////////////////////////////////

    // Collision modules
    
    //////////////////////////////////////////////////
    cout << "End of iteration, np = " << np << endl << endl;

    ////////////////////////////////////////////////////
    
    // Output info
    
    ////////////////////////////////////////////////////

    if ((iter+1)%25 == 0) {
      for (int i = 0; i < nn; ++i) {
        FieldFile << iter << " " << rho[i] << " " << phi[i] << " " << E_field[i];
	FieldFile << endl;
      }
      for (int i = 0; i < np; ++i) {
        ParticleFile << iter << " " << part_x[i] << " " << part_vx[i] << " ";
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
  delete(coll_CS);

  return 0;
}

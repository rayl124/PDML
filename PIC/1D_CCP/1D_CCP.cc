#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <cstdlib>

#include "species.h"
#include "1D_Poisson.h"
#include "collisionModules.h"

using namespace std;

void get_Weights(double part_x,
		 double *weights,
		 int *node_index,
		 double dx) {

*node_index = floor(part_x/dx);
double hx = (part_x -  (*node_index) * hx)/dx;
weights[0] = 1 - hx;
weights[1] = hx;
}

void legendrePoly(double x_target, double *x, 
		double *l_coeff, int order) {

  for (int j = 0; j < order; j++) {
    l_coeff[j] = 1.0;
    for (int k = 0; k < order; ++k) {
      if (j!=k) {
	l_coeff[j] *= (x_target-x[k])/(x[j]-x[k]);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
//                                                                    //
//                              Main Loop                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////

int main(void) {
  // Random Seed
  srand(time(NULL));

  // Constants
  const double epsilon0 = 8.854e-12; // Permittivity of free space
  const double e = 1.602e-19; // Elementary charge [C]
  const double k_B = 1.381e-23; // Boltzman constant [J/K]
  const double AMU = 1.661e-27; // Atomic Mass Unit [kg]

  // Input settings
  double phi0 = 0.0; // Reference potential
  double T_n = 300; // Neutral temp, [K]
  double n_n = 0;

  // Problem discretization
  int nn = 191; // # of x1 nodes
  int ts = 200; // # of time steps
  double dx= 5.0e-4; // length of each cell
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
  double L_inner = 0.085-0.01;

  //////////////////////////////////////////////////////////////
  //
  // Load data for collisions
  //
  ///////////////////////////////////////////////////////////
  
  // See cstrs.dat.txt for what each column is
  int *N_coll = new int[4]; // Types of collisions for each reaction
  N_coll[0] = 4;
  N_coll[1] = 2;
  N_coll[2] = 2;
  N_coll[3] = 1;
  
  int data_set_length = 5996;
  double *CS_energy = new double[data_set_length]; // energy values
  double *e_n_CS = new double[N_coll[0]*data_set_length]; // electron-neutral
  double *e_ex_CS = new double[N_coll[1]*data_set_length]; // electron-excited
  double *i_n_CS = new double [N_coll[2]*data_set_length]; // ion-neutral

  ifstream coll_data("crtrs.dat.txt");
  // Ignore the first 2 lines
  coll_data.ignore(1e5, '\n');
  coll_data.ignore(1e5, '\n');
  
  for (int i = 0; i < data_set_length; ++i) {
    coll_data >> CS_energy[i];
    for (int j = 0; j < N_coll[0]; ++j) {
      coll_data >> e_n_CS[j*data_set_length + i];
    }
    for (int j = 0; j < N_coll[1]; ++j) {
      coll_data >> e_ex_CS[j*data_set_length + i];
    }
    for (int j = 0; j < N_coll[2]; ++j) {
      coll_data >> i_n_CS[j*data_set_length + i];
    }
  }

  coll_data.close();


  ///////////////////////////////////////////////////////////////
  //
  // Particle variables
  //
  // //////////////////////////////////////////////////////////
  int max_part = 1e6; // max_particles if none leave during time history

  // Electron particle data
  species electron;
  electron.initialize(max_part);
  electron.m = 9.109e-31; //[kg]
  electron.q = -1.0*e; //[C]
  electron.np = 0;

  // Ion particle data
  species ion;
  ion.initialize(max_part);
  ion.m = 39.948*AMU; //[kg]
  ion.q = e; //[c]
  ion.np = 0;

  double m_n = ion.m;

  // Calculate collisions setup
  double max_index, nu_max, P_max;
  int N_c = 0;
  int rand_index, type;
  double epsilon_exc, epsilon_ion;


  /////////////////////////////////////////////////////////
  // 
  // Field variables
  //
  ////////////////////////////////////////////////////////////
  
  // Charge density
  double *weights = new double[2];
  double *rho = new double[nn]; // Density at each node

  // Electric potential
  int nn_inner = elec_range[2] - elec_range[1] - 1; // Interior nodes
  double *RHS = new double[nn];
  double *a = new double[nn_inner];
  double *b = new double[nn_inner];
  double *c = new double[nn_inner];

  // Electric field
  double *E_field = new double[nn];

  // Move Particles
  double part_E;
  
  // Set up phi boundaries and initial values
  double *phi = new double[nn];
  //double *phi_cc = new double[nn-1];
  //double *l_coeff = new double[3];
  double t;
  /* Cell-centered discretization
  double *x_cc_left = new double[3];
  double *x_cc_right = new double[3];
  x_cc_left[0] = -0.5*dx;
  x_cc_left[1] = 0.5*dx;
  x_cc_left[2] = 1.5*dx;
  x_cc_right[0] = (0.075)-1.5*dx;
  x_cc_right[1] = 0.075-0.5*dx;
  x_cc_right[2] = 0.075+0.5*dx;
  */ 

  // Set up exact solution;
  double *phi_exact = new double[nn];
  //double *phi_cc_exact = new double[nn-1];
  for (int i = 0; i < nn; ++i) {
     phi[i] = 0.0;
     phi_exact[i] = 0.0;
     /*
     if (i < nn-1) {
       phi_cc[i] = 0.0;
       phi_cc_exact[i] = 0.0;
     } */
  }
  //
 
  //////////////////////////////////////////////////////////
  //
  // Output files
  //
  //////////////////////////////////////////////////////////
  
  ofstream FieldFile("Results/ESFieldData.txt");
  //ofstream ElectronFile("Results/ElectronInfo.txt");
  //ofstream IonFile("Results/NumParticles.txt");

  FieldFile << "Iteration / Node x / Electric Potential /";
  FieldFile << "Electric Field" << endl;

  //ParticleFile << "Iteration / x / v / spwt / q";
  //ParticleFile << endl;

  //NumFile << "Iteration / Number of Macroparticles" << endl;
  
  //////////////////////////////////////////////////////////
  //
  // Main Loop
  //
  /////////////////////////////////////////////////////////
  
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
      
      if (i >=elec_range[0] && i <= elec_range[1]) {
        phi_exact[i] = phi_left;
      } else if (i >= elec_range[2] && i <= elec_range[3]) {
	phi_exact[i] = phi_right;
      } else {
	phi_exact[i] = (phi_right-phi_left)/(L_inner)*(dx*i-(elec_range[1]*dx)) 
		+ phi_left;
      }
      /*

      if (i >=elec_range[0] && i <= elec_range[1]-1) {
        phi_cc_exact[i] = phi_left;
      } else if (i >= elec_range[2] && i <= elec_range[3]) {
	phi_cc_exact[i] = phi_right;
      } else {
	phi_cc_exact[i] = (phi_right-phi_left)/(L_inner)*(dx*(i+0.5)-(elec_range[1]*dx)) 
		+ phi_left;
      }*/
    }
    
    cout << "Iter: " << iter << endl;
    
    /////////////////////////////////////////////////////////
    //
    // Compute ion charge density 
    //
    ////////////////////////////////////////////////////////
    
    cout << "Computing ion charge density..." << endl;

    // Get charge from particles
    // Electrons
    for (int p = 0; p < electron.np; ++p) {
      get_Weights(electron.x[p], weights, &electron.node_index[p], dx);

      rho[electron.node_index[p]] += electron.spwt[p]*
	      			     (electron.q)*weights[0];
      rho[electron.node_index[p]+1] += electron.spwt[p]*
	      			       (electron.q)*weights[1];
    }
    // Ions
    for (int p = 0; p < ion.np; ++p) {
      get_Weights(ion.x[p], weights, &ion.node_index[p], dx);

      rho[ion.node_index[p]] += ion.spwt[p]*ion.q*weights[0];
      rho[ion.node_index[p]+1] += ion.spwt[p]*ion.q*weights[1];
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
    
    //////////////////////////////////////////////////
    //
    // Compute electric potential
    //
    /////////////////////////////////////////////////

    cout << "Computing electric potential..." << endl;
    for (int i = elec_range[0]; i <= elec_range[1]; ++i) {
      phi[i] = phi_left;
      /*if (i < elec_range[1]) {
        /phi_cc[i] = phi_left;
      }*/
    }
    for (int i = elec_range[2]; i <= elec_range[3]; ++i) {
      phi[i] = phi_right;
      /*if (i < elec_range[3]) {
        /phi_cc[i] = phi_right;
      }*/
    }
    
    RHS[elec_range[1] + 1] -= phi_left/(dx*dx);
    RHS[elec_range[2] - 1] -= phi_right/(dx*dx);

    for (int i = 0; i < nn_inner; ++i) {
      a[i] = 1.0;
      b[i] = -2.0;
      c[i] = 1.0;
      RHS[i] *= dx*dx;
    }
    a[0] = 0.0;
    c[nn_inner-1] = 0.0;
    
    triDiagSolver(&phi[elec_range[1] + 1], a, b, c, 
		    &RHS[elec_range[1] + 1], nn_inner);

    // Cell Centered
    /*
    c[nn_inner-1] = 1.0;
    c[nn_inner] = 0.0;

    legendrePoly(0.0, x_cc_left, l_coeff, 3);
    b[0] -= l_coeff[1]/l_coeff[0];
    c[0] -= l_coeff[2]/l_coeff[0];

    for (int i = 0; i < nn; ++i) {
      RHS[i] = 0.0;
    }

    RHS[elec_range[1]] -= phi_left/(dx*dx*l_coeff[0]);

    legendrePoly(0.075, x_cc_right, l_coeff, 3);
    b[nn_inner] -= l_coeff[1]/l_coeff[2];
    a[nn_inner] -= l_coeff[0]/l_coeff[2];
    RHS[elec_range[2] - 1] -= phi_right/(dx*dx*l_coeff[2]);

    for (int i = 0; i < nn; ++i) {
      RHS[i] *= dx*dx;
    }

    triDiagSolver(&phi_cc[elec_range[1]], a, b, c, 
		    &RHS[elec_range[1]], nn_inner + 1);
    */


    
     /////////////////////////////////////////////////
    //
    // Compute electric field
    //
    //////////////////////////////////////////////////
    
    cout << "Computing electric field..." << endl;

    // Finite difference accuracy 2 order
    for (int i = 0; i < nn; ++i) {
      // Left
      if (i == 0) {
	E_field[i] = -(-0.5*phi[2] + 2.0*phi[1] - 1.5*phi[0])/dx;
      } // Right
      else if (i == nn-1) {
	E_field[i] = -(0.5*phi[i-2] -2.0*phi[i-1] + 1.5*phi[i])/dx;
      }
      else {
	E_field[i] = -(phi[i+1] - phi[i-1])/(2.0*dx);
      }
    }

    //////////////////////////////////////////////////
    //
    // Move particles
    //
    //////////////////////////////////////////////////

    cout << "Moving Particles..." << endl;

    // Uses same weights as charge density
    for (int p = 0; p < electron.np; ++p) {      
      get_Weights(electron.x[p], weights, &electron.node_index[p], dx);

      part_E = E_field[electron.node_index[p]]*weights[0] +
	          E_field[electron.node_index[p] + 1]*weights[1];

      electron.vx[p] = electron.vx[p] + electron.q/(electron.m) *
	      	       part_E*dt;
      electron.x[p] = electron.x[p] + electron.vx[p]*dt;

      electron.epsilon[p] = 0.5*electron.m*pow(getv(electron.vx[p], 
			      electron.vy[p], electron.vz[p]),2.0)/e;

      if (electron.epsilon[p] > electron.max_epsilon) {
        electron.max_epsilon = electron.epsilon[p];
      }
    }
    for (int p = 0; p < ion.np; ++p) {      
      get_Weights(ion.x[p], weights, &ion.node_index[p], dx);

      part_E = E_field[ion.node_index[p]]*weights[0] +
	          E_field[ion.node_index[p] + 1]*weights[1];

      ion.vx[p] = ion.vx[p] + ion.q/(ion.m) *
	      	       part_E*dt;
      ion.x[p] = ion.x[p] + ion.vx[p]*dt;

      ion.epsilon[p] = 0.5*ion.m*pow(getv(ion.vx[p], 
			      ion.vy[p], ion.vz[p]),2.0)/e;

      if (ion.epsilon[p] > ion.max_epsilon) {
        ion.max_epsilon = ion.epsilon[p];
      }
    }
    // Boundaries    

    //////////////////////////////////////////////////
    //
    // Collision modules
    //
    //////////////////////////////////////////////////
    
    cout << "Calculating collisions..." << endl;

    // Get number of particles for electron - neutral collisions
    if (electron.np > 0) {
    getNullCollPart(CS_energy, e_n_CS, electron.max_epsilon, 
		  &nu_max, &P_max, &N_c,
		  electron.m, n_n, dt, electron.np, 
		  N_coll[0], data_set_length);
    }

    for (int i = 0; i < N_c; ++i) {
      rand_index = round((double(rand())/RAND_MAX)*(electron.np-1));
      type = getCollType(CS_energy, e_n_CS, electron.epsilon[rand_index],
		    nu_max, electron.m, n_n, N_coll[0], data_set_length);
    
      // switch-case for electron-neutral collisions
      // 0 - elastic
      // 1 - excitation 1
      // 2 - excitation 2
      // 3 - ionization
      // 4 - null
      switch(type) {
        case 0:
	  e_elastic(&electron.vx[rand_index], &electron.vy[rand_index],
		  &electron.vz[rand_index], electron.epsilon[rand_index],
		  electron.m, m_n);
          electron.epsilon[rand_index] = 0.5*electron.m*
		  pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
		  electron.vz[rand_index]),2.0)/e;
	  continue;
        case 1:
	  epsilon_exc = 1.160330e1; // From crtrs.dat.txt
	  e_excitation(&electron.vx[rand_index], &electron.vy[rand_index],
		  &electron.vz[rand_index], electron.epsilon[rand_index],
		  epsilon_exc);
          electron.epsilon[rand_index] = 0.5*electron.m*
		  pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
		  electron.vz[rand_index]),2.0)/e;

	  continue;
        case 2:
	  epsilon_exc = 1.31041e1; // From crtrs.dat.txt
	  e_excitation(&electron.vx[rand_index], &electron.vy[rand_index],
		  &electron.vz[rand_index], electron.epsilon[rand_index],
		  epsilon_exc);
          electron.epsilon[rand_index] = 0.5*electron.m*
		  pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
		  electron.vz[rand_index]),2.0)/e;
	  continue;
        case 3:
	  electron.np += 1;
	  electron.x[electron.np-1] = 0.0;
	  electron.x[electron.np-1] = 0.0;
	  electron.x[electron.np-1] = 0.0;
	  epsilon_ion = 1.60055e1; // From crtrs.dat.txt
	  e_ionization(&electron.vx[rand_index], &electron.vy[rand_index],
		  &electron.vz[rand_index], &electron.vx[electron.np-1], 
		  &electron.vy[electron.np-1], &electron.vz[electron.np-1], 
		  electron.epsilon[rand_index], epsilon_ion);
	  electron.epsilon[rand_index] = 0.5*electron.m*pow(getv(electron.vx[rand_index], 
			electron.vy[rand_index], electron.vz[rand_index]),2.0)/e; //[eV]
	  electron.epsilon[electron.np-1] = 0.5*electron.m*
		  pow(getv(electron.vx[electron.np-1],
		  electron.vy[electron.np-1], 
	 	  electron.vz[electron.np-1]),2.0)/e; //[eV]
	  continue;	
      case 4:
	  continue;
      }
    }
    
    // Ion-neutral collisions
    if (ion.np > 0) {
    getNullCollPart(CS_energy, i_n_CS, ion.max_epsilon, &nu_max, &P_max, &N_c,
		  ion.m, n_n, dt, ion.np, N_coll[1], data_set_length);
    }

    for (int i = 0; i < N_c; ++i) {
      rand_index = round((double(rand())/RAND_MAX)*(ion.np-1));
      type = getCollType(CS_energy, i_n_CS, ion.epsilon[rand_index],
		    nu_max, ion.m, n_n, N_coll[1], data_set_length);
    
      // switch-case for electron-neutral collisions
      // 0 - charge exchange
      // 1 - back scattering
      // 2 - null
      switch(type) {
        case 0:
	  thermalVelSample(&ion.vx[rand_index], &ion.vy[rand_index],
			  &ion.vz[rand_index], T_n, m_n);
	  ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
			ion.vy[rand_index], ion.vz[rand_index]),2.0)/e; //[eV]
	  continue;
        case 1:
	  i_scattering(&ion.vx[rand_index], &ion.vy[rand_index],
		       &ion.vz[rand_index], ion.epsilon[rand_index],
		       ion.m, ion.m, T_n);
	  ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
			ion.vy[rand_index], ion.vz[rand_index]),2.0)/e; //[eV]

	  continue;
        case 2:
	  continue;
      }
    }


    //cout << "End of iteration, np = " << np << endl << endl;

    ////////////////////////////////////////////////////
    //
    // Output info
    //
    ////////////////////////////////////////////////////

    if ((iter+1)%5 == 0) {
      for (int i = 0; i < nn; ++i) {
        FieldFile << iter << " " << dx*i << " " << phi[i] << " ";
	FieldFile << E_field[i] << endl;
      }
      /*
      for (int i = 0; i < electron.np; ++i) {
        ParticleFile << iter << " " << part_x[i] << " " << part_vx[i] << " ";
	ParticleFile << part_spwt[i] << " " << part_q[i] << endl;
      }
      for (int i = 0; i < ion.np; ++i) {

      }
      NumFile << iter << " " << np << endl;
      */
    }
  }
  delete(elec_range);
  delete(CS_energy);
  delete(e_n_CS);
  delete(e_ex_CS);
  delete(i_n_CS);
  delete(rho);
  delete(weights);
  delete(RHS);
  delete(E_field);
  delete(phi);
  delete(phi_exact);
  
  electron.clean();
  ion.clean();

  return 0;
}

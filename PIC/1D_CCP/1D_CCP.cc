#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <cstdlib>

#include "species.h"
#include "solverModules.h"
#include "collisionModules.h"

using namespace std;

// Node centered weights
void get_WeightsNC(double part_x,
		 double *weights,
		 int *node_index,
		 double dx) {

*node_index = floor(part_x/dx);
double hx = (part_x -  (*node_index) * dx)/dx;
weights[0] = 1 - hx;
weights[1] = hx;
}

// Cell centered weights
void get_WeightsCC(double part_x,
		 double *weights,
		 int *cell_index,
		 int *elec_range,
		 double dx) {

  *cell_index = floor((part_x-0.5*dx)/dx);
  if (*cell_index < elec_range[1]) {
    *cell_index = elec_range[1];
    weights[0] = 1.0;
    weights[1] = 0.0;
  } else if (part_x > elec_range[2]*dx - 0.5*dx) {
    weights[0] = 1.0;
    weights[1] = 0.0;
  } else {
  double hx = (part_x -  (*cell_index + 0.5) * dx)/dx;
  weights[0] = 1 - hx;
  weights[1] = hx;
  }
}


void lagrangePoly(double x_target, double *x, 
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

  // Input Parameters
  double P = 0.3;

  // Problem discretization
  int nn = 191; // # of x1 nodes
  int n_cell = nn-1; // # of cells
  double x_start = 0.0;
  double x_end = 0.095;

  double L = x_end-x_start; // Domain length m, including electrodes
  int ts = 200*100; // # of time steps
  double dx= L/((double) nn-1); // length of each cell
  double dt = 1.0/(27.12e6*200.0); 

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
  int max_part = 4e6; // max_particles if none leave during time history
  double real_part = 7.5e-6*0.3*L_inner/(k_B*300.0);

  // Electron particle data
  particles electron;
  electron.initialize(max_part, n_cell);
  electron.m = 9.109e-31; //[kg]
  electron.q = -1.0*e; //[C]
  electron.np = 1e4;
  electron.T = 300.0; //[K]
  electron.spwt = real_part/electron.np;
  electron.gamma[0] = 0.2;

  // Ion particle data
  particles ion;
  ion.initialize(max_part, n_cell);
  ion.m = 39.948*AMU; //[kg]
  ion.q = e; //[c]
  ion.np = 1e4;
  ion.T = 300.0; //[K]
  ion.spwt = real_part/electron.np;
  ion.gamma[1] = 0.15;

  // Uniformly distribute initial positions, assign initial velocity from
  // thermal distribution
  for (int i = 0; i < electron.np; ++i) {
    electron.x[i] = double(rand())/RAND_MAX*(elec_range[2]-elec_range[1])*dx + 
	    elec_range[1]*dx;
    electron.thermalVelocity(i);
  }

  for (int i = 0; i < ion.np; ++i) {
    ion.x[i] = double(rand())/RAND_MAX*(elec_range[2]-elec_range[1])*dx + 
	    elec_range[1]*dx;
    ion.thermalVelocity(i);
  }

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
  
  fluid neutral;
  neutral.m = 39.948*AMU;
  neutral.initialize(n_cell);

  fluid excited;
  excited.m = 39.948*AMU;
  excited.initialize(n_cell);

  // Fluid density
  double epsilon_LJ = 93.30; // [K]
  double d_LJ = 3.542; // [Angstrom]

  for (int i = 0; i < n_cell; ++i) {
    neutral.n_dot[i] = 0.0;
    neutral.T[i] = 300.0;


    excited.n_dot[i] = 0.0;
    excited.T[i] = 300.0;

    if (i >= elec_range[1] && i < elec_range[2]) {
      neutral.n[i] = P/(k_B*neutral.T[i]);
      excited.n[i] = P*5.0e-5/(k_B*excited.T[i]);
    } else {
      neutral.n[i] = 0.0;
      excited.n[i] = 0.0;
    }
  }

  // Charge density
  double *weights = new double[2];
  // double *n_n = new double[nn];
  double *rho = new double[n_cell]; // Charge density at each cell center

  // Electric potential
  double *phi = new double[n_cell];
  double t;  // time
  double *RHS = new double[n_cell]; // Right hand side
  double *a = new double[n_cell];  // Tridiagonal coeffs
  double *b = new double[n_cell];
  double *c = new double[n_cell];

  double *lagrangeLeft = new double[2];
  double *lagrangeRight = new double[2];
  double *l_coeffL = new double[2];
  double *l_coeffR = new double[2];
  double ghostL, ghostR;

  for (int i = 0; i < 2; ++i) {
    lagrangeLeft[i] = elec_range[1]*dx + (i - 0.5)*dx;
    lagrangeRight[i] = elec_range[2]*dx + (i - 1.5)*dx;
  }

  // Electric field
  double *E_field = new double[nn]; // electric field at the node

  // Move Particles
  double part_E;
  double T_n;

  //////////////////////////////////////////////////////////
  //
  // Output files
  //
  //////////////////////////////////////////////////////////
  
  int write_iter = 25; // Write every x number of iterations
  string simNum ("004");
  
  ofstream InputFile("Results/InputData/Input"+simNum+".txt");
  ofstream FieldCCFile("Results/FieldData/FieldCCData"+simNum+".txt");
  ofstream FieldNCFile("Results/FieldData/FieldNCData"+simNum+".txt");

  ofstream NumFile("Results/NumberPartData/NumberPart"+simNum+".txt");

  InputFile << "Misc comments: Fixed E_field bug" << endl;
  InputFile << "Pressure [Pa] / electron.np / ion.np / electron.spwt / ion.spwt / ";
  InputFile << "V_hf / V_lf / f_hf / f_lf / Total steps / dt / NumNodes" << endl;
  InputFile << P << " " << electron.np << " " << ion.np << " " << electron.spwt;
  InputFile << " " << ion.spwt << " " << V_hf << " " << V_lf << " " << f_hf;
  InputFile << " " << f_lf << " " << ts << " " << dt << " " << nn << endl;


  FieldCCFile << "Iteration / Time / Cell x / Charge Density / ";
  FieldCCFile << "Electric Potential / "<< endl;

  FieldNCFile << "Iteration / Time / Node x / Electric Field" << endl;

  //ParticleFile << "Iteration / x / v / spwt / q";
  //ParticleFile << endl;

  NumFile << "Iteration / Electron np / Ion np" << endl;
  
  //////////////////////////////////////////////////////////
  //
  // Main Loop
  //
  /////////////////////////////////////////////////////////
  
  for (int iter = 0; iter < ts; ++iter) {
    // Electrode Potential boundaries
    t = iter*dt;
    phi_left = V_hf*sin(2*M_PI*f_hf*t) +  
	    V_lf*sin(2*M_PI*f_lf*t);
    // Reset variables
    electron.max_epsilon = 0.0;
    ion.max_epsilon = 0.0;
    
    for (int i = 0; i< n_cell; ++i) {
      phi[i] = 0.0;
      rho[i] = 0.0;
      electron.n[i] = 0.0;
      ion.n[i] = 0.0;
      E_field[i] = 0.0;
      RHS[i] = 0.0;
    }
    E_field[nn-1] = 0.0;
    
    cout << "Iter: " << iter << endl;

    /////////////////////////////////////////////////////////
    //
    // Compute charge density 
    //
    ////////////////////////////////////////////////////////
    
    cout << "Computing charge density..." << endl;

    // Get number densities
    // Electrons
    for (int p = 0; p < electron.np; ++p) {
      get_WeightsCC(electron.x[p], weights, &electron.cell_index[p],
		     elec_range, dx);

      electron.n[electron.cell_index[p]] += 
	      				electron.spwt*weights[0];
      electron.n[electron.cell_index[p]+1] += 
	      				electron.spwt*weights[1];
     
    }
    // Ions
    for (int p = 0; p < ion.np; ++p) {
      get_WeightsCC(ion.x[p], weights, &ion.cell_index[p], 
		      elec_range, dx);

      ion.n[ion.cell_index[p]] += ion.spwt*weights[0];
      ion.n[ion.cell_index[p]+1] += ion.spwt*weights[1];
    }

    // Account for volume
    for (int i = 0; i < n_cell; ++i) {
      electron.n[i] /= dx;
      ion.n[i] /= dx;

      // Charge density
      rho[i] = e*(ion.n[i]-electron.n[i]);

      // Right hand side of Poisson equation
      RHS[i] = -rho[i]/epsilon0;
    }

    /////////////////////////////////////////////////////////
    //
    // Compute fluid density
    //
    ////////////////////////////////////////////////////////
    
    cout << "Computing fluid density" << endl;
  
    if (iter > 0) {

      for (int i = elec_range[1]; i < elec_range[2]; ++i) {
        neutral.D[i] = 3.0*sqrt(k_B*neutral.T[i]) /
	      (8.0*sqrt(M_PI*neutral.m)*d_LJ*d_LJ*1.0e-20);
        neutral.D[i] /= neutral.n[i] + excited.n[i] + ion.n[i];

        excited.D[i] = 3.0*sqrt(k_B*excited.T[i]) /
	      (8.0*sqrt(M_PI*excited.m)*d_LJ*d_LJ*1.0e-20);
        excited.D[i] /= neutral.n[i] + excited.n[i] + ion.n[i];
      }

      double n_ghostL = neutral.flux_L*dx +
	                neutral.D[elec_range[1]]*neutral.n[elec_range[1]];
      n_ghostL /= neutral.D[elec_range[1]];
      double n_ghostR = neutral.flux_R*dx + 
	                neutral.D[elec_range[2]-1]*neutral.n[elec_range[2]-1];
      n_ghostL /= neutral.D[elec_range[2]-1];

      driftDiffusionFVExplicit2(&neutral.n[elec_range[1]],
		      &neutral.n_dot[elec_range[1]], 
		      &neutral.D[elec_range[1]],
		      n_ghostL, 
		      n_ghostR, 
		      dx, dt, elec_range[2]-elec_range[1]);
      driftDiffusionFVExplicit2(&excited.n[elec_range[1]], 
		      &excited.n_dot[elec_range[1]],
		      &excited.D[elec_range[1]],
		      0.0, 0.0, dx, dt, elec_range[2]-elec_range[1]);
    }


    //////////////////////////////////////////////////
    //
    // Compute electric potential
    //
    /////////////////////////////////////////////////

    cout << "Computing electric potential..." << endl;

    // Interpolation coefficients for the boundaries
    lagrangePoly(elec_range[1]*dx, lagrangeLeft, l_coeffL, 2);
    lagrangePoly(elec_range[2]*dx, lagrangeRight, l_coeffR, 2);
    RHS[elec_range[1]] -= phi_left/(dx*dx*l_coeffL[0]);
    RHS[elec_range[2] - 1] -= phi_right/(dx*dx*l_coeffR[1]);

    // Tridiagonal coefficients
    for (int i = 0; i < n_cell; ++i) {
      if (i >= elec_range[1] && i < elec_range[2]) {	    
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
        RHS[i] *= dx*dx;
      }
      else {
        a[i] = 0.0;
	b[i] = 1.0;
	c[i] = 0.0;
	if (i < elec_range[1]) {
	  RHS[i] = phi_left;
	} else {
	  RHS[i] = phi_right;
	}
      }
    }

    // Adjust RHS to account for BC
    a[elec_range[1]] = 0.0;
    b[elec_range[1]] -= l_coeffL[1]/l_coeffL[0];
    //c[elec_range[1]] -= l_coeffL[2]/l_coeffL[1];

    //a[elec_range[2] - 1] -= l_coeffR[0]/l_coeffR[2];
    b[elec_range[2] - 1] -= l_coeffR[0]/l_coeffR[1];
    c[elec_range[2] - 1] = 0.0;
     
    // Solve
    triDiagSolver(phi, a, b, c, RHS, n_cell);
    
    // Get ghost node values
    ghostL = phi_left-l_coeffL[1]*phi[elec_range[1]];
//	     -l_coeffL[2]*phi[elec_range[1]+1];
    ghostL /= l_coeffL[0];

    ghostR = phi_right-l_coeffR[0]*phi[elec_range[2]-1];
//	     -l_coeffR[1]*phi[elec_range[2]];
    ghostR /= l_coeffR[1];


    /////////////////////////////////////////////////
    //
    // Compute electric field
    //
    //////////////////////////////////////////////////
    
    cout << "Computing electric field..." << endl;

    // Finite difference
    for (int i = elec_range[1]; i <= elec_range[2]; ++i) {
      // Left
      if (i == elec_range[1]) {
	E_field[i] = -(phi[i] - ghostL)/dx;
      } // Right
      else if (i == elec_range[2]) {
	E_field[i] = -(ghostR - phi[i-1])/dx;
      }
      else {
	E_field[i] = -(phi[i] - phi[i-1])/(dx);
      }
    }

    //////////////////////////////////////////////////
    //
    // Move particles
    //
    //////////////////////////////////////////////////

    cout << "Moving Particles..." << endl;

    for (int p = 0; p < electron.np; ++p) {      
      // Interpolate electric field at each particle
      get_WeightsNC(electron.x[p], weights, &electron.node_index[p], dx);
      part_E = E_field[electron.node_index[p]]*weights[0] +
	          E_field[electron.node_index[p] + 1]*weights[1];

      // Push electrons
      electron.vx[p] = electron.vx[p] + electron.q/(electron.m) *
	      	       part_E*dt;
      electron.x[p] = electron.x[p] + electron.vx[p]*dt;

      // Update energy in eV
      electron.epsilon[p] = 0.5*electron.m*pow(getv(electron.vx[p], 
			      electron.vy[p], electron.vz[p]),2.0)/e;

      // At boundary, if elastic collision, reflect particle
      // if not, absorb particle
      if (electron.x[p] < elec_range[1]*dx ||
	  electron.x[p] > elec_range[2]*dx) {
	
	// Elastic      
	if(double(rand())/RAND_MAX < electron.gamma[0]) {
	  electron.vx[p] = -electron.vx[p];
	  if (electron.x[p] < elec_range[1]*dx) {
	    electron.x[p] = elec_range[1]*dx+
		    (elec_range[1]*dx-electron.x[p]);
	  } else {
	    electron.x[p] = elec_range[2]*dx - 
		    (electron.x[p]-elec_range[2]*dx);
	  }
	}
       
	// Absorption
	else {
	  electron.remove_part(p);
	--p;
	} 
      }
     
      // If not at boundary, update max_epsilon if needed	
      else if (electron.epsilon[p] > electron.max_epsilon) {
        electron.max_epsilon = electron.epsilon[p];
      }
    }

    // Boundaries
    for (int p = 0; p < ion.np; ++p) {      
      // Interpolate electric field
      get_WeightsNC(ion.x[p], weights, &ion.node_index[p], dx);
      part_E = E_field[ion.node_index[p]]*weights[0] +
	          E_field[ion.node_index[p] + 1]*weights[1];

      // Push ions
      ion.vx[p] = ion.vx[p] + ion.q/(ion.m) *
	      	       part_E*dt;
      ion.x[p] = ion.x[p] + ion.vx[p]*dt;


      ion.epsilon[p] = 0.5*ion.m*pow(getv(ion.vx[p], 
			      ion.vy[p], ion.vz[p]),2.0)/e;

      // Electron emit or absorption
      if (ion.x[p] < elec_range[1]*dx ||
	  ion.x[p] > elec_range[2]*dx) {
	// Inject electron at thermal velocity
	if(double(rand())/RAND_MAX < ion.gamma[1]) {
	  electron.np += 1;
	  electron.x[electron.np-1] = ion.x[p]-ion.vx[p]*dt; // Previous ion location
	  electron.thermalVelocity(electron.np-1);
	  electron.epsilon[electron.np-1] = 0.5*electron.m*
		  		pow(getv(electron.vx[electron.np-1],
			        electron.vy[electron.np-1],
			       	electron.vz[electron.np-1]),2.0)/e;
	}

	if (ion.x[p] < elec_range[1]*dx) {neutral.flux_L += ion.spwt/(dt);}
	else {neutral.flux_R += ion.spwt/(dt);}
	ion.remove_part(p);
	--p;
      }
      
      // Update max epsilon if needed
      else if (ion.epsilon[p] > ion.max_epsilon) {
        ion.max_epsilon = ion.epsilon[p];
      }
    }


    //////////////////////////////////////////////////
    //
    // Collision modules
    //
    //////////////////////////////////////////////////
    
    cout << "Calculating collisions..." << endl;

    // Reset n_dot
    for (int i = 0; i < n_cell; ++i) {
      neutral.n_dot[i] = 0.0;
      excited.n_dot[i] = 0.0;
    }

    // Get number of particles for electron - neutral collisions
    if (electron.np > 0) {
    getNullCollPart(CS_energy, e_n_CS, electron.max_epsilon, 
		  &nu_max, &P_max, &N_c,
		  electron.m, neutral.getMax_n(n_cell), dt, electron.np, 
		  N_coll[0], data_set_length);
    } else {
      nu_max = 0.0;
      N_c = 0;
    }

    for (int i = 0; i < N_c; ++i) {
      rand_index = round((double(rand())/RAND_MAX)*(electron.np-1));
      type = getCollType(CS_energy, e_n_CS, electron.epsilon[rand_index],
		    nu_max, electron.m, neutral.getMax_n(n_cell),
		    N_coll[0], data_set_length);
    
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
		  electron.m, neutral.m);
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

	  // Update n_dot
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);
	  neutral.n_dot[electron.cell_index[rand_index]] -= 
		  electron.spwt*weights[0]/(dx*dt);
	  neutral.n_dot[electron.cell_index[rand_index] + 1] -=
		  electron.spwt*weights[1]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index]] += 
		  electron.spwt*weights[0]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index] + 1] +=
		  electron.spwt*weights[1]/(dx*dt);
	  continue;
        case 2:
	  epsilon_exc = 1.31041e1; // From crtrs.dat.txt
	  e_excitation(&electron.vx[rand_index], &electron.vy[rand_index],
		  &electron.vz[rand_index], electron.epsilon[rand_index],
		  epsilon_exc);
          electron.epsilon[rand_index] = 0.5*electron.m*
		  pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
		  electron.vz[rand_index]),2.0)/e;

	  // Update n_dot
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);
	  neutral.n_dot[electron.cell_index[rand_index]] -= 
		  electron.spwt*weights[0]/(dx*dt);
	  neutral.n_dot[electron.cell_index[rand_index] + 1] -=
		  electron.spwt*weights[1]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index]] += 
		  electron.spwt*weights[0]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index] + 1] +=
		  electron.spwt*weights[1]/(dx*dt);

	  continue;
        case 3:
	  electron.np += 1;
	  electron.vx[electron.np-1] = 0.0;
	  electron.vy[electron.np-1] = 0.0;
	  electron.vz[electron.np-1] = 0.0;
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

	  // Update n_dot
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);
	  neutral.n_dot[electron.cell_index[rand_index]] -= 
		  electron.spwt*weights[0]/(dt*dx);
	  neutral.n_dot[electron.cell_index[rand_index] + 1] -=
		  electron.spwt*weights[1]/(dt*dx);
	  continue;	
      case 4:
	  continue;
      }
    }
    
    // Ion-neutral collisions
    if (ion.np > 0) {
    getNullCollPart(CS_energy, i_n_CS, ion.max_epsilon, &nu_max, &P_max, &N_c,
		  ion.m, neutral.getMax_n(n_cell), dt, ion.np, N_coll[2], data_set_length);
    }  else {
      nu_max = 0.0;
      N_c = 0;
    }

    for (int i = 0; i < N_c; ++i) {
      rand_index = round((double(rand())/RAND_MAX)*(ion.np-1));
      type = getCollType(CS_energy, i_n_CS, ion.epsilon[rand_index],
		    nu_max, ion.m, neutral.getMax_n(n_cell), N_coll[2], data_set_length);
    
      // switch-case for electron-neutral collisions
      // 0 - isotropic
      // 1 - charge exchange
      // 2 - null
      switch(type) {
        case 0:
	  T_n = neutral.T[electron.cell_index[rand_index]]*weights[0] +
		     neutral.T[electron.cell_index[rand_index] + 1]*weights[1];

	  i_scattering(&ion.vx[rand_index], &ion.vy[rand_index],
		       &ion.vz[rand_index], ion.epsilon[rand_index],
		       ion.m, ion.m, T_n);
	  ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
			ion.vy[rand_index], ion.vz[rand_index]),2.0)/e; //[eV]
	  continue;
        case 1:
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);

	  T_n = neutral.T[electron.cell_index[rand_index]]*weights[0] +
		     neutral.T[electron.cell_index[rand_index] + 1]*weights[1];

	  thermalVelSample(&ion.vx[rand_index], &ion.vy[rand_index],
			  &ion.vz[rand_index], T_n, neutral.m);
	  ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
			ion.vy[rand_index], ion.vz[rand_index]),2.0)/e; //[eV]
	  continue;
        case 2:
	  continue;
      }
    }

    // Electron-excited
    if (electron.np > 0) {
    getNullCollPart(CS_energy, e_ex_CS, electron.max_epsilon, 
		    &nu_max, &P_max, &N_c,
		    electron.m, excited.getMax_n(n_cell), 
		    dt, electron.np,
		    N_coll[1], data_set_length);
    }  else {
      nu_max = 0.0;
      N_c = 0;
    }

    for (int i = 0; i < N_c; ++i) {
      rand_index = round((double(rand())/RAND_MAX)*(ion.np-1));
      type = getCollType(CS_energy, e_ex_CS, electron.epsilon[rand_index],
		    nu_max, electron.m, excited.getMax_n(n_cell), N_coll[1], data_set_length);
    
      // switch-case for electron-neutral collisions
      // 0 - Step Ionization
      // 1 - Superelastic (De-excitation)
      // 2 - null
      switch(type) {
        case 0:
	  // Update n_dot
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);
	  excited.n_dot[electron.cell_index[rand_index]] -= 
		  electron.spwt*weights[0]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index] + 1] -=
		  electron.spwt*weights[1]/(dx*dt);

	  continue;
        case 1:
	  // Update n_dot
	  get_WeightsCC(electron.x[rand_index], weights, 
			  &electron.cell_index[rand_index], elec_range, dx);
	  excited.n_dot[electron.cell_index[rand_index]] -= 
		  electron.spwt*weights[0]/(dx*dt);
	  excited.n_dot[electron.cell_index[rand_index] + 1] -=
		  electron.spwt*weights[1]/(dx*dt);
	  neutral.n_dot[electron.cell_index[rand_index]] += 
		  electron.spwt*weights[0]/(dx*dt);
	  neutral.n_dot[electron.cell_index[rand_index] + 1] +=
		  electron.spwt*weights[1]/(dx*dt);

	  continue;
        case 2:
	  continue;
      }
    }



    cout << "End of iteration" << endl;
    cout << "electron.np = " << electron.np << endl;
    cout << "ion.np = " << ion.np << endl;

    ////////////////////////////////////////////////////
    //
    // Output info
    //
    ////////////////////////////////////////////////////

    if ((iter+1)%write_iter == 0) {
      for (int i = 0; i < n_cell; ++i) {
        FieldCCFile << iter << " " << t << " " << dx*(i+0.5) << " ";
	FieldCCFile << rho[i] << " " << phi[i] << " " << neutral.n[i] << endl;
      }
      for (int i = 0; i < nn; ++i) {
	FieldNCFile << iter << " " << t << " " << dx*i << " ";
	FieldNCFile << E_field[i]  << endl;
      }

      NumFile << iter << " " << electron.np << " " << ion.np << endl;
     
    }
  }

  // Clean up
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
  delete(l_coeffL);
  delete(l_coeffR);
  delete(lagrangeLeft);
  delete(lagrangeRight);
  
  electron.clean();
  ion.clean();
  neutral.clean();
  excited.clean();

  return 0;
}

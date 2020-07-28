#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include "mpi.h"
#include <stdio.h>
#include <string>

#include "species.h"
#include "solverModules.h"
#include "collisionModules.h"

using namespace std;
using namespace std::chrono;


//////////////////////////////////////////////////////////////////////
//
//	1D_CCP
//	Runs 1D PIC/MCC low temperature simulation
//	
//	To run in comet: 
//	module purge
//	module load gnu openmpi_ib
//	./1D_CCP (Pressure)
//
//	For parallel, have to run a batch script
//	Made by: Raymond Lau for PDML
//
////////////////////////////////////////////////////////////////////

// Exclude zero and one, used for particle generation
double ranf(void)
{
	double dtmp = 1.0;
	while(dtmp>=1.0 || dtmp <= 0.0) {
		dtmp = ((double)rand())/((double)RAND_MAX);
	}
	return dtmp;
}

// Exclude 1
// Used for particle sample
double ranf0(void)
{ 
	double dtmp = 1.0;
	while(dtmp>=1.0 || dtmp < 0.0) {
		dtmp = ((double)rand())/((double)RAND_MAX);
	}
	return dtmp;
}

// Node centered weights
void get_WeightsNC(double part_x,
		 double *weights,
		 int *node_index,
		 double dx) {

*node_index = floor(part_x/dx);
double hx = (part_x -	(*node_index) * dx)/dx;
weights[0] = 1.0 - hx;
weights[1] = hx;
}

// Cell centered weights
void get_WeightsCC(double part_x,
		 double *weights,
		 int *cell_index,
		 int *elec_range,
		 double dx) {

	*cell_index = floor((part_x-0.5*dx)/dx);
	if (part_x < elec_range[1]*dx + 0.5*dx) {
		*cell_index = elec_range[1];
		weights[0] = 1.0;
		weights[1] = 0.0;
	} else if (part_x >= elec_range[2]*dx - 0.5*dx) {
		weights[0] = 1.0;
		weights[1] = 0.0;
	} else {
	double hx = (part_x -	(*cell_index + 0.5) * dx)/dx;
	weights[0] = 1.0 - hx;
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
//																		//
//							Main Loop									//
//																			//
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) 
{	
	
	MPI_Init(&argc, &argv);
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
 
	// Random Seed
	srand(mpi_rank*time(NULL));

	// Constants
	const double epsilon0 = 8.854e-12; // Permittivity of free space
	const double ech = 1.602e-19; // Elementary charge [C]
	const double k_B = 1.381e-23; // Boltzman constant [J/K]
	const double AMU = 1.661e-27; // Atomic Mass Unit [kg]
	const double N_A = 6.022e23; // Avogadros Number

	// Input Parameters
	//double P = 0.3; //[Pa]
	double P = atof(argv[2]);
	double T_gas = 300.0; // [K]
	double f_ground = 1.0;
	double f_excite = 5.0e-5;
	//double f_excite = 2.0e-4;
	double f_ion = 7.5e-6;
	double f_tot = f_ground + f_excite + f_ion;
	f_ground /= f_tot;
	f_excite /= f_tot;
	f_ion /= f_tot;

	// Problem discretization
	int nn = 191; // # of x1 nodes
	int n_cell = nn-1; // # of cells
	double x_start = 0.0;
	double x_end = 0.095;
	double x_bias = 0.01;
	double x_ground = 0.085;

	double L = x_end-x_start; // Domain length m, including electrodes
	double dx= L/((double) n_cell); // length of each cell
	double esmall = 1.0e-9;

	double L_inner = x_ground-x_bias;
	int		N_inner = (int) round(L_inner/dx + esmall);

	// Electrode info
	// Bias electrode left, grounded electrode right
	int *elec_range = new int[4]; 
	elec_range[0] = 0; // Left bias electrode
	elec_range[1] = round((x_bias-x_start)/dx + esmall);
	elec_range[2] = round((x_ground-x_start)/dx + esmall); // left side ground
	elec_range[3] = round ((x_end-x_start)/dx + esmall); // right side ground
	double f_hf = 27.12e6; // Hz
	double f_lf = 271.2e3; // Hz
	double V_hf = 110.0; // V
	double V_lf = 281.9; // V
	double phi_left;
	double phi_right = 0.0;

	// Time
	// Timesteps per HF * Cycles HF per LF * # LF
	int ts_hf = 200;
	int n_cycle = 120;
	int n_cycle_ss = 5;
	int hf_2_lf = 100; // Cycles HF per LF
	int ts = ts_hf*hf_2_lf*n_cycle; // # of time steps for main simulation
	int ts_ss = ts_hf*hf_2_lf*n_cycle_ss;			// # of steady state time steps for bulk analysis
	double dt = 1.0/(f_hf*((double)ts_hf)); 

	//////////////////////////////////////////////////////////////
	//
	// Load data for collisions
	//
	///////////////////////////////////////////////////////////
	if(mpi_rank==0){cout<<"Loading collision data"<<endl;}
	// See cstrs.dat.txt for what each column is
	int *N_coll = new int[3]; // Types of collisions for each reaction
	N_coll[0] = 4;
	N_coll[1] = 2;
	N_coll[2] = 2;

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
	if(mpi_rank==0){cout<<"Loading particle data"<<endl;}
	int max_part = 6e6; // max_particles if none leave during time history
	int total_init_part = atoi(argv[1]);
	int init_part = floor(total_init_part/((double)mpi_size)); // initial # particles per processor
	double real_part = f_ion*P*L_inner/(k_B*T_gas);

	int mod_part = total_init_part % mpi_size;
	if(mpi_rank < mod_part) {init_part += 1;}

	// Electron particle data
	particles electron;
	electron.initialize(max_part, n_cell);
	electron.m = 9.109e-31; //[kg]
	electron.q = -1.0*ech; //[C]
	electron.np = init_part;
	electron.T = 300.0; //[K]
	electron.spwt = real_part/((double)total_init_part);
	electron.gamma[0] = 0.0;
	electron.max_epsilon = 0.0;

	// Ion particle data
	particles ion;
	ion.initialize(max_part, n_cell);
	ion.m = 39.948*AMU; //[kg]
	ion.q = ech; //[c]
	ion.np = init_part;
	ion.T = 300.0; //[K]
	ion.spwt = real_part/((double)total_init_part);
	ion.gamma[1] = 0.15;
	ion.max_epsilon = 0.0;

	// Uniformly distribute initial positions, assign initial velocity from
	// thermal distribution
	for (int i = 0; i < electron.np; ++i) {
		electron.x[i] = ranf()*L_inner + x_bias;
		//cout << "electron.x[i] = " << electron.x[i] << endl;
		if (electron.x[i] < x_bias) { 
			cout << "error" << endl;}
		electron.thermalVelocity(i);
		electron.epsilon[i] = 0.5*electron.m*pow(getv(electron.vx[i], electron.vy[i],
							electron.vz[i]),2.0)/ech;

		if (electron.epsilon[i] > electron.max_epsilon) 
		{electron.max_epsilon = electron.epsilon[i];}
	}

	for (int i = 0; i < ion.np; ++i) {
		ion.x[i] = ranf()*L_inner + x_bias;
		ion.thermalVelocity(i);
		ion.epsilon[i] = 0.5*ion.m*pow(getv(ion.vx[i], ion.vy[i],
							ion.vz[i]),2.0)/ech;

		if (ion.epsilon[i] > ion.max_epsilon) 
		{ion.max_epsilon = ion.epsilon[i];}

	}

	// Calculate collisions setup
	double max_index, nu_max_en, nu_max_eex, nu_max_in, P_max;
	int N_c_en = 0;
	int N_c_in = 0;
	int N_c_eex = 0;
	int rand_index, type, null_count, elastic_count;
	double epsilon_exc, epsilon_ion;


	/////////////////////////////////////////////////////////
	// 
	// Field variables
	//
	////////////////////////////////////////////////////////////
	if(mpi_rank == 0){cout<<"Initializing field data" <<endl;}
	double R; // General random number variable
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
		neutral.D[i] = 0.0;
		neutral.n_dot[i] = 0.0;
		neutral.T[i] = 300.0;

		excited.D[i] = 0.0;
		excited.n_dot[i] = 0.0;
		excited.T[i] = 300.0;

		if (i >= elec_range[1] && i < elec_range[2]) {
			neutral.n[i] = f_ground*P/(k_B*neutral.T[i]);
			excited.n[i] = f_excite*P/(k_B*excited.T[i]);
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
	double t;	// time
	double *RHS = new double[n_cell]; // Right hand side
	double *a = new double[n_cell];	// Tridiagonal coeffs
	double *b = new double[n_cell];
	double *c = new double[n_cell];
	double ghostL, ghostR;

	// Electric field
	double *E_field = new double[nn]; // electric field at the node

	// Move Particles
	double part_E, T_n;
 
	// Collision modules
	int search_index;
	int ionized_num = 0;
	double sigma, penning_rate;
	int stepwise_coll = 0;
	//e-ex vars
	double R10, R11, Rtmp1, atmp1, vtmp1, theta1;

	auto start = steady_clock::now();
	auto stop = steady_clock::now();
	auto duration = duration_cast<microseconds>(stop-start);

	auto start_total = steady_clock::now();

	//////////////////////////////////////////////////////////
	//
	// Output files
	//
	//////////////////////////////////////////////////////////
	if(mpi_rank == 0) {cout << "Setting up output files" << endl;}
	int write_iter = 200; 	// Write every x number of iterations
	int write_restart = ts/5; 	// Write restart file every 20% of the simulation
						// Make it so this is a multiple of write_iter 
	int restart_opt = 1; 		// 0 to not write restart, 1 to write restart file
	//int write_ss_iter = 25;	// Write steady state infor every x # of iterations
	// New storage convention
	ofstream InputFile("../Input.txt");
	ofstream FieldCCFile("../FieldCCData.txt");
	ofstream FieldNCFile("../FieldNCData.txt");
	ofstream FieldAverageFile("../FieldAverageData.txt");
	ofstream NumFile("../NumberPart.txt");

	// MAKE ELECTRON AND ION FOLDERS IN SLURM FILE
	ofstream ElectronFile("../Electron/Electrons"+to_string(mpi_rank)+".txt");
	ofstream ElectronDiagFile("../Electron/ElectronDiagnostics.txt");
	ofstream IonFile("../Ion/Ions"+to_string(mpi_rank)+".txt");
	ofstream IonDiagFile("../Ion/IonDiagnostics.txt");
	ofstream IEDFFile("../Ion/IEDF"+to_string(mpi_rank)+".txt");
	ofstream TimerFile("../Timer.txt");
	 
	InputFile << "Misc comments: Slurm job, 0.3 Pa, 2e6 Particles" << endl;
	InputFile << "Pressure [Pa] / electron.np / ion.np / electron.spwt / ion.spwt / ";
	InputFile << "V_hf / V_lf / f_hf / f_lf / Total steps / Steady State steps / dt / NumNodes / NumRanks" << endl;
	InputFile << P << " " << total_init_part << " " << total_init_part << " " << electron.spwt;
	InputFile << " " << ion.spwt << " " << V_hf << " " << V_lf << " " << f_hf;
	InputFile << " " << f_lf << " " << ts << " " << ts_ss << " " <<	dt << " " << nn << " " << mpi_size << endl;

	FieldAverageFile << "Iteration / Neutral Density / Excited Density / ";
	FieldAverageFile << "Ion Density / Electron Density" << endl;

	FieldCCFile << "Iteration / Time / Cell x / Charge Density / ";
	FieldCCFile << "Electric Potential / Neutral Density / Excited Density / ";
	FieldCCFile << "Ion Density / Electron Density"<< endl;

	FieldNCFile << "Iteration / Time / Node x / Electric Field" << endl;
 
	if (restart_opt == 1) {
		ElectronFile << "Iteration / x / vx / vy / vz / Epsilon " << endl;
		IonFile << "Iteration / x / vx / vy / vz / Epsilon " << endl;
	}

	ElectronDiagFile << "Cell x / n / ux / uy / uz / mean epsilon / T" << endl;
	IonDiagFile << "Cell x / n / ux / uy / uz / mean epsilon / T" << endl;
	IEDFFile << "Iter / x / vx / vy / vz / Epsilon" << endl;

	NumFile << "Iteration / Electron np / Ion np / Electron Inner np / Ion Inner np / Total density (mass conservation)" << endl;

	TimerFile << "Time in microseconds. Iter / Coll / Push1 / Fluid / Rho / Phi / E / Push2";
	TimerFile << " / Total min elapsed "<< endl;

	//////////////////////////////////////////////////////////////////////////////
	//
	//	
	// Initial Conditions
	//
	// Get E^0
	// n^0 and v^0 already determined beforehand
	//
	//
	//////////////////////////////////////////////////////////////////////////////

	// Electrode Potential boundaries
	t = 0.0;
	phi_left = V_hf*sin(2*M_PI*f_hf*t) + V_lf*sin(2*M_PI*f_lf*t);
	// Reset variables
	null_count = 0;
	elastic_count = 0;

	for (int i = 0; i< n_cell; ++i) {
		phi[i] = 0.0;
		rho[i] = 0.0;
		electron.n[i] = 0.0;
		ion.n[i] = 0.0;
		E_field[i] = 0.0;
	}
	E_field[nn-1] = 0.0;

	/////////////////////////////////////////////////////////
	//
	// Compute charge density 
	//
	////////////////////////////////////////////////////////
		
	if(mpi_rank == 0){cout << "Computing initial charge density..." << endl;}

	// Get number densities
	// Electrons
	for (int p = 0; p < electron.np; ++p) {
		get_WeightsCC(electron.x[p], weights, &electron.cell_index[p],
				 elec_range, dx);

		electron.n[electron.cell_index[p]] += electron.spwt*weights[0];
		electron.n[electron.cell_index[p]+1] += electron.spwt*weights[1];
		 
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
		rho[i] = ech*(ion.n[i]-electron.n[i]);
	}

	// Collect total charge density between ranks
	MPI_Allreduce(MPI_IN_PLACE, electron.n, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ion.n, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	MPI_Allreduce(MPI_IN_PLACE, rho, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//////////////////////////////////////////////////
	//
	// Compute electric potential
	//
	/////////////////////////////////////////////////

	if(mpi_rank==0){cout << "Computing initial electric potential..." << endl;}

	// Tridiagonal coefficients
	for (int i = 0; i < n_cell; ++i) {
		if (i >= elec_range[1] && i < elec_range[2]) {			
			a[i] = 1.0;
			b[i] = -2.0;
			c[i] = 1.0;
			RHS[i] = -rho[i]/epsilon0*dx*dx;
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
	b[elec_range[1]] -= 1.0; 
	RHS[elec_range[1]] -= 2.0*phi_left;
 
	b[elec_range[2] - 1] -= 1.0;
	c[elec_range[2] - 1] = 0.0;
	RHS[elec_range[2]-1] -= 2.0*phi_right;

	// Solve
	triDiagSolver(phi, a, b, c, RHS, n_cell);
		
	// Get ghost node values
	ghostL = 2.0*phi_left - phi[elec_range[1]];
	ghostR = 2.0*phi_right - phi[elec_range[2]-1];

	/////////////////////////////////////////////////
	//
	// Compute electric field
	// phi^{n+1} -> E^{n+1}
	//
	//////////////////////////////////////////////////
		
	if(mpi_rank==0){cout << "Computing initial electric field..." << endl;}

	// Finite difference cell_interface
	for (int i = elec_range[1]; i <= elec_range[2]; ++i) {
		// Left
		if (i == elec_range[1]) {
			E_field[i] = -(phi[i] - ghostL)/dx;
		} // Right
		else if (i == elec_range[2]) {
			E_field[i] = -(ghostR - phi[i-1])/dx;
		}
		else if (i > elec_range[1] && i < elec_range[2]) {
			E_field[i] = -(phi[i] - phi[i-1])/(dx);
		}
		else {
			E_field[i] = 0.0;			 
		}
	}
	for (int p = 0; p < electron.np; ++p) {			
		// Interpolate electric field at each particle
		get_WeightsNC(electron.x[p], weights, &electron.node_index[p], dx);
		electron.E[p] = E_field[electron.node_index[p]]*weights[0] +
					 E_field[electron.node_index[p] + 1]*weights[1];
	}
	for (int p = 0; p < ion.np; ++p) {			
		// Interpolate electric field
		get_WeightsNC(ion.x[p], weights, &ion.node_index[p], dx);
		ion.E[p] = E_field[ion.node_index[p]]*weights[0] +
				 E_field[ion.node_index[p] + 1]*weights[1];
	}

	////////////////////////////////////////////////////
	//
	// Output info
	//
	////////////////////////////////////////////////////

	neutral.n_bar = 0.0;
	excited.n_bar = 0.0;
	ion.n_bar = 0.0;
	electron.n_bar = 0.0;

	for (int i = 0; i < electron.np; ++i) {
		if (electron.x[i] > (0.25*L_inner+x_bias) && electron.x[i] < (x_ground-0.25*L_inner)) {
			electron.inner_np += 1;
		}
	}

	for (int i = 0; i < ion.np; ++i) {
		if (ion.x[i] > (0.25*L_inner+x_bias) && ion.x[i] < (x_ground-0.25*L_inner)) {
			ion.inner_np += 1;
		}
	}


	// Energy density and flux density
	for (int i = 0; i < n_cell; ++i) {
		if (mpi_rank == 0) {
			neutral.n_bar += neutral.n[i];
			excited.n_bar += excited.n[i];
			ion.n_bar += ion.n[i];
			electron.n_bar += electron.n[i];

			FieldCCFile << 0 << " " << t << " " << dx*(i+0.5) << " ";
			FieldCCFile << rho[i] << " " << phi[i] << " " << neutral.n[i] << " ";
			FieldCCFile << excited.n[i] << " " << ion.n[i] << " ";
			FieldCCFile << electron.n[i] << endl;
		}
	}

	if(mpi_rank == 0) {
		for (int i = 0; i < nn; ++i) {
			FieldNCFile << 0 << " " << t << " " << dx*i << " ";
			FieldNCFile << E_field[i]	<< endl;
		}
	}
			
	neutral.n_bar /= elec_range[2] - elec_range[1];
	excited.n_bar /= elec_range[2] - elec_range[1];
	ion.n_bar /= elec_range[2] - elec_range[1];
	electron.n_bar /= elec_range[2] - elec_range[1];
			
	MPI_Reduce(&electron.np, &electron.np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&ion.np, &ion.np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&electron.inner_np, &electron.inner_np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&ion.inner_np, &ion.inner_np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	double dn_tot = 0.0;
	for(int i = 0; i < n_cell; ++i) {
		dn_tot += ion.n[i] + neutral.n[i] + excited.n[i];
	} 
	dn_tot /= (double) (elec_range[2]-elec_range[1]);

	if(mpi_rank == 0){
		FieldAverageFile << 0	<< " " << neutral.n_bar << " " << excited.n_bar;
		FieldAverageFile << " " << ion.n_bar << " " << electron.n_bar << endl;
		NumFile << 0 << " " << electron.np_total << " " << ion.np_total << " ";
		NumFile << electron.inner_np_total << " " << ion.inner_np_total << " " << dn_tot << endl;
	}
/////////////////////////////////////////////////////////////////////////////////////////////////	
	//////////////////////////////////////////////////////////
	//
	// Main Loop
	// -- v^n -> v^{n+1/2}
	/// -- x^n -> x^{n_1}
	// -- E^{n+1} based on x^{n+1}
	// -- v^{n+1/2} -> v^{n+1}
	//
	/////////////////////////////////////////////////////////
	
	for (int iter = 1; iter <= ts + ts_ss; ++iter) {
		//////////////////////////////////////////////////
		//
		// Collision modules
		//
		//////////////////////////////////////////////////

		if(mpi_rank==0){cout << "Iter: " << iter << endl;}

		start = steady_clock::now();

		//for (int i = 0; i < n_cell; ++i) {
		// Get number of particles for electron - neutral collisions
		if (electron.np > 0) {
			getNullCollPart(CS_energy, e_n_CS, electron.max_epsilon, 0.0,
							&nu_max_en, &P_max, &N_c_en,
							electron.m, neutral.getMax_n(n_cell), 
							dt, electron.np, 
							N_coll[0], data_set_length);
		} else {
			nu_max_en = 0.0;
			N_c_en = 0;
		}
		
		// Electron-excited
		if (electron.np > 0) {
			getNullCollPart(CS_energy, e_ex_CS, electron.max_epsilon, 0.0,
							&nu_max_eex, &P_max, &N_c_eex,
							electron.m, excited.getMax_n(n_cell), 
							dt, electron.np,
							N_coll[1], data_set_length);
		} else {
			nu_max_eex = 0.0;
			N_c_eex = 0;
		}

		// Ion-neutral collisions
		// Assume neutral at thermal velocity is comparable to ions,
		// max_relative_velocity = max_ion_velocity + v_th_ion
		if (ion.np > 0) {
			getNullCollPart(CS_energy, i_n_CS, ion.max_epsilon, 
					-sqrt(2.0*k_B*ion.T/ion.m),
					&nu_max_in, &P_max, &N_c_in,
					ion.m, neutral.getMax_n(n_cell), 
					dt, ion.np, 
					N_coll[2], data_set_length);
		} else {
			nu_max_in = 0.0;
			N_c_in = 0;
		} 

		if(mpi_rank==0){cout << "Calculating electron-neutral collisions..." << endl;}
		//cout << N_c_en << " " << electron.max_epsilon << endl;

		for (int i = 0; i < N_c_en; ++i) {
			rand_index = floor(ranf0()*(electron.np));
			if(rand_index == electron.np) {rand_index--;}
			
			while (isnan(electron.epsilon[rand_index])) {
					cout << "Error at 622" << endl;
					return -1;
			}
			// Use the neutral density of the cell that the randomly chose particle is in
			type = getCollType(CS_energy, e_n_CS, electron.epsilon[rand_index], 
				nu_max_en, electron.m, neutral.n[electron.node_index[rand_index]],
				N_coll[0], data_set_length);
			
			double elec_prev = electron.epsilon[rand_index];
			double elec_vx_prev = electron.vx[rand_index];
			double elec_vy_prev = electron.vy[rand_index];
			double elec_vz_prev = electron.vz[rand_index];
			
			// switch-case for electron-neutral collisions
			// 0 - elastic
			// 1 - excitation 1
			// 2 - excitation 2
			// 3 - ionization
			// 4 - null
			switch(type) {
				case 0:
					elastic_count += 1;
					e_elastic(&electron.vx[rand_index], &electron.vy[rand_index],
							&electron.vz[rand_index], electron.epsilon[rand_index],
							electron.m, neutral.m);
					electron.epsilon[rand_index] = 0.5*electron.m*
							pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
							electron.vz[rand_index]),2.0)/ech;

					while (isnan(electron.epsilon[rand_index])) {
						cout << "Error at 665" << endl;
						cout << rand_index << " " << electron.m << " " << neutral.m << " " << electron.np << endl;
						cout << elec_prev << " " << elec_vx_prev << " " << elec_vy_prev << " " <<	elec_vz_prev << endl;
						cout << electron.vx[rand_index] << " " << electron.vy[rand_index]<< " " << electron.vz[rand_index] << endl;
						return -1;
					}			

					continue;
				case 1:
					epsilon_exc = 1.159330e1; // From crtrs.dat.txt
					e_excitation(&electron.vx[rand_index], &electron.vy[rand_index],
							&electron.vz[rand_index], electron.epsilon[rand_index],
							epsilon_exc);
					electron.epsilon[rand_index] = 0.5*electron.m*
							pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
							electron.vz[rand_index]),2.0)/ech;

					while (isnan(electron.epsilon[rand_index])) {
							cout << "Error at 679" << endl;
							return -1;
					}	 


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
					epsilon_exc = 1.30940e1; // From crtrs.dat.txt
					e_excitation(&electron.vx[rand_index], &electron.vy[rand_index],
							&electron.vz[rand_index], electron.epsilon[rand_index],
							epsilon_exc);
					electron.epsilon[rand_index] = 0.5*electron.m*
							pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
							electron.vz[rand_index]),2.0)/ech;

					while (isnan(electron.epsilon[rand_index])) {
						cout << "Error at 705" << endl;
						return -1;
					}		


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
					ionized_num += 1;
					electron.np += 1;
					electron.x[electron.np-1] = electron.x[rand_index];
					electron.vx[electron.np-1] = 0.0;
					electron.vy[electron.np-1] = 0.0;
					electron.vz[electron.np-1] = 0.0;
					epsilon_ion = 1.599550e1; // From crtrs.dat.txt
					e_ionization(&electron.vx[rand_index], &electron.vy[rand_index],
							&electron.vz[rand_index], &electron.vx[electron.np-1], 
							&electron.vy[electron.np-1], &electron.vz[electron.np-1], 
					electron.epsilon[rand_index], epsilon_ion);
					electron.epsilon[rand_index] = 0.5*electron.m*pow(getv(electron.vx[rand_index], 
					electron.vy[rand_index], electron.vz[rand_index]),2.0)/ech; //[eV]
					electron.epsilon[electron.np-1] = 0.5*electron.m*
							pow(getv(electron.vx[electron.np-1],
					electron.vy[electron.np-1], 
	 				electron.vz[electron.np-1]),2.0)/ech; //[eV]

					// Create ion
					ion.np += 1;
					ion.x[ion.np-1] = electron.x[rand_index];
					ion.thermalVelocity(ion.np-1);
					ion.epsilon[ion.np-1] = 0.5*ion.m*pow(getv(ion.vx[ion.np-1], ion.vy[ion.np-1],
											ion.vz[ion.np-1]),2.0)/ech;


					while (isnan(electron.epsilon[rand_index])) {
						cout << "Error at 740" << endl;
						return -1;
					}

					while (isnan(electron.epsilon[electron.np-1])) {
						cout << "Error at 744" << endl;
						return -1;
					}			

					// Update n_dot
					get_WeightsCC(electron.x[rand_index], weights, 
						&electron.cell_index[rand_index], elec_range, dx);
					neutral.n_dot[electron.cell_index[rand_index]] -= 
						electron.spwt*weights[0]/(dt*dx);
					neutral.n_dot[electron.cell_index[rand_index] + 1] -=
						electron.spwt*weights[1]/(dt*dx);
					continue;	
				case 4:
					null_count += 1;
					continue;
			}
		}
		//cout << elastic_count << " " << null_count << endl;

		if(mpi_rank==0){cout << "Calculating electron-excited collisions..." << endl;}
		
		//cout << N_c_eex << endl;

		elastic_count = 0;
		null_count = 0;
		for (int i = 0; i < N_c_eex; ++i) {
			rand_index = floor(ranf0()*(electron.np));
			if(rand_index == electron.np) {rand_index--;}
			
			type = getCollType(CS_energy, e_ex_CS, electron.epsilon[rand_index],
				nu_max_eex, electron.m, excited.n[electron.node_index[rand_index]],
				N_coll[1], data_set_length);

			// switch-case for electron-excited collisions
			// 0 - Step Ionization
			/// 1 - Superelastic (De-excitation)
			// 2 - null
			switch(type) {
				case 0:
					electron.np += 1;
					electron.x[electron.np-1] = electron.x[rand_index];
					electron.vx[electron.np-1] = 0.0;
					electron.vy[electron.np-1] = 0.0;
					electron.vz[electron.np-1] = 0.0;
					epsilon_ion = 4.425; // From crtrs.dat.txt
					e_ionization(&electron.vx[rand_index], &electron.vy[rand_index],
							&electron.vz[rand_index], &electron.vx[electron.np-1], 
							&electron.vy[electron.np-1], &electron.vz[electron.np-1], 
							electron.epsilon[rand_index], epsilon_ion);
					electron.epsilon[rand_index] = 0.5*electron.m*pow(getv(electron.vx[rand_index], 
					electron.vy[rand_index], electron.vz[rand_index]),2.0)/ech; //[eV]
					electron.epsilon[electron.np-1] = 0.5*electron.m*
							pow(getv(electron.vx[electron.np-1],
					electron.vy[electron.np-1], 
			 		electron.vz[electron.np-1]),2.0)/ech; //[eV]

					while (isnan(electron.epsilon[rand_index])) {
						cout << "Error at 854" << endl;
						return -1;
					}
					while (isnan(electron.epsilon[electron.np-1])) {
						cout << "Error at 854" << endl;
						return -1;
					}

					// Create ion	
					ion.np += 1;
					ion.x[ion.np-1] = electron.x[rand_index];
					ion.thermalVelocity(ion.np-1);
					ion.epsilon[ion.np-1] = 0.5*ion.m*pow(getv(ion.vx[ion.np-1], ion.vy[ion.np-1],
											ion.vz[ion.np-1]),2.0)/ech;
					// Update n_dot
					get_WeightsCC(electron.x[rand_index], weights, 
							&electron.cell_index[rand_index], elec_range, dx);
					excited.n_dot[electron.cell_index[rand_index]] -= 
						electron.spwt*weights[0]/(dx*dt);
					excited.n_dot[electron.cell_index[rand_index] + 1] -=
						electron.spwt*weights[1]/(dx*dt);
					stepwise_coll += 1;
					continue;
				case 1:
					// Add excitation energy to the electron and isotropically scatter
					electron.epsilon[rand_index] += 1.30940e1;
					R10 = ranf();
					R11 = ranf();
					theta1 = 2.0*M_PI*R10;
					Rtmp1 = -1.0 + 2.0*R11;
					atmp1 = sqrt(1.0-Rtmp1*Rtmp1);
					vtmp1 = sqrt(2.0*ech*electron.epsilon[rand_index]/electron.m);
			
					electron.vx[rand_index] = vtmp1*atmp1*cos(theta1);
					electron.vy[rand_index] = vtmp1*atmp1*sin(theta1);
					electron.vz[rand_index] = vtmp1*Rtmp1;
					
					electron.epsilon[rand_index] = 0.5*electron.m*
							pow(getv(electron.vx[rand_index], electron.vy[rand_index], 
							electron.vz[rand_index]),2.0)/ech;
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
					null_count += 1;
					continue;
			}
		}

		//cout << elastic_count << " " << null_count << endl;
		if(mpi_rank==0){cout << "Calculating ion-neutral collisions..." << endl;}
		//cout << N_c_in << endl;

		elastic_count = 0;
		null_count = 0;

		for (int i = 0; i < N_c_in; ++i) {
			rand_index = floor(ranf0()*(ion.np));
			if(rand_index == ion.np) {rand_index--;}
			type = getCollType(CS_energy, i_n_CS, ion.epsilon[rand_index],
				nu_max_in, ion.m, neutral.n[ion.node_index[rand_index]],
				N_coll[2], data_set_length);
			// switch-case for electron-neutral collisions
			// 0 - isotropic
			// 1 - charge exchange
			// 2 - null
			switch(type) {
				case 0:
					elastic_count += 1;
					get_WeightsCC(electron.x[rand_index], weights, 
							&electron.cell_index[rand_index], elec_range, dx);
	
					T_n = neutral.T[electron.cell_index[rand_index]]*weights[0] +
							 neutral.T[electron.cell_index[rand_index] + 1]*weights[1];

					i_scattering(&ion.vx[rand_index], &ion.vy[rand_index],
							 &ion.vz[rand_index], ion.epsilon[rand_index],
							 ion.m, ion.m, T_n);
					ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
					ion.vy[rand_index], ion.vz[rand_index]),2.0)/ech; //[eV]
					continue;
				case 1:
					get_WeightsCC(electron.x[rand_index], weights, 
						&electron.cell_index[rand_index], elec_range, dx);

					T_n = neutral.T[electron.cell_index[rand_index]]*weights[0] +
						 neutral.T[electron.cell_index[rand_index] + 1]*weights[1];

					thermalVelSample(&ion.vx[rand_index], &ion.vy[rand_index],
						&ion.vz[rand_index], T_n, neutral.m);
					ion.epsilon[rand_index] = 0.5*ion.m*pow(getv(ion.vx[rand_index], 
					ion.vy[rand_index], ion.vz[rand_index]),2.0)/ech; //[eV]
					continue;
				case 2:
					null_count += 1;
					continue;
			}
		}
		//cout << elastic_count << " " << null_count << endl;

		// Penning Ionization
		for (int i = elec_range[1]; i < elec_range[2]; ++i) {
			penning_rate = 5.0e-17*pow(excited.n[i], 2.0); // AR* consumed [real part/m^3s]
			penning_rate /= ((double) mpi_size);
			double new_ion_double = penning_rate*dx*dt/ion.spwt;
			int new_ion = floor(penning_rate*dx*dt/ion.spwt); // Ions formed [macro part]
			double d_frac_ion = new_ion_double - new_ion;

			R = ranf();
			if (R < d_frac_ion) {new_ion++;}

			for (int p = 0; p < new_ion; ++p) {
				ion.np += 1;
				ion.x[ion.np-1] = ranf()*dx + dx*i;
				ion.thermalVelocity(ion.np-1);
				ion.epsilon[ion.np-1] = 0.5*ion.m*pow(getv(ion.vx[ion.np-1], ion.vy[ion.np-1],
										ion.vz[ion.np-1]),2.0)/ech;

				electron.np += 1;
				electron.x[electron.np-1] = ranf()*dx + dx*i;

				// AR* + AR* -> AR+ + AR
				// 2*excitation - ionization
				// excitation = average of 2 excitation energies
				electron.epsilon[electron.np-1] = (1.30940e1+1.159330e1)-1.599550e1;
				R10 = ranf();
				R11 = ranf();
				theta1 = 2.0*M_PI*R10;
				Rtmp1 = -1.0 + 2.0*R11;
				atmp1 = sqrt(1.0-Rtmp1*Rtmp1);
				vtmp1 = sqrt(2.0*ech*electron.epsilon[rand_index]/electron.m);
	
				electron.vx[electron.np-1] = vtmp1*atmp1*cos(theta1);
				electron.vy[electron.np-1] = vtmp1*atmp1*sin(theta1);
				electron.vz[electron.np-1] = vtmp1*Rtmp1;
	
				// Get actual kinetic energy
				electron.epsilon[electron.np-1] = 0.5*electron.m*pow(getv(electron.vx[electron.np-1], 
						electron.vy[electron.np-1], electron.vz[electron.np-1]),2.0)/ech;
			}
			neutral.n_dot[i] += penning_rate; // Increase in density
			excited.n_dot[i] -= 2.0*penning_rate; // Decrease in density
		}
		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_coll = duration.count();
				
		/////////////////////////////////////////////////
		//
		//	Reset parameters
		//
		///////////////////////////////////////////////////

		// Electrode Potential boundaries
		t = iter*dt;
		phi_left = V_hf*sin(2*M_PI*f_hf*t) + V_lf*sin(2*M_PI*f_lf*t);
		// Reset variables
		electron.max_epsilon = 0.0;
		ion.max_epsilon = 0.0;
		null_count = 0;
		elastic_count = 0;

		//////////////////////////////////////////////////
		//
		// Accelerate/Move particles
		// -- n to n+1
		//
		// v^n -> v^{n+1/2}
		// x^n -> x^{n+1}
		//
		//////////////////////////////////////////////////

		if(mpi_rank==0){cout << "Moving Particles..." << endl;}
		start = steady_clock::now();

		for (int p = 0; p < electron.np; ++p) {			
			get_WeightsNC(electron.x[p], weights, &electron.node_index[p], dx);
			electron.E[p] = E_field[electron.node_index[p]]*weights[0] +
				 E_field[electron.node_index[p] + 1]*weights[1];
		 
			// Push electrons
			electron.vx[p] = electron.vx[p] + electron.q/(electron.m) *
								 electron.E[p]*(0.5*dt);
			electron.x[p] = electron.x[p] + electron.vx[p]*dt;


			// At boundary, if elastic collision, reflect particle
			// if not, absorb particle
			if (electron.x[p] <= x_bias ||
				electron.x[p] >= x_ground) {
				// Absorption
				if (electron.x[p] <= x_bias) {electron.flux_L -= electron.spwt/dt;}
				else {electron.flux_R += electron.spwt/dt;}
				electron.remove_part(p);
			 	--p; 
			}
		}



		for (int p = 0; p < ion.np; ++p) {			
			// Interpolate electric field
			get_WeightsNC(ion.x[p], weights, &ion.node_index[p], dx);
			ion.E[p] = E_field[ion.node_index[p]]*weights[0] +
				 E_field[ion.node_index[p] + 1]*weights[1];

			// Push ions
			ion.vx[p] = ion.vx[p] + ion.q/(ion.m)*ion.E[p]*(0.5*dt);
			ion.x[p] = ion.x[p] + ion.vx[p]*dt;


			ion.epsilon[p] = 0.5*ion.m*pow(getv(ion.vx[p], 
					 ion.vy[p], ion.vz[p]),2.0)/ech;
			
			double x_old;
			// Electron emit or absorption
			if (ion.x[p] <= x_bias ||
				ion.x[p] >= x_ground) {

				R = ranf();

				// Inject electron at thermal velocity
				if (R < ion.gamma[1]) {
					electron.np += 1;
					electron.injectionVelocity(electron.np-1);

					if (ion.x[p] <= x_bias) {
						// Time after ion collides at wall
						time_coll = (ion.x[p] - x_bias)/ion.vx[p];
						electron.x[electron.np-1] = x_bias + electron.vx[electron.np-1]*time_coll;
					}
	
					if (ion.x[p] >= x_ground) {
						time_coll = (ion.x[p] - x_ground)/ion.vx[p];
						electron.vx[electron.np-1] *= -1.0;
						electron.x[electron.np-1] = x_ground + electron.vx[electron.np-1]*time_coll;
					}

					electron.epsilon[electron.np-1] = 0.5*electron.m*
							pow(getv(electron.vx[electron.np-1],
							electron.vy[electron.np-1],
						 	electron.vz[electron.np-1]),2.0)/ech;
				}		

				// AR+ flux, it is subtracted in the fluid solver
				if (ion.x[p] <= x_bias) {
					ion.flux_L -= ion.spwt/(dt);
					if (iter > ts) {
						x_old = ion.x[p] - dt*ion.vx[p];
						IEDFFile << iter << " " << x_old << " " << ion.vx[p] << " " << ion.vy[p] << " ";
						IEDFFile << ion.vz[p] << " " << ion.epsilon[p] << endl;
					}
				} else {
					ion.flux_R += ion.spwt/(dt);
					if (iter > ts) {
						x_old = ion.x[p] - dt*ion.vx[p];
						IEDFFile << iter << " " << x_old << " " << ion.vx[p] << " " << ion.vy[p] << " ";
						IEDFFile << ion.vz[p] << " " << ion.epsilon[p] << endl;
					}
				}

				ion.remove_part(p);
				--p;
			}			 
			// Update max epsilon if needed
			else if (ion.epsilon[p] > ion.max_epsilon) {
				ion.max_epsilon = ion.epsilon[p];
			}
		}
		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		auto time_push1 = duration.count();

		/////////////////////////////////////////////////////////
		//
		// Compute fluid density
		// excited.n and neutral.n from n to n+1
		// ion.n, neutral.n excited.n at nth time step
		// n_dots also at nth time step
		//
		////////////////////////////////////////////////////////
		
		if(mpi_rank==0){cout << "Computing fluid density" << endl;}
		start = steady_clock::now();

		double n_tot;
		double von_neumann;
		int C_factor = 2; // Run fluid solver ever C PIC time steps
		double fluid_dt = C_factor*dt;


		if (iter%(C_factor) == 0) {

			// Get densities for fluids
			// ion.n[i] should already be reduced
			// Gather n_dot from collisions
			MPI_Allreduce(MPI_IN_PLACE, excited.n_dot, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, neutral.n_dot, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			for (int i = 0; i < n_cell; ++i) {
				if (i >= elec_range[1] && i < elec_range[2]) {
					n_tot =	neutral.n[i] + excited.n[i] + ion.n[i];
					neutral.D[i] = 3.0*sqrt(k_B*neutral.T[i]) /
						(8.0*sqrt(M_PI*neutral.m)*d_LJ*d_LJ*1.0e-20*n_tot);
					excited.D[i] = 3.0*sqrt(k_B*excited.T[i]) /
						(8.0*sqrt(M_PI*excited.m)*d_LJ*d_LJ*1.0e-20*n_tot);
					von_neumann = neutral.D[i]*fluid_dt/(dx*dx);
					cout << von_neumann << endl;
				} else {
					neutral.D[i] = 0.0;
					excited.D[i] = 0.0;
				}
			}
			driftDiffusionFVExplicit2(&excited.n[elec_range[1]], 
									&excited.n_dot[elec_range[1]],
		 							&excited.D[elec_range[1]],
									0.0, 0.0, dx, fluid_dt, elec_range[2]-elec_range[1],
									&excited.flux_L, &excited.flux_R);

			double n_ghostL = neutral.n[elec_range[1]];
			double n_ghostR = neutral.n[elec_range[2]-1];

			driftDiffusionFVExplicit2(&neutral.n[elec_range[1]],
							&neutral.n_dot[elec_range[1]], 
							&neutral.D[elec_range[1]],
							n_ghostL, n_ghostR, 
							dx, fluid_dt, elec_range[2]-elec_range[1],
							&neutral.flux_L, &neutral.flux_R);

			neutral.n[elec_range[1]] += -dt/dx*excited.flux_L;
			neutral.n[elec_range[2]-1] += dt/dx*excited.flux_R;

			// Gather the fluxes from all the ranks
			MPI_Allreduce(MPI_IN_PLACE, &electron.flux_L, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, &electron.flux_R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, &ion.flux_L, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, &ion.flux_R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			neutral.n[elec_range[1]] += -dt/dx*ion.flux_L;
			neutral.n[elec_range[2]-1] += dt/dx*ion.flux_R;

			while(abs(neutral.flux_L) > 1e-5 || abs(neutral.flux_R) > 1e-5) {
				cout << "Error: We have set Neumann" << endl;
				return -1; 
			}
			while (excited.flux_L > esmall || excited.flux_R < -esmall) {
				cout << "Error: Excited is generated" << endl;
				return -1;
			}
			while (ion.flux_L > esmall || ion.flux_R < -esmall) {
				cout << "Error: Ion is generated" << endl;
				return -1;
			}

			// Reset boundaries
			ion.flux_L = 0.0;
			ion.flux_R = 0.0;		
			electron.flux_L = 0.0;
			electron.flux_R = 0.0;
			
			// Reset n_dot
			for (int i = 0; i < n_cell; ++i) {
				neutral.n_dot[i] = 0.0;
				excited.n_dot[i] = 0.0;
			}
		}

		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_fluid = duration.count();


		/////////////////////////////////////////////////////////
		//
		// Compute charge density 
		//
		////////////////////////////////////////////////////////
		
		if(mpi_rank==0){cout << "Computing charge density..." << endl;}

		start = steady_clock::now();

		for (int i = 0; i< n_cell; ++i) {
			electron.n[i] = 0.0;
			ion.n[i] = 0.0;
		}

		// Get number densities
		// Electrons
		for (int p = 0; p < electron.np; ++p) {
			get_WeightsCC(electron.x[p], weights, &electron.cell_index[p],
				 elec_range, dx);

			electron.n[electron.cell_index[p]] += electron.spwt*weights[0];
			electron.n[electron.cell_index[p]+1] += electron.spwt*weights[1];
		 
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
			rho[i] = ech*(ion.n[i]-electron.n[i]);
		}

		MPI_Allreduce(MPI_IN_PLACE, electron.n, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, ion.n, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, rho, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_rho = duration.count();

		//////////////////////////////////////////////////
		//
		// Compute electric potential
		//
		/////////////////////////////////////////////////

		if(mpi_rank==0){cout << "Computing electric potential..." << endl;}

		start = steady_clock::now();

		// Tridiagonal coefficients
		for (int i = 0; i < n_cell; ++i) {
			if (i >= elec_range[1] && i < elec_range[2]) {			
				a[i] = 1.0;
				b[i] = -2.0;
				c[i] = 1.0;
				RHS[i] = -rho[i]/epsilon0*dx*dx;
			}
			else {
				a[i] = 0.0;
				b[i] = 1.0;
				c[i] = 0.0;
				if (i < elec_range[1]) {RHS[i] = phi_left;}
				else {RHS[i] = phi_right;}
			}
		}

		// Adjust RHS to account for BC
		a[elec_range[1]] = 0.0;
		b[elec_range[1]] -= 1.0;
		RHS[elec_range[1]] -= 2.0*phi_left;

		b[elec_range[2] - 1] -= 1.0;
		c[elec_range[2] - 1] = 0.0;
		RHS[elec_range[2]-1] -= 2.0*phi_right;
		 
		// Solve
		triDiagSolver(phi, a, b, c, RHS, n_cell);
		
		// Get ghost node values
		ghostL = 2.0*phi_left - phi[elec_range[1]];
		ghostR = 2.0*phi_right - phi[elec_range[2]-1];

		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_phi = duration.count();


		/////////////////////////////////////////////////
		//
		// Compute electric field
		// phi^{n+1} -> E^{n+1}
		//
		//////////////////////////////////////////////////
		
		if(mpi_rank==0){cout << "Computing electric field..." << endl;}

		start = steady_clock::now();

		// Finite difference cell_interface
		for (int i = elec_range[1]; i <= elec_range[2]; ++i) {
			// Left
			if (i == elec_range[1]) {
				E_field[i] = -(phi[i] - ghostL)/dx;
			} // Right
			else if (i == elec_range[2]) {
				E_field[i] = -(ghostR - phi[i-1])/dx;
			}
			else if (i > elec_range[1] && i < elec_range[2]) {
				E_field[i] = -(phi[i] - phi[i-1])/(dx);
			}
			else {
				E_field[i] = 0.0;			 
			}
		}

		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_E = duration.count();

		//////////////////////////////////////////////////
		//
		// Accelerate particles
		// Part 2
		// v^{n+1/2} -> v^{n+1}
		//
		//////////////////////////////////////////////////

		if(mpi_rank == 0){cout << "Moving Particles..." << endl;}

		start = steady_clock::now();

		for (int p = 0; p < electron.np; ++p) {			
			// Interpolate electric field at each particle
			get_WeightsNC(electron.x[p], weights, &electron.node_index[p], dx);
			electron.E[p] = E_field[electron.node_index[p]]*weights[0] +
						 E_field[electron.node_index[p] + 1]*weights[1];

			// Push electrons
			electron.vx[p] = electron.vx[p] + electron.q/(electron.m) *
								 electron.E[p]*(0.5*dt);

			// Update energy in eV
			electron.epsilon[p] = 0.5*electron.m*pow(getv(electron.vx[p], 
						electron.vy[p], electron.vz[p]),2.0)/ech;

		 
			// Update max_epsilon if needed	
			if (electron.epsilon[p] > electron.max_epsilon) {
				electron.max_epsilon = electron.epsilon[p];
			}
		}

		for (int p = 0; p < ion.np; ++p) {			
			// Interpolate electric field
			get_WeightsNC(ion.x[p], weights, &ion.node_index[p], dx);
			ion.E[p] = E_field[ion.node_index[p]]*weights[0] +
					 E_field[ion.node_index[p] + 1]*weights[1];

			// Push ions
			ion.vx[p] = ion.vx[p] + ion.q/(ion.m) * ion.E[p]*(0.5*dt);


			ion.epsilon[p] = 0.5*ion.m*pow(getv(ion.vx[p], ion.vy[p], ion.vz[p]),2.0)/ech;
			
		if (ion.epsilon[p] > ion.max_epsilon) {
				ion.max_epsilon = ion.epsilon[p];
			}
		}

		stop = steady_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		auto time_push2 = duration.count();

		
		// Check mass conservation
		dn_tot = 0.0;
		for(int k = elec_range[1]; k < elec_range[2]; ++k) {
			dn_tot += ion.n[k]+neutral.n[k]+excited.n[k];
		}
		dn_tot /= ((double)(elec_range[2]-elec_range[1]));

		/*
		// Check charge conservation
		double dn_tote = 0.0;
		for(int k = elec_range[1]; k < elec_range[2]; ++k) {
			dn_tote += (ion.n[k]+dt/dx*(ion.flux_R-ion.flux_L)) -
					 (electron.n[k] + dt/dx*(electron.flux_R-electron.flux_L));
		}*/
	
		// Get totals for output

		////////////////////////////////////////////////////
		//
		// Output info
		//
		////////////////////////////////////////////////////
		// Display end of iteration
		auto stop_total = steady_clock::now();
		auto duration_total = duration_cast<minutes>(stop_total-start_total);
		if(mpi_rank == 0) {
			cout << "End of iteration" << endl;
			cout << "Mass conservation check: vol average density " << dn_tot << endl;
			cout << "Total Stepwise Ionization Count: " << stepwise_coll << endl;
	//cout << "Charge conservation check: total charge - lost to wall " << dn_tote << endl;
		}

	//////////////// Main Writing output ///////////////////////////////////////////
		if ((iter)%write_iter == 0 && iter <= ts) {
			// Total particles
			MPI_Reduce(&electron.np, &electron.np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&ion.np, &ion.np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			
			if(mpi_rank == 0) {
			cout << "electron.np = " << electron.np_total << endl;
			cout << "ion.np = " << ion.np_total << endl;
			}

			neutral.n_bar = 0.0;
			excited.n_bar = 0.0;
			ion.n_bar = 0.0;
			electron.n_bar = 0.0;

			ion.inner_np = 0;
			electron.inner_np = 0;

			//	Particle Output
			for (int i = 0; i < electron.np; ++i) {
				if (restart_opt == 1 && (iter)%write_restart == 0) {
					ElectronFile << iter << " " << electron.x[i] << " " << electron.vx[i] << " ";
					ElectronFile << electron.vy[i] << " " << electron.vz[i] << " " << electron.epsilon[i] << endl;
				}

				if (electron.x[i] > (0.25*L_inner+x_bias) && electron.x[i] < (x_ground-0.25*L_inner)) {
					electron.inner_np += 1;
				}
			}
			for (int i = 0; i < ion.np; ++i) {
				if (restart_opt == 1 && (iter)%write_restart == 0) {
					IonFile << iter << " " << ion.x[i] << " " << ion.vx[i] << " ";
					IonFile << ion.vy[i] << " " << ion.vz[i] << " " << ion.epsilon[i] << endl;
				}
				if (ion.x[i] > (0.25*L_inner+x_bias) && ion.x[i] < (x_ground-0.25*L_inner)) {
					ion.inner_np += 1;
				}
			}

			// Field Output
			for (int i = 0; i < n_cell; ++i) {
				if (mpi_rank == 0) {
					neutral.n_bar += neutral.n[i];
					excited.n_bar += excited.n[i];
					ion.n_bar += ion.n[i];
					electron.n_bar += electron.n[i];
	
					FieldCCFile << iter << " " << t << " " << dx*(i+0.5) << " ";
					FieldCCFile << rho[i] << " " << phi[i] << " " << neutral.n[i] << " ";
					FieldCCFile << excited.n[i] << " " << ion.n[i] << " ";
					FieldCCFile << electron.n[i] << endl;
				}
			}
			if (mpi_rank == 0) {
				for (int i = 0; i < nn; ++i) {
					FieldNCFile << iter << " " << t << " " << dx*i << " ";
					FieldNCFile << E_field[i]	<< endl;
				}
			}
			
			neutral.n_bar /= elec_range[2] - elec_range[1];
			excited.n_bar /= elec_range[2] - elec_range[1];
			ion.n_bar /= elec_range[2] - elec_range[1];
			electron.n_bar /= elec_range[2] - elec_range[1];
			
			// Inner number of particles
			MPI_Reduce(&electron.inner_np, &electron.inner_np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&ion.inner_np, &ion.inner_np_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

			// Rest of Output
			if (mpi_rank == 0) {
				FieldAverageFile << iter << " " << neutral.n_bar << " " << excited.n_bar;
				FieldAverageFile << " " << ion.n_bar << " " << electron.n_bar << " " << endl;
				NumFile << iter << " " << electron.np_total << " " << ion.np_total << " ";
				NumFile << electron.inner_np_total << " " << ion.inner_np_total << " " << dn_tot << endl;
				TimerFile << iter << " " << time_coll << " " << time_push1 << " ";
				TimerFile << time_fluid << " " << time_rho << " " << time_phi << " ";
				TimerFile << time_E << " " << time_push2 << " " << duration_total.count() << endl;
				cout << "Total minutes: " << duration_total.count() << endl;
			}
		}
		///////////////////// Steady State output //////////////////////////////////////////////
		if (iter > ts) {
			if(mpi_rank == 0) {
			cout << "Diagnostic calculations" << endl;
			}

			//	Particle Output
			/*
			for (int i = 0; i < electron.np; ++i) {
	if (restart_opt == 1 && (iter)%write_restart == 0) {
		ElectronFile << iter << " " << electron.x[i] << " " << electron.vx[i] << " ";
		ElectronFile << electron.vy[i] << " " << electron.vz[i] << " " << electron.epsilon[i] << endl;
	}
			}
			for (int i = 0; i < ion.np; ++i) {
	if (restart_opt == 1 && (iter)%write_restart == 0) {
		IonFile << iter << " " << ion.x[i] << " " << ion.vx[i] << " ";
		IonFile << ion.vy[i] << " " << ion.vz[i] << " " << ion.epsilon[i] << endl;
	}
			}*/

			// Get n, ux, uy, uz, mean energy
			for (int p = 0; p < electron.np; ++p) {
	get_WeightsCC(electron.x[p], weights, &electron.cell_index[p], elec_range, dx);
				electron.getDiagnosticsLocal(weights, p, dx);
			}		 
			for(int p = 0; p < ion.np; ++p) {
	get_WeightsCC(ion.x[p], weights, &ion.cell_index[p], elec_range, dx);
				ion.getDiagnosticsLocal(weights, p, dx);
			}		 
		}
	}
	
	// Gather ranks
	MPI_Allreduce(MPI_IN_PLACE, electron.n_ss, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, electron.mean_epsilon, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	MPI_Allreduce(MPI_IN_PLACE, electron.ux, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, electron.uy, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, electron.uz, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, ion.n_ss, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	MPI_Allreduce(MPI_IN_PLACE, ion.mean_epsilon, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ion.ux, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ion.uy, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ion.uz, n_cell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// Normalize by density and number of time steps averaged from 
	for (int i = 0; i < n_cell; ++i) {
		if (electron.n_ss[i] > 1.0) {
			electron.n_ss[i] /= ts_ss;
			electron.ux[i] /= electron.n_ss[i]*ts_ss;
			electron.uy[i] /= electron.n_ss[i]*ts_ss;
			electron.uz[i] /= electron.n_ss[i]*ts_ss;
			electron.mean_epsilon[i] /= electron.n_ss[i]*ts_ss;
		}
		if (ion.n_ss[i] > 1.0) {
			ion.n_ss[i] /= ts_ss;
			ion.ux[i] /= ion.n_ss[i]*ts_ss;
			ion.uy[i] /= ion.n_ss[i]*ts_ss;
			ion.uz[i] /= ion.n_ss[i]*ts_ss;
			ion.mean_epsilon[i] /= ion.n_ss[i]*ts_ss;
		}
		electron.fieldT[i] = (2.0/3.0)*(electron.mean_epsilon[i] - 0.5*electron.m/(ech)*(
			 pow(getv(electron.ux[i], electron.uy[i], electron.uz[i]), 2.0)));
		ion.fieldT[i] = (2.0/3.0)*(ion.mean_epsilon[i] - 0.5*ion.m/(ech)*
				pow(getv(ion.ux[i], ion.uy[i], ion.uz[i]), 2.0));

		if(mpi_rank == 0) {
			ElectronDiagFile << (i+0.5)*dx << " " << electron.n_ss[i] << " " << electron.ux[i] << " ";
			ElectronDiagFile << electron.uy[i] << " " << electron.uz[i] << " " << electron.mean_epsilon[i] << " " << electron.fieldT[i] << endl;
			IonDiagFile << (i+0.5)*dx << " " << ion.n_ss[i] << " " << ion.ux[i] << " ";
			IonDiagFile << ion.uy[i] << " " << ion.uz[i] << " " << ion.mean_epsilon[i] << " " << ion.fieldT[i] << endl;
		}
	}
	
	// Clean up
	delete[]elec_range;
	delete[]CS_energy;
	delete[]e_n_CS;
	delete[]e_ex_CS;
	delete[]i_n_CS;
	delete[]rho;
	delete[]weights;
	delete[]RHS;
	delete[]E_field;
	delete[]phi;

	electron.clean();
	ion.clean();
	neutral.clean();
	excited.clean();
	
	MPI_Finalize();
	return 0;
}

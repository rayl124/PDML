#ifndef _COLLISIONMODULES_H
#define _COLLISIONMODULES_H

#include <iostream>
#include <ctime>
#include <math.h>

using namespace std;

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
/* searchIndex
 * Given a set of data and a target value,
 * this will return the index of the value that is
 * the closest but under the target value
 *
 * target: target value
 * data_set: the array to search through
 * data_set_length: number of elements in array
 *
 */
 int searchIndex(double target,
		double *data_set,
		int data_set_length) 
{
	int guess = round((data_set_length - 1)/2.0);
	int a = 0;
	int b = data_set_length - 1;
 
	if (target == data_set[0]) {
		return 0;
	}
	if (target > data_set[data_set_length-1]) {
		return data_set_length-1;
	}
	while(1) {
		if (target < data_set[0]) {
			cout << "ERROR!!! Target is "<< target << endl;
		}
		if (data_set[guess] == target) {
			return guess;
		} else if (data_set[guess - 1] == target) {
			return guess - 1;
		} else if (data_set[guess + 1] == target) {
			return guess + 1;
		}
		if (data_set[guess] < target) {
			if (data_set[guess + 1] > target) {
				return guess;
			} 
			else {
				a = guess;
			}
		}
		else {
			if (data_set[guess - 1] < target) {
				return guess - 1;
			} 
			else {
				b = guess;
			}
		}
		guess = round((b+a)/2.0);
	}
}
/* linInterp
 * Linearly interpolates between two values
 * Returns interpolated y value
 *
 * target_x: value interpolated at
 * x1, x2: left and right x bounds
 * y1, y2: left and right y bounds
 *
 */
double linInterp(double target_x,
		double x1, double x2,
		double y1, double y2) 
{
	return (y2-y1)*(target_x - x1)/(x2-x1) + y1;
}

/* getP
 * Calculates the probability of a collision
 *
 * epsilon: energy of particle [eV]
 * sigma: cross section [m^2]
 * m: mass [kg]
 * n_target: number density [m^-3]
 * dt: timestep [s]
 */
double getP(double g,
			double sigma,
			double n_target,
			double dt) 
{
	double nu = g*sigma*n_target;
	return 1.0-exp(-dt*nu);
}

/* getv
 * Get velocity magnitude from three components
 */
double getv(double vx, double vy, double vz) 
{
	return sqrt(vx*vx + vy*vy + vz*vz);
}

double getSigmaTgMax(double *CS_energy, double *CS_data, double mass, int N_coll, int data_set_length) 
{
	const double e = 1.602e-19;
	double sigmag_max = 0.0;
	double sigma_Tg;
	for (int i = 0; i < data_set_length; ++i) {
		sigma_Tg = 0.0;
		for (int j = 0; j < N_coll; ++j) {
			sigma_Tg += CS_data[j*data_set_length + i]*sqrt(2.0*e*CS_energy[i]/mass);
		}
		if (sigma_Tg > sigmag_max) {
			sigmag_max = sigma_Tg;
		}
	}
	return sigmag_max;
}


double getSigmaTMax(double *CS_data, int N_coll, int data_set_length) 
{
	double sigma_max = 0.0;
	double sigma_T;
	for (int i = 0; i < data_set_length; ++i) {
		sigma_T = 0.0;
		for (int j = 0; j < N_coll; ++j) {
			sigma_T += CS_data[j*data_set_length + i];
		}
		if (sigma_T > sigma_max) {
			sigma_max = sigma_T;
		}
	}
	return sigma_max;
}

/* getNullCollPart
 * Returns number of particles that partake in collisions
 * Also returns max frequency for collision type calculations
 *
 * CS_energy: array of energy vallues of length data_set_length
 * CS_data: array of cross section values for each type of collision
 * 			of length N_coll*data_set_length
 * epsilon_max: max energy of all particles [eV]
 * nu_max: frequency to be returned
 * N_c: number of particles to partake in collisions
 * m_source: mass of source particle [kg]
 * n_target: density of target particle [#/m^3]
 * dt: time step [s]
 * np: number of total particles of source species
 * e: elementary charge [C]
 * N_coll: number of different collision types
 * data_set_length: length of arrays
*/ 	
void getNullCollPart(double sigma_total,
				 double *CS_data,
				 double epsilon_max, double v_target,
				 double *nu_max, double *P_max, int *N_c,
				 double m_source, double n_target, 
				 double dt, double np,
				 int N_coll, int data_set_length) 
{
	//int search_index;
	const double e = 1.602e-19;

	/*
	// Get the closest search index under the target value
	search_index = searchIndex(epsilon_max, CS_energy,
				data_set_length);
	
	// Find sigma_total
	for (int i = 0; i < N_coll; ++i) {
		if (search_index == data_set_length - 1) {
		sigma_total += linInterp(epsilon_max, CS_energy[search_index - 1],
					 CS_energy[search_index],
					 CS_data[i*data_set_length + search_index - 1],
					 CS_data[i*data_set_length + search_index]);
		}
		else {
			sigma_total += linInterp(epsilon_max, CS_energy[search_index],
					 CS_energy[search_index + 1],
					 CS_data[i*data_set_length + search_index],
					 CS_data[i*data_set_length + search_index + 1]);
		}
	}*/


	// Convert energy to joules then calculate velocity
	double v = sqrt(2.0*e*epsilon_max/m_source); 

	// Relative valocity
	double g = v - v_target;
	// Returned values
	*nu_max = g*sigma_total*n_target;

	//*nu_max *= 2.0
	// Probability
	*P_max = 1.0-exp(-dt*(*nu_max));

	// Number of particles in collision
	double n_total = (double) np*(*P_max);
	*N_c = floor(n_total);
	double R = (double) rand()/((double) RAND_MAX);
	if (R < (n_total - *N_c)) {(*N_c)++;}
}


/* Return type of collision from 0 through N_coll, where N_coll
 * is the null collision
 *
 * CS_energy: array of energy vallues of length data_set_length
 * CS_data: array of cross section values for each type of collision
 * 			of length N_coll*data_set_length
 * epsilon: energy of particle [eV]
 * nu_max: max frequency, found in getNullCollPart
 * m_source: mass of source particle [kg]
 * n_target: density of target particle [#/m^3]
 * e: elementary charge [C]
 * N_coll: number of different collision types
 * data_set_length: length of arrays
*/ 	

int getCollType(double *CS_energy, double *CS_data,
				double vx, double vy, double vz,
				double vx_target, double vy_target, double vz_target,	
				double nu_max, double m_source, double n_target,
				int N_coll, int data_set_length)
{
	double *P_vec = new double[N_coll + 1]; // Include P0 = 0 
	double *sigma = new double[N_coll];	
	double *nu = new double[N_coll];
	double sigma_total;
	int search_index;
	const double e = 1.602e-19;
	
	double g = getv(vx-vx_target, vy-vy_target, vz-vz_target);
	double epsilon = 0.5*m_source*g*g/e;

	// Find cross section at the relative velocity
	search_index = searchIndex(epsilon, CS_energy, data_set_length);
	for (int i = 0; i < N_coll; ++i) {
			if(epsilon < 1.0e-7) {sigma[i] = CS_data[i*data_set_length + 1];}
			else { 
				sigma[i] = linInterp(epsilon, CS_energy[search_index],
				 	CS_energy[search_index + 1], 
				 	CS_data[i*data_set_length + search_index],
				 	CS_data[i*data_set_length + search_index + 1]);
			}
	}
		
	//double v = sqrt(2.0*e*epsilon/m_source);

	// P0 = 0, P1 = nu_1/nu_max, p2 ....
	P_vec[0] = 0.0;
	for (int i = 0; i < N_coll; ++i) {
		P_vec[i+1] = 0.0;
		nu[i] = g*sigma[i]*n_target;
		for (int j = 0; j <= i; ++j) {
			P_vec[i+1] += nu[j];
		}
		P_vec[i+1] = P_vec[i+1]/nu_max;
	}
	// Error catching: for rare case of ion.epsilon < 1e-7, renormalize so P=1
	if (P_vec[N_coll] > 1) {
		cout << "getCollType: P > 1 "<< nu_max << endl;
		for(int i = 0; i < N_coll; ++i){
			cout << "nu["<<i<<"]="<<nu[i]<<", sigma["<<i<<"]= "<<sigma[i]<<endl;	
		}
		for(int i = 0; i < N_coll + 1; ++i) {
			P_vec[i] = P_vec[i]*nu_max/P_vec[N_coll];	
		}
	}

	double R = (double)rand()/((double)RAND_MAX);
	// if P_i < R < P_i+1, collision is type i, where tyoe N_coll is null
	search_index = searchIndex(R, P_vec, N_coll + 1);
	
	delete[]P_vec;
	delete[]sigma;
	delete[]nu;

	return search_index;
}


/*	thermalVelSample
 *	Gives a particle a thermal velocity based on
 *	Maxwellian distribution
 *
 *	*part_vx/y/z: Velocities to be returned
 *	T_part: Particle temperature [K]
 *	m_part: Particle mass [kg]
 */
void thermalVelSample(double *part_vx,
			 	double *part_vy,
				double *part_vz,
				double T_part,
				double m_part) 
{
	const double k_B = 1.381e-23;
	// Thermal velocity
	double v_th = sqrt(k_B*T_part/m_part);
	double R1 = (double) rand()/((double)RAND_MAX);
	double R2 = (double) rand()/((double)RAND_MAX);
	double vtmp1 = sqrt(2.0)*v_th*sqrt(-log(R1));
	double vtmp2 = sqrt(2.0)*v_th*sqrt(-log(R2));

	double phi = ((double)(rand()/((double)RAND_MAX)))*2.0*M_PI;
	double theta = ((double)(rand()/((double)RAND_MAX)))*2.0*M_PI;

	// Replace velocity with sampled velocity
	*part_vx = vtmp1*cos(phi);
	*part_vy = vtmp2*sin(theta);
	*part_vz = vtmp2*cos(theta);
}

void thermalBiasVelSample
(double *part_vx, double *part_vy, double *part_vz,
 double T_part, double m_part) 
{
	const double k_B = 1.381e-23;

	double v_th = sqrt(k_B*T_part/m_part);
	double vtmp = sqrt(2.0)*v_th*sqrt(-log(ranf()));
	double theta = ((double)rand()/((double)RAND_MAX))*2.0*M_PI;

	*part_vx = sqrt(2.0)*v_th*sqrt(-log(ranf()));
	*part_vy = vtmp*sin(theta);
	*part_vz = vtmp*cos(theta);
}

void isotropic_e_elastic
(double *part_vx, double *part_vy, double *part_vz)
{
	double v = getv(*part_vx, *part_vy, *part_vz);
	double R1 = (double)rand()/((double)RAND_MAX);
	R1 = -1.0 + 2.0*R1;
	double R3 = (double)rand()/((double)RAND_MAX);

	*part_vx = v*cos(2.0*M_PI*R3)*sqrt(1.0-R1*R1);	
	*part_vx = v*sin(2.0*M_PI*R3)*sqrt(1.0-R1*R1);
	*part_vz = v*R1;
}

/*	e_elastic
 *	Simulates electron-neutral elastic collision
 *	
 *	*partvx/y/z: velocities to be returned
 *	epsilon: particle energy [eV]
 *	m_e: mass of electron [kg]
 *	m_n: mass of neutral [kg]
 */
void e_elastic(double *part_vx,
		double *part_vy,
		double *part_vz,
		double epsilon,
		double m_e,
		double m_n) 
{
	double R1 = (double)rand()/((double)RAND_MAX);
	double R2 = (double)rand()/((double)RAND_MAX);

	// Scattering angles
	// Fix undefines due to rounding
	//double inner = (2.0+epsilon - 2.0*pow(1.0+epsilon, R1))/epsilon;
	//if (inner > 
	double chi = acos((2.0+epsilon-2.0*pow(1.0+epsilon, R1))/epsilon);
	while(isnan(chi)) {chi = acos((2.0+epsilon-2.0*pow(1.0+epsilon,((double)rand()/((double)RAND_MAX))))/epsilon);}

	double phi = R2*2.0*M_PI;

	// Scatter the velocity
	double v_newx;
	double v_newy;
	double v_newz;
	double v = getv(*part_vx, *part_vy, *part_vz);

	v_newx = (*part_vx)*cos(chi) +
		 (*part_vy)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vx)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newy = (*part_vy)*cos(chi) -
		 (*part_vx)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vy)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newz = (*part_vz)*cos(chi) -
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) *
		 sin(chi)*cos(phi);

	// Energy loss correction
	double delta_epsilon = 2.0*m_e*(1.0-cos(chi))*epsilon/m_n;
	double alpha = sqrt(1.0-delta_epsilon/epsilon);
	if (delta_epsilon > epsilon) {cout << "Error in elastic: epislon was " << epsilon << endl;}
	*part_vx = alpha*v_newx;
	*part_vy = alpha*v_newy;
	*part_vz = alpha*v_newz;
	if (isnan(*part_vx)) {
		cout << "Error in elastic" << endl;
		cout << alpha << " " << epsilon << " " << chi << " " << R1 << " " << R2 << endl;
		cout << v_newx << " " << v_newy << " " << v_newz << endl;
	}
}


/*	e_excitation
 *	Electron-neutral excitation collision
 *
 *	*part_vx/y/z: Velocity particle to be returned
 *	epsilon_inc: energy of incident electron [eV]
 *	epsilon_exc: excitation energy [eV]
 *
 */
void e_excitation(double *part_vx,
		double *part_vy,
		double *part_vz,
		double epsilon_inc,
		double epsilon_exc) 
{
	double R1 = (double)rand()/((double)RAND_MAX);
	double R2 = (double)rand()/((double)RAND_MAX);

	double epsilon = epsilon_inc - epsilon_exc;
	// Scattering angles
	double chi = acos((2.0+epsilon-2.0*pow(1.0+epsilon, R1))/epsilon);
	while(isnan(chi)) {chi = acos((2.0+epsilon-2.0*pow(1.0+epsilon,((double)rand()/((double)RAND_MAX))))/epsilon);}
	double phi = R2*2.0*M_PI;

	// Scatter the velocity
	double v_newx;
	double v_newy;
	double v_newz;
	double v = getv(*part_vx, *part_vy, *part_vz);

	v_newx = (*part_vx)*cos(chi) +
		 (*part_vy)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vx)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newy = (*part_vy)*cos(chi) -
		 (*part_vx)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vy)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newz = (*part_vz)*cos(chi) -
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) *
		 sin(chi)*cos(phi);

	// Energy correction factor
	double alpha = sqrt(1.0-epsilon_exc/(epsilon_inc));

	*part_vx = alpha*v_newx;
	*part_vy = alpha*v_newy;
	*part_vz = alpha*v_newz;

	// MAKE SURE TO GIVE THE NEWLY CREATED EXCITED MOLECULE A NEW VELOCITY
	// FROM THE THERMAL VELOCITY AFTER THIS IS CALLED AND THE SAME X
	// COORDINGATE AS THE ELECTRON
}

/*	e_ionization
 *	electron-neutral ionization
 *
 *	*part_vx/y/z: Incident electron velocity to be updated
 *	*part_ej_vx/y/z: Ejected electron velocity to be returned
 *	epsilon_inc: Incident electron energy [eV]
 *	epsilon_ion: Ionization energy [eV]
 */
void e_ionization(double *part_vx,
		double *part_vy,
		double *part_vz,
		double *part_ej_vx,
		double *part_ej_vy,
		double *part_ej_vz,
		double epsilon_inc,
		double epsilon_ion)
{
	double B = 10.0; // for Argon, Vehedi (1994), [eV]
				 // Might make this an input so other species can use this?

	double R1 = ((double)rand())/((double)RAND_MAX);
	double R2 = ((double)rand())/((double)RAND_MAX);
	double R3 = ((double)rand())/((double)RAND_MAX);
	double R4 = ((double)rand())/((double)RAND_MAX);
	double R5 = ((double)rand())/((double)RAND_MAX);

	double epsilon_ej = B*tan(R1*atan((epsilon_inc-epsilon_ion)/(2.0*B)));
	double epsilon_sc = epsilon_inc - epsilon_ion - epsilon_ej;

	double chi_ej = acos((2.0+epsilon_ej-2.0*pow(1.0+epsilon_ej, R2))/epsilon_ej);
	double chi_sc = acos((2.0+epsilon_sc-2.0*pow(1.0+epsilon_sc, R3))/epsilon_sc);

	double phi_ej = R4*2.0*M_PI;
	double phi_sc = R5*2.0*M_PI;

	double v_newx;
	double v_newy;
	double v_newz;
	double v = getv(*part_vx, *part_vy, *part_vz);

	// Ejected electron
	v_newx = (*part_vx)*cos(chi_ej) +
		 (*part_vy)*v*sin(chi_ej)*sin(phi_ej)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vx)*(*part_vz)*sin(chi_ej)*cos(phi_ej)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newy = (*part_vy)*cos(chi_ej) -
		 (*part_vx)*v*sin(chi_ej)*sin(phi_ej)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vy)*(*part_vz)*sin(chi_ej)*cos(phi_ej)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newz = (*part_vz)*cos(chi_ej) -
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) *
		 sin(chi_ej)*cos(phi_ej);

	double alpha = sqrt(epsilon_ej/epsilon_inc);
	*part_ej_vx = alpha*v_newx;
	*part_ej_vy = alpha*v_newy;
	*part_ej_vz = alpha*v_newz;


	// Scattered original electron
	v_newx = (*part_vx)*cos(chi_sc) +
		 (*part_vy)*v*sin(chi_sc)*sin(phi_sc)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vx)*(*part_vz)*sin(chi_sc)*cos(phi_sc)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newy = (*part_vy)*cos(chi_sc) -
		 (*part_vx)*v*sin(chi_sc)*sin(phi_sc)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vy)*(*part_vz)*sin(chi_sc)*cos(phi_sc)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newz = (*part_vz)*cos(chi_sc) -
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) *
		 sin(chi_sc)*cos(phi_sc);

	alpha = sqrt(epsilon_sc/epsilon_inc);
	*part_vx = alpha*v_newx;
	*part_vy = alpha*v_newy;
	*part_vz = alpha*v_newz;

}

/*	i_scattering
 *	Ion-neutral scattering
 *
 *	*part_vx/y/z: Particle velocity to be returned
 *	epsilon_inc: Incident particle energy [eV]
 *	m_source: Incident particle mass [kg]
 *	m_target: Target particle mass [kg]
 *	T_target: Temperature of the target [K] *
 */
void i_scattering(double *part_vx, double *part_vy,	double *part_vz,
		double part_t_vx, double part_t_vy, double part_t_vz,
		double epsilon_inc,	double m_source, double m_target,
		double T_target)
{
	/*
	double part_t_vx = 0.0;
	double part_t_vy = 0.0;
	double part_t_vz = 0.0;

	thermalVelSample(&part_t_vx, &part_t_vy, &part_t_vz,
			 T_target, m_target);
	*/
	
	// Change to reference frame where target is at rest
	*part_vx -= part_t_vx;
	*part_vy -= part_t_vy;
	*part_vz -= part_t_vz;
	

	// Assumes that source and target are of the same mass
	// such as Ar and Ar+
	// DO NOT USE for scattering of different masses
	double R = (double) rand()/((double)RAND_MAX);

	// Neutral frame scattering angle
	double chi = acos(sqrt(1.0-R)); 
	
	double alpha_L = sin(chi)*sin(chi);
	double phi = 2.0*M_PI*((double)rand()/((double)RAND_MAX));

	
	double v_newx;
	double v_newy;
	double v_newz;
	double v = getv(*part_vx, *part_vy, *part_vz);

	// Scatter in neutral frame
	v_newx = (*part_vx)*cos(chi) +
		 (*part_vy)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vx)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newy = (*part_vy)*cos(chi) -
		 (*part_vx)*v*sin(chi)*sin(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) +	
		 (*part_vy)*(*part_vz)*sin(chi)*cos(phi)/
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0));

	v_newz = (*part_vz)*cos(chi) -
		 sqrt(pow(*part_vx,2.0) + pow(*part_vy,2.0)) *
		 sin(chi)*cos(phi);
	
	// Change of energy in neutral frame
	double delta_epsilon = alpha_L*epsilon_inc;
	double alpha = sqrt(1.0 - delta_epsilon/epsilon_inc);

	// Apply energy correction in neutral frame
	*part_vx = alpha*v_newx;
	*part_vy = alpha*v_newy;
	*part_vz = alpha*v_newz;

	// Change back to lab frame
	*part_vx += part_t_vx;
	*part_vy += part_t_vy;
	*part_vz += part_t_vz;
}
#endif

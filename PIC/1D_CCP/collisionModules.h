#ifndef _COLLISIONMODULES_H
#define _COLLISIONMODULES_H

#include <iostream>
#include <ctime>
#include <math.h>

using namespace std;

/* Given a set of data and a target value,
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
		int data_set_length) {
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
/* Linearly interpolates between two values
 * Returns interpolated y value
 * target_x: value interpolated at
 * x1, x2: left and right x bounds
 * y1, y2: left and right y bounds
 *
 */
double linInterp(double target_x,
		double x1, double x2,
		double y1, double y2) {
  return (y2-y1)*(target_x - x1)/(x2-x1) + y1;
}


double getP(double epsilon,
	    double sigma,
	    double m,
	    double n_target,
	    double dt) {
  double e = 1.602e-19;
  double v = sqrt(2*e*epsilon/m);
  double nu = v*sigma*n_target;
  return 1.0-exp(-dt*nu);
}

double getv(double vx, double vy, double vz) {
  return sqrt(vx*vx + vy*vy + vz*vz);
}

/* Returns number of particles that partake in collisions
 * Also returns max frequency for collision type calculations
 *
 * CS_energy: array of energy vallues of length data_set_length
 * CS_data: array of cross section values for each type of collision
 * 	    of length N_coll*data_set_length
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
void getNullCollPart(double *CS_energy,
		     double *CS_data,
		     double epsilon_max,
		     double *nu_max, double *P_max, int *N_c,
		     double m_source, double n_target, 
		     double dt, double np,
		     int N_coll, int data_set_length) {

  int search_index;
  double sigma_total = 0.0;
  const double e = 1.602e-19;

  search_index = searchIndex(epsilon_max, CS_energy,
		    data_set_length);
  
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
  }

  // Convert energy to joules then calculate velocity
  double v = sqrt(2*e*epsilon_max/m_source); 
  *P_max = getP(epsilon_max, sigma_total, m_source,
			n_target, dt);
  // Returned values
  *nu_max = v*sigma_total*n_target;

  // Add padding to P_max and nu_max to ensure you get the max
  // probability
  *P_max *= 2.0;
  *nu_max *= 2.0;

  *N_c = round(np*(*P_max));
}


/* Return type of collision from 0 through N_coll, where N_coll
 * is the null collision
 *
 * CS_energy: array of energy vallues of length data_set_length
 * CS_data: array of cross section values for each type of collision
 * 	    of length N_coll*data_set_length
 * epsilon: energy of particle [eV]
 * nu_max: max frequency, found in getNullCollPart
 * m_source: mass of source particle [kg]
 * n_target: density of target particle [#/m^3]
 * e: elementary charge [C]
 * N_coll: number of different collision types
 * data_set_length: length of arrays
*/ 	

int getCollType(double *CS_energy,
		   double *CS_data,
		   double epsilon,
		   double nu_max,
		   double m_source, double n_target,
		   int N_coll, int data_set_length) {

  double *P_vec = new double[N_coll + 1]; // Include P0 = 0
  double *sigma = new double[N_coll];  
  double *nu = new double[N_coll];
  double sigma_total;
  int search_index;
  const double e = 1.602e-19;
  
  search_index = searchIndex(epsilon, CS_energy, data_set_length);
  for (int i = 0; i < N_coll; ++i) {
       sigma[i] = linInterp(epsilon, CS_energy[search_index],
	       CS_energy[search_index + 1], 
	       CS_data[i*data_set_length + search_index],
	       CS_data[i*data_set_length + search_index + 1]);
  }
    
  double v = sqrt(2*e*epsilon/m_source);
  // P0 = 0, P1 = nu_1/nu_max, p2 ....
  P_vec[0] = 0.0;
  for (int i = 0; i < N_coll; ++i) {
    P_vec[i+1] = 0.0;
    nu[i] = v*sigma[i]*n_target;
    for (int j = 0; j <= i; ++j) {
      P_vec[i+1] += nu[j];
    }
    P_vec[i+1] = P_vec[i+1]/nu_max;
  }

  double R = double(rand())/RAND_MAX;
  // if P_i < R < P_i+1, collision is type i, where tyoe N_coll is null
  search_index = searchIndex(R, P_vec, N_coll + 1);

  delete(P_vec);
  delete(sigma);
  delete(nu);

  return search_index;
}


// Take a maxwellian sample of the thermal velocity
// Also used to simulate ion-neutral collision
void thermalVelSample(double *part_vx,
			 	double *part_vy,
				double *part_vz,
				double T_part,
				double m_part) {

  const double k_B = 1.381e-23;
  // Thermal velocity
  double v_th = sqrt(2*k_B*T_part/m_part);
  
  int M = 3;
  double f_M = 0.0;

  
  // Maxwellian frequency
  for (int j = 1; j <= M; ++j) {
  f_M += double(rand())/RAND_MAX;
  }
  f_M = (f_M - M/2.0)*(1.0/(sqrt(M/12.0)));

  double phi = (double(rand())/RAND_MAX)*2*M_PI;
  double theta = (double(rand())/RAND_MAX)*2*M_PI;

  // Replace velocity with sampled velocity
  *part_vx = f_M*v_th*cos(phi)*cos(theta);
  *part_vy = f_M*v_th*cos(phi)*sin(theta);
  *part_vz = f_M*v_th*sin(phi);
}

void e_elastic(double *part_vx,
		double *part_vy,
		double *part_vz,
		double epsilon,
		double m_e,
		double m_n) {

  double R1 = double(rand())/RAND_MAX;
  double R2 = double(rand())/RAND_MAX;

  double chi = acos((2+epsilon-2*pow(1+epsilon, R1))/epsilon);
  double phi = R2*2*M_PI;

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

  double delta_epsilon = 2*m_e*(1-cos(chi))*epsilon/m_n;
  double alpha = sqrt(1-delta_epsilon/epsilon);

  *part_vx = alpha*v_newx;
  *part_vy = alpha*v_newy;
  *part_vz = alpha*v_newz;
}


// Excitation, keep epsilon_exc as positive,
// Deexcitation, make epsilon_exc negative
void e_excitation(double *part_vx,
		double *part_vy,
		double *part_vz,
		double epsilon_inc,
		double epsilon_exc) 
{
  double R1 = double(rand())/RAND_MAX;
  double R2 = double(rand())/RAND_MAX;

  double epsilon = epsilon_inc - epsilon_exc;
  double chi = acos((2+epsilon-2*pow(1+epsilon, R1))/epsilon);
  double phi = R2*2*M_PI;

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

  double alpha = sqrt(1-epsilon_exc/(epsilon_inc));

  *part_vx = alpha*v_newx;
  *part_vy = alpha*v_newy;
  *part_vz = alpha*v_newz;

  // MAKE SURE TO GIVE THE NEWLY CREATED EXCITED MOLECULE A NEW VELOCITY
  // FROM THE THERMAL VELOCITY AFTER THIS IS CALLED AND THE SAME X
  // COORDINGATE AS THE ELECTRON
}

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

  double R1 = double(rand())/RAND_MAX;
  double R2 = double(rand())/RAND_MAX;
  double R3 = double(rand())/RAND_MAX;
  double R4 = double(rand())/RAND_MAX;
  double R5 = double(rand())/RAND_MAX;

  double epsilon_ej = B*tan(R1*atan((epsilon_inc-epsilon_ion)/(2*B)));
  double epsilon_sc = epsilon_inc - epsilon_ion - epsilon_ej;

  double chi_ej = acos((2+epsilon_ej-2*pow(1+epsilon_ej, R2))/epsilon_ej);
  double chi_sc = acos((2+epsilon_sc-2*pow(1+epsilon_sc, R3))/epsilon_sc);

  double phi_ej = R4*2*M_PI;
  double phi_sc = R5*2*M_PI;

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
  while (isnan(alpha)) {
    cout << "alpha 1 is wrong" << endl;
  }
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
  while (isnan(alpha)) {
    cout << "alpha 2 is wrong" << endl;
  }

  *part_vx = alpha*v_newx;
  *part_vy = alpha*v_newy;
  *part_vz = alpha*v_newz;

}

void i_scattering(double *part_vx,
		double *part_vy,
		double *part_vz,
		double epsilon_inc,
		double m_source,
		double m_target,
		double T_target)
{

  double part_t_vx = 0.0;
  double part_t_vy = 0.0;
  double part_t_vz = 0.0;

  thermalVelSample(&part_t_vx, &part_t_vy, &part_t_vz,
		   T_target, m_target);

  // Change to reference frame where target is at rest
  *part_vx -= part_t_vx;
  *part_vy -= part_t_vy;
  *part_vz -= part_t_vz;

  // Assumes that source and target are of the same mass
  // such as Ar and Ar+
  // DO NOT USE for scattering of different masses
    double chi = acos(sqrt(1.0-(double(rand())/RAND_MAX)));
  /* Add this when you figure out what chi is in the case
   * that m_target =/= m_source
  double alpha = 2*m_target*m_source*(1.0-cos(theta))/
	  pow(m_source+m_target,2.0);
*/
  double alpha = sin(chi)*sin(chi);
  double phi = 2.0*M_PI*(double (rand())/RAND_MAX);


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

  *part_vx = alpha*v_newx;
  *part_vy = alpha*v_newy;
  *part_vz = alpha*v_newz;

  // Change back to lab frame
  *part_vx += part_t_vx;
  *part_vy += part_t_vy;
  *part_vz += part_t_vz;
}
#endif

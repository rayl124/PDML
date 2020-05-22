#ifndef _COLLISIONMODULES_H
#define _COLLISIONMODULES_H

#include <iostream>
#include <ctime>
#include <math.h>

double getP(double epsilon,
	    double sigma,
	    double m,
	    double n_target,
	    double dt) {
  double v = sqrt(2*epsilon/m);
  double nu = v*sigma*n_target;
  return 1.0-exp(-dt*nu);
}

int searchIndex(double target,
		double *data_set,
		int data_set_length) {
  int guess = round((data_set_length - 1)/2.0);
  int a = 0;
  int b = data_set_length - 1;

  while(1) {
    if (data_set[guess] == target) {
      return guess;
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
        return guess;
      } 
      else {
        b = guess;
      }
    }
    guess = round((b+a)/2.0);
  }
  
}

double linInterp(double target_x,
		double x1, double x2,
		double y1, double y2) {
  return (y2-y1)*(target_x - x1)/(x2-x1) + y1;
}

int getCollType(double *CS_energy,
		   double *CS_data,
		   double epsilon,
		   double epsilon_max,
		   double n_target,
		   double m_source,
		   double np, double dt,
		   int N_coll, int data_set_length) {

  double *P_vec = new double[N_coll + 1]; // Include P0 = 0
  double *sigma = new double[N_coll];

  double v_max = sqrt(2*epsilon_max/m_source);
  
  double sigma_total;
  int search_index;

  for (int i = 0; i < N_coll; ++i) {
    search_index = searchIndex(epsilon, CS_data, data_set_length);
    sigma[y] = linInterp(epsilon, CS_energy[search_index],
	       CS_energy[search_index + 1], 
	       CS_data[i*data_set_length + search_index],
	       CS_data[i*data_set_length + search_index + 1]);
    sigma_total += sigma[i];
  }
  
  double nu_max = v_max*sigma_total*n_target;
  
  // P0 = 0, P1 = nu_1/nu_max, p2 ....
  P_vec[0] = 0.0;
  for (int i = 0; i < N_coll; ++i) {
    nu[i] = v*sigma[i]*n_target;
    for (int j = 0; j <= i; ++j) {
      P_vec[i+1] += nu[i];
    }
    P_vec[i] = P_vec[i]/nu_max;
  }

  double R = double(rand())/RAND_MAX;
  // if P_i < R < P_i+1, collision is type i+1
  search_index = searchIndex(R, P_vec, N_coll + 1);

  delete(P_vec);
  delete(sigma);

  return search_index;
}


// Models isotropic charge exchange OR
// is used to give newly created excited atoms
// a thermal velocity
void thermalVelSample(double *part_vx,
			 	double *part_vy,
				double *part_vz,
				double T_n,
				double k_B,
				double m_n) {
  // Neutral thermal velocity
  double v_nth = sqrt(2*k_B*T_n/m_n);

  int M = 3;
  double f_M = new double[3];


  // Maxwellian frequency
  for (int i = 0; i < 4; ++i) {
    for (int j = 1; j <= M; ++j) {
      f_M[i] += double(rand())/RAND_MAX;
    }
    f_M[i] = (f_M - M/2.0)*(1/(sqrt(M/12)));
  }

  // Replace ion velocity with neutral velocity
  *part_vx = f_M[0]*v_nth;
  *part_vy = f_M[1]*v_nth;
  *part_vz = f_M[2]*v_nth;

  delete(f_M);
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
  double v = sqrt(pow(*part_vx,2.0) +
		  pow(*part_vy,2.0) +
		  pow(*part_vz,2.0));

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
		double epsilon_exc,
		double m_e,
		double m_n) {

  double R1 = double(rand())/RAND_MAX;
  double R2 = double(rand())/RAND_MAX;

  double epsilon = epsilon_inc - epsilon_exc;
  double chi = acos((2+epsilon-2*pow(1+epsilon, R1))/epsilon);
  double phi = R2*2*M_PI;

  double v_newx;
  double v_newy;
  double v_newz;
  double v = sqrt(pow(*part_vx,2.0) +
		  pow(*part_vy,2.0) +
		  pow(*part_vz,2.0));

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
		double *part_ej,vz,
		double epsilon_inc,
		double epsilon_ion,
		double m_e,
		double m_n) {
  double B = 10.0; // for Argon, Vehedi (1994), [eV]
  		   // Might make this an input so other species can use this?

  double R1 = double(rand())/RAND_MAX;
  double R2 = double(rand())/RAND_MAX;
  double R3 = double(rand())/RAND_MAX;
  double R4 = double(rand())/RAND_MAX;
  double R5 = double(rand())/RAND_MAX;

  double epsilon_ej = B*tan(R1*atan((epsilon_inc-epsilon_ion)/2*B));
  double epsilon_sc = epsilon_inc - epsilon_ion - epsilon_ej;

  double chi_ej = acos((2+epsilon_ej-2*pow(1+epsilon_ej, R2))/epsilon_ej);
  double chi_sc = acos((2+epsilon_sc-2*pow(1+epsilon_sc, R3))/epsilon_sc);

  double phi_ej = R4*2*M_PI;
  double phi_sc = R5*2*M_PI;

  double v_newx;
  double v_newy;
  double v_newz;
  double v = sqrt(pow(*part_vx,2.0) +
		  pow(*part_vy,2.0) +
		  pow(*part_vz,2.0));

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
#endif

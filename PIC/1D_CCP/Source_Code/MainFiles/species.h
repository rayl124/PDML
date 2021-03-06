#include <iostream>
#include "collisionModules.h"

class particles {
  int max_part;  // Max particles

  public:
    double spwt; // Specific weight (can be changed to be variable)
    double *x;   // Position in m
    double *vx;  // x velocity in m/s
    double *vy;  // y velocity in m/s
    double *vz;  // z velocity in m/s
    double *epsilon;  // energy in eV
    double *gamma; // scattering coefficients,
    		   // 0 is reflection, 1 is secondary electron emission
    double *E;
	double *v_max; // Max possible velocity
    // Keeps track of global nodes the particle is between
    int *node_index;
    // Keeps track of cell center the particle is between;
    int *cell_index;

    // Field arrays
    double *n;   	   // Number density in m^-3, local
    double *n_ss;	   // Number density for steady state to be averaged at the end
    double *mean_epsilon;  // Mean energy
    double *ux;		   // Bulk vx
    double *uy;		   // Bulk vy
    double *uz;		   // Bulk vz
    double *fieldT;
    int *fieldnp;    

    double T;  // Temperature in K
    double m;  // Mass in kg
    double q;  // Charge in C
    double flux_L = 0.0;
    double flux_R = 0.0;
    double n_bar = 0.0; // Average density

    double max_epsilon = 0.0;  // Max energy
    int np = 0;  // Number of current particles
    int inner_np = 0; // Particles in the middle 50% of the domain

    // Global total of particles
    int np_total;
    int inner_np_total;


    void initialize(int max_part, int n_cell);
    void clean(void);
    void remove_part(int index);
    void thermalVelocity(int index);
    void injectionVelocity(int index);
    void getDiagnosticsLocal(double *weights, int p, double dx);
};

// Initializes arrays, required every time a new
// species is created
void particles::initialize(int max_part, int n_cell) {
  //spwt = new double[max_part];
  x = new double[max_part];
  vx = new double[max_part];
  vy = new double [max_part];
  vz = new double [max_part];
  epsilon = new double[max_part];
  // Keeps track of node closest to
  node_index = new int[max_part];
  // Keeps track of cell center closest to
  cell_index = new int[max_part];
  E = new double[max_part];


  n = new double[n_cell];
  n_ss = new double[n_cell];
  mean_epsilon = new double[n_cell];
  ux = new double[n_cell];
  uy = new double[n_cell];
  uz = new double[n_cell];
  fieldT = new double[n_cell];
  fieldnp = new int[n_cell];

  for (int i = 0; i < n_cell; ++i) {
    n_ss[i] = 0.0;
    mean_epsilon[i] = 0.0;
    ux[i] = 0.0;
    uy[i] = 0.0;
    uz[i] = 0.0;
    fieldT[i] = 0.0;
    fieldnp[i] = 0.0;
  }
  gamma = new double[2]; 
}


// Frees memory when simulation is done
void particles::clean(void) {
  delete[]x;
  delete[]vx;
  delete[]vy;
  delete[]vz;
  delete[]epsilon;
  delete[]node_index;
  delete[]cell_index;
  delete[]n;
  delete[]n_ss;
  delete[]mean_epsilon;
  delete[]ux;
  delete[]uy;
  delete[]uz;
  delete[]fieldT;
  delete[]gamma;
  delete[]E;
}

// Removed a particle when it leaves the domain
void particles::remove_part(int index) {
  //spwt[index] = spwt[np-1];
  x[index] = x[np-1];
  vx[index] = vx[np-1];
  vy[index] = vy[np-1];
  vz[index] = vz[np-1];
  epsilon[index] = epsilon[np-1];
  node_index[index] = node_index[np-1];
  cell_index[index] = cell_index[np-1];
  E[index] = E[np-1];
  np -= 1;
}

// Gives a particle a velocity from its thermal velocity
void particles::thermalVelocity(int index) {
  thermalVelSample(&vx[index], &vy[index], &vz[index], T, m);
}

void particles::injectionVelocity(int index) {
  double injectionT = 3.0; // eV
  vx[index] = sqrt(2.0*injectionT*1.602e-19/m);
  vy[index] = 0.0;
  vz[index] = 0.0;
  //convert injectionT to T: T*ech/k_B
  //injectionT = 1.602e-19*injectionT/(1.381e-23);
  //thermalBiasVelSample(&vx[index], &vy[index], &vz[index], injectionT, m);
}

void particles::getDiagnosticsLocal(double *weights, int index, double dx)
{
  n_ss[cell_index[index]] += spwt*weights[0]/dx;
  n_ss[cell_index[index]+1] += spwt*weights[1]/dx;

  ux[cell_index[index]] += spwt*weights[0]*vx[index]/dx;
  ux[cell_index[index] + 1] += spwt*weights[1]*vx[index]/dx;
  uy[cell_index[index]] += spwt*weights[0]*vy[index]/dx;
  uy[cell_index[index] + 1] += spwt*weights[1]*vy[index]/dx;
  uz[cell_index[index]] += spwt*weights[0]*vz[index]/dx;
  uz[cell_index[index] + 1] += spwt*weights[1]*vz[index]/dx;

  mean_epsilon[cell_index[index]] += spwt*weights[0]*epsilon[index]/dx;
  mean_epsilon[cell_index[index] + 1] += spwt*weights[1]*epsilon[index]/dx; 

}

class fluid {
  public:
    double *n;
    double *n_ss;
    double *n_dot;
    double *f_coeff;
    double *T;
    double *v;
    double *beta;
    double *D;

    double m;
    double flux_L;
    double flux_R;
    double n_bar = 0.0;

    void initialize(int n_cell);
    void clean(void);
    double getMomentumCS(int index);
    double getMax_n(int n_cell);
};

void fluid::initialize(int n_cell) {
  n = new double[n_cell];
  n_ss = new double[n_cell];
  n_dot = new double[n_cell];
  f_coeff = new double[n_cell];
  T = new double[n_cell];
  v = new double[n_cell];
  beta = new double[n_cell];
  D = new double[n_cell];
  flux_L = 0.0;
  flux_R = 0.0;

  for (int i = 0; i < n_cell; ++i) {
    n_ss[i] = 0.0;
	n_dot[i] = 0.0;
  }
}

void fluid::clean() {
  delete[]n;
  delete[]n_dot;
  delete[]f_coeff;
  delete[]T;
  delete[]v;
  delete[]beta;
  delete[]D;
  delete[]n_ss;
}

double fluid::getMomentumCS(int index) {
  double *c_HEKhrapak = new double[4]; // Coefficients from Khrapak(2014)
  double *c_LEKhrapak = new double[4];

  c_HEKhrapak[0] = -0.692;
  c_HEKhrapak[1] = 9.594;
  c_HEKhrapak[2] = -8.284;
  c_HEKhrapak[3] = -2.355;
  c_LEKhrapak[0] = -0.019;
  c_LEKhrapak[1] = 0.038;
  c_LEKhrapak[2] = -0.049;
  c_LEKhrapak[3] = 0.015;

  double f = 1.0;

  double sigma;

  double beta_sample = beta[index];

  if (beta_sample < 0.506) {
    for (int i = 0; i < 4; ++i) {
      f += c_HEKhrapak[i]*pow(beta_sample,(i+1));
    }
    sigma = 4.507*pow(beta_sample,1.0/6.0)*f;
  } else {
    for (int i = 0; i < 4; ++i) {
      f += c_LEKhrapak[i]*pow(beta_sample,-(i+1));
    }
    sigma = 9.866*pow(beta_sample,1.0/3.0)*f;
  }

  return sigma;
}

double fluid::getMax_n(int n_cell) {
  double max = 0.0;

  for (int i = 0; i < n_cell; ++i) {
    if (n[i] > max) {
      max = n[i];
    }
  }

  return max;
}

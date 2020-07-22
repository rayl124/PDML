#include <iostream>
#include "collisionModules_old.h"

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

    // Keeps track of global nodes the particle is between
    int *node_index;
    // Keeps track of cell center the particle is between;
    int *cell_index;

    // Field arrays
    double *n;   // Number density in m^-3, local
    double *n_master; // Global number density
    double *epsilon_bulk;
    double *vx_bulk;

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
  n_master = new double[n_cell];
  epsilon_bulk = new double[n_cell];
  vx_bulk = new double[n_cell];
  gamma = new double[2]; 
}

// Frees memory when simulation is done
void particles::clean(void) {
  delete(x);
  delete(vx);
  delete(vy);
  delete(vz);
  delete(epsilon);
  delete(node_index);
  delete(cell_index);
  delete(n);
  delete(n_master);
  delete(epsilon_bulk);
  delete(vx_bulk);
  delete(gamma);
  delete(E);
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
  np -= 1;
}

// Gives a particle a velocity from its thermal velocity
void particles::thermalVelocity(int index) {
  thermalVelSample(&vx[index], &vy[index], &vz[index],
		  T, m);
}

void particles::injectionVelocity(int index) {
  thermalBiasVelSample(&vx[index], &vy[index], &vz[index],
		  T, m);
}

class fluid {
  public:
    double *n;
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
  n_dot = new double[n_cell];
  f_coeff = new double[n_cell];
  T = new double[n_cell];
  v = new double[n_cell];
  beta = new double[n_cell];
  D = new double[n_cell];
  flux_L = 0.0;
  flux_R = 0.0;
}

void fluid::clean() {
  delete(n);
  delete(n_dot);
  delete(f_coeff);
  delete(T);
  delete(v);
  delete(beta);
  delete(D);
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

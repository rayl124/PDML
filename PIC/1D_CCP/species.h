#include <iostream>
#include "collisionModules.h"

class species {
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
    double *n;   // Number density in m^-3

    // Keeps track of global nodes the particle is between
    int *node_index;
    // Keeps track of cell center the particle is between;
    int *cell_index;


    double T;  // Temperature in K
    double m;  // Mass in kg
    double q;  // Charge in C

    double max_epsilon = 0.0;  // Max energy
    int np = 0;  // Number of current particles

    void initialize(int max_part, int n_cell);
    void clean(void);
    void remove_part(int index);
    void thermalVelocity(int index);
};

// Initializes arrays, required every time a new
// species is created
void species::initialize(int max_part, int n_cell) {
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

  n = new double[n_cell];
  gamma = new double[2]; 
}

// Frees memory when simulation is done
void species::clean(void) {
  delete(x);
  delete(vx);
  delete(vy);
  delete(vz);
  delete(epsilon);
  delete(node_index);
  delete(cell_index);
  delete(n);
  delete(gamma);
}

// Removed a particle when it leaves the domain
void species::remove_part(int index) {
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
void species::thermalVelocity(int index) {
  thermalVelSample(&vx[index], &vy[index], &vz[index],
		  T, m);
}

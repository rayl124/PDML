#include <iostream>
#include "collisionModules.h"

class species {
  int max_part;

  public:
    double spwt;
    double *x;
    double *vx;
    double *vy;    
    double *vz;   
    double *epsilon;
    double *gamma; // 0 is reflection, 1 is see    
    // Keeps track of node closest to
    int *node_index;

    double T;
    double m;
    double q;
    double max_epsilon = 0.0;
    int np = 0;

    void initialize(int max_part);
    void clean(void);
    void remove_part(int index);
    void thermalVelocity(int index);
};

void species::initialize(int max_part) {
  //spwt = new double[max_part];
  x = new double[max_part];
  vx = new double[max_part];
  vy = new double [max_part];
  vz = new double [max_part];
  epsilon = new double[max_part];
  // Keeps track of node closest to
  node_index = new int[max_part];
  gamma = new double[2]; 
}

void species::clean(void) {
  //delete(spwt);
  delete(x);
  delete(vx);
  delete(vy);
  delete(vz);
  delete(epsilon);
  delete(gamma);
}

void species::remove_part(int index) {
  //spwt[index] = spwt[np-1];
  x[index] = x[np-1];
  vx[index] = vx[np-1];
  vy[index] = vy[np-1];
  vz[index] = vz[np-1];
  epsilon[index] = epsilon[np-1];
  np -= 1;
}


void species::thermalVelocity(int index) {
  thermalVelSample(&vx[index], &vy[index], &vz[index],
		  T, m);
}

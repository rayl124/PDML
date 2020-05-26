#include <iostream>

class species {
  int max_part;

  public:
    species (int);
    double *spwt = new double[max_part];
    double *x = new double[max_part];
    double *vx = new double[max_part];
    double *vy = new double [max_part];
    double *vz = new double [max_part];
    double *epsilon = new double[max_part];
    // Keeps track of node closest to
    int *node_index = new int[max_part]; 
    double m;
    double q;
    double max_epsilon = 0.0;
    int np = 0;

    void clean(void);
    void remove_part(int index);
};

species::species(int a) {
  max_part = a;
}

void species::clean(void) {
  delete(spwt);
  delete(x);
  delete(vx);
  delete(vy);
  delete(vz);
  delete(epsilon);
}

void species::remove_part(int index) {
  spwt[index] = spwt[np-1];
  x[index] = x[np-1];
  vx[index] = vx[np-1];
  vy[index] = vy[np-1];
  vz[index] = vz[np-1];
  epsilon[index] = epsilon[np-1];
  np -= 1;
}

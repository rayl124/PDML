#include "collisionModules.h"


int main(void) {
  int N_c = 3754;
  double max_epsilon = 266.265;
  double nu_max = 6.58146e7;

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


}

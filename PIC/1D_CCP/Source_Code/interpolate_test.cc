#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

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


int main(void) {

  int *N_coll = new int[4]; // Types of collisions for each reaction
  N_coll[0] = 4;
  N_coll[1] = 2;
  N_coll[2] = 2;
  N_coll[3] = 1;

  int data_length = 5996;
  double *CS_energy = new double[5996]; // energy values
  double *el_n_CS = new double[N_coll[0]*5996]; // electron-neutral
  double *el_ex_CS = new double[N_coll[1]*5996]; // electron-excited
  double *i_n_CS = new double [N_coll[2]*5996]; // ion-neutral

  ifstream coll_data("crtrs.dat.txt");
  // Ignore the first 2 lines
  coll_data.ignore(1e5, '\n');
  coll_data.ignore(1e5, '\n');
  
  for (int i = 0; i < 5996; ++i) {
    coll_data >> CS_energy[i];
    for (int j = 0; j < N_coll[0]; ++j) {
      coll_data >> el_n_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[1]; ++j) {
      coll_data >> el_ex_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[2]; ++j) {
      coll_data >> i_n_CS[j*5996 + i];
    }
  }
  
  coll_data.close();

  

  double test_epsilon = 3.925e-2;

  int index;
  index  = searchIndex(test_epsilon, CS_energy, data_length);

  double interp_value;
  interp_value = ((el_n_CS[index+1]-el_n_CS[index])/
		 (CS_energy[index+1]-CS_energy[index]))*
	  	 (test_epsilon - CS_energy[index]) + el_n_CS[index];

  cout <<"Test epsilon:  " << test_epsilon << " , Cross Section Value: ";
  cout << interp_value << endl;

}


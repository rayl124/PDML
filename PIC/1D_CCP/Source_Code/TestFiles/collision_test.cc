#include "collisionModules.h"
#include <fstream>

using namespace std;

int main(void) {
  srand(time(NULL));
  // Load data set

  int data_set_length = 5996;
  int *N_coll = new int[4];
  N_coll[0] = 4;
  N_coll[1] = 2;
  N_coll[2] = 2;
  N_coll[3] = 1;

  double *CS_energy = new double[5996]; // energy values
  double *e_n_CS = new double[N_coll[0]*5996]; // electron-neutral
  double *e_ex_CS = new double[N_coll[1]*5996]; // electron-excited
  double *i_n_CS = new double [N_coll[2]*5996]; // ion-neutral
  double *penning_CS = new double[N_coll[3]*5996]; // penning_excitation
  ifstream coll_data("crtrs.dat.txt");
  // Ignore the first 2 lines
  coll_data.ignore(1e5, '\n');
  coll_data.ignore(1e5, '\n');
  
  for (int i = 0; i < 5996; ++i) {
    coll_data >> CS_energy[i];
    for (int j = 0; j < N_coll[0]; ++j) {
      coll_data >> e_n_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[1]; ++j) {
      coll_data >> e_ex_CS[j*5996 + i];
    }
    for (int j = 0; j < N_coll[2]; ++j) {
      coll_data >> i_n_CS[j*5996 + i];
    }
  }

  coll_data.close();  
  
  // Load test particles
   // Constants
  const double epsilon0 = 8.854e-12; // Permittivity of free space
  const double e = 1.602e-19; // Elementary charge [C]
  const double k_B = 1.381e-23; // Boltzman constant [J/K]
  const double AMU = 1.661e-27; // Atomic Mass Unit [kg]
  const double m_n = 39.948*AMU; // Ion mass of Ar [kg]
  const double m_i = m_n;
  const double m_e = 9.109e-31; // Electron mass [kg]

  // Input settings
  double n0 = 1.0e12; // Electron density [#/m^3]
  double phi0 = 0.0; // Reference potential
  double T_n = 300; // Neutral temp, [K]
  double T_i = T_n;
  double T_e = 30000; // Electron temp

  // Problem discretization
  int nn = 191; // # of x1 nodes
  int ts = 4e6; // # of time steps
  double dx = 5.0e-4; // length of each cell
  double dt = 1.0/(27.12e6*200.0); 
  double L = 0.095; // Domain length m

  // Particle info
  int max_part = 1e6; // max_particles if none leave during time history
  
  // Electron particle data
  double *part_e_x = new double[max_part]; // Position
  double *part_e_vx = new double[max_part]; // Velocity
  double *part_e_vy = new double[max_part];
  double *part_e_vz = new double[max_part];
  // To convert J to eV divide by e (elementary charge)
  double *part_e_epsilon = new double[max_part]; // Energy [J]

  /*
  // Ion particle data
  double *part_i_spwt = new double[max_part]; // Particle weight
  double *part_i_x = new double[max_part]; // Position
  double *part_i_vx = new double[max_part]; // Velocity
  double *part_i_vy = new double[max_part];
  double *part_i_vz = new double[max_part];
  // To convert J to eV divide by e (elementary charge)
  double *part_i_epsilon = new double[max_part]; // Energy [J]
  */

  double n_e = 1e13;
  double n_n = 1e19;

  double part_q;
  double max_epsilon = 0.0;

  double vth_i = sqrt(2*k_B*T_i/m_i);
  double vth_e = sqrt(2*k_B*T_e/m_e);
  
  double lambda_d = sqrt(epsilon0*k_B*T_e/(e*e*n_e));
  int np_e = round(lambda_d*lambda_d*lambda_d*n_e);
  
  int max_index;
  cout << "np_e = " << np_e << endl;


  for (int i = 0; i < np_e; ++i) {
    thermalVelSample(&part_e_vx[i], &part_e_vy[i], &part_e_vz[i],
		    T_e, m_e);
    //part_e_vx[i] = part_e_vx[i]*6.0*(double (rand())/RAND_MAX);
    //part_e_vy[i] = part_e_vy[i]*6.0*(double (rand())/RAND_MAX);
    //part_e_vz[i] = part_e_vz[i]*6.0*(double (rand())/RAND_MAX);


    part_e_epsilon[i] = 0.5*m_e*pow(getv(part_e_vx[i], 
			part_e_vy[i], part_e_vz[i]),2.0)/e; //[eV]
    if (part_e_epsilon[i] > max_epsilon) {
      max_epsilon = part_e_epsilon[i];
      max_index = i;
    }
  }
  cout << "max_epsilon = " << max_epsilon << endl;
  double nu_max, P_max, N_c;
  // Get number of particles for electron - neutral collisions
  getNullCollPart(CS_energy, e_n_CS, max_epsilon, &nu_max, &P_max, &N_c,
		  m_e, n_n, dt, np_e, N_coll[0], data_set_length);

  cout << "nu_max = " << nu_max << endl;
  cout << "N_c = " << N_c << endl;
  cout << "P_max = " << P_max << endl; 

  int rand_index;
  int type;
  double epsilon_exc, epsilon_ion;

  for (int i = 0; i < N_c; ++i) {
    rand_index = round((double(rand())/RAND_MAX)*(np_e-1));
    type = getCollType(CS_energy, e_n_CS, part_e_epsilon[rand_index],
		    nu_max, m_e, n_n, N_coll[0], data_set_length);
    cout << "Particle  = " << rand_index << ", Type = " << type << endl;

    // switch-case for ion-neutral collisions
    // 0 - elastic
    // 1 - excitation 1
    // 2 - excitation 2
    // 3 - ionization
    // 4 - null
    switch(type) {
      case 0:
	cout << "Previous velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;
	e_elastic(&part_e_vx[rand_index], &part_e_vy[rand_index],
		  &part_e_vz[rand_index], part_e_epsilon[rand_index],
		  m_e, m_n);
        part_e_epsilon[rand_index] = 0.5*m_e*pow(getv(part_e_vx[rand_index], 
			part_e_vy[rand_index], part_e_vz[rand_index]),2.0)/e; //[eV]

	cout << "New velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;
	continue;

      case 1:
        cout << "Previous velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;

	epsilon_exc = 1.160330e1; // From crtrs.dat.txt
	e_excitation(&part_e_vx[rand_index], &part_e_vy[rand_index],
		  &part_e_vz[rand_index], part_e_epsilon[rand_index], 
		  epsilon_exc);
	part_e_epsilon[rand_index] = 0.5*m_e*pow(getv(part_e_vx[rand_index], 
			part_e_vy[rand_index], part_e_vz[rand_index]),2.0)/e; //[eV]
	cout << "New velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;
	continue;


      case 2:
	cout << "Previous velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;

	epsilon_exc = 1.31041e1; // From crtrs.dat.txt
	e_excitation(&part_e_vx[rand_index], &part_e_vy[rand_index],
		  &part_e_vz[rand_index], part_e_epsilon[rand_index], 
		  epsilon_exc);
	part_e_epsilon[rand_index] = 0.5*m_e*pow(getv(part_e_vx[rand_index], 
			part_e_vy[rand_index], part_e_vz[rand_index]),2.0)/e; //[eV]
	cout << "New velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;
	continue;


      case 3:
	cout << "Previous velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;np_e += 1;
	part_e_vx[np_e-1] = 0.0;
	part_e_vy[np_e-1] = 0.0;
	part_e_vz[np_e-1] = 0.0;
	epsilon_ion = 1.60055e1; // From crtrs.dat.txt
	e_ionization(&part_e_vx[rand_index], &part_e_vy[rand_index],
		  &part_e_vz[rand_index], &part_e_vx[np_e-1], &part_e_vy[np_e-1],
		  &part_e_vz[np_e-1], part_e_epsilon[rand_index], 
		  epsilon_ion);
	part_e_epsilon[rand_index] = 0.5*m_e*pow(getv(part_e_vx[rand_index], 
			part_e_vy[rand_index], part_e_vz[rand_index]),2.0)/e; //[eV]
	part_e_epsilon[np_e-1] = 0.5*m_e*pow(getv(part_e_vx[np_e-1], 
			part_e_vy[np_e-1], part_e_vz[np_e-1]),2.0)/e; //[eV]
	cout << "New velocity: " << part_e_vx[rand_index];
	cout << ", " << part_e_vy[rand_index] << ", ";
	cout << part_e_vz[rand_index] << endl;
	continue;

	
      case 4:
	continue;
    }
  }
   
  delete(part_e_x);
  delete(part_e_vx);
  delete(part_e_vy);
  delete(part_e_vz);
  delete(N_coll);
  
  return 0;
}

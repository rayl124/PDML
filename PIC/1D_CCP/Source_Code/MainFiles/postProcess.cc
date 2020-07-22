#include <math.h>
#include <stdio.h>
#include "collisionModules.h"
#include <iostream>
#include <fstream>
using namespace std;

// Cell centered weights
void get_WeightsCC(double part_x,
		 double *weights,
		 int *cell_index,
		 int *elec_range,
		 double dx) {

  *cell_index = floor((part_x-0.5*dx)/dx);
  if (part_x < elec_range[1]*dx + 0.5*dx) {
    *cell_index = elec_range[1];
    weights[0] = 1.0;
    weights[1] = 0.0;
  } else if (part_x >= elec_range[2]*dx - 0.5*dx) {
    weights[0] = 1.0;
    weights[1] = 0.0;
  } else {
  double hx = (part_x -  (*cell_index + 0.5) * dx)/dx;
  weights[0] = 1.0 - hx;
  weights[1] = hx;
  }
}

int main(int argc, char **argv) 
{
  const double e_mass = 9.109e-31;
  const double i_mass = 39.948*1.66e-27;
  const double ech = 1.602e-19;
  
  // Problem discretization
  int nn = 191; // # of x1 nodes
  int n_cell = nn-1; // # of cells
  double x_start = 0.0;
  double x_end = 0.095;
  double x_bias = 0.01;
  double x_ground = 0.085;

  double L = x_end-x_start; // Domain length m, including electrodes
  double dx= L/((double) n_cell); // length of each cell
  double esmall = 1.0e-9;

  double L_inner = x_ground-x_bias;
  int    N_inner = (int) round(L_inner/dx + esmall);

  // Electrode info
  // Bias electrode left, grounded electrode right
  int *elec_range = new int[4]; 
  elec_range[0] = 0; // Left bias electrode
  elec_range[1] = round((x_bias-x_start)/dx + esmall);
  elec_range[2] = round((x_ground-x_start)/dx + esmall); // left side ground
  elec_range[3] = round ((x_end-x_start)/dx + esmall); // right side ground
  
  int iter_total = 2.4e6;
  int iter_write = 200;
  int num_iter = iter_total/iter_write + 1; // Number of entries in numberpart.txt
  
  string run_num = argv[1];
  string resultPath = "../../Comet_Results/1D_CCP_Results.run"+run_num+"/";
  
  ifstream input_data(resultPath+"Input.txt");
  input_data.ignore(1e5,'\n');
  input_data.ignore(1e5,'\n');
  double electron_spwt;
  for (int i = 0; i < 4; ++i) {input_data >> electron_spwt;}

  ifstream num_data(resultPath+"NumberPart.txt");
  num_data.ignore(1e5, '\n');
  
  int *iters = new int[num_iter];
  int *e_total = new int[num_iter];
  int *i_total = new int[num_iter];
  int *e_inner = new int[num_iter];
  int *i_inner = new int[num_iter];
  double mass;
  for (int i = 0; i < num_iter; ++i) {
    num_data >> iters[i];
    //cout << iters[i] << endl;
    num_data >> e_total[i];
    num_data >> i_total[i];
    num_data >> e_inner[i];
    num_data >> i_inner[i];
    num_data >> mass;    
  }
  num_data.close();

  double *weights = new double[2];
  double *n 	  = new double[5*n_cell];
  double *ux 	  = new double[5*n_cell];
  double *uy 	  = new double[5*n_cell];
  double *uz      = new double[5*n_cell];
  double *Te 	  = new double[5*n_cell];
  double *average_epsilon = new double[5*n_cell];
  int cell_index;
 
  for (int i = 0; i < 5*n_cell; ++i) {n[i] = 0.0; Te[i] = 0.0; ux[i] = 0.0; uy[i] = 0.0; uz[i] = 0.0; average_epsilon[i] = 0.0;}
    
  // Populate electron arrays
  int total_part = 0;
  
  int iter_interval = iter_total/5;
  for (int i = 0; i < 5; ++i) {
      int index = (i+1)*(num_iter-1)/(5); 
      total_part += e_total[index];
  }
  int    *e_iter = new int[total_part];
  double *x      = new double[total_part];
  double *vx     = new double[total_part];
  double *vy     = new double[total_part];
  double *vz     = new double[total_part];
  double epsilon;
  
  cout << "Loading electron arrays" << endl;
  int count = 0;
  for (int rank = 0; rank < 24; ++rank) {
    cout << rank << endl;
    string rank_string = to_string(rank);
    ifstream electron_data(resultPath+"Electron/Electrons"+rank_string+".txt");
    electron_data.ignore(1e5, '\n');    
 
    electron_data >> e_iter[count];
    electron_data >> x[count];
    electron_data >> vx[count];
    electron_data >> vy[count];
    electron_data >> vz[count];
    electron_data >> epsilon;
    count++;
    while (!electron_data.eof()) {
      //cout << count << endl;
      electron_data >> e_iter[count];
      electron_data >> x[count];
      electron_data >> vx[count];
      electron_data >> vy[count];
      electron_data >> vz[count];
      electron_data >> epsilon; 
      if (x[count]!=0) {count++;}
    }
    //cout << x[count-1]  << endl;
    //cout << "Count: " << count << endl;
    electron_data.close();
  }



  int curr_iter, offset;
  int three_count = 0;
  cout << "Calculating diagnostics" << endl; 
  for (int i = 0; i < total_part; ++i) {
    curr_iter = e_iter[i]; 
    offset = curr_iter/iter_interval - 1;
    if (offset==3) {three_count++;}
    get_WeightsCC(x[i], weights, &cell_index, elec_range, dx);
    
    n[offset*n_cell + cell_index] += electron_spwt*weights[0]/dx;
    n[offset*n_cell + cell_index + 1] += electron_spwt*weights[1]/dx;
    ux[offset*n_cell + cell_index] += electron_spwt*weights[0]*vx[i]/dx;
    ux[offset*n_cell + cell_index + 1] += electron_spwt*weights[1]*vx[i]/dx;
    uy[offset*n_cell + cell_index] += electron_spwt*weights[0]*vy[i]/dx;
    uy[offset*n_cell + cell_index + 1] += electron_spwt*weights[1]*vy[i]/dx;
    uz[offset*n_cell + cell_index] += electron_spwt*weights[0]*vz[i]/dx;
    uz[offset*n_cell + cell_index + 1] += electron_spwt*weights[1]*vz[i]/dx;

    average_epsilon[offset*n_cell + cell_index] += 0.5*e_mass*electron_spwt*weights[0]/(ech*dx) *
				      		   (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
 
    average_epsilon[offset*n_cell + cell_index + 1] += 0.5*e_mass*electron_spwt*weights[1]/(ech*dx) *
				      		       (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); 
  }
  
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < n_cell; ++j) {
      if (n[i*n_cell+j] > 1.0) {
        ux[i*n_cell + j] /= n[i*n_cell+j];
        uy[i*n_cell + j] /= n[i*n_cell+j];
        uz[i*n_cell + j] /= n[i*n_cell+j];
        average_epsilon[i*n_cell+j] /= n[i*n_cell+j];
     }
      Te[i*n_cell+j] = (2.0/3.0)*(average_epsilon[i*n_cell+j] -
		       0.5*e_mass/(ech)*(pow(ux[i*n_cell+j],2.0) + pow(uy[i*n_cell+j],2.0) + pow(uz[i*n_cell+j],2.0)));		       
    }     
  }
  /*
  for (int i = 0; i < total_part; ++i) { 
    curr_iter = e_iter[i]; 
    offset = curr_iter/iter_interval - 1;
    
    get_WeightsCC(x[i], weights, &cell_index, elec_range, dx);
    
    Te[offset*n_cell + cell_index] += 0.5*e_mass*electron_spwt*weights[0]/(ech*dx) *
				      (pow(vx[i] - ux[offset*n_cell + cell_index], 2.0) +
				      pow(vy[i] - uy[offset*n_cell + cell_index], 2.0) +
				      pow(vz[i] - uz[offset*n_cell + cell_index], 2.0)); 
    
    Te[offset*n_cell + cell_index + 1] += 0.5*e_mass*electron_spwt*weights[1]/(ech*dx) *
				          (pow(vx[i] - ux[offset*n_cell + cell_index], 2.0) +
				          pow(vy[i] - uy[offset*n_cell + cell_index], 2.0) +
				          pow(vz[i] - uz[offset*n_cell + cell_index], 2.0)); 
  }
  for (int i = 0; i < 5*n_cell; ++i) {
    if (n[i] > 0.0) {
      Te[i] /= n[i];
    }
  }*/
  
  ofstream diagFile(resultPath+"Electron/ElectronDiagnostics.txt");
  diagFile << "Iter / x / n / ux / uy / uz / average epsilon/ Te" << endl;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < n_cell; ++j) {
      diagFile << (i+1)*iter_interval << " " << (j+0.5)*dx << " " << n[i*n_cell + j] << " ";
      diagFile << ux[i*n_cell + j] << " " << uy[i*n_cell+j] << " " << uz[i*n_cell+j] << " ";
      diagFile << average_epsilon[i*n_cell + j] << " " << Te[i*n_cell + j] << endl;
    }
  }
  delete(elec_range);
  delete(iters);
  delete(e_total);
  delete(i_total);
  delete(e_inner);
  delete(i_inner);
  delete(e_iter);
  delete(x);
  delete(vx);
  delete(vy);
  delete(vz);
  delete(weights);
  delete(Te);
  delete(n);
  delete(ux);
  delete(uy);
  delete(uz);
  delete(average_epsilon);
}

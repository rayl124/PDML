#include "fluidModules.h"
#include <iostream>

using namespace std;

int main(void)
{
  double k_B = 1.381e-23;
  double AMU = 1.661e-27;
  double P = 101325.0; // atm
  double T = 353.2*k_B;

  double epsilon = 124.0*k_B; // K
  double d = 3.418e-10; // Angstrom
  double m = 39.948e-3; // [g/mol]
  double D;
  D = diffusionCoeff(T, epsilon, d, P, m);

  
 
  cout << D << endl;
}

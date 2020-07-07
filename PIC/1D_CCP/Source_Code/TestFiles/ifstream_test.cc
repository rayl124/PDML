#include <iostream>
#include <fstream>

using namespace std;

int main(void) {
  ifstream test_data("arp_ar_backscattering_02.xsn");
  double digit;
  double digit2;
  
  // Ignore the first 12 lines
  for (int i = 0; i < 12; i++) {
    test_data.ignore(1e5, '\n');
  }

  while(test_data >> digit >> digit2) {
    cout << digit2 << endl;
  }
}

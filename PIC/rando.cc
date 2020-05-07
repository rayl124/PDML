#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;
int main() {
srand(time(nullptr));
  for (int i = 0; i < 10 ; ++i) {
	  
	  cout<< double(rand())/(RAND_MAX) + double(rand())/(RAND_MAX)<< std::endl;
  }

}

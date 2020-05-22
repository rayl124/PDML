#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;
double inside_fun(void) {
	return double(rand())/(RAND_MAX);

}

int main(void) {
srand(time(nullptr));
cout << "running rando" << endl;
  for (int i = 0; i < 10 ; ++i) {

     
      cout << inside_fun() << endl;
   }
return 0;
}

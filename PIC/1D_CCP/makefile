CC = g++

ifstream_test: ifstream_test.cc
	$(CC) -g -o ifstream_test ifstream_test.cc

1D_Poisson_test: 1D_Poisson_test.cc 1D_Poisson.h
	$(CC) -g -o 1D_Poisson_test 1D_Poisson_test.cc -lm

interpolate_test: interpolate_test.cc
	$(CC) -g -o interpolate_test interpolate_test.cc

collision_test: collision_test.cc collisionModules.h
	$(CC) -g -o collision_test collision_test.cc -lm

1D_CCP: 1D_CCP.cc
	$(CC) -g -o 1D_CCP 1D_CCP.cc -lm

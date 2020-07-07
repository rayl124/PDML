CC = g++
MPI_CC = mpicxx

ifstream_test: ifstream_test.cc
	$(CC) -g -o ifstream_test ifstream_test.cc

1D_Poisson_test: 1D_Poisson_test.cc solverModules.h
	$(CC) -g -o 1D_Poisson_test 1D_Poisson_test.cc -lm

interpolate_test: interpolate_test.cc
	$(CC) -g -o interpolate_test interpolate_test.cc

collision_test: collision_test.cc collisionModules.h
	$(CC) -g -o collision_test collision_test.cc -lm

driftDiffusionTest: driftDiffusionTest.cc
	$(CC) -g -o driftDiffusionTest driftDiffusionTest.cc -lm

simpsonsTest: simpsonsTest.cc
	$(CC) -g -o simpsonsTest simpsonsTest.cc -lm

getDTest: getDTest.cc
	$(CC) -g -o getDTest getDTest.cc -lm

1D_CCP: 1D_CCP.cc
	$(MPI_CC) -Ofast -o 1D_CCP 1D_CCP.cc -lm

uploadGIT: 
	cd ..; \
	git add -A 1D_CCP/
	git commit -m "New simulation"
	git push origin master

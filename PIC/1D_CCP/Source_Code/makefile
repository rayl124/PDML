CC = g++
MPI_CC = mpiCC

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
	$(MPI_CC) -g -o 1D_CCP.out 1D_CCP.cc -lm

uploadGIT: 
	cd ..; \
	git add -A Source_Code/
	git commit -m "Updated Source Code"
	git push origin master

clean:
	rm ifstream_test 1D_Poisson_test interpolate_test collision_test driftDiffusionTest simpsonsTest getDTest 1D_CCP

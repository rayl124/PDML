#!/bin/bash

#SBATCH --job-name=1D_CCP_Results	#Set job name
#SBATCH --time=00:01:00			#Set time limits
#SBATCH --nodes=1			#Set nodes
#SBATCH --ntasks-per-node=32		#32 cores
#SBATCH --mem=4G			#4 GB per node
#SBATCH --output=1D_CCPOut.%j.out	#Output File

##Optional job specifications
#SBATCH --mail-type=BEGIN,END,FAIL      #email events
#SBATCH --mail-user=rl49@stanford.edu	#email


NJOBS=2
JOBNAME=$SLURM_JOB_NAME

echo -n 'Started job at : ' ; date
echo " "

mkdir ../Results/$JOBNAME.run1
cp -r MainFiles ../Results/$JOBNAME.run1 
cd ../Results/$JOBNAME.run1/MainFiles
mpiicc -Ofast -o 1D_CCP 1D_CCP.cc -lm
mpiexec -n 16 -v ./1D_CCP 0.3 & 

cd ..

mkdir ../$JOBNAME.run2
cp -r MainFiles ../$JOBNAME.run2 
cd ../$JOBNAME.run2/MainFiles
mpiicc -Ofast -o 1D_CCP 1D_CCP.cc -lm
mpiexec -n 16 -v ./1D_CCP 1.4 &

wait

echo -n 'Ended job at : ' ; date
echo " "
exit


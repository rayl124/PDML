#!/bin/bash
#SBATCH --job-name=1D_CCP_Results	# Job name
#SBATCH --output="1D_CCP.%j.out"	# Output file
#SBATCH --partition=compute		# Which queue to use
#SBATCH --nodes=1			# Number of nodes
#SBATCH --ntasks-per-node=16		# Cores per node
#SBATCH --export=ALL			# Unsure what this is
#SBATCH -t 00:03:30			# Time limit
					# At minimum, time*cores = 1 hour

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rl49@stanford.edu

JOBNAME=$SLURM_JOB_NAME

# Modules
module purge
module load gnu openmpi_ib

# Create Folder
mkdir ../Results/$JOBNAME.run1
cp -r MainFiles ../Results/$JOBNAME.run1
cd ../Results/$JOBNAME.run1/MainFiles

# Compilation code: 
mpicxx -Ofast -o 1D_CCP 1D_CCP.cc -lm
ibrun -v ./1D_CCP 0.3



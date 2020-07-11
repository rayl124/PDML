#/bin/bash
#SBATCH --job-name=1D_CCP_Results	# Job name
#SBATCH --output="1D_CCP.%j.out"	# Output file
#SBATCH --partition=compute		# Which queue to use
#SBATCH --nodes=1			# Number of nodes
#SBATCH --ntasks-per-node=24		# Cores per node
#SBATCH --export=ALL			
#SBATCH -t 04:00:00			# Time limit
					# At minimum, time*cores = 1 hour

#SBATCH --mail-type=BEGIN,END,FAIL	# What events will be mailed to you
#SBATCH --mail-user=rl49@stanford.edu	# Email these events will be sent to

JOBNAME=$SLURM_JOB_NAME

# Modules
module purge				# Erases any current modules
module load gnu openmpi_ib		# Load module for compiling

# Create folder to store results in
mkdir ../Results/$JOBNAME.run9
cp -r MainFiles ../Results/$JOBNAME.run9
cd ../Results/$JOBNAME.run9/MainFiles
mkdir ../Ion
mkdir ../Electron

# Compilation code: 
mpicxx -Ofast -o 1D_CCP 1D_CCP.cc -lm 	# Compiles code
ibrun -n 24 --binding-level=core ./1D_CCP 500000 2.6 # Runs code with an input of 0.3 Pascals

# Uncomment for extra job
#cd ../../../Source_Code
#mkdir ../Results/$JOBNAME.run7
#cp -r MainFiles ../Results/$JOBNAME.run7
#cd ../Results/$JOBNAME.run7/MainFiles
#mkdir ../Ion
#mkdir ../Electron

#mpicxx -Ofast -o 1D_CCP 1D_CCP.cc -lm 	# Compiles code
#ibrun -n 8 -o 16 --binding-level=core ./1D_CCP 100000 0.3 &	# Runs code with an input of 0.3 Pascals

wait

exit

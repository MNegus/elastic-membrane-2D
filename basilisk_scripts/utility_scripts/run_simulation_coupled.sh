#!/bin/bash

# run_simulation.sh
# Bash script to run a simulation. It takes two inputs:
# Input 1: Name of the C file (without the .C extension)
# Input 2: Number of threads to run the simulation on (default 1)
# It removes any previous outputs, runs the code and then moves the output into 
# a directory one level up called "raw_data"

# Saves script name, which will also be the name of the directory that the 
#output gets saved for (crucially this does not contain the .c extension)
script_name=$1



# Make a directory
mkdir ${script_name}

# Copy the relevant files into the directory
cp {${script_name}.c,parameters.h,wave-equation.c,wave-equation.h} ${script_name}

# Move into the directory
cd ${script_name}/

# Makes the wave equation object
gcc -O2 -c wave-equation.c 

# Makes the main object
qcc -O2 -c ${script_name}.c 

# Compiles
qcc wave-equation.o ${script_name}.o -o a.out -lm -llapacke -llapack -O2 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -I$BASILISK/Makefile.defs -Wall

# Sets the number of OpenMP threads
export OMP_NUM_THREADS=$2

# Runs
./a.out 

# # Moves back upwards
# cd ../

# # Remove any existing raw data in the parent directory
# rm -r ../raw_data 

# # Moves the new data to the parent directory
# mv ${script_name} ../raw_data


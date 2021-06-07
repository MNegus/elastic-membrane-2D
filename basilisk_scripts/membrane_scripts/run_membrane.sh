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

# Removes previous membrane output folder
rm -r ../membrane_outputs
mkdir ../membrane_outputs

# Removes last code
rm a.out

# Compiles program
gcc $script_name wave-equation.c -llapacke -llapack -lm

# Runs the executable
./a.out     



#!/bin/bash

# Script to run the validation procedure, which compares the finite difference 
# code to the exact solutions for varying computational parameters

# Make a valdidation data directory
mkdir validation/forced_validation_data

# Go into the C_code directory
cd C_code/

# Remove old output directories and create new ones
rm -r exact_outputs
mkdir exact_outputs

rm -r mitchell_outputs
mkdir mitchell_outputs



##############################
# Varying spatial grid size
##############################

# Pick a specific timestep value
DT=1e-4
sed -i "/const double DELTA_T/c\const double DELTA_T = $DT; // Timestep size" parameters.h

# Loop over grid sizes
for N in 16 32 64 128 256 512 1024 2048 4096
do
    # Edit the value of N_MEMBRANE in the parameters.h file
    sed -i "/const int N_MEMBRANE/c\const int N_MEMBRANE = $N; // Number of grid points on the membrane" parameters.h

    # Run the exact solution
    ./run_code.sh exact_forced_membrane.c 

    # Run the implicit solution
    ./run_code.sh forced_mitchell_fd_wave.c 

    # Make a new directory in the validation data directory
    rm -r ../validation/forced_validation_data/N_MEMBRANE_$N 
    mkdir ../validation/forced_validation_data/N_MEMBRANE_$N 

    # Copy parameters file into the data directory
    cp parameters.h ../validation/forced_validation_data/N_MEMBRANE_$N 
    cp anim_all.gp ../validation/forced_validation_data/N_MEMBRANE_$N 

    # Move outputs into the data directory
    mv exact_outputs ../validation/forced_validation_data/N_MEMBRANE_$N/
    mv mitchell_outputs ../validation/forced_validation_data/N_MEMBRANE_$N 

    # Create new directories for next time
    mkdir exact_outputs
    mkdir mitchell_outputs

    # Output N 
    echo Finished N = $N
done

##############################
# Varying timestep size
##############################

# Pick a specific membrane size
N=4096
sed -i "/const int N_MEMBRANE/c\const int N_MEMBRANE = $N; // Number of grid points on the membrane" parameters.h

# Loop over timestep values
# for DTPOWER in 1 2 3 3.5 4 4.5 5
DTPOWER=1
for DT in 1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5
do
    # DT=1e-($(DTPOWER))
    # Edit the value of N_MEMBRANE in the parameters.h file
    sed -i "/const double DELTA_T/c\const double DELTA_T = $DT; // Timestep size" parameters.h

    # Run the exact solution
    ./run_code.sh exact_forced_membrane.c 

    # # Run the implicit solution
    # ./run_code.sh mitchell_fd_wave.c 
    # Run the mitchell solution
    ./run_code.sh forced_mitchell_fd_wave.c 

    # Make a new directory in the validation data directory
    rm -r ../validation/forced_validation_data/DT_$DTPOWER
    mkdir ../validation/forced_validation_data/DT_$DTPOWER

    # Copy parameters file into the data directory
    cp parameters.h ../validation/forced_validation_data/DT_$DTPOWER
    cp anim_all.gp ../validation/forced_validation_data/DT_$DTPOWER

    # Move outputs into the data directory
    mv exact_outputs ../validation/forced_validation_data/DT_$DTPOWER/
    mv mitchell_outputs ../validation/forced_validation_data/DT_$DTPOWER 

    # Create new directories for next time
    mkdir exact_outputs
    mkdir mitchell_outputs

    # Output DT
    echo Finished DT = $DT

    DTPOWER=$((DTPOWER+1))
done
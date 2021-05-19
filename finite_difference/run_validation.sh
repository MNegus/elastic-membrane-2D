#!/bin/bash

# Script to run the validation procedure, which compares the finite difference 
# code to the exact solutions for varying computational parameters

# Go into the C_code directory
cd C_code/

# Remove old output directories and create new ones
rm -r exact_outputs
mkdir exact_outputs

# rm -r implicit_outputs
# mkdir implicit_outputs
rm -r mitchell_outputs
mkdir mitchell_outputs



##############################
# Varying spatial grid size
##############################

# # Pick a specific timestep value
# DT=1e-3
# sed -i "/const double DELTA_T/c\const double DELTA_T = $DT; // Timestep size" parameters.h

# # Loop over grid sizes
# for N in 16 32 64 128 256 512 1024  
# do
#     # Edit the value of N_MEMBRANE in the parameters.h file
#     sed -i "/const int N_MEMBRANE/c\const int N_MEMBRANE = $N; // Number of grid points on the membrane" parameters.h

#     # Run the exact solution
#     ./run_code.sh exact_membrane.c 

#     # Run the implicit solution
#     ./run_code.sh FIXED_ORIGIN_mitchell_fd_wave.c 

#     # Make a new directory in the validation data directory
#     rm -r ../validation/validation_data/N_MEMBRANE_$N 
#     mkdir ../validation/validation_data/N_MEMBRANE_$N 

#     # Copy parameters file into the data directory
#     cp parameters.h ../validation/validation_data/N_MEMBRANE_$N 
#     cp anim_all.gp ../validation/validation_data/N_MEMBRANE_$N 

#     # Move outputs into the data directory
#     mv exact_outputs ../validation/validation_data/N_MEMBRANE_$N/
#     mv mitchell_outputs ../validation/validation_data/N_MEMBRANE_$N 

#     # Create new directories for next time
#     mkdir exact_outputs
#     mkdir mitchell_outputs

#     # Output N 
#     echo Finished N = $N
# done

##############################
# Varying timestep size
##############################

# # Pick a specific membrane size
N=1024
sed -i "/const int N_MEMBRANE/c\const int N_MEMBRANE = $N; // Number of grid points on the membrane" parameters.h

# Loop over timestep values
for DTPOWER in 1 2 3 4
do
    DT=1e-$DTPOWER
    # Edit the value of N_MEMBRANE in the parameters.h file
    sed -i "/const double DELTA_T/c\const double DELTA_T = $DT; // Timestep size" parameters.h

    # Run the exact solution
    ./run_code.sh exact_membrane.c 

    # # Run the implicit solution
    # ./run_code.sh mitchell_fd_wave.c 
    # Run the mitchell solution
    ./run_code.sh FIXED_ORIGIN_mitchell_fd_wave.c 

    # Make a new directory in the validation data directory
    rm -r ../validation/validation_data/DT_$DTPOWER
    mkdir ../validation/validation_data/DT_$DTPOWER

    # Copy parameters file into the data directory
    cp parameters.h ../validation/validation_data/DT_$DTPOWER
    cp anim_all.gp ../validation/validation_data/DT_$DTPOWER

    # Move outputs into the data directory
    mv exact_outputs ../validation/validation_data/DT_$DTPOWER/
    mv mitchell_outputs ../validation/validation_data/DT_$DTPOWER 

    # Create new directories for next time
    mkdir exact_outputs
    mkdir mitchell_outputs

    # Output DT
    echo Finished DT = $DT
done
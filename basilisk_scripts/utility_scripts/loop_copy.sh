#!/bin/bash

code_dir=$1
dest_dir=$2

for GAMMA in 0.01 0.1 1 10 100
do
    # Changes the parameters file
    # sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h
    sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" parameters.h

    # Copies over the code 
    ./code_copy.sh $code_dir $dest_dir gamma_$GAMMA
done
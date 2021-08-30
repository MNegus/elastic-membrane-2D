#!/bin/bash

code_dir=$1
dest_dir=$2

ALPHA=1;

for GAMMA in 0.8 1.6 3.2 6.4 12.8
do
    for BETA in 0
    do 
        # Changes the parameters file
        sed -i "/BETA/c\const double BETA = $BETA; // Tension term" parameters.h
        sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" parameters.h

        # Copies over the code 
        ./code_copy.sh $code_dir $dest_dir alpha_$ALPHA-beta_$BETA-gamma_$GAMMA
    done
done

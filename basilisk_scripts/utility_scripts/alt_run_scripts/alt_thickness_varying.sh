#!/bin/bash

# Code for running the thickness varying scripts

dest_dir=$1
code_name=$2

parent_dir=$dest_dir/thickness_varying

BETA=0;
for ALPHA in 3 4 6 8
do
    # Determines GAMMA = 2 * ALPHA^3
    GAMMA=$(bc <<< "scale=2;2*$ALPHA*$ALPHA*$ALPHA")
    echo Thickness varying: ALPHA=$ALPHA, BETA=$BETA, GAMMA=$GAMMA

    # Moves into the code directory
    cd $parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code

    # Runs the script
    nohup ./run_simulation.sh $code_name 8 &> nohup_alt_thickness.out
done
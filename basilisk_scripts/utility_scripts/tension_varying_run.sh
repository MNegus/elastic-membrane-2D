#!/bin/bash

# Code for running the thickness varying scripts

dest_dir=$1
code_name=$2

parent_dir=$dest_dir/tension_varying

ALPHA=1;
GAMMA=2;

for BETA in 10 40 160 640 2560 10240
do
    echo Tension varying: ALPHA=$ALPHA, BETA=$BETA, GAMMA=$GAMMA

    # Moves into the code directory
    cd $parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code

    # Runs the script
    ./run_simulation.sh $code_name 8
done
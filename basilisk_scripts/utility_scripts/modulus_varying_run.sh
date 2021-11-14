#!/bin/bash

# Code for running the Young's modulus varying scripts

dest_dir=$1
code_name=$2

parent_dir=$dest_dir/modulus_varying

ALPHA=1;
BETA=0;

for GAMMA in 8 32 128 512 2048 8192
do
    echo Young\'s modulus varying: ALPHA=$ALPHA, BETA=$BETA, GAMMA=$GAMMA

    # Moves into the code directory
    cd $parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code

    # Runs the script
    ./run_simulation.sh $code_name 8
done
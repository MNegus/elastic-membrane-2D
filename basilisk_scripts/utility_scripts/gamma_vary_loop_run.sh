#!/bin/bash

code_dir=$1
dest_dir=$2

ALPHA=1;

for GAMMA in 0.8 1.6 3.2 6.4 12.8
do
    for BETA in 0
    do 
        echo BETA = $BETA, GAMMA = $GAMMA
        
        cd $parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code

        ./run_simulation.sh $code_name 8
    done
done
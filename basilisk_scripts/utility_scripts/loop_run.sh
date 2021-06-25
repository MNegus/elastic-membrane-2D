#!/bin/bash

parent_dir=$1
code_name=$2


for GAMMA in 100 10 1 0.1 0.01
do
    # echo Max level = $MAXLEVEL
    # cd $parent_dir/max_level_$MAXLEVEL-coupled_$COUPLED/code
    
    echo GAMMA = $GAMMA
    cd $parent_dir/gamma_$GAMMA/code

    ./run_simulation.sh $code_name 4

done
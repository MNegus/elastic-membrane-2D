#!/bin/bash

parent_dir=$1
code_name=$2


for BETA in 100 10 1 0.1 0.01
do
    # echo Max level = $MAXLEVEL
    # cd $parent_dir/max_level_$MAXLEVEL-coupled_$COUPLED/code
    
    echo BETA = $BETA
    cd $parent_dir/beta_$BETA/code

    ./run_simulation.sh $code_name 8

done
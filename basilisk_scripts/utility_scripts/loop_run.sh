#!/bin/bash

parent_dir=$1
code_name=$2

for COUPLED in 1 0 
do
    for MAXLEVEL in 8 9 10 11 12
    do
        echo Max level = $MAXLEVEL
        cd $parent_dir/max_level_$MAXLEVEL-coupled_$COUPLED/code
        # ./run_simulation.sh $code_name 4
        ./run_manual.sh $code_name 4
    done
done
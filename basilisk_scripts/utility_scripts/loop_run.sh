#!/bin/bash

parent_dir=$1
code_name=$2
cores=$3

for MAXLEVEL in 13 14 
do
    echo Max level = $MAXLEVEL
    cd $parent_dir/max_level_$MAXLEVEL/code

    ./run_simulation.sh $code_name $cores
done

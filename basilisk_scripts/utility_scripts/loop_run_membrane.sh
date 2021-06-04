#!/bin/bash

parent_dir=$1
code_name=$2

for MAXLEVEL in 8 9 10 11 12
do
    echo Max level = $MAXLEVEL
    cd $parent_dir/max_level_$MAXLEVEL/code
    ./run_membrane.sh $code_name
done
#!/bin/bash

code_dir=$1
dest_dir=$2

for COUPLED in 1 0
do
    # Changes the parameters file
    sed -i "/COUPLED/c\const int COUPLED = $COUPLED; // Set to 1 if coupled" parameters.h
    for MAXLEVEL in 8 9 10 11 12
    do
        # Changes the parameters file
        sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h

        # Copies over the code 
        ./code_copy.sh $code_dir $dest_dir max_level_$MAXLEVEL-coupled_$COUPLED
    done
done
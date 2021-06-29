#!/bin/bash

code_dir=$1
dest_dir=$2

# for GAMMA in 0.01 0.1 1 10 100
for XMIN in 0 1 2 3 4
do
    sed -i "/x_min_height/c\const int x_min_height = $XMIN;" parameters.h
    for YMIN in 0 1 2 3
    do
        # Changes the parameters file
        # sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h
        sed -i "/y_min_height/c\const int y_min_height = $YMIN;" parameters.h

        # Copies over the code 
        ./code_copy.sh $code_dir $dest_dir x_min_$XMIN-y_min_$YMIN

    done
    
done
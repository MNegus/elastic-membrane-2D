#!/bin/bash

code_dir=$1
dest_dir=$2


for IMPOSEDCOEFF in 0 0.05
do
    # Change the axisymmetric value in parameters file
    sed -i "/IMPOSED_COEFF/c\const double IMPOSED_COEFF = $IMPOSEDCOEFF;" parameters.h

    # Make the respective directory
    mkdir $dest_dir/imposed_coeff_$IMPOSEDCOEFF
    
    for MAXLEVEL in 10 11 12 13 14
    do
        # Changes the parameters file
        sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h

        # Copies over the code 
        ./code_copy.sh $code_dir $dest_dir/imposed_coeff_$IMPOSEDCOEFF max_level_$MAXLEVEL
    done
done
#!/bin/bash

code_dir=$1
dest_dir=$2

for TRANSPOSED in 0 1
do 
    # Changes the parameters file
    sed -i "/TRANSPOSED/c\#define TRANSPOSED $TRANSPOSED // Transposes so the membrane is along y" test_parameters.h

    for WALL in 0 1
    do 
        # Changes the parameters file
        sed -i "/WALL/c\#define WALL $WALL // Droplet along the wall" test_parameters.h

        for MEMBRANE in 0 1
        do
            # Changes the parameters file
            sed -i "/#define MEMBRANE/c\#define MEMBRANE $MEMBRANE // Impose a membrane to deform the droplet" test_parameters.h

            for MOVING in 0 1
            do
                # Skip the case where we have the moving frame but no 
                # deformation, as this is pointless
                if [[ $MEMBRANE = 0 ]] && [[ $MOVING = 1 ]]; then
                    continue
                fi

                # Changes the parameters file
                sed -i "/MOVING/c\#define MOVING $MOVING // Moving frame adjustment" test_parameters.h

                # Copy the code into directory
                ./test_copy.sh $1 $2 MOVING_$MOVING-MEMBRANE_$MEMBRANE-WALL_$WALL-TRANSPOSED_$TRANSPOSED

                echo MOVING_$MOVING-MEMBRANE_$MEMBRANE-WALL_$WALL-TRANSPOSED_$TRANSPOSED
            done
        done
    done
done
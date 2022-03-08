#!/bin/bash

dest_dir=$1
code_name=$2

for TRANSPOSED in 1 0
do 
    for WALL in 0 1
    do 
        for MEMBRANE in 0 1
        do
            for MOVING in 0 1
            do
                # Skip the case where we have the moving frame but no 
                # deformation, as this is pointless
                if [[ $MEMBRANE = 0 ]] && [[ $MOVING = 1 ]]; then
                    continue
                fi

                # Move into the code directory
                cd ${dest_dir}/MOVING_$MOVING-MEMBRANE_$MEMBRANE-WALL_$WALL-TRANSPOSED_$TRANSPOSED/code
                echo Entered $PWD
                
                # Run the code
                echo Running code
                ./run_test.sh ${code_name} 1
                echo Finished code
                
            done
        done
    done
done
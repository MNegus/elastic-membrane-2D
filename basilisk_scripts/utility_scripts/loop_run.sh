#!/bin/bash

parent_dir=$1
code_name=$2


# for GAMMA in 100 10 1 0.1 0.01
# do
#     # echo Max level = $MAXLEVEL
#     # cd $parent_dir/max_level_$MAXLEVEL-coupled_$COUPLED/code
    
#     echo GAMMA = $GAMMA
#     cd $parent_dir/gamma_$GAMMA/code

#     ./run_simulation.sh $code_name 4

# done

for XMIN in 0 1 2 3 4
do
    for YMIN in 0 1 2 3
    do
        echo XMIN = $XMIN, YMIN = $YMIN
        
        cd $parent_dir/x_min_$XMIN-y_min_$YMIN/code

        ./run_simulation.sh $code_name 4

    done
    
done
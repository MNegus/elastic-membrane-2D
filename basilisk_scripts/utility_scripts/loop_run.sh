#!/bin/bash

parent_dir=$1
code_name=$2
move_dest=$3

for GAMMA in 0.1 0.2 0.4 0.8 1.6 3.2 6.4 12.8
do
    echo GAMMA = $GAMMA

    cd $parent_dir/gamma_$GAMMA/code

    ./run_simulation.sh $code_name 4

    cd $parent_dir

    mv gamma_$GAMMA $move_dest
done

# for COARSEN in 0 1 2 3
# do
#     echo COARSEN = $COARSEN

#     cd $parent_dir/coarsen_$COARSEN/code

#     ./run_simulation.sh $code_name 8
# done

# for MAXLEVEL in 8 9 10 11 12 13
# do
#     echo Max level = $MAXLEVEL
#     cd $parent_dir/max_level_$MAXLEVEL/code

#     ./run_simulation.sh $code_name 8
# done

# for XMIN in 0 1 2 3 4
# do
#     for YMIN in 0 1 2 3
#     do
#         echo XMIN = $XMIN, YMIN = $YMIN
        
#         cd $parent_dir/x_min_$XMIN-y_min_$YMIN/code

#         ./run_simulation.sh $code_name 4

#     done
    
# done
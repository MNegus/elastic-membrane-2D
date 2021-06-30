#!/bin/bash

code_dir=$1
dest_dir=$2


for MAXLEVEL in 8 9 10 11 12 13
do
    # Changes the parameters file
    sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h

    # Copies over the code 
    ./code_copy.sh $code_dir $dest_dir max_level_$MAXLEVEL
done

# for XMIN in 0 1 2 3 4
# do
#     sed -i "/x_min_height/c\const int x_min_height = $XMIN;" parameters.h
#     for YMIN in 0 1 2 3
#     do
#         # Changes the parameters file
#         # sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h
#         sed -i "/y_min_height/c\const int y_min_height = $YMIN;" parameters.h

#         # Copies over the code 
#         ./code_copy.sh $code_dir $dest_dir x_min_$XMIN-y_min_$YMIN

#     done
    
# done
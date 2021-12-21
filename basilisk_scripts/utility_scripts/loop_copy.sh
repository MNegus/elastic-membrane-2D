#!/bin/bash

code_dir=$1
dest_dir=$2

ALPHA=1;

# for GAMMA in 0.4 0.01
# do
#     for BETA in 0.1 0.2 0.4 0.8 1.6 3.2 6.4 12.8
#     do 
#         # Changes the parameters file
#         sed -i "/BETA/c\const double BETA = $BETA; // Tension term" parameters.h
#         sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" parameters.h

#         # Copies over the code 
#         ./code_copy.sh $code_dir $dest_dir alpha_$ALPHA-beta_$BETA-gamma_$GAMMA
#     done
# done




# for GAMMA in 0.1 0.2 0.4 0.8 1.6 3.2 6.4 12.8
# do 
#     # Changes the parameters file
#     sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" parameters.h

#     # Copies over the code 
#     ./code_copy.sh $code_dir $dest_dir gamma_$GAMMA

# done

# for COARSEN in 0 1 2 3
# do
#     # Changes the parameters file
#     sed -i "/FD_COARSEN_LEVEL/c\const int FD_COARSEN_LEVEL = $COARSEN;" parameters.h

#     # Copies over the code 
#     ./code_copy.sh $code_dir $dest_dir coarsen_$COARSEN
# done

for MAXLEVEL in 12 13
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
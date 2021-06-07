#!/bin/bash

# loop_clean.sh
# Repeatedly runs the output_clean.sh script to clean a selection of outputs

PARENTDIR=$1 # Parent directory i.e. containing the max_level_* directories

for MAXLEVEL in 8 9 10 11 12
do
    echo Max level = $MAXLEVEL
    ./output_clean.sh $PARENTDIR/max_level_$MAXLEVEL
done
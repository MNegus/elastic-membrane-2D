#!/bin/bash

parent_dir=$1

for DIR_NAME in $parent_dir/*/
do
    echo $DIR_NAME

    ./output_clean.sh $DIR_NAME
done

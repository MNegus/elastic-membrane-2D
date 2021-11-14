#!/bin/bash

# Code for copying scripts for the parameter sweeping

code_dir=$1
dest_dir=$2

# Thickness varying
parent_dir=$dest_dir/thickness_varying
mkdir $parent_dir

BETA=0;
for ALPHA in 1 1.5 2 3 4 6 8
do
    # Determines GAMMA = 2 * ALPHA^3
    GAMMA=$(bc <<< "scale=2;2*$ALPHA*$ALPHA*$ALPHA")
    echo Thickness varying: $ALPHA $BETA $GAMMA

    # Copies over the code 
    ./code_copy.sh $code_dir $parent_dir alpha_$ALPHA-beta_$BETA-gamma_$GAMMA

    # Changes the parameters file in the destination
    param_file=$parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code/parameters.h
    sed -i "/ALPHA/c\const double ALPHA = $ALPHA; // Mass ratio" $param_file
    sed -i "/BETA/c\const double BETA = $BETA; // Tension term" $param_file
    sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" $param_file
done

# Tension varying
parent_dir=$dest_dir/tension_varying
mkdir $parent_dir

ALPHA=1;
GAMMA=2;

for BETA in 10 40 160 640 2560 10240
do
    echo Tension varying: $ALPHA $BETA $GAMMA

    # Copies over the code 
    ./code_copy.sh $code_dir $parent_dir alpha_$ALPHA-beta_$BETA-gamma_$GAMMA

    # Changes the parameters file in the destination
    param_file=$parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code/parameters.h
    sed -i "/ALPHA/c\const double ALPHA = $ALPHA; // Mass ratio" $param_file
    sed -i "/BETA/c\const double BETA = $BETA; // Tension term" $param_file
    sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" $param_file
done

# Young's modulus varying
parent_dir=$dest_dir/modulus_varying
mkdir $parent_dir

ALPHA=1;
BETA=0;

for GAMMA in 8 32 128 512 2048 8192
do
    echo Young\'s modulus varying: $ALPHA $BETA $GAMMA

    # Copies over the code 
    ./code_copy.sh $code_dir $parent_dir alpha_$ALPHA-beta_$BETA-gamma_$GAMMA

    # Changes the parameters file in the destination
    param_file=$parent_dir/alpha_$ALPHA-beta_$BETA-gamma_$GAMMA/code/parameters.h
    sed -i "/ALPHA/c\const double ALPHA = $ALPHA; // Mass ratio" $param_file
    sed -i "/BETA/c\const double BETA = $BETA; // Tension term" $param_file
    sed -i "/GAMMA/c\const double GAMMA = $GAMMA; // Bending term" $param_file
done



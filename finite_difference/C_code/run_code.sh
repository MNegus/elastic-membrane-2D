#!/bin/bash

# Compiles and runs the C code given by CODENAME

CODENAME=$1

# Removes any previous output files
rm outputs/*

# Compiles program
gcc $CODENAME -llapacke -llapack -lm

# Runs the executable
./a.out     

#!/bin/bash

# code_copy.sh
# Script to copy the code into the directory it is to be run in (i.e. the 
# /scratch directory)
# The first input is the local directory of the code (e.g. solid_wall, 
# coupled_plate etc), the second directory is the parent destination directory
# and the third input is the sub-directory name (i.e. the name of the directory
# the code will be stored in within the scratch)

LOCAL_DIR=$1 # Local directory of the code
DEST_DIR=$2 # Destination of the parent directory
SUB_DIR_NAME=$3 # Name of the sub-directory

# Makes the sub-directory
mkdir ${DEST_DIR}/${SUB_DIR_NAME}
mkdir ${DEST_DIR}/${SUB_DIR_NAME}/code

# Copies the code over to the destination
cp -a ${LOCAL_DIR}/. ${DEST_DIR}/${SUB_DIR_NAME}/code

# Copies the run script, Makefile and parameters over to the destination
# cp {run_simulation.sh,run_simulation_coupled.sh,run_manual.sh,Makefile,parameters.h} ${DEST_DIR}/${SUB_DIR_NAME}/code
cp run_test.sh ${DEST_DIR}/${SUB_DIR_NAME}/code
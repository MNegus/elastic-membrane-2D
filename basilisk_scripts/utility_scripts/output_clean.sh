#!/bin/bash

# output_clean.sh
# Cleans the data in a results directory. This includes going into the raw_data/ 
# directory and sorting the pressure output files into ascending order

PARENTDIR=$1 # Parent directory i.e. containing code/ and raw_data/

cd $PARENTDIR

rm -r boundary_outputs
mkdir boundary_outputs # Directory to store boundary output files

rm -r membrane_outputs
mkdir membrane_outputs # Directory to store membrane outputs

rm -r videos
mkdir videos

# Works out how many boundary_output_*.txt files there are in raw_data/
NOFILES=$(ls raw_data/ | grep 'boundary_output_.*\.txt' | wc -l)

# # Loops over all boundary_output files, sorting them by x and outputting them
# # into the pressure_outputs directory
# for filenum in $(seq 0 $(($NOFILES - 1)))
# do
#     sort -g -k1 raw_data/boundary_output_$filenum.txt > boundary_outputs/boundary_output_$filenum.txt
#     echo Sorted boundary_output_$filenum.txt
# done

# Move files
mv raw_data/boundary_output_*.txt boundary_outputs
mv raw_data/w_*.txt membrane_outputs
mv raw_data/p_*.txt membrane_outputs
mv raw_data/*.mp4 videos

# Copies the log file into the upper directory 
cp raw_data/log .


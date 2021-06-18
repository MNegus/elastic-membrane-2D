# run_simulation.sh
# Runs the Basilisk code script_name.c (remember to not include the .c when
# calling the script).
 
script_name=$1 # Calls in script name from first input

# Removes and makes a new directory with name script-name
rm -r ${script_name} 
mkdir ${script_name}

# Not sure what this does or if it's needed, but here for now
$BASILISK/qcc -MD -o ${script_name}.s.d ${script_name}.c

# Compiles script-name.c using the desired flags with qcc
qcc -O2 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -c ${script_name}.c 

# Compiles wave-equation.c using gcc
gcc -O2 -c wave-equation.c 

# Links wave-equation.o and script_name.o, using the same flags as before PLUS lapacke and lapack
qcc wave-equation.o ${script_name}.o -O2 -L$BASILISK/gl -llapacke -llapack -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -o ${script_name}/${script_name} -lm

# Move into the script_name directory
cd ${script_name}

# Set the number of threads to be the second input
export OMP_NUM_THREADS=$2

# Runs the script_name script
./${script_name}

# Moves back upwards
cd ../

# Remove any existing raw data in the parent directory
rm -r ../raw_data 

# Moves the new data to the parent directory
mv ${script_name} ../raw_data
# run_simulation.sh
# Runs the Basilisk code basilisk_script.c (remember to not include the .c when
# calling the script).
 
basilisk_script=$1 # Calls in script name of Basilisk code from first input
membrane_script=$2 # Membrane script, either wave-equation or membrane-equation

# Removes and makes a new directory with name script-name
rm -r ${basilisk_script} 
mkdir ${basilisk_script}

# Not sure what this does or if it's needed, but here for now
$BASILISK/qcc -MD -o ${basilisk_script}.s.d ${basilisk_script}.c

# Compiles basilisk_script.c using the desired flags with qcc
qcc -O2 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -c ${basilisk_script}.c 

# Compiles membrane_script.c using gcc
gcc -O2 -c ${membrane_script}.c 

# Links membrane_script.o and basilisk_script.o, using the same flags as before PLUS lapacke and lapack
qcc ${membrane_script}.o ${basilisk_script}.o -O2 -L$BASILISK/gl -llapacke -llapack -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -o ${basilisk_script}/${basilisk_script} -lm

# Move into the basilisk_script directory
cd ${basilisk_script}

# Set the number of threads to be the third input
export OMP_NUM_THREADS=$3

# Runs the basilisk_script script
./${basilisk_script}

# Moves back upwards
cd ../

# Remove any existing raw data in the parent directory
rm -r ../raw_data 

# Moves the new data to the parent directory
mv ${basilisk_script} ../raw_data
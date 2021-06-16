script_name=$1

rm -r ${script_name}
mkdir ${script_name}

/home/michael/basilisk_0/src/qcc -MD -o ${script_name}.s.d ${script_name}.c

# qcc -O2 -fopenmp -c ${script_name}.c 
qcc -O2 -L/home/michael/basilisk_0/src/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -c ${script_name}.c 

gcc -O2 -c wave-equation.c 

qcc wave-equation.o ${script_name}.o -O2 -L/home/michael/basilisk_0/src/gl -llapacke -llapack -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL -o ${script_name}/${script_name} -lm

cd ${script_name}

export OMP_NUM_THREADS=$2

./${script_name}

# Moves back upwards
cd ../

# Remove any existing raw data in the parent directory
rm -r ../raw_data 

# Moves the new data to the parent directory
mv ${script_name} ../raw_data
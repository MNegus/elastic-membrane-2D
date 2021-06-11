/* main_prescribed_pressure.c
Solves the wave equation given a set of pressure outputs from a stationary 
membrane droplet impact code. 
*/

#include "parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave-equation.h"

// Constants
double Deltax; //  Spatial grid size
int M; // Number of grid nodes
char pressure_dirname[40] = "../pressure_outputs"; // Directory for pressures

// Arrays
double *w_previous, *w, *w_next, *w_deriv;
double *p_previous, *p, *p_next;

// Variables
double t; // Time value
int k; // Timestep number
int kmax; // Maximum value of k


// Function declarations
void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr);
void read_pressure(double *p_arr, int k);


int main (int argc, const char * argv[]) {
    t = 0;
    k = 0;

    M = (int) floor(pow(2, MAXLEVEL) * MEMBRANE_RADIUS / BOX_WIDTH);
    Deltax = BOX_WIDTH / M;

    // Initialise membrane arrays
    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_next = malloc(M * sizeof(double)); // w at next timestep
    w_deriv = malloc(M * sizeof(double)); // Time derivative of w
    p_previous = malloc(M * sizeof(double)); // p at previous timestep
    p = malloc(M * sizeof(double)); // p at current timestep
    p_next = malloc(M * sizeof(double)); // p at next timestep

    /* For the first few timesteps, output a stationary membrane, waiting for 
    the Basilisk pressure output to calm down */
    for (int i = 0; i < M; i++) {
            w[i] = 0.0;
            w_deriv[i] = 0.0;
            p[i] = 0.0;
    }

    for (k = 0; k < 10; k++) {
        output_arrays(w, w_deriv, p);
        t += DELTA_T;
    }
    
    // Initialise pressures at k = 2
    t = 10 * DELTA_T;
    k = 10;
    read_pressure(p_previous, k - 1);
    read_pressure(p, k);
    read_pressure(p_next, k + 1);

    // Initialise membrane
    initialise_membrane(w_previous, w, w_deriv, p_previous, p, p_next, M + 1, DELTA_T, MEMBRANE_RADIUS, ALPHA, BETA);

    // Output w_previous
    output_arrays(w_previous, w_deriv, p_previous);

    // Update time and pressures
    t += DELTA_T;
    k++;
    double *temp = p_previous;
    p_previous = p;
    p = p_next;
    p_next = temp;
    read_pressure(p_next, k + 1);

    while (t <= MAX_TIME) {
        printf("t = %g\n", t);

        // Compute a timestep to find w_next
        membrane_timestep(w_previous, w, w_next, w_deriv, \
            p_previous, p, p_next, DELTA_T);

        // Output w, w_deriv and p
        output_arrays(w, w_deriv, p);

        // Incremement time
        t += DELTA_T;
        k++;

        // Updates w and p arrays if we are not at the end
        if (t + DELTA_T <= MAX_TIME) {
            double *temp1 = w_previous;
            w_previous = w;
            w = w_next;
            w_next = temp1;

            double *temp2 = p_previous;
            p_previous = p;
            p = p_next;
            p_next = temp2;
            read_pressure(p_next, k + 1);
        }
    }
}


void read_pressure(double *p_arr, int k) {
/* read_pressure
Reads the values of the pressure at timestep k into the given array from the 
pressure outputs directory
*/
    // Line reading values
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    // Reads the intial condition for w_previous from the initial_condition.txt file
    char pressure_filename[80];
    sprintf(pressure_filename, "%s/p_%d.txt", pressure_dirname, k);
    FILE *pressure_file = fopen(pressure_filename, "r");

    int i = 0;
    while ((read = getline(&line, &len, pressure_file)) != -1) {
        // Read x and w_val from the file
        double x, p_val;
        sscanf(line, "%lf %lf", &x, &p_val);
        
        // Saves w_val into w_previous
        if (i < M) {
            p_arr[i] = p_val;
        }
       
        // Increment i
        i++;
    }
    fclose(pressure_file);

}


void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "../membrane_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "../membrane_outputs/w_deriv_%d.txt", k);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    char p_filename[40];
    sprintf(p_filename, "../membrane_outputs/p_%d.txt", k);
    FILE *p_file = fopen(p_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[i]);
        fprintf(w_deriv_file, "%.10f, %.10f\n", x, w_deriv_arr[i]);
        fprintf(p_file, "%.10f, %.10f\n", x, p_arr[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f\n", x, 0.0);
    fprintf(p_file, "%.10f, %.10f\n", x, 0.0);
    fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    fclose(p_file);
    fclose(w_deriv_file);

}



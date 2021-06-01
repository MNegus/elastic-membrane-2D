
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave-equation.h"

#include "parameters.h"

// Globals
double Deltax;
int M;
double *w_previous, *w, *w_next, *w_deriv;
double *p_previous, *p, *p_next;

double t;
int k;


// Function declarations
void analytical_pressure(double *p_arr, double t);
void output_membrane(double *w_arr);


int main (int argc, const char * argv[]) {
    t = 0;
    k = 0;

    Deltax = L / (N_MEMBRANE - 1);
    M = N_MEMBRANE - 1;

    // Initialise membrane arrays
    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_next = malloc(M * sizeof(double)); // w at next timestep
    p_previous = malloc(M * sizeof(double));
    p = malloc(M * sizeof(double));
    p_next = malloc(M * sizeof(double));
    
    // Initialise pressures
    analytical_pressure(p_previous, -DELTA_T);
    analytical_pressure(p, 0.0);
    analytical_pressure(p_next, DELTA_T);

    // Initialise membrane
    initialise_membrane(w_previous, w, p_previous, p, p_next, N_MEMBRANE, DELTA_T, L, ALPHA, BETA);

    // Output w_previous
    output_membrane(w_previous);
    

    // Output w
    t += DELTA_T;
    k++;
    output_membrane(w);

    // Update pressures
    double *temp = p_previous;
    p_previous = p;
    p = p_next;
    p_next = temp;
    analytical_pressure(p_next, t + DELTA_T);

    // Loops over all timesteps
    while (t < T_MAX) {

        // Find w_next
        membrane_timestep(w_previous, w, w_next, p_previous, p, p_next);

        // Output w
        t += DELTA_T;
        k++;

        double *temp1 = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp1;
        output_membrane(w);

        // Update pressures
        double *temp2 = p_previous;
        p_previous = p;
        p = p_next;
        p_next = temp2;
        analytical_pressure(p_next, t + DELTA_T);

    }

}


void analytical_pressure(double *p_arr, double t) {
/* analytical_pressure
Outputs an analytical pressure solution into the length M array p_arr, which is
given by the function
    p(x, t) = sum(A_n * cos(k_n t) * cos(lambda_n x))
for given coefficients A_n, k_n and lambda_n.
*/

    double As[3] = {10, 5, 2.5};
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        double p_val = 0;
        for (int n = 1; n <= 3; n++) {
            double lambda = M_PI * (2 * n - 1) / (2 * L);
            double l = sqrt(BETA * pow(lambda, 2) + GAMMA * pow(lambda, 4)) / sqrt(ALPHA);
            double k = 10 * l;
            p_val += As[n - 1] * cos(k * t) * cos(lambda * x);
        }
        p_arr[i] = p_val;
    }
}


void output_membrane(double *w_arr) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "test_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char p_filename[40];
    sprintf(p_filename, "test_outputs/p_%d.txt", k);
    FILE *p_file = fopen(p_filename, "w");

    // char w_deriv_filename[40];
    // sprintf(w_deriv_filename, "mitchell_outputs/w_deriv_%d.txt", k);
    // FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[i]);
        fprintf(p_file, "%.10f, %.10f\n", x, p[i]);
        // fprintf(w_deriv_file, "%.10f, %.10f\n", x, q[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f\n", x, 0.0);
    fprintf(p_file, "%.10f, %.10f\n", x, 0.0);
    // fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    fclose(p_file);
    // fclose(w_deriv_file);

}



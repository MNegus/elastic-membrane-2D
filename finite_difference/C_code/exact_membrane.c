/* implicit_fd_membrane.c 
Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
with boundary conditions
    w_x = w_xxx at x = 0, w = w_xx = 0 at x = L,
using an exact solution given by
    w(x, t) = sum_{n=1}^\infinity A_n cos(l_n t) cos(lambda_n x),
with  
    lambda_n = (2n -1) pi / (2L) and 
    l_n^2 = (BETA lambda_n^2 + GAMMA lambda_n^4) / ALPHA,
and A_n given by integrating the initial condition
    A_n = (2/L) int_0^L w(x, 0) cos(lambda_n x) dx.
We consider simple initial conditions of the form 
    w(x, 0) = sum_{n=1}^{N0} A_n cos(lambda_n x),
for a finite value of N0, so that we just know A_n exactly at the start.

. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size Deltax = L / (N_MEMBRANE - 1), and
we discretise in time with a timestep of DELTA_T.

Due to the w = 0 term at x = L, we do not need to calculate the 
n = N_MEMBRANE - 1 term, and hence we are left with a matrix equation
A w^(k+1) = 2 * w^k - w^(k-1),
with A being an MxM banded matrix, with M = N_MEMBRANE - 1. To this end, we use 
the dgsbv subroutine from LAPACKE, the C-wrapper to LAPACK, to solve the banded 
matrix equation. 

Author: Michael Negus
*/

#include <stdio.h> // For text output
#include <math.h> // For trig terms
#include "parameters.h" // Parameter file

int main (int argc, const char * argv[]) {

    /* Parameters */
    int N0 = 3; // Number of terms to take in the sum
    double As[3] = {1, 0.5, 0.25}; // Coefficients
    double Deltax = L / (N_MEMBRANE - 1); // Spatial grid size

    /* Outputs initial condition */
    char w_filename[40] = "initial_condition.txt";
    FILE *w_file = fopen(w_filename, "w");
    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = i * Deltax;
        double w = 0;
        for (int n = 1; n <= N0; n++) {
            double lambda = M_PI * (2 * n - 1) / (2 * L);
            w += As[n - 1] * cos(lambda * x);
        }
        fprintf(w_file, "%.10f, %.10f\n", x, w);
    }
    fclose(w_file);

    /* Loops over all time and saves the exact solution */
    double t = 0;
    int k = 0;
    while (t < T_MAX) {
        char w_filename[40];
        sprintf(w_filename, "exact_outputs/w_%d.txt", k);
        FILE *w_file = fopen(w_filename, "w");

        for (int i = 0; i < N_MEMBRANE; i++) {
            double x = i * Deltax;
            double w = 0;

            // Loops over terms
            for (int n = 1; n <= N0; n++) {
                double lambda = M_PI * (2 * n - 1) / (2 * L);
                double l = sqrt(BETA * pow(lambda, 2) + GAMMA * pow(lambda, 4)) / sqrt(ALPHA);
                w += As[n - 1] * cos(l * t) * cos(lambda * x);
            }

            // Outputs solution
            fprintf(w_file, "%.10f, %.10f\n", x, w);
        }
        fclose(w_file);

        // Incremements t and k
        t += DELTA_T;
        k++;
    }

}
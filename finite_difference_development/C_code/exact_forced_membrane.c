/* exact_forced_membrane.c 
Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = p(x, t),
with boundary conditions
    w_x = w_xxx at x = 0, w = w_xx = 0 at x = L,
and pressure given by
    p(x, t) = sum_{n=1}^{N0} A_n cos(k_n t) cos(lambda_n x),
using an exact solution given by
    w(x, t) = sum_{n=1}^\infinity a_n(t) cos(lambda_n x),
with  
    lambda_n = (2n -1) pi / (2L) and 
    l_n^2 = (BETA lambda_n^2 + GAMMA lambda_n^4) / ALPHA,
    k_n =/= l_n,
and 
    a_n(t) = A_n / (ALPHA (k_n^2 - l_n^2)) * (cos(l_n t) - cos(k_n t)),
and A_n are given parameters, N0 is a finite number.

Author: Michael Negus
*/

#include <stdio.h> // For text output
#include <math.h> // For trig terms
#include "parameters.h" // Parameter file

int main (int argc, const char * argv[]) {

    /* Parameters */
    int N0 = 3; // Number of terms to take in the sum
    double As[3] = {10, 5, 2.5}; // Coefficients
    double Deltax = L / (N_MEMBRANE - 1); // Spatial grid size

    /* Outputs initial condition */
    char w_filename[40] = "initial_condition.txt";
    FILE *w_file = fopen(w_filename, "w");
    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = i * Deltax;
        double w = 0;
        fprintf(w_file, "%.10f, %.10f\n", x, w);
    }
    fclose(w_file);

    /* Loops over all time and saves the exact solution */
    double t = 0;
    int k = 0;
    while (t < T_MAX) {
        printf("t = %g\n", t);
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
                double k = 10 * l;
                double a = As[n - 1] * (cos(l * t) - cos(k * t)) / (ALPHA * (k * k - l * l));
                w += a * cos(lambda * x);
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
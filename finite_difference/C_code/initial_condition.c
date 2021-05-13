/* initial_condition.c
Script to output the initial condition used for the simulations, based on the 
parameters given in parameters.h. It outputs a text file initial_condition.txt,
with the first column giving the values of x and second the values of w. In 
order to be compatible with the exact solution, the initial condition is of the
form

w(x, 0) = sum_{n=1}^{N0} A_n cos(lambda_n x),

where N0 is a finite number and we pick suitable values of A_n. 

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

    /* Outputs solution */
    char w_filename[40] = "initial_condition.txt";
    FILE *w_file = fopen(w_filename, "w");
    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = i * Deltax;
        double w = 0;
        for (int n = 1; n <= N0; n++) {
            double lambda = M_PI * (2 * n - 1) / (2 * L);
            w += As[n - 1] * cos(lambda * x);
        }
        fprintf(w_file, "%.8f, %.8f\n", x, w);
    }
    fclose(w_file);

}
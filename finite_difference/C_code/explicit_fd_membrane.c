/* explicit_fd_membrane.c 
Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
with boundary conditions
    w_x = w_xxx at x = 0, w = w_xx = 0 at x = L,
using an implicit finite difference method. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size Deltax = L / (N_MEMBRANE - 1), and
we discretise in time with a timestep of DELTA_T.


Author: Michael Negus
*/

#include <stdio.h>
#include <stdlib.h> // For malloc
#include <math.h>
#include "parameters.h"

/* Global variables */
// Finite-difference parameters
double Deltax, Cbeta2, Cgamma2; 


// Array definitions
double *w_previous, *w, *w_next, *w_deriv;

// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

/* Function definitions */
void init();
void initialise_coefficient_matrix();
void initialise_membrane();
void boundary_conditions(double *w_arr);
void output_membrane();
void run();

// Function definitions



int main (int argc, const char * argv[]) {
    /* Intialise problem */
    init();

    int N_max = 1 + L * sqrt(2 * ALPHA / (DELTA_T * (BETA * DELTA_T + sqrt(16 * ALPHA * GAMMA + BETA * BETA * DELTA_T))));;
    printf("N_max = %d\n", N_max);

    /* Loops over all timesteps */
    run();

    /* Finish */
    printf("Finished with t = %g\n", t);
}


void init() {
    /* Derived parameters */
    Deltax = L / (N_MEMBRANE - 1); // May need to make it L / N
    Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax);
    Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax); 

    /* Initialise w arrays */
    w_previous = malloc(N_MEMBRANE * sizeof(double)); // w at previous timestep
    w = malloc(N_MEMBRANE * sizeof(double)); // w at current timestep
    w_next = malloc(N_MEMBRANE * sizeof(double));
    w_deriv = malloc(N_MEMBRANE * sizeof(double)); // Time derivative of w

    // Initialise w_previous and w using second-order initial conditions
    initialise_membrane();

    /* Outputs w_previous */
    output_membrane();

    // Increments times
    t += DELTA_T;
    k++;
}


void initialise_membrane() {
/* initialise_membrane 
Initialises the membrane position arrays w_previous and w using a second-order
in time scheme. w_previous is set to a known position, and then w at k = 1 is 
set by solving the membrane equation under the conditions that w_t = 0 at t = 0.
*/

    // Reads the intial condition for w from the initial_condition.txt file

    // Line reading values
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    // IC file
    FILE *initial_condition_file = fopen("initial_condition.txt", "r");

    // Loops over all lines in file
    int i = 0;
    while ((read = getline(&line, &len, initial_condition_file)) != -1) {
        // printf("%s", line);

        // Read x and w_val from the file
        double x, w_val;
        sscanf(line, "%lf, %lf", &x, &w_val);
        
        // Save w_val into w_previous
        w_previous[i] = w_val;
        
        // Initial condition for w_deriv
        w_deriv[i] = 0.;

        // Incremement i
        i++;
    }
    fclose(initial_condition_file);

    // // Initialise w_previous and w_deriv
    // for (int i = 0; i < N_MEMBRANE; i++) {
    //     double x = Deltax * i;

    //     // Prescribed initial conditions for w and its first derivative
    //     w_previous[i] = cos(M_PI * x / (2 * L));
    //     w_deriv[i] = 0;
    // }

    // Initialise w at the k = 1 timestep for the bulk nodes
    for (int i = 2; i < N_MEMBRANE - 2; i++) {
        w[i] = 0.5 * (- Cgamma2 * w_previous[i - 2] \
            + (Cbeta2 + 4. * Cgamma2) * w_previous[i - 1] \
            + 2 * (1 - Cbeta2 - 3. * Cgamma2) * w_previous[i] \
            + (Cbeta2 + 4. * Cgamma2) * w_previous[i + 1] \
            - Cgamma2 * w_previous[i + 2]);
    }
    // Applies boundary conditions
    boundary_conditions(w);
}

void boundary_conditions(double *w_arr) {
/* Applies the boundary conditions to a w array given by w_arr */
    w_arr[0] = (39. / 17.) * w_arr[2] \
        - (28. / 17.) * w_arr[3] \
        + (6. / 17.) * w_arr[4];
    w_arr[1] = (67. / 34.) * w_arr[2] \
        - (21. / 17.) * w_arr[3] \
        + (9. / 34.) * w_arr[4];
    w_arr[N_MEMBRANE - 2] = - (1. / 5.) * w_arr[N_MEMBRANE - 4] + (4. / 5.) * w_arr[N_MEMBRANE - 3];
    w_arr[N_MEMBRANE - 1] = 0;
}


void output_membrane() {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "explicit_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "explicit_outputs/w_deriv_%d.txt", k);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w_previous[i]);
        fprintf(w_deriv_file, "%.10f, %.10f\n", x, w_deriv[i]);
    }
    fclose(w_file);
    fclose(w_deriv_file);

}


void run() {
/* Loops over all time values solving the equations */
    while (t <= T_MAX) {
        // printf("t = %g\n", t);
        
        /* Determines bulk nodes */
        for (int i = 2; i < N_MEMBRANE - 2; i++) {
            w_next[i] = \
                - Cgamma2 * w[i - 2] \
                + (Cbeta2 + 4. * Cgamma2) * w[i - 1] \
                + 2 * (1 - Cbeta2 - 3. * Cgamma2) * w[i] \
                + (Cbeta2 + 4. * Cgamma2) * w[i + 1] \
                - Cgamma2 * w[i + 2] \
                - w_previous[i];
        } 

        /* Determines boundary nodes */
        boundary_conditions(w_next);
        
        /* Determines w_deriv */
        for (int i = 0; i < N_MEMBRANE; i++) {
            w_deriv[i] = (w_next[i] - w_previous[i]) / DELTA_T;
        }

        // Shifts arrays, updating w
        double *temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;

        // Outputs w_previous and w_deriv
        output_membrane();

        // Increments time
        t += DELTA_T;
        k++;
    }
}


// int main() {

//     // Determines M based on the stability condition
//     // double M_max = log2(L) \
//     //     - 0.5 * log2(Deltat * (BETA * Deltat + sqrt(16 * ALPHA * GAMMA + pow(BETA * Deltat, 2))) / (2 * ALPHA));
//     // M = (int) floor(M_max);

//     // Derived constants
//     N_membrane = N_MEMBRANE;
//     Deltax = MEMBRANE_RADIUS / (N_membrane - 1); // Grid size
//     CBETA2 = BETA * Deltat * Deltat / (ALPHA * pow(Deltax, 2));
//     CGAMMA2 = GAMMA * Deltat * Deltat / (ALPHA * pow(Deltax, 4));
//     kmax = floor(HARD_MAX_TIME / Deltat); // Determines max number of timesteps

//     // Arrays for solutions
//     w_previous = malloc(sizeof(double) * N_membrane);
//     w = malloc(sizeof(double) * N_membrane);
//     w_next = malloc(sizeof(double) * N_membrane);
//     w_t = malloc(sizeof(double) * N_membrane);

//     /* Initialisation */

//     // Set w and w_next to zero at t = 0 and output
//     char output_filename[80] = "membrane_position_0.txt";
//     FILE *output_file = fopen(output_filename, "w");

//     char pressure_filename[80] = "pressure_output_0.txt";
//     FILE *pressure_file = fopen(pressure_filename, "w");

//     char w_t_filename[80] = "membrane_derivative_0.txt";
//     FILE *w_t_file = fopen(w_t_filename, "w");

//     for (int n = 0; n < N_membrane; n++) {
//         w[n] = 0.;
//         w_next[n] = 0.;
//         w_t[n] = 0.;

//         double x = Deltax * n;
//         fprintf(output_file, "%g %g\n", x, w[n]);
//         fprintf(pressure_file, "%g %g\n", x, prescribed_pressure(x, t));
//         fprintf(w_t_file, "%g %g\n", x, -w_t[n]);
//     }
//     fclose(output_file);
//     fclose(pressure_file);
//     fclose(w_t_file);


//     // Shifts arrays
//     // w_previous = w;
//     // w = w_next;
//     // w_next = w_previous;

//     /* Loops over all timesteps */
//     for (int k = 1; k <= kmax; k++) {
//         double t = k * Deltat; // Updates time

//         // Outputs solution
//         char output_filename[80];
//         sprintf(output_filename, "membrane_position_%d.txt", k);
//         FILE *output_file = fopen(output_filename, "w");

//         char pressure_filename[80];
//         sprintf(pressure_filename, "pressure_output_%d.txt", k);
//         FILE *pressure_file = fopen(pressure_filename, "w");

//         char w_t_filename[80];
//         sprintf(w_t_filename, "membrane_derivative_%d.txt", k);
//         FILE *w_t_file = fopen(w_t_filename, "w");

//         for (int n = 0; n < N_membrane; n++) {
//             double x = Deltax * n;
//             fprintf(output_file, "%g %g\n", x, w[n]);
//             fprintf(pressure_file, "%g %g\n", x, prescribed_pressure(x, t));
//             fprintf(w_t_file, "%g %g\n", x, -w_t[n]);
//         }
//         fclose(output_file);
//         fclose(pressure_file);
//         fclose(w_t_file);

//         /* Updates numerical solution */
//         update_position(w_next, w, w_previous, w_t, t);

//         // Shifts arrays
//         double *temp = w_previous;
//         w_previous = w;
//         w = w_next;
//         w_next = temp;

//     }

//     // Frees arrays
//     free(w_previous); 
//     free(w); 

// }


// void update_position(double *w_next, double *w, double *w_previous, \
//     double *w_t, double t) {
//     // Bulk nodes
//     double x;
//     // fprintf(stderr, "Before: w_next = %g, w = %g, w_previous = %g\n", w_next[10], w[10],  w_previous[10]);
//     for (int n = 2; n < N_membrane - 2; n++) {
//         x = n * Deltax;
        
//         w_next[n] = \
//             - CGAMMA2 * w[n - 2] \
//             + (CBETA2 + 4. * CGAMMA2) * w[n - 1] \
//             + 2 * (1 - CBETA2 - 3. * CGAMMA2) * w[n] \
//             + (CBETA2 + 4. * CGAMMA2) * w[n + 1] \
//             - CGAMMA2 * w[n + 2] \
//             - w_previous[n] + (Deltat * Deltat / ALPHA) * prescribed_pressure(x, t);
//         // w_t[n] = (w_next[n] - w[n]) / dt;
//         w_t[n] = (w_next[n] - w_previous[n]) / (2 * dt);
//     } 
//     // fprintf(stderr, "After: w_next = %g, w = %g, w_previous = %g\n", w_next[10], w[10],  w_previous[10]);


//     // Boundary nodes
//     boundary_conditions(w_next);
//     w_t[0] = (w_next[0] - w_previous[0]) / (2 * dt);
//     w_t[1] = (w_next[1] - w_previous[1]) / (2 * dt);
//     w_t[N_membrane - 2] = (w_next[N_membrane - 2] - w_previous[N_membrane - 2]) / (2 * dt);
//     w_t[N_membrane - 1] = (w_next[N_membrane - 1] - w_previous[N_membrane - 1]) / (2 * dt);
// }

// void boundary_conditions(double *w_next) {
//     w_next[0] = (39. / 17.) * w_next[2] \
//         - (28. / 17.) * w_next[3] \
//         + (6. / 17.) * w_next[4];
//     w_next[1] = (67. / 34.) * w_next[2] \
//         - (21. / 17.) * w_next[3] \
//         + (9. / 34.) * w_next[4];
//     w_next[N_membrane - 2] = - (1. / 5.) * w_next[N_membrane - 4] + (4. / 5.) * w_next[N_membrane - 3];
//     w_next[N_membrane - 1] = 0;
// }

// double prescribed_pressure(double x, double t) {
//     /* Outer pressure */
//     double scaled_t = t - INITIAL_DROP_HEIGHT;
//     if (scaled_t <= 0) {
//         return 0;
//     } else if (x < 2 * sqrt(scaled_t)) {
//         return 2. / sqrt(4 * scaled_t- x * x); 
//     } else {
//         return 0;
//     }

//     /* Time dependent sum */
//     // int Np = 5;
//     // double scaled_t = t - INITIAL_DROP_HEIGHT;
//     // double pval = 0;
//     // if (scaled_t > 0) {
//     //     for (int m = 1; m <= Np; m++) {
//     //         pval += scaled_t * cos(x * (2 * m - 1) * pi / (2 * MEMBRANE_RADIUS)) / sqrt(MEMBRANE_RADIUS);
//     //     }
//     // }
//     // return pval;
// }   

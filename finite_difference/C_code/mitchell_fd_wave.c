/* mitchell_fd_wave.c 
Solves the wave equation 
    w_tt = c^2 w_xx,
with boundary conditions
    w_x = 0 at x = 0, w = 0 at x = L,
using the Mitchell method. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size Deltax = L / (N_MEMBRANE - 1), and
we discretise in time with a timestep of DELTA_T.

In general, the Mitchell method involves discretising a second-order in time PDE
    w_tt = L(x, t, w, w_x, w_xx, w_xxx, w_xxxx),
in the following way
    (w_n^(k-1) - 2 w_n^k + w_n^(k+1))/(DELTA_T^2) \
        = L_n^(k-1) / 4 + L_n^k / 2 + L_n^(k+1) / 4,
which is in theory an implicit, second-order in time and space discretisation. 
For our application of the wave equation, where L = c^2 w_xx, we are left with a
matrix equation
    A W^(k+1) = B W^(k) - A W^(k-1),
where A and B are tri-diagonal matrices, which we solve using LAPACK's dgsv 
subroutine for solving banded matrix systems.

Author: Michael Negus
*/

#include <stdio.h> // For text output
#include <lapacke.h> // For solving linear algebra
#include <math.h> // For pi
#include <string.h> // For memcpy
#include "parameters.h" // Parameter file

/* Global variables */
// Finite-difference parameters
double Deltax, Cbeta2; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *w_previous, *w, *rhs, *A_static, *A;

// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

/* Function definitions */
void init();
void initialise_coefficient_matrix();
void initialise_membrane();
void output_membrane(double *w_arr);
void run();


int main (int argc, const char * argv[]) {

    /* Intialise problem */
    init();

    /* Loops over all timesteps */
    run();

    /* Finish */
    printf("Finished with t = %g\n", t);
}


void init() {
    /* Derived parameters */
    Deltax = L / (N_MEMBRANE - 1); // Spatial grid size
    Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax);
    M = N_MEMBRANE - 1; // Size of matrix once the last row has been removed

    /* Matrix equations set-up */

    // LAPACKE constants
    info; // Output for the dgbsv LAPACKE subroutine
    ipiv = malloc(M * sizeof(int)); // Pivot array for dgbsv
    kl = 1; // Number of sub-diagonals
    ku = 1; // Number of super-diagonals
    nrhs = 1; // Number of right-hand side vectors
    ldab = M; // Leading dimension of A
    ldb = 1; // Leading dimension of the right-hand side vector
    noRows = 2 * kl + ku + 1; // Number of rows in coefficient matrix

    // Creates coefficient matrix
    A_static = malloc(M * noRows * sizeof(double));

    // Fills in the coefficient matrix
    initialise_coefficient_matrix();

    /* Initialise w arrays */
    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    rhs = malloc(M * sizeof(double)); // Right-hand-side vector
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

    // Initialise w_previous and w using second-order initial conditions
    initialise_membrane();

}


void initialise_coefficient_matrix() {
/* initialise_coefficient_matrix 
Function to fill in the Mx4 coefficient matrix A, where each row is one of the 
diagonals. Cbeta2 is the scaled term equal to:
Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax).
*/

    // Stores upper diagonal in second row
    for (int colNum = 1; colNum < M; colNum++) {
        double value;
        if (colNum == 1) {
            value = - Cbeta2 / 2.;
        } else {
            value = - Cbeta2 / 4.;
        }
        A_static[M + colNum] = value;
    }

    // Stores main diagonal in third row
    for (int colNum = 0; colNum < M; colNum++) {
        A_static[2 * M + colNum] = 1 + Cbeta2 / 2;
    }

    // Stores lower diagonal in fourth row
    for (int colNum = 0; colNum < M - 1; colNum++) {
        A_static[3 * M + colNum] = - Cbeta2 / 4;
    }
}


void initialise_membrane() {
/* initialise_membrane 
Initialises the membrane position arrays w and q using a second-order
in time scheme. w is set to a known position, and q is set to 0
*/

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));
    
    /* Initialise w_previous from file */
    // Line reading values
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    // Reads the intial condition for w_previous from the initial_condition.txt file
    FILE *initial_condition_file = fopen("initial_condition.txt", "r");

    // Loops over all lines in file
    int i = 0;
    while ((read = getline(&line, &len, initial_condition_file)) != -1) {
        // Read x and w_val from the file
        double x, w_val;
        sscanf(line, "%lf, %lf", &x, &w_val);
        
        // Saves w_val into w_previous
        if (i < M) {
            w_previous[i] = w_val;
        }
       
        // Increment i
        i++;
    }
    fclose(initial_condition_file);

    // Outputs w_previous
    output_membrane(w_previous);

    // Increments the time t and timestep k
    t += DELTA_T;
    k++;

    /* Initialise w by solving the matrix equation */
    // Configures right-hand-side vector
    rhs[0] = 0.5 * ((2 - Cbeta2) * w_previous[0] + Cbeta2 * w_previous[1]);
    for (int i = 1; i < M - 1; i++) {
        rhs[i] = 0.5 * (Cbeta2 * w_previous[i - 1] / 2 + (2 - Cbeta2) * w_previous[i] + Cbeta2 * w_previous[i + 1] / 2);
    }
    rhs[M - 1] = 0.5 * (Cbeta2 * w_previous[M - 2] / 2 + (2 - Cbeta2) * w_previous[M - 1]);

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in rhs
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

    // Swaps rhs and w
    double *temp = w; 
    w = rhs;
    rhs = w;

    // Outputs w
    output_membrane(w);
}


void output_membrane(double *w_arr) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "mitchell_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    // char w_deriv_filename[40];
    // sprintf(w_deriv_filename, "mitchell_outputs/w_deriv_%d.txt", k);
    // FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[i]);
        // fprintf(w_deriv_file, "%.10f, %.10f\n", x, q[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f", x, 0.0);
    // fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    // fclose(w_deriv_file);

}


void run() {
/* Loops over all time values solving the equations */
    while (t <= T_MAX) {
        printf("t = %g\n", t);
        
        /* Configures right-hand-side vector */
        rhs[0] = (2 - Cbeta2) * w[0] + Cbeta2 * w[1] \
            - ((1 + Cbeta2 / 2) * w_previous[0] - Cbeta2 * w_previous[1] / 2);
        
        for (int i = 1; i < M - 1; i++) {
            rhs[i] = Cbeta2 * w[i - 1] / 2 + (2 - Cbeta2) * w[i] + Cbeta2 * w[i + 1] / 2 \
                + Cbeta2 * w_previous[i - 1] / 4 - (1 + Cbeta2 / 2) * w_previous[i] +  Cbeta2 * w_previous[i + 1] / 4;
        }
        rhs[M - 1] = Cbeta2 * w[M - 2] / 2 + (2 - Cbeta2) * w[M - 1] \
                + Cbeta2 * w_previous[M - 2] / 4 - (1 + Cbeta2 / 2) * w_previous[M - 1];

        /* Solves matrix  equation */
        // Copies over elements to A
        memcpy(A, A_static, M * noRows * sizeof(double));

        // Uses LAPACK to solve the matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        /* Outputting and swapping */
        // Swaps arrays, updating w and w_previous
        double *temp = w_previous;
        w_previous = w;
        w = rhs;
        rhs = temp;

        // Increments time
        t += DELTA_T;
        k++;

        // Outputs the new value of w
        output_membrane(w);
        
    }
}
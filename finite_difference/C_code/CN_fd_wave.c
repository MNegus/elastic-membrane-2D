/* CN_fd_wave.c 
Solves the wave equation 
    w_tt = c^2 w_xx,
with boundary conditions
    w_x = 0 at x = 0, w = 0 at x = L,
using the Crank-Nicolson (CN) method. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size Deltax = L / (N_MEMBRANE - 1), and
we discretise in time with a timestep of DELTA_T.

We introduce a variable q, such that q = w_t, in order to conduct the CN method.
Due to the w = 0 term at x = L, we do not need to calculate the 
n = N_MEMBRANE - 1 term, and hence we are left with a matrix equation
    A W^(k+1) = B W^k + DELTA_T Q^k,
where A and B are MxM banded matrices, with M = N_MEMBRANE - 1, and W and Q are
length M vectors representing w and q for n = 0 to M - 1. To this end, we use 
the dgsbv subroutine from LAPACKE, the C-wrapper to LAPACK, to solve the banded 
matrix equation. We then have an explicit equation for Q^(k+1) given by
    Q^(k+1) = (2 / DELTA_T) (W^(k+1) - W^k) - Q^k.

Author: Michael Negus
*/

#include <stdio.h> // For text output
#include <lapacke.h> // For solving linear algebra
#include <math.h> // For pi
#include <string.h> // For memcpy
#include "parameters.h" // Parameter file

/* Global variables */
// Finite-difference parameters
double Deltax, Cbeta2, Cgamma2; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *w, *q, *q_next, *rhs, *A_static, *A;

// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

/* Function definitions */
void init();
void initialise_coefficient_matrix();
void initialise_membrane(); and w_deriv
void output_membrane();
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
    Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax); 
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
    w = malloc(M * sizeof(double)); // w at current timestep
    q = malloc(M * sizeof(double)); // q at current timestep
    q_next = malloc(M * sizeof(double)); // q at next timestep
    rhs = malloc(M * sizeof(double)); // Right-hand-side vector
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

    // Initialise w_previous and w using second-order initial conditions
    initialise_membrane();

    /* Outputs w_previous */
    output_membrane();

}


void initialise_coefficient_matrix() {
/* initialise_coefficient_matrix 
Function to fill in the Nx3 coefficient matrix A, where each row is one of the 
diagonals. Cbeta2 and Cgamma2 are the scaled terms equal to:
Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax),
Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax)
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
    
    /* Initialise w and q */

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
        
        // Saves w_val 
        if (i < M) {
            w[i] = w_val;
        
            // Initial condition for q
            q[i] = 0.;
        }
       
        // Incremement i
        i++;
    }
    fclose(initial_condition_file);

}


void output_membrane() {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "CN_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "CN_outputs/w_deriv_%d.txt", k);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w[i]);
        fprintf(w_deriv_file, "%.10f, %.10f\n", x, q[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f", x, 0.0);
    fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    fclose(w_deriv_file);

}


void run() {
/* Loops over all time values solving the equations */
    while (t <= T_MAX) {
        printf("t = %g\n", t);
        
        // Configures right-hand-side vector
        rhs[0] = (1 - Cbeta2 / 2) * w[0] + Cbeta2 * w[1] / 2 + DELTA_T * q[0];
        for (int i = 1; i < M; i++) {
            rhs[i] =  Cbeta2 * w[i - 1] / 4 \
                    + (1 - Cbeta2 / 2) * w[i] \
                    + Cbeta2 * w[i + 1] / 4 \
                    + DELTA_T * q[i];
        }

        // Copies over elements to A
        memcpy(A, A_static, M * noRows * sizeof(double));

        // Solves matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        // Determines q_next, using the fact that rhs = w_next
        for (int i = 0; i < M; i++) {
            q_next[i] = 2 * (rhs[i] - w[i]) / DELTA_T - q[i];
        }

        // Swaps arrays, updating w and q
        double *temp1 = w;
        w = rhs;
        rhs = temp1;

        double *temp2 = q;
        q = q_next;
        q_next = temp2;

        // Outputs w_previous and w_deriv
        output_membrane();

        // Increments time
        t += DELTA_T;
        k++;
    }
}
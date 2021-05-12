/* implicit_fd_membrane.c 
Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
using an implicit finite difference method.
This ends up becoming a matrix 
equation
A w^(k+1) = 2 * w^k - w^(k-1),
with A being a banded matrix. To this end, we use the dgsbv subroutine from
LAPACKE, the C-wrapper to LAPACK, to solve the banded matrix equation. 

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

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *w_previous, *w, *w_deriv, *rhs, *A_static, *A;

// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

/* Function definitions */
void init();
void initialise_coefficient_matrix();
void initialise_membrane();
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
    Deltax = L / (N_MEMBRANE - 1); // May need to make it L / N
    Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax);
    Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax); 

    /* Matrix equations set-up */

    // LAPACKE constants
    info; // Output for the dgbsv LAPACKE subroutine
    ipiv = malloc(N_MEMBRANE * sizeof(int)); // Pivot array for dgbsv
    kl = 2; // Number of sub-diagonals
    ku = 2; // Number of super-diagonals
    nrhs = 1; // Number of right-hand side vectors
    ldab = N_MEMBRANE; // Leading dimension of A
    ldb = 1; // Leading dimension of the right-hand side vector
    noRows = 2 * kl + ku + 1; // Number of rows in coefficient matrix

    // Creates coefficient matrix
    A_static = malloc(N_MEMBRANE * noRows * sizeof(double));

    // Fills in the coefficient matrix
    initialise_coefficient_matrix();

    /* Initialise w arrays */
    w_previous = malloc(N_MEMBRANE * sizeof(double)); // w at previous timestep
    w = malloc(N_MEMBRANE * sizeof(double)); // w at current timestep
    w_deriv = malloc(N_MEMBRANE * sizeof(double)); // Time derivative of w
    rhs = malloc(N_MEMBRANE * sizeof(double)); // Right-hand-side vector
    A = malloc(N_MEMBRANE * noRows * sizeof(double)); // Coefficient matrix

    // Initialise w_previous and w using second-order initial conditions
    initialise_membrane();

    /* Outputs w_previous */
    output_membrane();

    // Increments times
    t += DELTA_T;
    k++;
}


void initialise_coefficient_matrix() {
/* initialise_coefficient_matrix 
Function to fill in the Nx5 coefficient matrix A, where each row is one of the 
diagonals. Cbeta2 and Cgamma2 are the scaled terms equal to:
Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax),
Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax)
*/


    // Stores upper-upper-diagonal in third row
    for (int colNum = 2; colNum < N_MEMBRANE; colNum++) {
        double value;
        if (colNum == 2) {
            value = 2 * Cgamma2;
        } else {
            value = Cgamma2;
        }
        A_static[2 * N_MEMBRANE + colNum] = value;
    }

    // Stores upper diagonal in fourth row
    for (int colNum = 1; colNum < N_MEMBRANE; colNum++) {
        double value;
        if (colNum == 1) {
            value = - 2 * (Cbeta2 + 4 * Cgamma2);
        } else {
            value = -(Cbeta2 + 4 * Cgamma2);
        }
        A_static[3 * N_MEMBRANE + colNum] = value;
    }

    // Stores main diagonal in fifth row
    for (int colNum = 0; colNum < N_MEMBRANE; colNum++) {
        double value;
        if (colNum == 1) {
            value = 1 + 2 * Cbeta2 + 7 * Cgamma2;
        } else {
            value = 1 + 2 * Cbeta2 + 6 * Cgamma2;
        }
        A_static[4 * N_MEMBRANE + colNum] = value;
    }

    // Stores lower diagonal in sixth row
    for (int colNum = 0; colNum < N_MEMBRANE - 1; colNum++) {
        double value;
        if (colNum == N_MEMBRANE - 2) {
            value = -(Cbeta2 + 3 * Cgamma2);
        } else {
            value = -(Cbeta2 + 4 * Cgamma2);
        }
        A_static[5 * N_MEMBRANE + colNum] = value;
    }

    // Stores lower-lower diagonal in seventh row
    for (int colNum = 0; colNum < N_MEMBRANE - 2; colNum++) {
        A_static[6 * N_MEMBRANE + colNum] = Cgamma2;
    }

}


void initialise_membrane() {
/* initialise_membrane 
Initialises the membrane position arrays w_previous and w using a second-order
in time scheme. w_previous is set to a known position, and then w at k = 1 is 
set by solving the membrane equation under the conditions that w_t = 0 at t = 0.
*/

    // Copies over elements to A
    memcpy(A, A_static, N_MEMBRANE * noRows * sizeof(double));
    
    // Initialise w_previous and w_deriv
    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = Deltax * i;
        w_previous[i] = cos(M_PI * x / L) + 1;
        w_deriv[i] = 0;
    }

    // Adjusts the main diagonal of A (i.e. mapping A -> A + I)
    for (int colNum = 0; colNum < N_MEMBRANE; colNum++) {
        A[4 * N_MEMBRANE + colNum] += 1;
    }

    // Configures right-hand-side vector
    for (int i = 0; i < N_MEMBRANE; i++) {
        rhs[i] = 2 * w_previous[i];
    }

    // Solves matrix equation, saving result in rhs
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, N_MEMBRANE, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);
    
    // Swaps arrays, updating w
    double *temp = w;
    w = rhs;
    rhs = temp;
}


void output_membrane() {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "outputs/w_deriv_%d.txt", k);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    for (int i = 0; i < N_MEMBRANE; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%g, %g\n", x, w_previous[i]);
        fprintf(w_deriv_file, "%g, %g\n", x, w_deriv[i]);
    }
    fclose(w_file);
    fclose(w_deriv_file);

}


void run() {
/* Loops over all time values solving the equations */
    while (t < T_MAX) {
        printf("t = %g\n", t);
        
        // Configures right-hand-side vector
        for (int i = 0; i < N_MEMBRANE; i++) {
            rhs[i] = 2 * w[i] - w_previous[i];
        }

        // Copies over elements to A
        memcpy(A, A_static, N_MEMBRANE * noRows * sizeof(double));

        // Solves matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, N_MEMBRANE, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        // Determines w_deriv at current timestep 
        for (int i = 0; i < N_MEMBRANE; i++) {
            w_deriv[i] = (rhs[i] - w_previous[i]) / (2 * DELTA_T);
        }

        // Outputs w and w_deriv
        output_membrane();

        // Shifts arrays, updating w
        double *temp = w_previous;
        w_previous = w;
        w = rhs;
        rhs = temp;

        // Increments time
        t += DELTA_T;
        k++;
    }
}
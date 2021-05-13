/* implicit_fd_membrane.c 
Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
with boundary conditions
    w_x = w_xxx at x = 0, w = w_xx = 0 at x = L,
using an implicit finite difference method. We discretise the spatial domain 
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
    Deltax = L / (N_MEMBRANE - 1); // Spatial grid size
    Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax);
    Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax * Deltax * Deltax); 
    M = N_MEMBRANE - 1; // Size of matrix once the last row has been removed

    /* Matrix equations set-up */

    // LAPACKE constants
    info; // Output for the dgbsv LAPACKE subroutine
    ipiv = malloc(M * sizeof(int)); // Pivot array for dgbsv
    kl = 2; // Number of sub-diagonals
    ku = 2; // Number of super-diagonals
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
    w_deriv = malloc(M * sizeof(double)); // Time derivative of w
    rhs = malloc(M * sizeof(double)); // Right-hand-side vector
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

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
    for (int colNum = 2; colNum < M; colNum++) {
        double value;
        if (colNum == 2) {
            value = 2 * Cgamma2;
        } else {
            value = Cgamma2;
        }
        A_static[2 * M + colNum] = value;
    }

    // Stores upper diagonal in fourth row
    for (int colNum = 1; colNum < M; colNum++) {
        double value;
        if (colNum == 1) {
            value = - 2 * (Cbeta2 + 4 * Cgamma2);
        } else {
            value = -(Cbeta2 + 4 * Cgamma2);
        }
        A_static[3 * M + colNum] = value;
    }

    // Stores main diagonal in fifth row
    for (int colNum = 0; colNum < M; colNum++) {
        double value;
        if (colNum == 1) {
            value = 1 + 2 * Cbeta2 + 7 * Cgamma2;
        } else if (colNum == M - 1) {
            value = 1 + 2 * Cbeta2 + 5 * Cgamma2;
        } else {
            value = 1 + 2 * Cbeta2 + 6 * Cgamma2;
        }
        A_static[4 * M + colNum] = value;
    }

    // Stores lower diagonal in sixth row
    for (int colNum = 0; colNum < M - 1; colNum++) {
        double value;
        // if (colNum == M - 2) {
        //     value = -(Cbeta2 + 5 * Cgamma2);
        // } else {
        //     value = -(Cbeta2 + 4 * Cgamma2);
        // }
        value = -(Cbeta2 + 4 * Cgamma2);
        A_static[5 * M + colNum] = value;
    }

    // Stores lower-lower diagonal in seventh row
    for (int colNum = 0; colNum < M - 2; colNum++) {
        A_static[6 * M + colNum] = Cgamma2;
    }

}


void initialise_membrane() {
/* initialise_membrane 
Initialises the membrane position arrays w_previous and w using a second-order
in time scheme. w_previous is set to a known position, and then w at k = 1 is 
set by solving the membrane equation under the conditions that w_t = 0 at t = 0.
*/

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));
    
    
    /* Initialise w_previous and w_deriv */

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
        if (i < M) {
            w_previous[i] = w_val;
        
            // Initial condition for w_deriv
            w_deriv[i] = 0.;
        }
       
        // Incremement i
        i++;
    }
    fclose(initial_condition_file);

    // for (int i = 0; i < M; i++) {
    //     double x = Deltax * i;
    //     w_previous[i] = cos(M_PI * x / (2 * L));
    //     w_deriv[i] = 0;
    // }

    // Adjusts the main diagonal of A (i.e. mapping A -> A + I)
    for (int colNum = 0; colNum < M; colNum++) {
        A[4 * M + colNum] += 1;
    }

    // Configures right-hand-side vector
    for (int i = 0; i < M; i++) {
        rhs[i] = 2 * w_previous[i];
    }

    // Solves matrix equation, saving result in rhs
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);
    
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
    sprintf(w_filename, "implicit_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "implicit_outputs/w_deriv_%d.txt", k);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.8f, %.8f\n", x, w_previous[i]);
        fprintf(w_deriv_file, "%.8f, %.8f\n", x, w_deriv[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.8f, %.8f", x, 0.0);
    fprintf(w_deriv_file, "%.8f, %.8f", x, 0.0);

    fclose(w_file);
    fclose(w_deriv_file);

}


void run() {
/* Loops over all time values solving the equations */
    while (t <= T_MAX) {
        // printf("t = %g\n", t);
        
        // Configures right-hand-side vector
        for (int i = 0; i < M; i++) {
            rhs[i] = 2 * w[i] - w_previous[i];
        }

        // Copies over elements to A
        memcpy(A, A_static, M * noRows * sizeof(double));

        // Solves matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        // Determines w_deriv at current timestep 
        for (int i = 0; i < M; i++) {
            w_deriv[i] = (rhs[i] - w_previous[i]) / (2 * DELTA_T);
        }

        // Shifts arrays, updating w
        double *temp = w_previous;
        w_previous = w;
        w = rhs;
        rhs = temp;

        // Outputs w_previous and w_deriv
        output_membrane();

        // Increments time
        t += DELTA_T;
        k++;
    }
}
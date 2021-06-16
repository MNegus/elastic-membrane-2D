/* wave-equation.c
Solves the wave equation 
    ALPHA w_tt - BETA w_xx = p(x, t),
with boundary conditions
    w_x = 0 at x = 0, w = 0 at x = L,
using the Mitchell method. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size DELTA_X = L / (N_MEMBRANE - 1), and
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
where A and B are tri-diagonal matrices, which we solve using LAPACK's dgbsv
subroutine for solving banded matrix systems.

Author: Michael Negus
*/

#include "wave-equation.h"
#include <stdio.h> // For text output
#include <stdlib.h>
#include <lapacke.h> // For solving linear algebra
#include <string.h> // For memcpy

/* Global variables */
// Finite-difference parameters
double DELTA_X, Dbeta2, Dpressure; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *A_static, *B_static, *A;


void initialise_membrane(double *w_previous, double *w, double *w_deriv, 
    double *p_previous, double *p, double *p_next, int N_MEMBRANE, 
    double DELTA_T, double L, double ALPHA, double BETA) {
/* initialise_membrane
Function initialises the problem, setting w_previous and w, and requiring the 
relevant parameters and pressure arrays as input. 
*/

    /* Derived parameters */
    DELTA_X = L / (N_MEMBRANE - 1); // Spatial grid size
    Dbeta2 = (ALPHA * DELTA_X * DELTA_X) / (BETA * DELTA_T * DELTA_T);
    Dpressure = DELTA_X * DELTA_X / BETA; // Scaled term in front of pressure
    M = N_MEMBRANE - 1; // Size of matrix once the last row has been removed

    /* LAPACKE constants */
    info; // Output for the dgbsv LAPACKE subroutine
    ipiv = malloc(M * sizeof(int)); // Pivot array for dgbsv
    kl = 1; // Number of sub-diagonals
    ku = 1; // Number of super-diagonals
    nrhs = 1; // Number of right-hand side vectors
    ldab = M; // Leading dimension of A
    ldb = 1; // Leading dimension of the right-hand side vector
    noRows = 2 * kl + ku + 1; // Number of rows in coefficient matrix

    /* Initialise coefficient matrices */
    A_static = malloc(M * noRows * sizeof(double)); // Used for dgbsv
    B_static = malloc(M * noRows * sizeof(double)); // Used for matrix multiplication

    // Set all entries to zero to start
    for (int k = 0; k < M * noRows; k++) {
        A_static[k] = 0;
        B_static[k] = 0;
    }

    // Stores upper diagonal in second row
    int colNum = 1;
    A_static[M + colNum] = -0.5;
    B_static[M + colNum] = 1;
    for (colNum = 2; colNum < M; colNum++) {
        A_static[M + colNum] = -0.25;
        B_static[M + colNum] = 0.5;
    }

    // Stores main diagonal in third row
    for (colNum = 0; colNum < M; colNum++) {
        A_static[2 * M + colNum] = Dbeta2 + 0.5;
        B_static[2 * M + colNum] = 2 * Dbeta2 - 1;
    }

    // Stores lower diagonal in fourth row
    for (colNum = 0; colNum < M - 1; colNum++) {
        A_static[3 * M + colNum] = -0.25;
        B_static[3 * M + colNum] = 0.5;
    }

    /* Initialise arrays */
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

    /* Sets w_previous and w_deriv to 0 everywhere */
    for (int k = 0; k < M; k++) {
        w_previous[k] = 0.0;
        w_deriv[k] = 0.0;
    }

    /* Determines w at the first timestep using a second-order accurate initial 
    condition */
    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Sets rhs = w to be equal to 0.5 * B * w_previous
    multiply_matrix(w, B_static, w_previous, 0.5, 0);

    // Adds pressure terms onto rhs = w
    for (int k = 0; k < M; k++) {
        w[k] += 0.5 * Dpressure \
            * (0.25 * p_previous[k] + 0.5 * p[k] + 0.25 * p_next[k]);
    }

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in w
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, w, ldb);

    // Exits if dgbsv fails
    if (info) {
        fprintf(stderr, "LAPACKE_dgbsv failed with info = %d\n", info);
        exit(1);
    }
}


void membrane_timestep(double *w_previous, double *w, double *w_next, 
    double *w_deriv, double *p_previous, double *p, double *p_next, 
    double DELTA_T) {
/* membrane_timestep 
Function performs once timestep of the Mitchell method algorithm, determining 
the value of w_next. 
*/

    /* Configures right-hand-side vector */
    // Sets w_next = B * w
    multiply_matrix(w_next, B_static, w, 1, 0);

    // Sets w_next = w_next - A * w_previous
    multiply_matrix(w_next, A_static, w_previous, -1, 1);

    // Sets w_next = w_next + pressure term
    for (int k = 0; k < M; k++) {
        w_next[k] += Dpressure * 0.25 * (p_previous[k] + 2 * p[k] + p_next[k]);
    }

    /* Solves matrix equation */
    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in w_next
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, w_next, ldb);

    // Exits if dgbsv fails
    if (info) {
        fprintf(stderr, "LAPACKE_dgbsv failed with info = %d\n", info);
        exit(1);
    }
} 


void multiply_matrix(double *y_arr, double *matrix_arr, double *x_arr, \
    double scale, int ADD) {
/* multiply_matrix
Function to compute the result of the matrix muliplication 
    y_arr = ADD * y_arr + scale * matrix_arr * x_arr, 
where scale is a real scaling factor, y_arr and x_arr are length M arrays and 
matrix_arr is a coefficient matrix in banded format, equal to either A_static or
B_static. If ADD = 0, then this sets y_arr = scale * matrix_arr * x_arr, whereas if
ADD = 1, then this adds scale * matrix_arr * x_arr to y
*/
    // First entry, ignoring the lower diagonal
    y_arr[0] = ADD * y_arr[0] \
        + scale * (matrix_arr[2 * M] * x_arr[0] \
                 + matrix_arr[M + 1] * x_arr[1]);

    // Main entries, with all diagonals
    for (int k = 1; k < M - 1; k++) {
        y_arr[k] = ADD * y_arr[k] \
            + scale * (matrix_arr[3 * M + k - 1] * x_arr[k - 1]\
                     + matrix_arr[2 * M + k] * x_arr[k] \
                     + matrix_arr[M + k + 1] * x_arr[k + 1]);
    }

    // Last entry, ignoring the upper diagonal
    y_arr[M - 1] = ADD * y_arr[M - 1] \
        + scale * (matrix_arr[3 * M + M - 2] * x_arr[M - 2] \
                 + matrix_arr[2 * M + M - 1] * x_arr[M - 1]);
}

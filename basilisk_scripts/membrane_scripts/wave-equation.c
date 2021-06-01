/* wave-equation.c
Solves the wave equation 
    ALPHA w_tt - BETA w_xx = p(x, t),
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
where A and B are tri-diagonal matrices, which we solve using LAPACK's dgbsv
subroutine for solving banded matrix systems.

Author: Michael Negus
*/

#include <stdio.h> // For text output
#include <stdlib.h>
#include <lapacke.h> // For solving linear algebra
#include <string.h> // For memcpy

/* Global variables */
// Finite-difference parameters
double Deltax, Dbeta2, Dpressure; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *A_static, *B_static, *A;


void initialise_membrane(double *w_previous, double *w, double *p_previous, \
    double *p, double *p_next, int N_MEMBRANE, double DELTA_T, double L, \
    double ALPHA, double BETA) {
/* initialise_membrane
Function initialises the problem, setting w_previous and w, and requiring the 
relevant parameters and pressure arrays as input. 
*/

    /* Derived parameters */
    Deltax = L / (N_MEMBRANE - 1); // Spatial grid size
    Dbeta2 = (ALPHA * Deltax * Deltax) / (BETA * DELTA_T * DELTA_T);
    Dpressure = Deltax * Deltax / BETA; // Scaled term in front of pressure
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

    /* Sets w_previous to 0 everywhere */
    for (int i = 0; i < M; i++) {
        w_previous[i] = 0.0;
    }

    /* Determines w at the first timestep using a second-order accurate initial 
    condition */
    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Sets rhs = w to be equal to 0.5 * B * w_previous
    matrix_multiply(w, B_static, w_previous, 0.5, 0);

    // Adds pressure terms onto rhs = w
    for (int i = 0; i < M; i++) {
        w[i] += 0.5 * Dpressure \
            * (0.25 * p_previous[i] + 0.5 * p[i] + 0.25 * p_next[i]);
    }

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in w
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, w, ldb);

}


void membrane_timestep(double *w_previous, double *w, double *w_next, \
    double *p_previous, double *p, double *p_next) {
/* membrane_timestep 
Function performs once timestep of the Mitchell method algorithm, determining 
the value of w_next. 
*/

    /* Configures right-hand-side vector */
    // Sets w_next = B * w
    matrix_multiply(w_next, B_static, w, 1, 0);

    // Sets w_next = w_next - A * w_previous
    matrix_multiply(w_next, A_static, w_previous, -1, 1);

    // Sets w_next = w_next + pressure term
    for (int i = 0; i < M; i++) {
        w_next[i] += Dpressure * (0.25 * p_previous[i] + 0.5 * p[i] + 0.25 * p_next[i]);
    }


    /* Solves matrix equation */
    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in w_next
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, w_next, ldb);

} 


void matrix_multiply(double *y_arr, double *matrix_arr, double *x_arr, \
    double scale, int ADD) {
/* matrix_multiply
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
    for (int i = 1; i < M - 1; i++) {
        y_arr[i] = ADD * y_arr[i] \
            + scale * (matrix_arr[3 * M + i - 1] * x_arr[i - 1]\
                     + matrix_arr[2 * M + i] * x_arr[i] \
                     + matrix_arr[M + i + 1] * x_arr[i + 1]);
    }

    // Last entry, ignoring the upper diagonal
    y_arr[M - 1] = ADD * y_arr[M - 1] \
        + scale * (matrix_arr[3 * M + M - 2] * x_arr[M - 2] \
                 + matrix_arr[2 * M + M - 1] * x_arr[M - 1]);
}

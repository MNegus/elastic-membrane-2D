/* membrane-equation.c
Solves the membrane equation 
    ALPHA w_tt - BETA w_xx + GAMMA w_xxxx = p(x, t),
with boundary conditions
    w_x, w_xxx = 0 at x = 0, w = = w_xx = 0 at x = L,
using the Mitchell method. We discretise the spatial domain 
with N_MEMBRANE points, with a grid size DELTA_X = L / (N_MEMBRANE - 1), and
we discretise in time with a timestep of DELTA_T.

In general, the Mitchell method involves discretising a second-order in time PDE
    w_tt = L(x, t, w, w_x, w_xx, w_xxx, w_xxxx),
in the following way
    (w_n^(k-1) - 2 w_n^k + w_n^(k+1))/(DELTA_T^2) \
        = L_n^(k-1) / 4 + L_n^k / 2 + L_n^(k+1) / 4,
which is in theory an implicit, second-order in time and space discretisation. 
For our application of the wave equation, where 
    L = (1 / ALPHA) * (BETA w_xx - GAMMA w_xxxx + p(x, t)), 
and we are left with a matrix equation
    A W^(k+1) = B W^(k) - A W^(k-1),
where A and B are banded matrices with 5 diagonals, which we solve using 
LAPACK's dgbsv subroutine for solving banded matrix systems.

Author: Michael Negus
*/

#include "membrane-equation.h"
#include <stdio.h> // For text output
#include <stdlib.h>
#include <lapacke.h> // For solving linear algebra
#include <string.h> // For memcpy
#include <math.h> // For pow, fmin and fmax

/* Global variables */
// Finite-difference parameters
double DELTA_X, Cpressure; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *A_static, *B_static, *A;


void initialise_membrane(double *w_previous, double *w, double *w_deriv, 
    double *p_previous, double *p, double *p_next, int N_MEMBRANE, 
    double DELTA_T, double L, double ALPHA, double BETA, double GAMMA) {
/* initialise_membrane
Function initialises the problem, setting w_previous and w, and requiring the 
relevant parameters and pressure arrays as input. 
*/

    // Checks if GAMMA == 0, in which case use wave-equation.c
    // if (GAMMA == 0) {
    //     fprintf(stderr, "Error: GAMMA == 0. Use wave-equation.c instead\n");
    //     exit(1);
    // }

    /* Derived parameters */
    DELTA_X = L / (N_MEMBRANE - 1); // Spatial grid size
    M = N_MEMBRANE - 1; // Size of matrix once the last row has been removed
    if (GAMMA == 0) {
        Cpressure = DELTA_X * DELTA_X / BETA; // Scaled term in front of pressure;
    } else {
        Cpressure = 4 * pow(DELTA_X, 4) / GAMMA; // Scaled term in front of pressure
    }
    
    /* LAPACKE constants */
    info; // Output for the dgbsv LAPACKE subroutine
    ipiv = malloc(M * sizeof(int)); // Pivot array for dgbsv
    nrhs = 1; // Number of right-hand side vectors
    ldab = M; // Leading dimension of A
    ldb = 1; // Leading dimension of the right-hand side vector
    if (GAMMA == 0) {
        kl = 1; // Number of sub-diagonals
        ku = 1; // Number of super-diagonals        
    } else {
        kl = 2; // Number of sub-diagonals
        ku = 2; // Number of super-diagonals
    }
    noRows = 2 * kl + ku + 1; // Number of rows in coefficient matrix

    /* Allocate memory for coefficient matrices */
    A_static = malloc(M * noRows * sizeof(double)); // Used for dgbsv
    B_static = malloc(M * noRows * sizeof(double)); // Used for matrix multiplication

    /* Initialise coefficient matrices depending on value of GAMMA */
    initialise_coefficient_matrices(DELTA_T, ALPHA, BETA, GAMMA);

    // for (int k = 0; k < M * noRows; k++) {
    //     fprintf(stderr, "A_static[%d] = %g\n", k, A_static[k]);
    //     fprintf(stderr, "B_static[%d] = %g\n", k, B_static[k]);
    // }

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
        w[k] += 0.5 * Cpressure \
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
        w_next[k] += Cpressure * 0.25 * (p_previous[k] + 2 * p[k] + p_next[k]);
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

    /* Determines w_deriv */
    for (int k = 0; k < M; k++) {
        w_deriv[k] = (w_next[k] - w_previous[k]) / (2 * DELTA_T);
    }
} 


void initialise_coefficient_matrices(double DELTA_T, \
    double ALPHA, double BETA, double GAMMA) {
/* initialise_coefficient_matrices
Function to initialise the coefficient matrices A_static and B_static, which 
are used when solving the membrane equation. If GAMMA == 0, then these will be
tri-diagonal matrices to solve the wave equation, else, they will be banded 
matrices with 2 sub- and super-diagonals to solve the full membrane equation. 
*/

    // Set all entries to zero to start
    for (int k = 0; k < M * noRows; k++) {
        A_static[k] = 0;
        B_static[k] = 0;
    }

    if (GAMMA == 0) {
        /* GAMMA == 0 case, with 1 sub- and super-diagonal */

        // Beta-dependent scaling parameters
        double Dbeta2 = (ALPHA * DELTA_X * DELTA_X) / (BETA * DELTA_T * DELTA_T);
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
    } else {
        /* GAMMA > 0 case, with 2 sub- and super-diagonals */
        // Gamma-dependent scaling parameters
        double Calpha = ALPHA * pow(DELTA_X, 4) / (GAMMA * pow(DELTA_T, 2));
        double Cbeta = BETA * pow(DELTA_X, 2) / GAMMA;

        // Stores upper-upper-diagonal in third row
        for (int colNum = 2; colNum < M; colNum++) {
            double A_value, B_value;
            if (colNum == 2) {
                A_value = 2;
                B_value = -4;
            } else {
                A_value = 1;
                B_value = -2;
            }
            A_static[2 * M + colNum] = A_value;
            B_static[2 * M + colNum] = B_value;
        }

        // Stores upper-diagonal in fourth row
        for (int colNum = 1; colNum < M; colNum++) {
            double A_value, B_value;
            if (colNum == 1) {
                A_value = -2 * Cbeta - 8;
                B_value = 4 * Cbeta + 16;
            } else {
                A_value = -Cbeta - 4;
                B_value = 2 * Cbeta + 8;
            }
            A_static[3 * M + colNum] = A_value;
            B_static[3 * M + colNum] = B_value;
        }

        // Stores main diagonal in fifth row
        for (int colNum = 0; colNum < M; colNum++) {
            double A_value, B_value;
            if (colNum == 1) {
                A_value = 4 * Calpha + 2 * Cbeta + 7;
                B_value = 8 * Calpha - 4 * Cbeta - 14; 
            } else if (colNum == M - 1) {
                A_value = 4 * Calpha + 2 * Cbeta + 5;
                B_value = 8 * Calpha - 4 * Cbeta - 10; 
            } else {
                A_value = 4 * Calpha + 2 * Cbeta + 6;
                B_value = 8 * Calpha - 4 * Cbeta - 12; 
            }
            A_static[4 * M + colNum] = A_value;
            B_static[4 * M + colNum] = B_value;
        }

        // Stores lower-diagonal in sixth row
        for (int colNum = 0; colNum < M - 1; colNum++) {
            A_static[5 * M + colNum] = - Cbeta - 4;
            B_static[5 * M + colNum] = 2 * Cbeta + 8;
        }

        // Stores lower-lower diagonal in seventh row
        for (int colNum = 0; colNum < M - 2; colNum++) {
            A_static[6 * M + colNum] = 1;
            B_static[6 * M + colNum] = -2;
        }        
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
    // Loops over all values in y_arr
    for (int k = 0; k < M; k++) {
        double y_val = 0;

        // Loop over the appropriate diagonals
        for (int j = fmax(k - kl, 0); j <= fmin(k + ku, M - 1); j++) {
            y_val += matrix_arr[(2 * kl + k - j) * M + j] * x_arr[j];
        }
        y_arr[k] = ADD * y_arr[k] + scale * y_val;
    }
}

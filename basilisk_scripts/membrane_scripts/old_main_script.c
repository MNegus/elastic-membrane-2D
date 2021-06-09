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
where A and B are tri-diagonal matrices, which we solve using LAPACK's dgbsv
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
double Deltax, Dbeta2; 
int M; // Size of matrix 

// LAPACK parameters
int info, *ipiv, kl, ku, nrhs, ldab, ldb, noRows;

// Array definitions
double *w_previous, *w, *w_next, *rhs, *A_static, *B_static, *A;
double *p_previous, *p, *p_next;
// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

/* Function definitions */
void analytical_pressure(double *p_arr, double t);
void init(double *w_previous, double *w, double *p_previous, double *p, double *p_next);
void matrix_multiply(double *y_arr, double *matrix_arr, double *x_arr, \
    double scale, int ADD);
void output_arrays(double *w_arr);
void membrane_timestep(double *w_previous, double *w, double *w_next, \
    double *p_previous, double *p, double *p_next);


int main (int argc, const char * argv[]) {

    /* Derived parameters */
    Deltax = L / (N_MEMBRANE - 1); // Spatial grid size
    Dbeta2 = (ALPHA * Deltax * Deltax) / (BETA * DELTA_T * DELTA_T);
    M = N_MEMBRANE - 1; // Size of matrix once the last row has been removed

    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_next = malloc(M * sizeof(double));

    p_previous = malloc(M * sizeof(double));
    p = malloc(M * sizeof(double));
    p_next = malloc(M * sizeof(double));
    analytical_pressure(p_previous, -DELTA_T);
    analytical_pressure(p, 0.0);
    analytical_pressure(p_next, DELTA_T);

    init(w_previous, w, p_previous, p, p_next);

    output_arrays(w_previous);

    k++;
    t += DELTA_T;

    output_arrays(w);

    // Swaps pressures
    double *temp = p_previous;
    p_previous = p;
    p = p_next;
    p_next = temp;
    analytical_pressure(p_next, t + DELTA_T);
    

    /* Loops over all timesteps */
    while (t < T_MAX) {
        membrane_timestep(w_previous, w, w_next, p_previous, p, p_next);
        
        /* Outputting and swapping */
        // Swaps arrays, updating w and w_previous
        for (int i = 0; i < M; i++) {
            printf("w_previous[%d] = %g, w[%d] = %g, w_next[%d] = %g\n", i, w_previous[i], i, w[i], i, w_next[i]);
        }
        double *temp1 = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp1;
        for (int i = 0; i < M; i++) {
            printf("w_previous[%d] = %g, w[%d] = %g, w_next[%d] = %g\n", i, w_previous[i], i, w[i], i, w_next[i]);
        }

        // Increments time
        t += DELTA_T;
        k++;

        // Swaps around pressures
        double *temp2 = p_previous;
        p_previous = p;
        p = p_next;
        p_next = temp2;
        analytical_pressure(p_next, t + DELTA_T);

        // Outputs the new value of w
        output_arrays(w);
    }
    

    /* Finish */
    printf("Finished with t = %g\n", t);
}


void analytical_pressure(double *p_arr, double t) {
/* analytical_pressure
Outputs an analytical pressure solution into the length M array p_arr, which is
given by the function
    p(x, t) = 2 / sqrt(d(t)^2 - x^2) for x < d(t), 
            = 0, otherwise
where
    d(t) = 2 * sqrt(t) for t > 0.
*/
    double As[3] = {10, 5, 2.5};
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        double p_val = 0;
        for (int n = 1; n <= 3; n++) {
            double lambda = M_PI * (2 * n - 1) / (2 * L);
            double l = sqrt(BETA * pow(lambda, 2) + GAMMA * pow(lambda, 4)) / sqrt(ALPHA);
            double k = 10 * l;
            p_val += As[n - 1] * cos(k * t) * cos(lambda * x);
        }
        p_arr[i] = p_val;
    }
}


void init(double *w_previous_arr, double *w_arr, double *p_previous_arr, double *p_arr, double *p_next_arr) {
    

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

    // Creates coefficient matrices
    A_static = malloc(M * noRows * sizeof(double)); // Used for dgbsv
    B_static = malloc(M * noRows * sizeof(double)); // Used for matrix multiplication

    /* Fills in A_static, used for dgbsv */
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
    
    rhs = malloc(M * sizeof(double)); // Right-hand-side vector
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

    /* Initialise w_previous from file */
    for (int i = 0; i < M; i++) {
        w_previous_arr[i] = 0.0;
        w_arr[i] = 0.0;
    }

    // /* Initialise w by solving the matrix equation */
    // // Copies over elements to A
    // memcpy(A, A_static, M * noRows * sizeof(double));

    // // Sets rhs to be equal to 0.5 * B * w_previous
    // matrix_multiply(rhs, B_static, w_previous, 0.5, 0);

    // // Adds pressure terms onto rhs
    // for (int i = 0; i < M; i++) {
    //     rhs[i] += 0.5 * (Deltax * Deltax / BETA) * (0.25 * p_previous[i] + 0.5 * p[i] + 0.25 * p_next[i]);
    // }

    // // Copies over elements to A
    // memcpy(A, A_static, M * noRows * sizeof(double));

    // // Uses LAPACK to solve the matrix equation, saving result in rhs
    // info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

    // // Swaps rhs and w
    // double *temp1 = w; 
    // w = rhs;
    // rhs = temp1;

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


void output_arrays(double *w_arr) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "test_outputs/w_%d.txt", k);
    FILE *w_file = fopen(w_filename, "w");

    // char p_filename[40];
    // sprintf(p_filename, "test_outputs/p_%d.txt", k);
    // FILE *p_file = fopen(p_filename, "w");

    // char w_deriv_filename[40];
    // sprintf(w_deriv_filename, "mitchell_outputs/w_deriv_%d.txt", k);
    // FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        // printf("x = %g, w = %g\n", x, w_arr[i]);
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[i]);
        // fprintf(p_file, "%.10f, %.10f\n", x, p[i]);
        // fprintf(w_deriv_file, "%.10f, %.10f\n", x, q[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f\n", x, 0.0);
    // fprintf(p_file, "%.10f, %.10f\n", x, 0.0);
    // fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    // fclose(p_file);
    // fclose(w_deriv_file);

}


void membrane_timestep(double *w_previous_arr, double *w_arr, double *w_next_arr, \
    double *p_previous_arr, double *p_arr, double *p_next_arr) {
/* Loops over all time values solving the equations */
    /* Configures right-hand-side vector, which for now is w_next_arr */

    // Sets w_next_arr = B * w
    matrix_multiply(w_next_arr, B_static, w_arr, 1, 0);
    for (int i = 0; i < M; i++) {
        printf("w_next_arr[%d] = %g\n", i, w_next_arr[i]);
    }

    // Sets w_next_arr = w_next_arr - A * w_previous
    matrix_multiply(w_next_arr, A_static, w_previous_arr, -1, 1);
    for (int i = 0; i < M; i++) {
        printf("w_next_arr[%d] = %g\n", i, w_next_arr[i]);
    }

    // Sets w_next_arr = w_next_arr + pressure term
    for (int i = 0; i < M; i++) {
        w_next_arr[i] += (Deltax * Deltax / BETA) * (0.25 * p_previous_arr[i] + 0.5 * p_arr[i] + 0.25 * p_next_arr[i]);
    }
    for (int i = 0; i < M; i++) {
        printf("w_next_arr[%d] = %g\n", i, w_next_arr[i]);
    }


    /* Solves matrix  equation */
    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in w_next_arr
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, w_next_arr, ldb);
    printf("info = %d\n", info);
    for (int i = 0; i < M; i++) {
        printf("w_next_arr[%d] = %g\n", i, w_next_arr[i]);
    }


    for (int i = 0; i < M; i++) {
        printf("w_next_arr[%d] = %g\n", i, w_next_arr[i]);
    }
}
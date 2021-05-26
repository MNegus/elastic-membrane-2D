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
double *w_previous, *w, *w_deriv, *rhs, *A_static, *B_static, *A;

// Time definitions
double t = 0; // Time variable 
int k = 0; // Timestep variable, such that t = DELTA_T * k

// Output values
char w_name[40] = "w";
char w_deriv_name[40] = "w_deriv";

/* Function definitions */
void init();
void initialise_coefficient_matrices();
void B_multiply(double *result, double *w_arr, double scale, int ADD);
void initialise_membrane();
void output_membrane(double *w_arr, char *name);
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
    Dbeta2 = (ALPHA * Deltax * Deltax) / (BETA * DELTA_T * DELTA_T);
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

    // Creates coefficient matrices
    A_static = malloc(M * noRows * sizeof(double)); // Used for dgbsv
    B_static = malloc(M * noRows * sizeof(double)); // Used for matrix multiplication

    // Fills in the coefficient matrix
    initialise_coefficient_matrices();

    /* Initialise w arrays */
    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_deriv = malloc(M * sizeof(double)); // Time derivative of w
    rhs = malloc(M * sizeof(double)); // Right-hand-side vector
    A = malloc(M * noRows * sizeof(double)); // Coefficient matrix

    // Initialise w_previous and w using second-order initial conditions
    initialise_membrane();

}


void initialise_coefficient_matrices() {
/* initialise_coefficient_matrices 
Function to fill in the Mx4 coefficient matrices A and B, where each row is one 
of the diagonals. Dbeta2 is the scaled term equal to:
Dbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * Deltax * Deltax).
*/

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

    // Loop over points in x_arr
    for (int i = 0; i < M; i++) {
        double y_val = 0;

        // Loop over the appropriate diagonals
        for (int j = fmax(i - 1, 0); j <= fmin(i + 1, M - 1); j++) {
            y_val += matrix_arr[(2 + i - j) * M + j] * x_arr[j];
        }
        y_arr[i] = ADD * y_arr[i] + scale * y_val;
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
        
        // Saves w_val into w_previous, and sets w_deriv to 0
        if (i < M) {
            w_previous[i] = w_val;
            w_deriv[i] = 0;
        }
       
        // Increment i
        i++;
    }
    fclose(initial_condition_file);

    // Outputs w_previous and w_deriv
    output_membrane(w_previous, w_name);
    output_membrane(w_deriv, w_deriv_name);

    // Increments the time t and timestep k
    t += DELTA_T;
    k++;

    /* Initialise w by solving the matrix equation */
    // Sets rhs to be equal to 0.5 * B * w_previous
    matrix_multiply(rhs, B_static, w_previous, 0.5, 0);

    // Copies over elements to A
    memcpy(A, A_static, M * noRows * sizeof(double));

    // Uses LAPACK to solve the matrix equation, saving result in rhs
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

    // Swaps rhs and w
    double *temp = w; 
    w = rhs;
    rhs = temp;

    /* First-order initial condition */
    // for (int i = 0; i < M; i++) {
    //     w[i] = w_previous[i];
    // }

    // Outputs w
    output_membrane(w, w_name);
}


void output_membrane(double *w_arr, char *name) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "mitchell_outputs/%s_%d.txt", name, k);
    FILE *w_file = fopen(w_filename, "w");

    // Outputs from x = 0 to L - dx
    for (int i = 0; i < M; i++) {
        double x = i * Deltax;
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[i]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * Deltax;
    fprintf(w_file, "%.10f, %.10f\n", x, 0.0);

    fclose(w_file);
}


void run() {
/* Loops over all time values solving the equations */
    while (t <= T_MAX) {
        printf("t = %g\n", t);
        
        /* Configures right-hand-side vector */

        // Sets rhs = B * w
        matrix_multiply(rhs, B_static, w, 1, 0);

        // Sets rhs = rhs - A * w_previous
        matrix_multiply(rhs, A_static, w_previous, -1, 1);

        /* Solves matrix  equation */
        // Copies over elements to A
        memcpy(A, A_static, M * noRows * sizeof(double));

        // Uses LAPACK to solve the matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, M, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        /* Saves w_deriv, given that rhs = w_next */
        for (int i = 0; i < M; i++) {
            w_deriv[i] = (rhs[i] - w_previous[i]) / (2 * DELTA_T);
        }
        output_membrane(w_deriv, w_deriv_name);

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
        output_membrane(w, w_name);
        
    }
}
#include <stdio.h>
#include <lapacke.h>
#include <math.h>
#include "parameters.h"

int main (int argc, const char * argv[])
{
    /* Solves the membrane equation 
    ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
    This ends up becoming a matrix 
    equation
    A w^(k+1) = 2 * w^k - w^(k-1),
    with some big complicated matrix A
    */

    /* Parameters */
    double pi = 3.1415926535;

    // Derived parameters
    double dx = L / (N_MEMBRANE - 1); // May need to make it L / N
    double Cbeta2 = BETA * DELTA_T * DELTA_T / (ALPHA * dx * dx);
    double Cgamma2 = GAMMA * DELTA_T * DELTA_T / (ALPHA * dx * dx * dx * dx); 
    int N = N_MEMBRANE;

    /* Matrix equations set-up. The first KL = 2 rows are empty, then the rows
    are the diagonals, appropriately shifted */
    int KL = 2; // Number of sub-diags
    int KU = 2; // Number of super-diags
    int noRows = 2 * KL + KU + 1; // Number of rows in coefficient matrix

    // Lapack constants
    lapack_int info, n, kl, ku, nrhs, ldab, ldb;
    lapack_int ipiv[N];
    n = N;
    kl = KL;
    ku = KU;
    nrhs = 1;
    ldab = N;
    ldb = 1;
    
    // Creates coefficient matrix
    double *A_static = malloc(N * noRows * sizeof(double));

    // Stores upper-upper-diagonal in third row
    for (int colNum = 2; colNum < N; colNum++) {
        double value;
        if (colNum == 2) {
            value = 2 * Cgamma2;
        } else {
            value = Cgamma2;
        }
        A_static[2 * N + colNum] = value;
    }

    // Stores upper diagonal in fourth row
    for (int colNum = 1; colNum < N; colNum++) {
        double value;
        if (colNum == 1) {
            value = - 2 * (Cbeta2 + 4 * Cgamma2);
        } else {
            value = -(Cbeta2 + 4 * Cgamma2);
        }
        A_static[3 * N + colNum] = value;
    }

    // Stores main diagonal in fifth row
    for (int colNum = 0; colNum < N; colNum++) {
        double value;
        if (colNum == 1) {
            value = 1 + 2 * Cbeta2 + 7 * Cgamma2;
        } else {
            value = 1 + 2 * Cbeta2 + 6 * Cgamma2;
        }
        A_static[4 * N + colNum] = value;
    }

    // Stores lower diagonal in sixth row
    for (int colNum = 0; colNum < N - 1; colNum++) {
        double value;
        if (colNum == N - 2) {
            value = -(Cbeta2 + 3 * Cgamma2);
        } else {
            value = -(Cbeta2 + 4 * Cgamma2);
        }
        A_static[5 * N + colNum] = value;
    }

    // Stores lower-lower diagonal in seventh row
    for (int colNum = 0; colNum < N - 2; colNum++) {
        A_static[6 * N + colNum] = Cgamma2;
    }

    // // Output A
    // printf("A = \n");
    // for(int rowNum=0; rowNum < noRows; rowNum++) {
    //     for(int colNum=0; colNum < N; colNum++) {
    //         printf("%g ", A_static[rowNum * N + colNum]);
    //     }
    //     printf("\n");
    // }

    /* Initialise w arrays */
    double *w_previous = malloc(N * sizeof(double));
    double *w = malloc(N * sizeof(double));
    double *rhs = malloc(N * sizeof(double));
    double *A = malloc(N * noRows * sizeof(double));

    // First-order initial conditions
    for (int i = 0; i < N; i++) {
        double x = dx * i;
        double w_val = cos(M_PI * x / L) + 1;
        w_previous[i] = w_val;
        w[i] = w_val;
    }

    /* Outputs w_previous */
    double t = 0;
    int k = 0;
    char output_filename[40];
    sprintf(output_filename, "outputs/output_%d.txt", k);
    FILE *output_file = fopen(output_filename, "w");
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        fprintf(output_file, "%g, %g\n", x, w_previous[i]);
    }
    fclose(output_file);

    // Increments times
    t += DELTA_T;
    k++;

    /* Loops over all timesteps */
    while (t < T_MAX) {
        printf("t = %g\n", t);
        // Outputs solution for w;
        char output_filename[40];
        sprintf(output_filename, "outputs/output_%d.txt", k);
        FILE *output_file = fopen(output_filename, "w");
        for (int i = 0; i < N; i++) {
            double x = i * dx;
            fprintf(output_file, "%g, %g\n", x, w[i]);
        }
        fclose(output_file);

        // Configures right-hand-side vector
        for (int i = 0; i < N; i++) {
            rhs[i] = 2 * w[i] - w_previous[i];
        }

        // Copies over elements to A
        for (int i = 0; i < N * noRows; i++) {
            A[i] = A_static[i];
        }
        // memcpy(A, A_static, N * noRows * sizeof(double));

        // Solves matrix equation, saving result in rhs
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, n, kl, ku, nrhs, A, ldab, ipiv, rhs, ldb);

        // Shifts arrays, updating w
        double *temp = w_previous;
        w_previous = w;
        w = rhs;
        rhs = temp;

        // Increments time
        t += DELTA_T;
        k++;
    }
    printf("Finished with t = %g\n", t);
}
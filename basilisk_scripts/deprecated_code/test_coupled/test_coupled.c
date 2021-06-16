
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "wave-equation.h" // For solving the wave equation


double *w_previous, *w, *w_next, *w_deriv; // Membrane position arrays
double *p_previous_arr, *p_arr, *p_next_arr; // Pressure arrays

int main() {

    int M = 15;
    double DELTA_T = 1e-4;
    double MEMBRANE_RADIUS = 4;
    double ALPHA = 1;
    double BETA = 1;

    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_next = malloc(M * sizeof(double)); // w at next timestep
    w_deriv = malloc(M * sizeof(double)); // Time derivative of w
    p_previous_arr = malloc(M * sizeof(double)); // p at previous timestep
    p_arr = malloc(M * sizeof(double)); // p at current timestep
    p_next_arr = malloc(M * sizeof(double)); // p at next timestep

    
    initialise_membrane(w_previous, w, w_deriv, p_previous_arr, p_arr, \
                p_next_arr, M + 1, DELTA_T, MEMBRANE_RADIUS, ALPHA, BETA);



    fprintf(stderr, "Starting\n");
    run();
    fprintf(stderr, "Finished\n");
}
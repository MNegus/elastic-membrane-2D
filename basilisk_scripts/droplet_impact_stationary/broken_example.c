/* broken_example.c
*/

#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "tension.h" // Surface tension of droplet


int main() {
/* Main function for running the simulation */

    /* Runs the simulation */
    fprintf(stderr, "Started!\n");
    run(); 
    fprintf(stderr, "Finished!\n");
}

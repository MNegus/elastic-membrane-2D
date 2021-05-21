/* parameters.h
 Provides a list of parameters to be used across difference C programs 
*/

/* Physical parameters */
const double ALPHA = 0.1; // Mass term
const double BETA = 10; // Tension term
const double GAMMA = 0; // Bending stiffness term
const double L = 4; // Width of domain in x


/* Computational parameters */
const int N_MEMBRANE = 4096; // Number of grid points on the membrane
const double T_MAX = 0.25; // Maximum value of time
const double DELTA_T = 1e-5; // Timestep size


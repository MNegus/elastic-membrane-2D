/* droplet_impact_coupled.c
    A 2D droplet falling towards an elastic membrane lying 
    along the boundary at y = 0. The solution for the membrane is given as a 
    function of pressure by the routines defined in wave-equation.h, and its 
    velocity is fed back into Basilisk by altering the boundary condition. 

    Runs until the turnover point approximately reaches the initial droplet 
    radius. 

    Author: Michael Negus
*/

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "parameters.h" // Includes all defined parameters
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "membrane-equation.h" // For solving the membrane equation
#include <omp.h> // For openMP parallel
#include <stdlib.h>


/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double DROP_REFINED_WIDTH; // Width of the refined area around the droplet
double MEMBRANE_REFINE_NO; // Number of grid cells above the membrane to refine
double MEMBRANE_REFINED_HEIGHT; // Width of the refined area above the membrane 
double DROP_CENTRE; // Initial centre position of the droplet
double IMPACT_TIME; // Theoretical time of impact

/* Membrane parameters and arrays */
double DELTA_X; // Spatial grid size
int M; // Number of grid nodes for membrane solution
double *w_previous, *w, *w_next, *w_deriv; // Membrane position arrays
double *p_previous_arr, *p_arr, *p_next_arr; // Pressure arrays

/* Global variables */
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 0; // Records how many GFS files have been outputted
int log_output_no = 0; // Records how many plate data files there have been
int membrane_output_no = 0; // Records how many membrane outputs there have been
int start_membrane = 0; // Boolean to indicate if membrane motion has started

/* Function definitions */
double membrane_bc(double x, double *w_deriv_arr);
void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr);

/* Boundary conditions */
// Symmetry on left boundary
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary

// Conditions on surface
uf.n[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);

// Conditions for entry from above
u.n[top] = neumann(0.); // Allows outflow through boundary
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the horizontal direction
u.n[right] = neumann(0.); // Allows outflow through boundary
p[right] = dirichlet(0.); // 0 pressure far from surface


int main() {
/* Main function for running the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time
    MEMBRANE_REFINE_NO = 8; // Number of cells above membrane to refine by
    MEMBRANE_REFINED_HEIGHT = MEMBRANE_REFINE_NO * MIN_CELL_SIZE; 
    M = (int) floor(pow(2, MAXLEVEL) * MEMBRANE_RADIUS / BOX_WIDTH);
    DELTA_X = MEMBRANE_RADIUS / M;

    /* Define membrane arrays */
    w_previous = malloc(M * sizeof(double)); // w at previous timestep
    w = malloc(M * sizeof(double)); // w at current timestep
    w_next = malloc(M * sizeof(double)); // w at next timestep
    w_deriv = malloc(M * sizeof(double)); // Time derivative of w
    p_previous_arr = malloc(M * sizeof(double)); // p at previous timestep
    p_arr = malloc(M * sizeof(double)); // p at current timestep
    p_next_arr = malloc(M * sizeof(double)); // p at next timestep

    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Runs the simulation */
    run(); 
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Refines around the droplet */
    refine((((sq(x) + sq(y - DROP_CENTRE) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) && (sq(x) + sq(y - DROP_CENTRE)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH))) || ((x < MEMBRANE_RADIUS) && (y <= MEMBRANE_REFINED_HEIGHT))) \
        && (level < MAXLEVEL));
    
    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x) - sq(y - DROP_CENTRE) + sq(DROP_RADIUS));

    /* Initialise the droplet velocity downwards */
    foreach() {
        u.y[] = DROP_VEL * f[];
    }
    boundary ((scalar *){u});

    /* Initialises membrane arrays */
    #pragma omp parallel for
    for (int k = 0; k < M; k++) {
        w_previous[k] = 0.0;
        w[k] = 0.0;
        w_deriv[k] = 0.0;
        p_previous_arr[k] = 0.0;
        p_arr[k] = 0.0;
    }
}


event refinement (i++) {
/* Adaptive grid refinement */

    // Adapts with respect to velocities and volume fraction 
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-4},
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    // Refines above the membrane
    refine((x < MEMBRANE_RADIUS) && (y <= MEMBRANE_REFINED_HEIGHT) \
        && level < MAXLEVEL);
}


event gravity (i++) {
/* Adds acceleration due to gravity in the vertical direction */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.y[] -= 1./sq(FR); // Adds acceleration due to gravity
}


event small_droplet_removal (i++) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    int remove_droplet_radius = min(20, (int)(0.2 / MIN_CELL_SIZE));
    // int remove_droplet_radius = (int)(0.25 / MIN_CELL_SIZE);
    remove_droplets(f, remove_droplet_radius);

    // Also remove air bubbles
    remove_droplets(f, remove_droplet_radius, 1e-4, true);
}

event update_membrane(t += DELTA_T) {
/* Updates the membrane arrays by solving the membrane equation, and outputs*/

    /* Update pressure arrays */
    // Swaps
    double *temp1 = p_previous_arr;
    p_previous_arr = p_arr;
    p_arr = p_next_arr;
    p_next_arr = temp1;

    // Fills pressure in from boundary nodes
    foreach_boundary(bottom) {
        if (x <= MEMBRANE_RADIUS) {
            int k = (int) (x / DELTA_X);
            p_next_arr[k] = p[];
        }
    }

    /* Updates membrane position after the start time */
    if (t >= MEMBRANE_START_TIME) {
        if (start_membrane == 0) {
            /* Initialise membrane motion */
            start_membrane = 1; // Indicates membrane motion has started

            // Initialises w and w_deriv
            initialise_membrane(w_previous, w, w_deriv, p_previous_arr, p_arr, \
                p_next_arr, M + 1, DELTA_T, MEMBRANE_RADIUS, ALPHA, BETA, GAMMA);
        } else {
            /* Solve membrane equation to determine w_next */
            membrane_timestep(w_previous, w, w_next, w_deriv, p_previous_arr, \
                p_arr, p_next_arr, DELTA_T);
        }
    }

    /* Updates boundary condition if COUPLED is set*/
    if (COUPLED) {
        uf.n[bottom] = dirichlet(membrane_bc(x, w_deriv)); 
    }

    /* Outputs membrane arrays */
    output_arrays(w, w_deriv, p_arr);

    // Swaps membrane arrays 
    double *temp2 = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp2;
}


event output_data (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */
    /* Outputs data to log file */
    fprintf(stderr, \
        "t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    fclose(logfile);

    /* Outputs info about cells along membrane */
    char output_filename[80];
    sprintf(output_filename, "boundary_output_%d.txt", log_output_no);
    FILE *output_file = fopen(output_filename, "w");

    foreach_boundary(bottom) {
        fprintf(output_file, "%g %g %g\n", x, p[], uf.y[]);
    }
    fclose(output_file);

    log_output_no++;
}

event gfs_output (t += GFS_OUTPUT_TIMESTEP) {
/* Saves a gfs file */
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}


event movies (t += 1e-3) {
/* Produces movies using bview */ 
    if (MOVIES) {
        // Creates a string with the time to put on the plots
        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        /* Zoomed out view */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 12.0, ty = -0.31, tx = -0.31);

        /* Movie of the volume fraction of the droplet */
        clear();
        draw_vof("f", lw = 2);
        squares("f", linear = true, spread = -1, linear = true, map = cool_warm); // RC - minor changes here and beyond
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Movie of the horiztonal velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.x", linear = false, spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("horizontal_vel.mp4");


        /* Movie of the vertical velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", linear = false, spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("vertical_vel.mp4");

        /* Movie of the pressure */
        clear();
        draw_vof("f", lw = 2);
        squares("p", linear = false, spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("pressure.mp4");

        /* Zoomed in view of pressure around entrapped bubble */
        /*
        // Set up bview box
        view (width = 1024, height = 1024, fov = 5.0, ty = -0.1, \
            quat = {0, 0, -0.707, 0.707});

        clear();
        draw_vof("f", lw = 2);
        squares("p", linear = false, spread = -1, linear = true, map = cool_warm);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("p", linear = false, spread = -1, linear = true, map = cool_warm);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("zoomed_pressure.mp4");
        */
    }
}


event end (t = MAX_TIME) {
/* Ends the simulation */ 

    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "Finished after %g seconds\n", end_wall_time - start_wall_time);
    fclose(logfile);

}

double membrane_bc(double x, double *w_deriv_arr) {
/* membrane_bc
Outputs the boundary condition for the vertical face velocity, uf.n, at the 
bottom boundary, which matches the velocity of the membrane 
*/
    if (x >= MEMBRANE_RADIUS) {
        return 0.;
    } else {
        int k = (int) (x / DELTA_X);
        return -w_deriv_arr[k];
    }
}

void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr) {
/* output_membrane
Outputs the x positions of the membrane into a text file
*/
    char w_filename[40];
    sprintf(w_filename, "w_%d.txt", membrane_output_no);
    FILE *w_file = fopen(w_filename, "w");

    char w_deriv_filename[40];
    sprintf(w_deriv_filename, "w_deriv_%d.txt", membrane_output_no);
    FILE *w_deriv_file = fopen(w_deriv_filename, "w");

    char p_filename[40];
    sprintf(p_filename, "p_%d.txt", membrane_output_no);
    FILE *p_file = fopen(p_filename, "w");

    // Outputs from x = 0 to L - dx
    #pragma omp parallel for
    for (int k = 0; k < M; k++) {
        double x = k * DELTA_X;
        fprintf(w_file, "%.10f, %.10f\n", x, w_arr[k]);
        fprintf(w_deriv_file, "%.10f, %.10f\n", x, w_deriv_arr[k]);
        fprintf(p_file, "%.10f, %.10f\n", x, p_arr[k]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * DELTA_X;
    fprintf(w_file, "%.10f, %.10f\n", x, 0.0);
    fprintf(p_file, "%.10f, %.10f\n", x, 0.0);
    fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

    fclose(w_file);
    fclose(p_file);
    fclose(w_deriv_file);

    membrane_output_no++;
}


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
#include "heights.h"
#include "contact.h" // For imposing contact angle on the surface
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
int num_threads; // Number of OpenMP threads

/* Membrane parameters and arrays */
double DELTA_X; // Spatial grid size
int M; // Number of grid nodes for membrane solution
double *w_previous, *w, *w_next, *w_deriv; // Membrane position arrays
double *p_previous_arr, *p_arr, *p_next_arr; // Pressure arrays

/* Turnover point search arrays */
// double *turnover_x, *turnover_y, *turnover_x_vel, *turnover_y_vel;

/* Global variables */
int impact = 0; // Sets to 1 once impact has happened
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 0; // Records how many GFS files have been outputted
int log_output_no = 0; // Records how many plate data files there have been
int interface_output_no = 0; // Records how many interface files there have been
int membrane_output_no = 0; // Records how many membrane outputs there have been
int start_membrane = 0; // Boolean to indicate if membrane motion has started
double drop_thresh = 1e-4; // Remove droplets threshold
double pinch_off_time = 0.; // Time of pinch-off


/* Function definitions */
double membrane_bc(double x, double *w_deriv_arr);
void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr);
void output_arrays_stationary(double *p_arr);
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);

/* Contact angle variables */ 
vector h[];  // Height function
double theta0 = 90;  // Contact angle in degrees

/* Boundary conditions */
// Symmetry on left boundary
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary

// Conditions on surface
uf.n[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);
h.t[bottom] = contact_angle (theta0*pi/180.);

// Conditions for entry from above
u.n[top] = neumann(0.); // Allows outflow through boundary
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the horizontal direction
u.n[right] = neumann(0.); // Allows outflow through boundary
p[right] = dirichlet(0.); // 0 pressure far from surface

/* Field definitions */
vector h_test[]; // Height function field

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
    f.height = h; // For contact angle calculation
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time
    MEMBRANE_REFINE_NO = 8; // Number of cells above membrane to refine by
    MEMBRANE_REFINED_HEIGHT = MEMBRANE_REFINE_NO * MIN_CELL_SIZE; 

    /* Allocates memory for turnover point search arrays */
    num_threads = atoi(getenv("OMP_NUM_THREADS"));
    // if (num_threads == 0) return(1);
    // turnover_x = malloc(num_threads * sizeof(double));
    // turnover_y = malloc(num_threads * sizeof(double));
    // turnover_x_vel = malloc(num_threads * sizeof(double));
    // turnover_y_vel = malloc(num_threads * sizeof(double));

    /* Membrane constants */
    // Number of grid points, depending on MAXLEVEL and FD_COARSEN_LEVEL
    // M = (int) floor(pow(2, MAXLEVEL) * MEMBRANE_RADIUS / (BOX_WIDTH * pow(2, FD_COARSEN_LEVEL)));
    // DELTA_X = MEMBRANE_RADIUS / M;
    DELTA_X = MIN_CELL_SIZE * pow(2, FD_COARSEN_LEVEL);
    M = (int) floor(MEMBRANE_RADIUS / DELTA_X) + 1;
    fprintf(stderr, "M = %d, MEMBRANE_RADIUS / DELTA_X = %g\n", M, MEMBRANE_RADIUS / DELTA_X);

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

    /* Creates turnover point file */
    FILE *turnover_point_file = fopen("turnover_points_basilisk.txt", "w");
    fclose(turnover_point_file);

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

    /* Runs the simulation */
    run(); 
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Refines around the droplet */
    refine((sq(x) + sq(y - DROP_CENTRE) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x) + sq(y - DROP_CENTRE)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
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
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-3, 1e-3, 1e-6},
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    /* Attempts to refine above the membrane, doubling the refine height until
    successful */
    double refine_height = MEMBRANE_REFINED_HEIGHT;
    int adequate_refinement = 0;

    while (adequate_refinement == 0) {

        // Attempts to refine
        refine((x < MEMBRANE_RADIUS) && (y <= refine_height) \
        && level < MAXLEVEL);

        // Refines a box near origin
        refine((x <= 0.2) && (y <= 0.1) && level < MAXLEVEL);

        // Loops and check if refinement was successful
        adequate_refinement = 1;
        foreach_boundary(bottom) {
            if ((x < MEMBRANE_RADIUS) && (level < MAXLEVEL)) {
                adequate_refinement = 0;
                break;
            }
        }

        // If refinement was unsuccessful, then double the refined height
        if (adequate_refinement == 0) refine_height = 2 * refine_height;
    }
    
}


event gravity (i++) {
/* Adds acceleration due to gravity in the vertical direction */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.y[] -= 1./sq(FR); // Adds acceleration due to gravity
}


event small_droplet_removal (t += 1e-4) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    int remove_droplet_radius = min(16, (int)(0.2 / MIN_CELL_SIZE));
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

    /* Checks if we have impacted yet, which is defined to be when the first 
    time a cell with f[] = 1 is found on the boundary. If so, then finds the 
    furthest-right cell to which f > 0, and then saves the value of the pressure
    one cell to the right of it. This will then be used as the value of pressure
    for any cells where 0 < f < 1. */

    if (impact == 0) {
        foreach_boundary(bottom) {
            // Set impact = 1 if this is the first time we have detected f[] == 1
            if (f[] == 1) {
                impact = 1;
                break;
            }
        }
    }

    /* Saves pressure along boundary. Sets pressure to 0 if either:
        * The boundary cell is mixed (i.e. 0 < f[] < 1)
        * The height in the y direction at the cell is less than y_min_height 
        * The height in the x direction at the cell is less than x_min_height
    */
    if (CUTOFF && impact) heights(f, h_test); // Associates h_test with the heights of f

    // Iterates over bottom boundary
    foreach_boundary(bottom) {

        // Skip if x is not above the membrane
        if (x > MEMBRANE_RADIUS) continue;

        // Determines index in array
        double fractpart, intpart;
        fractpart = modf((x - 0.5 * MIN_CELL_SIZE) / DELTA_X, &intpart);
        if (fractpart > 1e-6) continue;
        
        int k = (int) intpart;

        // Initially fills p_next_arr[k] with the current value of pressure
        p_next_arr[k] = p[];

        // If we are pre-impact or are not doing a cutoff, then continue in loop
        if ((CUTOFF == 0) || (impact == 0)) continue;

        /* Implements cutoff */
        if ((f[] > 0) && (f[] < 1)) {
            /* If a mixed cell, output 0 */
            // p_scale = 0;
            p_next_arr[k] = 0.0;
        } else {
            /* Else, f[] == 0 or 1, and we check the heights */

            // Check the vertical height, if defined
            if (h_test.y[] != nodata) {
                if (height(h_test.y[]) - 0.5 < y_min_height) {
                    // If the interface position is less than min_height
                    // p_scale = 0;
                    p_next_arr[k] = 0.0;
                }
            }

            // If too close to the edges, then skip the x cutoff
            if ((k < x_min_height) && (k > M - x_min_height)) continue;

            // x cutoff procedure
            #pragma omp parallel for
            for (int q = -x_min_height; q <= x_min_height; q++) {
                // Does not cutoff a non-mixed cell
                if ((f[q, 0] == 1) || (f[q, 0] == 0)) continue;
                
                // If a mixed cell is detected then set p_next_arr[k] = 0
                p_next_arr[k] = 0.0;
            }
        }
    }

    /* Updates membrane position after the start time */
    if ((t >= MEMBRANE_START_TIME) && (STATIONARY == 0)) {
        if (start_membrane == 0) {
            /* Initialise membrane motion */
            start_membrane = 1; // Indicates membrane motion has started

            // Initialises w and w_deriv
            initialise_membrane(w_previous, w, w_deriv, p_previous_arr, p_arr, \
                p_next_arr, M, DELTA_X, DELTA_T, MEMBRANE_RADIUS, ALPHA, BETA, GAMMA);
        } else {
            /* Solve membrane equation to determine w_next */
            membrane_timestep(w_previous, w, w_next, w_deriv, p_previous_arr, \
                p_arr, p_next_arr, M, DELTA_X, DELTA_T);
        }

        /* Updates boundary condition if COUPLED is set*/
        if (COUPLED) {
            uf.n[bottom] = dirichlet(membrane_bc(x, w_deriv)); 
        }
    }

    /* Outputs membrane arrays */
    if (STATIONARY) {
        output_arrays_stationary(p_arr);
    } else {
        output_arrays(w, w_deriv, p_arr);
    }

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


event output_interface (t += DELTA_T) {
/* Outputs the interface locations of the droplet */
    // Creates text file to save output to
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d.txt", interface_output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    // Outputs the interface locations and closes the file
    output_facets(f, interface_file);
    fclose(interface_file);

    interface_output_no++;
}


event output_turnover_point (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs the coordinates of the turnover point and its velocity */

    // Opens the turnover points file
    FILE *turnover_point_file = fopen("turnover_points_basilisk.txt", "a");
    
    if (impact == 0) {
        /* If we are pre-impact, output the turnover point to be at (0, 0) */
        fprintf(turnover_point_file, "%.4f, %.4f, %.4f, %.4f, %.4f\n", t, 0., 0., 0., 0.);
    } else {
        /* Else we find the turnover point to be the point along the bottom half
        of the droplet with the lowest value of x */
        // NOTE: WILL NOT WORK IF ENTRAPPED BUBBLE IS PRESENT

        // Maximum value of y to look is the half way point of the droplet
        // double max_y = INITIAL_DROP_HEIGHT + 0.5 * DROP_RADIUS - t;

        // Maximum value of y is taken to be 0.25, which the turnover point 
        // should not reach in the current timescale
        double max_y = 0.25;

        // // Initialises arrays
        // #pragma omp parallel for
        // for (int k = 0; k < num_threads; k++) {
        //     turnover_x[k] = 1e6;
        //     turnover_y[k] = 0;
        //     turnover_x_vel[k] = 0;
        //     turnover_y_vel[k] = 0; 
        // }

        double turnover_x = 1e6;

        // Variables for the y coordinate and velocities of the turnover point,
        // initialised to 0 in case a  turnover point is not found
        double turnover_y = 0;
        double turnover_x_vel = 0;
        double turnover_y_vel = 0;

        // Iterates over the interface to find the turnover point
        foreach(reduction(min:turnover_x)) {
            // Checks if current point is on the interface
            if ((f[] > 1e-6 && f[] < 1. - 1e-6)) {

                // Finds the segment of the interface
                // coord n = facet_normal (point, c, s);
                coord n = interface_normal(point, f);
                double alpha = plane_alpha (f[], n);
                coord segment[2];
                facets (n, alpha, segment);

                // Finds the end of the segment with the lowest value of x
                double current_x, current_y;
                if (segment[0].x < segment[1].x) {
                    current_x = x + segment[0].x * Delta;
                    current_y = y + segment[0].y * Delta;
                } else {
                    current_x = x + segment[1].x * Delta;
                    current_y = y + segment[1].y * Delta;
                }

                // If current_x and current_y are in the appropriate bounds, and
                // current_x < turnover_x, then save
                if ((current_x < turnover_x) \
                    && (current_x >= 0) && (current_y >= 0) \
                    && (current_y <= max_y)) {
                    turnover_x = current_x;
                    turnover_y = current_y;
                    turnover_x_vel = uf.x[];
                    turnover_y_vel = uf.y[];
                }

            } 
        }

        // Outputs the turnover point data
        fprintf(turnover_point_file, "%.4f, %.4f, %.4f, %.4f, %.4f\n", \
            t, turnover_x, turnover_y, turnover_x_vel, turnover_y_vel);
    }

    fclose(turnover_point_file);
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
        squares("f", linear = true, spread = -1, map = cool_warm); // RC - minor changes here and beyond
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Movie of the horiztonal velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.x", spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("horizontal_vel.mp4");


        /* Movie of the vertical velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("vertical_vel.mp4");

        /* Movie of the pressure */
        clear();
        draw_vof("f", lw = 2);
        squares("p", spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("pressure.mp4");

        /* Zoomed in view of pressure around entrapped bubble */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 2.0, ty = -0.05, tx = -0.05);
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("zoomed_vertical_vel.mp4");
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
        // fprintf(w_file, "%.10f, %.10f\n", x, w_arr[k]);
        // fprintf(w_deriv_file, "%.10f, %.10f\n", x, w_deriv_arr[k]);
        // fprintf(p_file, "%.10f, %.10f\n", x, p_arr[k]);
        fprintf(w_file, "%g, %g\n", x, w_arr[k]);
        fprintf(w_deriv_file, "%g, %g\n", x, w_deriv_arr[k]);
        fprintf(p_file, "%g, %g\n", x, p_arr[k]);
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


void output_arrays_stationary(double *p_arr) {
/* output_membrane_stationary
Outputs the x positions of the pressure in a text file
*/
    char p_filename[40];
    sprintf(p_filename, "p_%d.txt", membrane_output_no);
    FILE *p_file = fopen(p_filename, "w");

    // Outputs from x = 0 to L - dx
    #pragma omp parallel for
    for (int k = 0; k < M; k++) {
        double x = k * DELTA_X;
        fprintf(p_file, "%g, %g\n", x, p_arr[k]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * DELTA_X;
    fprintf(p_file, "%.10f, %.10f\n", x, 0.0);

    fclose(p_file);

    membrane_output_no++;
}

/* Alternative remove_droplets definitions */
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit) {
    scalar d[], f = p.f;
    double threshold = p.threshold ? p.threshold : 1e-4;
    foreach()
    d[] = (p.bubbles ? 1. - f[] : f[]) > threshold;
    int n = tag (d), size[n], keep_tags[n];

    for (int i = 0; i < n; i++) {
        size[i] = 0;
        keep_tags[i] = 1;
    }
    foreach_leaf() {
        if (d[] > 0) {
            int j = ((int) d[]) - 1;
            size[j]++;
            if ((x < ignore_region_x_limit) && (y < ignore_region_y_limit)) {
                keep_tags[j] = 0;
            }
        }
    }
    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif
    int minsize = pow (p.minsize ? p.minsize : 3, dimension);
    foreach() {
        int j = ((int) d[]) - 1;
        if (d[] > 0 && size[j] < minsize && keep_tags[j] == 1)
            f[] = p.bubbles;
    }
    boundary ({f});
}


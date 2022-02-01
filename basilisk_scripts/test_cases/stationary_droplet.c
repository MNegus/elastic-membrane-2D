/* stationary_circle.c



    Author: Michael Negus
*/

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#define MOVING 1

#include "parameters.h" // Includes all defined parameters
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview

#include "tag.h" // For removing small droplets
// #include "contact.h" // For imposing contact angle on the surface
// #include "membrane-equation.h" // For solving the membrane equation
#include <omp.h> // For openMP parallel
#include <stdlib.h>


/* Physical constants */
double REYNOLDS; // Reynolds number of liquid
double WEBER; // Weber number of liquid
double FROUDE; // Froude number of liquid
double RHO_R; // Density ratio
double MU_R; // Viscosity ratio

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
double *wt, *wtt, *wx, *wxx, *wxxx, *wtx, *wtxx; // Derivative arrays for membrane 
double mag = 0.25; // Magnitude of membrane displacement

/* Turnover point search arrays */
double jet_energy = 0; // Energy in jet
double *turnover_x_arr, *turnover_y_arr;
scalar positions_x[];
scalar positions_y[];

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
FILE * fp_stats; 

/* Function definitions */
double membrane_bc(double x, double *w_deriv_arr);
void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr);
void output_arrays_stationary(double *p_arr);
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);
double membrane_position(double x);

/* Contact angle variables */ 
// vector h[];  // Height function
// double theta0 = 90;  // Contact angle in degrees

/* Boundary conditions */
// Symmetry on left boundary
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary

// Conditions on surface
uf.n[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);
// h.t[bottom] = contact_angle (theta0*pi/180.);

// Conditions for entry from above
u.n[top] = neumann(0.); // Allows outflow through boundary
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the horizontal direction
u.n[right] = neumann(0.); // Allows outflow through boundary
p[right] = dirichlet(0.); // 0 pressure far from surface

/* Field definitions */
// vector h_test[]; // Height function field
scalar mean_curvature[];

/* Scalar fields for membrane displacement */
scalar W[]; // Membrane displacement
scalar Wt[]; // Time derivative of membrane displacement
scalar Wtt[]; // Second time derivative of membrane displacement
scalar Wx[]; // x derivative of membrane displacement
scalar Wxx[]; // Second x derivative of membrane displacement
scalar Wxxx[]; // Third x derivative of membrane displacement
scalar Wtx[]; // Mixed time and x derivative of membrane displacement
scalar Wtxx[]; // Mixed time and second x derivative of membrane displacement

#include "tension.h" // Surface tension of droplet
#include "heights.h"

int main() {
/* Main function for running the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Determine physical constants */
    REYNOLDS = RHO_L * V * R / MU_L; // Reynolds number of liquid
    WEBER = RHO_L * V * V * R / SIGMA; // Weber number of liquid
    FROUDE = V / sqrt(9.81 * R); // Froude number of liquid
    RHO_R = RHO_G / RHO_L; // Density ratio
    MU_R = MU_G / MU_L; // Viscosity ratio

    /* Set dimensionless constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface
    
    /* Contact angle determination */
    // f.height = h; // Associates height field with the tracer

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time
    MEMBRANE_REFINE_NO = 8; // Number of cells above membrane to refine by
    MEMBRANE_REFINED_HEIGHT = MEMBRANE_REFINE_NO * MIN_CELL_SIZE; 

    /* Allocates memory for turnover point search arrays */
    num_threads = atoi(getenv("OMP_NUM_THREADS"));
    if (num_threads == 0) return(1);
    turnover_x_arr = malloc(num_threads * sizeof(double));
    turnover_y_arr = malloc(num_threads * sizeof(double));

    /* Membrane constants */
    // Number of grid points, depending on MAXLEVEL and FD_COARSEN_LEVEL
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

    /* Derivatives of membrane arrays */
    wt = malloc(M * sizeof(double)); // Time derivative of membrane displacement
    wtt = malloc(M * sizeof(double)); // Second time derivative of membrane displacement
    wx = malloc(M * sizeof(double)); // x derivative of membrane displacement
    wxx = malloc(M * sizeof(double)); // Second x derivative of membrane displacement
    wxxx = malloc(M * sizeof(double)); // Third x derivative of membrane displacement
    wtx = malloc(M * sizeof(double)); // Mixed time and x derivative of membrane displacement
    wtxx = malloc(M * sizeof(double)); // Mixed time and second x derivative of membrane displacement

    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Creates turnover point file */
    // FILE *turnover_point_file = fopen("turnover_points_basilisk.txt", "w");
    // fclose(turnover_point_file);

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

    /* Open stats file */
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");

    /* Runs the simulation */
    run(); 

    // Close stats file
    fclose(fp_stats);
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Initialises membrane arrays */
    #pragma omp parallel for
    for (int k = 0; k < M; k++) {
        // Initialise membrane based on a single mode
        double x = k * DELTA_X;
        // w_previous[k] = mag * cos(pi * x / (2 * MEMBRANE_RADIUS));
        // w[k] = w_previous[k];
        // wx[k] = - (pi * mag / (2 * MEMBRANE_RADIUS)) * sin(pi * x / (2 * MEMBRANE_RADIUS));
        // wxx[k] = - pow(pi / (2 * MEMBRANE_RADIUS), 2) * mag * cos(pi * x / (2 * MEMBRANE_RADIUS));
        // wxxx[k] = pow(pi / (2 * MEMBRANE_RADIUS), 3) * mag * sin(pi * x / (2 * MEMBRANE_RADIUS));
        w_previous[k] = 0;
        w[k] = 0;
        wx[k] = 0;
        wxx[k] = 0;
        wxxx[k] = 0;

        wt[k] = 0;
        wtt[k] = 0;
        wtx[k] = 0;
        wtxx[k] = 0;

        w_deriv[k] = 0.0;
        p_previous_arr[k] = 0.0;
        p_arr[k] = 0.0;
    }


    /* Refines around the droplet */
    refine((sq(x) + sq(y - membrane_position(x) - DROP_CENTRE) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x) + sq(y - membrane_position(x) - DROP_CENTRE)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    
    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x) - sq(y - membrane_position(x) - DROP_CENTRE) + sq(DROP_RADIUS));

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
        // refine((x <= 0.2) && (y <= 0.1) && level < MAXLEVEL);

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

event set_fields (i++) {
/* Sets the scalar fields for the membrane displacement by interpolating from 
the finite differences solution */

    /* Adds the values to the cells along the bottom */
    foreach_boundary(bottom) {
        if (x <= MEMBRANE_RADIUS) {
            // Determines index in array
            double fractpart, intpart;
            fractpart = modf((x - 0.5 * MIN_CELL_SIZE) / DELTA_X, &intpart);
            if (fractpart > 1e-6) continue;
            
            int k = (int) intpart;

            // Set field values to be equal to the values in the FD arrays
            W[] = w[k];
            Wx[] = wx[k];
            Wxx[] = wxx[k];
            Wxxx[] = wxxx[k];
            Wt[] = wt[k];
            Wtt[] = wtt[k];
            Wtx[] = wtx[k];
            Wtxx[] = wtxx[k];
        } else {
            // Else we are not above the membrane and we set W = 0
            W[] = 0.; Wx[] = 0.; Wxx[] = 0.; Wxxx[] = 0.; 
            Wt[] = 0.; Wtt[] = 0.; Wtx[] = 0.; Wtxx[] = 0.;
        }
    }

    /* Loops over the remaining cells and interpolates them from the boundary */
    foreach() {
        // If above the boundary
        if (y > 0.5 * MIN_CELL_SIZE) {
            // Interpolates from the boundary cells
            W[] = interpolate(W, x, 0.5 * MIN_CELL_SIZE);
            Wx[] = interpolate(Wx, x, 0.5 * MIN_CELL_SIZE);
            Wxx[] = interpolate(Wxx, x, 0.5 * MIN_CELL_SIZE);
            Wxxx[] = interpolate(Wxxx, x, 0.5 * MIN_CELL_SIZE);
            Wt[] = interpolate(Wt, x, 0.5 * MIN_CELL_SIZE);
            Wtt[] = interpolate(Wtt, x, 0.5 * MIN_CELL_SIZE);
            Wtx[] = interpolate(Wtx, x, 0.5 * MIN_CELL_SIZE);
            Wtxx[] = interpolate(Wtxx, x, 0.5 * MIN_CELL_SIZE);
        }
    }
}

event additional_acceleration (i++) {
/* Impose additional acceleration due to the change of frame */
    face vector av = a; // Acceleration at each face
    
    foreach_face(x) {
        /* Determine the derivatives of the horizontal velocity, u */
        double uy = (u.x[0, 1] - u.x[0, -1]) / (2 * Delta); // y derivative of u
        double uyy = (u.x[0, 1] - 2 * u.x[] + u.x[0, -1]) \
            / (Delta * Delta); // Second y derivative of u
        double uxy = (u.x[1, 1] - u.x[1, -1] - u.x[-1, 1] + u.x[-1, -1]) \
            / (4 * Delta * Delta); // Mixed x and y derivative of u

        /* Determine derivatives of pressure, p */
        double py = (p[0, 1] - p[0, -1]) / (2 * Delta); // y derivative of p

        /* Increment horizontal acceleration */
        av.x[] += (1 / rho[]) * (- Wx[] * py \
            + mu.x[] * (Wxx[] * uy + 2 * Wx[] * uxy + Wx[] * Wx[] * uyy));
    }

    foreach_face(y) {
        /* Determine derivatives of horizontal velocity, u */
        double ut = av.x[]; // Time derivative i.e. the current velocity
        double ux = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta); // x derivative of u
        double uxx = (u.x[1, 0] - 2 * u.x[] + u.x[-1, 0]) \
            / (Delta * Delta); // Second x derivative of u
        double uy = (u.x[0, 1] - u.x[0, -1]) / (2 * Delta); // y derivative of u
        double uyy = (u.x[0, 1] - 2 * u.x[] + u.x[0, -1]) \
            / (Delta * Delta); // Second y derivative of u
        double uxy = (u.x[1, 1] - u.x[1, -1] - u.x[-1, 1] + u.x[-1, -1]) \
            / (4 * Delta * Delta); // Mixed x and y derivative of u
        
        /* Determine derivatives of vertical velocity, v */
        double vy = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta); // y derivative of v
        double vyy = (u.y[0, 1] - 2 * u.y[] + u.y[0, -1]) \ 
            / (Delta * Delta); // Second y derivative of v
        double vxy = (u.y[1, 1] - u.y[1, -1] - u.y[-1, 1] + u.y[-1, -1]) \
            / (4 * Delta * Delta); // Mixed x and y derivative of v

        /* Incremement vertical acceleration */
        av.y[] += Wtt[] + 2 * Wtx[] * u.x[] + Wx[] * ut \
            + Wxx[] * u.x[] * u.x[] + Wx[] * ux * u.x[] + Wx[] * uy * u.y[] \
            + (mu.x[] / rho[]) \
                * (Wxx[] * vy + 2 * Wx[] * vxy + Wx[] * Wx[] * vyy \
                    - (Wtxx[] + Wxxx[] * u.x[] + 2 * Wxxx[] * ux + Wx[] * uxx \
                        + Wx[] * uyy + 3 * Wx[] * Wxx[] * uy \
                        + 2 * Wx[] * Wx[] * uxy + pow(Wx[], 3) * uyy));
    }
}
// event gravity (i++) {
// /* Adds acceleration due to gravity in the vertical direction */
//     face vector av = a; // Acceleration at each face
//     foreach_face(x) av.y[] -= 1. / sq(FROUDE); // Adds acceleration due to gravity
// }


event small_droplet_removal (t += 1e-4) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    int remove_droplet_radius = min(16, (int)(0.2 / MIN_CELL_SIZE));
    remove_droplets(f, remove_droplet_radius);

    // Also remove air bubbles
    remove_droplets(f, remove_droplet_radius, 1e-4, true);
    
}


event output_data (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */

    /* Determine (lab frame) curvature */
    curvature(f, mean_curvature);

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
        // squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        squares("u.y", linear = true, spread = -1, map = cool_warm);
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


event logstats (t += 1e-4) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}


event end (t = MAX_TIME) {
/* Ends the simulation */ 

    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "Finished after %g seconds\n", end_wall_time - start_wall_time);
    fclose(logfile);

    /* Saves a gfs file */
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;

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

double membrane_position(double x) {
    if (x <= MEMBRANE_RADIUS) {
        return mag * cos(pi * x / (2 * MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}


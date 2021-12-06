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
#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "heights.h"
#include "contact.h" // For imposing contact angle on the surface
#include "contact-embed.h"
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
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);

/* Contact angle variables */ 
vector h[];  // Height function
double theta0 = 90;  // Contact angle in degrees

/* Embedded boundary */
vertex scalar phi[];

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
    origin (0, - 0.26); // Shift the bottom boundary

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
    f.height = h; // Associates height field with the tracer

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time
    MEMBRANE_REFINE_NO = 8; // Number of cells above membrane to refine by
    MEMBRANE_REFINED_HEIGHT = MEMBRANE_REFINE_NO * MIN_CELL_SIZE; 


    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

    /* Embedded contact angle */
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;

    /* Runs the simulation */
    run(); 
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Embedded boundary at y = 0 */
    foreach_vertex()
        phi[] = y;
    boundary ({phi});
    fractions (phi, cs, fs);

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
    foreach_face(x) av.y[] -= 1. / sq(FROUDE); // Adds acceleration due to gravity
}


// event small_droplet_removal (t += 1e-4) {
// /* Removes any small droplets or bubbles that have formed, that are smaller than
//     a specific size */
//     // Removes droplets of diameter 5 cells or less
//     int remove_droplet_radius = min(16, (int)(0.2 / MIN_CELL_SIZE));
//     remove_droplets(f, remove_droplet_radius);

//     // Also remove air bubbles
//     remove_droplets(f, remove_droplet_radius, 1e-4, true);
    
// }


event output_data (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */
    /* Outputs data to log file */
    fprintf(stderr, \
        "t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    fclose(logfile);

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
        isoline("phi", 0, 0, lw=2, lc = {1, 1, 1});
        squares("f", linear = true, spread = -1, map = cool_warm); // RC - minor changes here and beyond
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        // /* Movie of the horiztonal velocity */
        // clear();
        // draw_vof("f", lw = 2);
        // squares("u.x", spread = -1, linear = true, map = cool_warm);
        // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        // save ("horizontal_vel.mp4");


        // /* Movie of the vertical velocity */
        // clear();
        // draw_vof("f", lw = 2);
        // squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        // save ("vertical_vel.mp4");

        // /* Movie of the pressure */
        // clear();
        // draw_vof("f", lw = 2);
        // squares("p", spread = -1, linear = true, map = cool_warm);
        // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        // save ("pressure.mp4");

        // /* Zoomed in view of pressure around entrapped bubble */
        // // Set up bview box
        // view (width = 1024, height = 1024, fov = 2.0, ty = -0.05, tx = -0.05);
        // clear();
        // draw_vof("f", lw = 2);
        // squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        // save ("zoomed_vertical_vel.mp4");
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


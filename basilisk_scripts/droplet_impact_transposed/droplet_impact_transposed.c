/* droplet_impact_transposed.c
    2D droplet impact onto a membrane, such that x is the vertical coordinate 
    y is the horizontal. 
*/

// Filtering for large viscosity ratios
#define FILTERED

// Filtered viscosity field
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "parameters.h" // Includes all defined parameters
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "heights.h" // For height function routines
#include "contact.h" // For imposing contact angle on the surface
#include "membrane-equation.h" // For solving the membrane equation
#include <omp.h> // For openMP parallel
#include <stdlib.h> // For string manipulation

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
double DROP_CENTRE; // Initial centre of the droplet
double IMPACT_TIME; // Theoretical time of impact
int num_threads; // Number of OpenMP threads

/* Membrane parameters and arrays */
double DELTA_X; // Spatial grid size
int M; // Number of grid nodes for membrane solution
double *w_previous, *w, *w_next, *w_deriv; // Membrane position arrays
double *p_previous_arr, *p_arr, *p_next_arr; // Pressure arrays

/* Global variables */
int impact = 0; // Sets to 1 once impact has happened
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 0; // Records how many GFS files have been outputted
int log_output_no = 0; // Records how many log outputs files there have been
int membrane_output_no = 0; // Records how many membrane outputs there have been
int interface_output_no = 0; // Records how many interface files there have been
int start_membrane = 0; // Boolean to indicate if membrane motion has started
double pinch_off_time = 0.; // Time pinch-off of the entrapped bubble occurs
double drop_thresh = 1e-4; // Remove droplets threshold
double bubble_area = 0.; // Area of entrapped bubble
char interface_time_filename[80] \
    = "interface_times.txt"; // Stores the time the interface was outputted

/* Stats output */
FILE * fp_stats; 

/* Cutoff heights */
// vector h_test[]; // Height function field

/* Contact angle variables */ 
vector h[]; // Height function
double theta0 = 90; // Contact angle in degrees

/* Bubble counting */
scalar bubbles[]; // Tag field for bubbles

/* Boundary conditions */
// Conditions for entry from above
u.n[right] = neumann(0.); // Free flow condition
p[right] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the radial direction
u.n[top] = neumann(0.); // Allows outflow through boundary
u.t[top] = dirichlet(0.); // Stationary vertical flow
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions on surface
u.n[left] = dirichlet(0.); // No flow through surface
u.t[left] = dirichlet(0.); // No slip at surface
h.t[left] = contact_angle (theta0*pi/180.); // RC contact angle


// Function for removing droplets away from a specific region
double membrane_bc(double x, double *w_deriv_arr);
void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr);
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);
int main() {
/* Main function to set up the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Determine physical constants */
    REYNOLDS = RHO_L * V * R / MU_L; // Reynolds number of liquid
    WEBER = RHO_L * V * V * R / SIGMA; // Weber number of liquid
    FROUDE = V / sqrt(9.81 * R); // Froude number of liquid
    RHO_R = RHO_G / RHO_L; // Density ratio
    MU_R = MU_G / MU_L; // Viscosity ratio

    /* Set VOF constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    MEMBRANE_REFINED_HEIGHT = 0.01;
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS;
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL);
    DROP_REFINED_WIDTH = 0.04;

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


    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Open stats file */
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-5; // Possion solver tolerance (default 1e-3)

    // Run the simulation
    run();

    // Close stats file
    fclose(fp_stats);
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();


    /* Attempts to refine around the droplet and membrane */
    double refine_height = MIN_CELL_SIZE;
    int adequate_refinement = 0;

    while ((adequate_refinement == 0) && (refine_height < BOX_WIDTH)) {

        // Attempts to refine
        refine(((sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH) \
            && sq(x - DROP_CENTRE) + sq(y) > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
            || (y < MEMBRANE_RADIUS && x <= refine_height)) \
            && level < MAXLEVEL);

        // Loops and check if refinement was successful along boundary
        adequate_refinement = 1;
        foreach_boundary(left) {
            if ((y < MEMBRANE_RADIUS) && (level < MAXLEVEL)) {
                adequate_refinement = 0;
                break;
            }
        }

        // If refinement was unsuccessful, then double the refined height
        if (adequate_refinement == 0) refine_height += MIN_CELL_SIZE;
    }

    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS));

    /* Initialise the droplet velocity downwards */
    foreach() {
        u.x[] = DROP_VEL * f[];
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
    fprintf(stderr, "DROP_RADIUS = %g, DROP_REFINED_WIDTH = %g, MAXLEVEL = %d\n", DROP_RADIUS, DROP_REFINED_WIDTH, MAXLEVEL);
}


event refinement (i++) {
/* Refines the grid where appropriate */

    /* Adapts with respect to velocities and volume fraction */
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-6}, 
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refine above the membane */
    refine((y < MEMBRANE_RADIUS) && (x <= MEMBRANE_REFINED_HEIGHT) \
        && level < MAXLEVEL);
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
        foreach_boundary(left) {
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
    if (CUTOFF && impact) heights(f, h); // Associates h_test with the heights of f

    // Iterates over solid
    foreach_boundary(left) {

        // Skip if y is not above the membrane
        if (y > MEMBRANE_RADIUS) continue;

        // Determines index in array
        double fractpart, intpart;
        fractpart = modf((y - 0.5 * MIN_CELL_SIZE) / DELTA_X, &intpart);
        if (fractpart > 1e-6) continue;
        
        int k = (int) intpart;

        /* Initially fills p_next_arr[k] with the current value of the sum of
        the pressure and viscous stress */
        // Viscosity average in the cell above the plate
        double avg_mu = f[] * (mu1 - mu2) + mu2;

        // Viscous stress in the cell above the plate
        double viscous_stress = \
            - 2 * avg_mu * (u.x[1, 0] - u.x[]) / Delta;

        p_next_arr[k] = p[] - viscous_stress;

        // If we are pre-impact or are not doing a cutoff, then continue in loop
        if ((CUTOFF == 0) || (impact == 0)) continue;

        /* Implements cutoff */
        if ((f[] > 0) && (f[] < 1)) {
            /* If a mixed cell, output 0 */
            p_next_arr[k] = 0.0;
        } else {
            /* Else, f[] == 0 or 1, and we check the heights */

            // Check the vertical height, if defined
            if (h.x[] != nodata) {
                if (height(h.x[]) - 0.5 < x_min_height) {
                    // If the interface position is less than min_height
                    // p_scale = 0;
                    p_next_arr[k] = 0.0;
                }
            }

            // If too close to the edges, then skip the x cutoff
            if ((k < y_min_height) && (k > M - y_min_height)) continue;

            // x cutoff procedure
            #pragma omp parallel for
            for (int q = -y_min_height; q <= y_min_height; q++) {
                // Does not cutoff a non-mixed cell
                if ((f[0, q] == 1) || (f[0, q] == 0)) continue;
                
                // If a mixed cell is detected then set p_next_arr[k] = 0
                p_next_arr[k] = 0.0;
            }
        }
    }

    /* If constant acceleration, then use imposed solution to set w_next */
    if (CONST_ACC) {
        // Start motion after impact time
        if (t >= IMPACT_TIME) {
            // Loop over membrane
            #pragma omp parallel for
            for (int k = 0; k < M; k++) {
                double y = k * DELTA_X;

                /* Set w(x, t) = 0.5 * a * (t - impact_time)^2 * cos(pi / (2L)) */
                w_next[k] = 0.5 * MEMBRANE_ACC * sq(t + DELTA_T - IMPACT_TIME) \
                    * cos(pi * y / (2 * MEMBRANE_RADIUS)); 
                w_deriv[k] = MEMBRANE_ACC * (t - IMPACT_TIME) \
                    * cos(pi * y / (2 * MEMBRANE_RADIUS));
            }
        }

        // Set boundary condition on membrane
        u.n[left] = dirichlet(membrane_bc(y, w_deriv)); 
    } else if (IMPOSED) {
        // Start motion after impact time
        if (t >= IMPACT_TIME) {
            double tShift = t - IMPACT_TIME + DELTA_T;
            double lambda = 12.;

            double w_coeff = IMPOSED_COEFF * (1 - cos(lambda * tShift) + sq(tShift));
            double w_t_coeff = IMPOSED_COEFF * (2 * (tShift - DELTA_T) \
                + lambda * sin(lambda * (tShift - DELTA_T)));

            // Loop over membrane
            #pragma omp parallel for
            for (int k = 0; k < M; k++) {
                double y = k * DELTA_X;

                w_next[k] = w_coeff * cos(pi * y / (2 * MEMBRANE_RADIUS)); 
                w_deriv[k] = w_t_coeff * cos(pi * y / (2 * MEMBRANE_RADIUS));
            }
        }

        // Set boundary condition on membrane
        u.n[left] = dirichlet(membrane_bc(y, w_deriv)); 


    } else if (t >= MEMBRANE_START_TIME) {
        /* Else solve membrane equation PDE */

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
            u.n[left] = dirichlet(membrane_bc(y, w_deriv)); 
        }
    }

    /* Outputs membrane arrays */
    output_arrays(w, w_deriv, p_arr);

    // Swaps membrane arrays 
    double *temp2 = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp2;
}


event acceleration (i++) {
/* Adds acceleration due to gravity and the moving plate at each time step */
    face vector av = a; // Acceleration at each face

    /* Adds acceleration due to gravity and the plate */
    foreach_face(x){
        av.x[] += - 1./sq(FROUDE);
    }
}


event small_droplet_removal (t += 1e-3) { 
/* Removes any small droplets or bubbles that have formed, that are smaller than
 a specific size. Uses the remove_droplets_region code to leave the area near 
 the point of impact alone in order to properly resolve the entrapped bubble */

    // Minimum diameter (in cells) a droplet/bubble has to be, else it will be 
    // removed
    int drop_min_cell_width = 4;

    // Region to ignore
    double ignore_region_x_limit = 0.1; 
    double ignore_region_y_limit = 0.1; 
    
    // Counts the number of bubbles there are using the tag function
    foreach() {
        bubbles[] = 1. - f[] > drop_thresh;
    }
    int bubble_no = tag(bubbles);

    if (pinch_off_time == 0.) {
        /* The first time the bubble number is above 1, we define it to be the 
        pinch off time */
        if (bubble_no > 1) {
            pinch_off_time = t;
        }
    } else {
        /* Determine area of entrapped bubble */

        // We assume that entrapped bubbles are any bubble which is not the #
        // surrounding air, so we determine the tag of the surrounding air, and 
        // then add up the volume of all of the cells which do not have that 
        // tag

        // We initialise the tag to be 0, which is the tag the
        // bulk droplet will have, so that if the bubble is not found, the 
        // resulting area will be huge and easy to identify that we're wrong.
        int air_tag = 0; 

        // Serial foreach loop along the right boundary to find the tag of the 
        // air, which in theory should end after one iteration
        foreach_boundary(right, serial) {
            if (f[] == 0.) {
                air_tag = bubbles[];
                break;
            }
        }

        // Determine the area of all of the entrapped air cells, which will have
        // a tag not equal to 0 (which will be liquid) or air_tag, which is the
        // tag of the surrounding air
        bubble_area = 0.;
        foreach(reduction(+:bubble_area)) {
            if ((bubbles[] > 0) && (bubbles[] != air_tag)) {
                bubble_area += (1. - f[]) * dv();
            }
        }

        /* After the removal delay, remove drops and bubbles as necessary */    
        if (t >= pinch_off_time + REMOVAL_DELAY) {
            // Set up RemoveDroplets struct
            struct RemoveDroplets remove_struct;
            remove_struct.f = f;
            remove_struct.minsize = drop_min_cell_width;
            remove_struct.threshold = drop_thresh;
            remove_struct.bubbles = false;

            if (t < 0.3) {
                // Remove droplets outside of the specified region
                remove_droplets_region(remove_struct, ignore_region_x_limit, \
                    ignore_region_y_limit);

                // Remove bubbles outside of the specified region
                remove_struct.bubbles = true;
                remove_droplets_region(remove_struct, ignore_region_x_limit, \
                    ignore_region_y_limit);
            } else {
                remove_droplets_region(remove_struct, 0, 0);
                
                remove_struct.bubbles = true;
                remove_droplets_region(remove_struct, 0, 0);
            }

            // Remove the entrapped bubble if specified
            if (REMOVE_ENTRAPMENT) {
                foreach() { 
                    if (x < 0.01 && y < 2 * 0.05) {
                        f[] = 1.;
                    }
                }
            }
        }
    }
        
}


event output_data (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */

    /* Find the pressure at the origin */
    double p0 = 0;
    foreach_boundary(left) {
        if (y < MIN_CELL_SIZE) {
            p0 = p[];
            break;
        }
    }

    /* Outputs data to log file */
    fprintf(stderr, \
        "t = %.5f, v = %.8f, p0 = %8f, bubble_area = %.8f\n", t, \
            2 * pi * statsf(f).sum, p0, bubble_area);
    
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, v = %.8f, p0 = %8f, bubble_area = %.8f\n", t, \
        2 * pi * statsf(f).sum, p0, bubble_area);
    fclose(logfile);

    /* Outputs info about cells along membrane */
    char output_filename[80];
    sprintf(output_filename, "boundary_output_%d.txt", log_output_no);
    FILE *output_file = fopen(output_filename, "w");

    foreach_boundary(left) {
        /* Determine viscous stress */
        double avg_mu = f[] * (mu1 - mu2) + mu2;
        double viscous_stress = - 2 * avg_mu * (u.x[1, 0] - u.x[]) / Delta;

        fprintf(output_file, "%g %g %g %g\n", y, p[], viscous_stress, u.x[]);
    }
    fclose(output_file);

    log_output_no++;
}


event output_interface (t += PLATE_OUTPUT_TIMESTEP) {
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


#if MOVIES
event movies (t += 1e-3) {
/* Produces movies using bview */ 
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
    squares("u.y", spread = -1, linear = true, map = cool_warm);
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("horizontal_vel.mp4");


    /* Movie of the vertical velocity */
    clear();
    draw_vof("f", lw = 2);
    squares("u.x", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("vertical_vel.mp4");

    /* Movie of the pressure */
    clear();
    draw_vof("f", lw = 2);
    squares("p", spread = -1, linear = true, map = cool_warm);
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("pressure.mp4");
}
#endif


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

}

double membrane_bc(double y, double *w_deriv_arr) {
/* membrane_bc
Outputs the boundary condition for the vertical face velocity, uf.n, at the 
bottom boundary, which matches the velocity of the membrane 
*/
    if (y >= MEMBRANE_RADIUS) {
        return 0.;
    } else {
        int k = (int) (y / DELTA_X);
        return -w_deriv_arr[k];
    }
}

void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr) {
/* output_membrane
Outputs the position and time derivative of the membrane, as well as the 
pressure along the membrane, into a text file.
*/
   
    char arr_filename[40];
    sprintf(arr_filename, "membrane_arr_%d.txt", membrane_output_no);
    FILE *arr_file = fopen(arr_filename, "w");

    // Outputs from x = 0 to L - dx
    #pragma omp parallel for
    for (int k = 0; k < M; k++) {
        double x = k * DELTA_X;
        fprintf(arr_file, "%g, %g, %g, %g\n", x, w_arr[k], w_deriv_arr[k], p_arr[k]);
    }

    // Outputs x = L, where w and w_deriv = 0
    double x = M * DELTA_X;
    fprintf(arr_file, "%g, %g, %g, %g\n", x, 0.0, 0.0, 0.0);

    fclose(arr_file);

    membrane_output_no++;
}

/* Alternative remove_droplets definition */
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

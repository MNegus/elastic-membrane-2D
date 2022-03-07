/** spuriousTransform.c
 * A version of the spurious.c test from the Basilisk src, instead testing the 
 * transformed frame with a membrane. 
*/

#define MOVING 0 // Moving frame adjustment
#define AMR 1 // Adaptive mesh refinement
#define WALL 1 // Droplet along the wall

#define JACOBI 1

#include <vofi.h>
#include "test_parameters.h" // Includes all defined parameters

// Membrane scalar fields
scalar W[], Wx[], Wxx[];

#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"


/* Physical constants */
double LAPLACE; // Laplace number
double MU; // Viscosity 
double MEMBRANE_RADIUS; // Radius of the membrane
double mag; //  Magnitude of membrane 
double DIAMETER; // Diameter of droplet


/* Computational variables */
double TMAX; // Maximum time to run simulation
int gfs_output_no = 0; // Tracks how many gfs outputs there have been
double xCentre; // x position of centre of drop
double yCentre; // y position of centre of drop (in lab frame)
double DROP_REFINED_WIDTH = 0.02; // Width of refined region around drop
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
double DC = 1e-12; // Change in domain tolerance


/* Field definitions */
vector htest[]; // Height function for droplet interface
scalar c[], origc[], * interfaces = {c, origc}; // Volume fraction of droplet 
scalar kappa[], kappax[], kappay[]; // Curvature fields
scalar avX[], avY[]; // Acceleration in x and y direction

/* File names */
FILE * fp = NULL; // Output file 
FILE * fp_stats; // Stats file

/* Boundary conditions (NEUMANN) */
// #if WALL    
// u.n[left] = dirichlet(0.); // No flow in the x direction along boundary
// #else
// u.n[left] = neumann(0.);
// #endif

// // Zero Neumann conditions at far-field boundaries
// u.n[bottom] = neumann(0.);
// u.n[top] = neumann(0.);
// u.n[right] = neumann(0.);

/* Boundary conditions (DIRICHLET) */
u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);


double membrane_position(double x) {
/* Continuous function for the membrane position */
    #if MOVING
    if (x <= MEMBRANE_RADIUS) {
        // return mag * (1 - x * x / sq(MEMBRANE_RADIUS));
        return mag * x * x;
    } else {
        return 0.;
    }
    #else
    return 0.;
    #endif
}
}


double membrane_first_derivative(double x) {
/* Continuous function for the first derivative of the membrane position */
    #if MOVING
    if (x <= MEMBRANE_RADIUS) {
        // return mag * (-2 * x / sq(MEMBRANE_RADIUS));
        return 2 * mag * x;
    } else {
        return 0.;
    }
    #else
    return 0.;
    #endif
}


double membrane_second_derivative(double x) {
/* Continuous function for the second derivative of the membrane position */
    #if MOVING
    if (x <= MEMBRANE_RADIUS) {
        // return mag * (-2 / sq(MEMBRANE_RADIUS));
        return 2 * mag;
    } else {
        return 0.;
    }
    #else
    return 0.;
    #endif
}


double x_derivative(Point point, scalar q) {
/* Centred difference for x derivative */
    return (q[1, 0] - q[-1, 0]) / (2. * Delta);
}


double y_derivative(Point point, scalar q) {
/* Centred difference for y derivative */
    return (q[0, 1] - q[0, -1]) / (2. * Delta);
}


double xx_derivative(Point point, scalar q) {
/* Second difference for xx derivative */
    return (q[1, 0] - 2 * q[] + q[-1, 0]) / (Delta * Delta);
}


double yy_derivative(Point point, scalar q) {
/* Second difference for yy derivative */
    return (q[0, 1] - 2 * q[] + q[0, -1]) / (Delta * Delta);
}


double xy_derivative(Point point, scalar q) {
/* Centred difference for xy derivative */
    return (q[1, 1] - q[-1, 1] - q[1, -1] + q[-1, -1]) / (4. * Delta * Delta);
}



static double droplet_phi (creal xy[2]) {
/* Level-set function for the initial position of the droplet, where xy is an 
array with xy[0] = x and xy[1] = y */
  return sq(xy[0] - xCentre) + sq(xy[1] - membrane_position(xy[0]) - yCentre) \
    - sq(DROP_RADIUS);
}


static void vofi (scalar c, int levelmax) {
/* vofi function to define the free surface of the droplet */
  double fh = Get_fh (droplet_phi, NULL, BOX_WIDTH / (1 << levelmax), dimension, 0);
  foreach() {
    creal xy[2] = {x - Delta/2., y - Delta/2.};
    c[] = Get_cc (droplet_phi, xy, Delta, fh, dimension);
  }
}


int main() {

    /* Optional: Set stokes = true to just solve the Stokes equation */
    // stokes = true;  

    /* Sets quantities depending on if the droplet is at the wall or not */
    #if WALL
        xCentre = 0.0;
        yCentre = 3.0;
        MEMBRANE_RADIUS = 1.5;
        mag = 0.5;
    #else
        xCentre = 3.0;
        yCentre = 3.0;
        MEMBRANE_RADIUS = 5.;
        mag = 3.;
    #endif

    /* Determine the physical constants */
    LAPLACE =  2.5 * 12000;

    /* Set dimensionless constants */
    c.sigma = 1.;
    DIAMETER = 2.0 * DROP_RADIUS;
    MU = sqrt(DIAMETER/LAPLACE);
    TMAX = sq(DIAMETER) / MU;

    /* Poisson solver constants */
    TOLERANCE = 1e-6;
    
    /* Set up grid */
    #if AMR
        init_grid(1 << MINLEVEL);
    #else
        init_grid(1 << MAXLEVEL);
    #endif
    size(BOX_WIDTH); // Size of the domain

    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Open stats file */
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");

    /* Print variables to a file */
    char varName[80];
    sprintf(varName, "variables.txt");
    FILE * varFile = fopen(varName, "w");
    fprintf(varFile, "%g %g %g %g\n", MU, LAPLACE, DIAMETER, TMAX);
    fclose(varFile);


    /* Runs the simulation */
    run();
}

/**
We allocate a field to store the previous volume fraction field (to
check for stationary solutions). */

scalar cn[];

event init (i = 0) {

    /* Set constant viscosity field */
    const face vector muc[] = {MU,MU};
    mu = muc;

    /* Open a new file to store the evolution of the amplitude of
    spurious currents  */
    char name[80];
    sprintf (name, "amplitudes_%d.txt", MAXLEVEL);
    if (fp)
        fclose (fp);
    fp = fopen (name, "w");
  
    /* Define the volume fraction and the curvature */
    vofi (c, MINLEVEL);
    for (int l = MINLEVEL + 1; l <= MAXLEVEL; l++) {
        refine (c[] > 0. && c[] < 1. && level < l);
        vofi (c, l);
    }

    /* If using adaptive refinement, initialise with a region refined close to
    the droplet interface */
    #if AMR
    refine((sq(x - xCentre) + sq(y - membrane_position(x) - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - xCentre) + sq(y - membrane_position(x) - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    #endif

    /* Volume fraction for determining change in c */
    foreach()
        cn[] = c[];
    boundary ({cn});
}


event logfile (i++; t <= TMAX) {
    /* At every timestep, we check whether the volume fraction field has
    converged. */
    double dc = change (c, cn);
    if (i > 1 && dc < DC)
    return 1; /* stop */

    /* Output the evolution of the maximum velocity */
    scalar un[];
    foreach()
        un[] = norm(u);

    /* Print out to the amplitudes file */
    fprintf (fp, "%g %g %g %g %.8f\n",
        t, normf(un).max, normf(un).rms, dc, statsf(c).sum);
    fprintf (stderr, "%g %g %g %g %.8f\n",
        t, normf(un).max, normf(un).rms, dc, statsf(c).sum);
    
}


event logstats (i++) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}


event gfsOutput(t += 1) {
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);
    gfs_output_no++;
}


event error (t = end) {
  
  /**
  At the end of the simulation, we compute the equivalent radius of
  the droplet. */

  double vol = statsf(c).sum;
  double radius = sqrt(4.*vol/pi);

  /**
  We recompute the reference solution. */
  
  scalar cref[];
  fraction (cref, sq(DIAMETER/2) - sq(x) - sq(y));
  
  /**
  And compute the maximum error on the curvature *ekmax*, the norm of
  the velocity *un* and the shape error *ec*. */
  
  double ekmax = 0.;
  scalar un[], ec[], kappa[];
  curvature (c, kappa);
  foreach() {
    un[] = norm(u);
    ec[] = c[] - cref[];
    if (kappa[] != nodata) {
      double ek = fabs (kappa[] - (/*AXI*/ + 1.)/radius);
      if (ek > ekmax)
	ekmax = ek;
    }
  }
  
  /**
  We output these on standard error (i.e. the *log* file). */

  norm ne = normf (ec);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", 
	   MAXLEVEL, LAPLACE, 
	   normf(un).max*sqrt(DIAMETER), 
	   ne.avg, ne.rms, ne.max,
	   ekmax);
}


#if AMR
event refinement (i++) {
/* Adaptive grid refinement */
    // Adapts with respect to velocities and volume fraction 
    // adapt_wavelet ({u.x, u.y, f}, (double[]){1e-3, 1e-3, 1e-3},
    //     minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    refine((sq(x - xCentre) + sq(y - membrane_position(x) - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - xCentre) + sq(y - membrane_position(x) - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));


}
#endif

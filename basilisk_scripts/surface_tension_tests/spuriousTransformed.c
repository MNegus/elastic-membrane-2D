/** spuriousTransform.c
 * A version of the spurious.c test from the Basilisk src, instead testing the 
 * transformed frame with a membrane. 
*/

#define MOVING 0 // Moving frame adjustment
#define MEMBRANE 1 // Impose a membrane to deform the droplet
#define AMR 1 // Adaptive mesh refinement
#define WALL 0 // Droplet along the wall

#define TRANSPOSED 1 // Transposes so the membrane is along y
#define SINGLESTEP 1 // If 1, only performs one timestep 
#define REFINEMENTSTUDY 1 // Runs at multiple levels for a refinement study

#define JACOBI 1

#include <vofi.h>
#include "test_parameters.h" // Includes all defined parameters

// Membrane scalar fields
#if TRANSPOSED
scalar W[], Wy[], Wyy[];
#else
scalar W[], Wx[], Wxx[];
#endif

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
int gfs_output_no; // Tracks how many gfs outputs there have been
int interface_output_no; // Tracks how many interface outputs there have been
double xCentre; // x position of centre of drop
double yCentre; // y position of centre of drop (in lab frame)
double DROP_REFINED_WIDTH = 0.5; // Width of refined region around drop
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
double DC = 1e-12; // Change in domain tolerance
int refineLevel; // Variable level of refinement


/* Field definitions */
vector htest[]; // Height function for droplet interface
scalar c[], origc[], * interfaces = {c, origc}; // Volume fraction of droplet 
scalar cn[]; // Field to check for change of c
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
    #if MEMBRANE
    if (x <= MEMBRANE_RADIUS) {
        return mag * (1 - x * x / sq(MEMBRANE_RADIUS));
        // return mag * x * x;
    } else {
        return 0.;
    }
    #else
    return 0.;
    #endif
}


double membrane_first_derivative(double x) {
/* Continuous function for the first derivative of the membrane position */
    #if MEMBRANE
    if (x <= MEMBRANE_RADIUS) {
        return mag * (-2 * x / sq(MEMBRANE_RADIUS));
        // return 2 * mag * x;
    } else {
        return 0.;
    }
    #else
    return 0.;
    #endif
}


double membrane_second_derivative(double x) {
/* Continuous function for the second derivative of the membrane position */
    #if MEMBRANE
    if (x <= MEMBRANE_RADIUS) {
        return mag * (-2 / sq(MEMBRANE_RADIUS));
        // return 2 * mag;
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
    #if TRANSPOSED
    return sq(xy[0] - membrane_position(xy[1]) - xCentre) + sq(xy[1] - yCentre) \
        - sq(DROP_RADIUS);
    #else
    return sq(xy[0] - xCentre) + sq(xy[1] - membrane_position(xy[0]) - yCentre) \
    - sq(DROP_RADIUS);
    #endif
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
    #if TRANSPOSED
    xCentre = 0.5 * BOX_WIDTH;
    yCentre = 0.;
    #else
    xCentre = 0.0;
    yCentre = 0.5 * BOX_WIDTH;
    #endif // TRANSPOSED
    MEMBRANE_RADIUS = 1.5;
    mag = 0.5;
    #else
    xCentre = 3.0;
    yCentre = 3.0;
    MEMBRANE_RADIUS = 5.;
    mag = 0.5;
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
    #if REFINEMENTSTUDY
    for (refineLevel = MINLEVEL; refineLevel <= MAXLEVEL; refineLevel++)
        run();
    #else
    refineLevel = MAXLEVEL;
    run();
    #endif
}


event init (i = 0) {

    /* Set constant viscosity field */
    const face vector muc[] = {MU,MU};
    mu = muc;

    /* Open a new file to store the evolution of the amplitude of
    spurious currents  */
    char name[80];
    sprintf (name, "amplitudes_%d.txt", refineLevel);
    if (fp)
        fclose (fp);
    fp = fopen (name, "w");

    /* Initialise global variables */
    gfs_output_no = 0;
    interface_output_no = 0;
  
    /* Define the volume fraction and the curvature */
    vofi (c, MINLEVEL);
    for (int l = MINLEVEL + 1; l <= refineLevel; l++) {
        refine (c[] > 0. && c[] < 1. && level < l);
        vofi (c, l);
    }

    /* If using adaptive refinement, initialise with a region refined close to
    the droplet interface */
    #if AMR
    #if TRANSPOSED
    refine((sq(x - membrane_position(y) - xCentre) + sq(y  - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - membrane_position(y) - xCentre) + sq(y  - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < refineLevel));
    #else
    refine((sq(x - xCentre) + sq(y - membrane_position(x) - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - xCentre) + sq(y - membrane_position(x) - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < refineLevel));
    #endif
    #endif

    /* Set fields for membrane position */
    foreach() { 
        #if TRANSPOSED
        W[] = membrane_position(y);
        Wy[] = membrane_first_derivative(y);
        Wyy[] = membrane_second_derivative(y);
        #else
        W[] = membrane_position(x);
        Wx[] = membrane_first_derivative(x);
        Wxx[] = membrane_second_derivative(x);
        #endif
    }

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


event output_data(t += 1.0) {
    /* Set c to be f */
    // scalar c = f;

    /* Determine the curvature and heights */
    heights(c, htest);
    cstats s = curvature (c, kappa, sigma=1., add = false);
    foreach() { 
        if (c[] == 0 || c[] == 1) {
            kappax[] = nodata;
            kappay[] = nodata;
        } else {
            double kappaxVal = kappa_x(point, htest);
            double kappayVal = kappa_y(point, htest);

            if (fabs(kappaxVal) > 1e3) {
                kappax[] = nodata;
            } else {
                kappax[] = fabs(kappaxVal);
            }

            if (fabs(kappayVal) > 1e3) {
                kappay[] = nodata;
            } else {
                kappay[] = fabs(kappayVal);
            }
        }
    }

    /* Output the interface and curvature along it */
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d_%d.txt", interface_output_no, refineLevel);
    FILE *interface_file = fopen(interface_filename, "w");

    foreach() {
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {

            // Height function derivatives
            double hx = (htest.y[1, 0] - htest.y[-1, 0])/2.;
            double hxx = (htest.y[1, 0] + htest.y[-1, 0] - 2.*htest.y[])/Delta;
            double hy = (htest.x[0, 1] - htest.x[0, -1])/2.;
            double hyy = (htest.x[0, 1] + htest.x[0, -1] - 2.*htest.x[])/Delta;

            // W values
            #if TRANSPOSED
            double Wval = W[];
            double Wxval = Wy[];
            double Wxxval = Wyy[];
            #else
            double Wval = W[];
            double Wxval = Wx[];
            double Wxxval = Wxx[];
            #endif

            // Segment info
            coord n = interface_normal(point, c);
            double alpha = plane_alpha(c[], n);
            coord segment[2];
            if (facets(n, alpha, segment) == 2) {
                fprintf(interface_file, \
                    "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
                    x + segment[0].x * Delta, y + segment[0].y * Delta,
                    x + segment[1].x * Delta, y + segment[1].y * Delta,
                    kappa[], kappax[], kappay[],
                    htest.y[], hx, hxx, 
                    htest.x[], hy, hyy,
                    Wval, Wxval, Wxxval, 
                    x, y);
            }
        }
    }
    fclose(interface_file);

    interface_output_no++;

}


event gfsOutput(t += 1.0) {
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d_%d.gfs", gfs_output_no, refineLevel);
    output_gfs(file = gfs_filename);
    gfs_output_no++;
}


#if SINGLESTEP
event end (i = 0) {
    fprintf(stderr, "Finished refineLevel %d\n", refineLevel);
    fflush(stderr);
    return 1;
}
#endif

#if AMR
event refinement (i++) {
/* Adaptive grid refinement */
    // Adapts with respect to velocities and volume fraction 
    // adapt_wavelet ({u.x, u.y, c}, (double[]){1e-4, 1e-4, 0},
    //     minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    #if TRANSPOSE
    refine((sq(x - membrane_position(y) - xCentre) + sq(y  - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - membrane_position(y) - xCentre) + sq(y  - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    #else
    refine((sq(x - xCentre) + sq(y - membrane_position(x) - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - xCentre) + sq(y - membrane_position(x) - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    #endif


}
#endif

#if MOVING
event accAdjustment(i++) {

    face vector av = a;
    
    // y acceleration
    foreach_face(y) {
        
        double ut = av.x[];

        double ux = x_derivative(point, u.x);
        double uy = y_derivative(point, u.x);
        double uxx = xx_derivative(point, u.x);
        double uyy = yy_derivative(point, u.x);
        double uxy = xy_derivative(point, u.x);

        double v = u.y[];
        double vy = y_derivative(point, u.y);
        double vyy = yy_derivative(point, u.y);
        double vxy = xy_derivative(point, u.y);

        double Wxf = interpolate(Wx, x, y);
        double Wxxf = interpolate(Wxx, x, y);

        av.y[] += Wxf * ut + Wxxf * sq(u.x[]) + Wxf * ux * u.x[] + Wxf * uy * v \
            + (mu.y[] / rho[]) * (Wxxf * vy + 2. * Wxf * vxy + sq(Wxf) * vyy \
                - (2. * Wxxf * ux + Wxf * uxx + Wxf * uyy \
                    + 3. * Wxf * Wxxf * uy + 2. * sq(Wxf) * uxy \
                    + pow(Wxf, 3.) * uyy));
        avY[] = av.y[];
    }

    // x acceleration
    foreach_face(x) {
        double py = y_derivative(point, p);
        
        double uy = y_derivative(point, u.x);
        double uyy = yy_derivative(point, u.x);
        double uxy = xy_derivative(point, u.x);

        double Wxf = interpolate(Wx, x, y);
        double Wxxf = interpolate(Wxx, x, y);

        av.x[] += (1. / rho[]) * (-Wxf * py \
            + mu.x[] * (Wxxf * uy + 2. * Wxf * uxy + Wxf * Wxf * uyy));

        avX[] = av.x[];
    }

}
#endif
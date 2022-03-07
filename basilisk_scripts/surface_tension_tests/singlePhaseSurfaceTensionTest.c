/* surfaceTensionTest.c
    Membrane at the bottom, gravity pointing down, droplet center at x=0, y>0, 
    without AMR (adaptive mesh refinement), without parallelisation
*/

#define MOVING 0 // Moving frame adjustment
#define AMR 0 // Adaptive mesh refinement
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


/* Computational variables */
double TMAX; // Maximum time to run simulation
int gfs_output_no = 0; // Tracks how many gfs outputs there have been
double x_drop_centre; // x position of centre of drop
double y_drop_centre; // y position of centre of drop (in lab frame)
double DROP_REFINED_WIDTH = 0.02; // Width of refined region around drop
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished


/* Field definitions */
face vector av; // Acceleration at each face
vector htest[]; // Height function for droplet interface
scalar c[], origc[], * interfaces = {c, origc}; // Volume fraction of droplet 
// scalar origc[]; // Volume fraction of droplet at t = 0
scalar kappa[], kappax[], kappay[]; // Curvature fields
scalar avX[], avY[]; // Acceleration in x and y direction


/* Boundary conditions */
// Left-hand boundary
#if WALL
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary
#else
u.n[left] = neumann(0.); // Neumann condition if the droplet is not at the wall
#endif

// Zero Neumann conditions at far-field boundaries
u.n[bottom] = neumann(0.);
u.n[top] = neumann(0.);
u.n[right] = neumann(0.);

/* File names */
FILE * fp_stats; // Stats file


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
  return sq(xy[0] - x_drop_centre) + sq(xy[1] - membrane_position(xy[0]) - y_drop_centre) \
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
    stokes = true;
    
    /* Sets quantities depending on if the droplet is at the wall or not */
    #if WALL
        x_drop_centre = 0.;
        y_drop_centre = 3.0;
        MEMBRANE_RADIUS = 1.5;
        // mag = 2.;
        mag = 0.5;
    #else
        x_drop_centre = 3.;
        y_drop_centre = 1.5;
        MEMBRANE_RADIUS = 5.;
        mag = 3.;
    #endif

    /* Associates the acceleration with the face vector field av */
    a = av;

    /* Determine physical constants */
    LAPLACE = 3000; // Laplace number

    /* Set dimensionless constants */
    c.sigma = 1.;
    double diameter = 2.;
    MU = sqrt(diameter / LAPLACE);
    TMAX = sq(diameter) / MU;
    fprintf(stderr, "MU = %g, TMAX = %g\n", MU, TMAX);

    #if AMR
        init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    #else
        init_grid(1 << MAXLEVEL); // Create uniform grid at maximum level
    #endif 

    size(BOX_WIDTH); // Size of the domain

    // c.refine = c.prolongation = fraction_refine;
    // c.height = htest;
    
    /* Poisson solver constants */
    // DT = 1.0e-4; // Minimum timestep
    // NITERMIN = 1; // Min number of iterations (default 1)
    // NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Open stats file */
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");


    /* Runs the simulation */
    run();
}

event init (i = 0) {

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Define the viscosity */
    const face vector muc[] = {MU,MU};
    mu = muc;

    /* Define the volume fraction and the curvature */
    vofi (c, MINLEVEL);
    for (int l = MINLEVEL + 1; l <= MAXLEVEL; l++) {
        refine (c[] > 0. && c[] < 1. && level < l);
        vofi (c, l);
    }

    #if AMR
    refine((sq(x - x_drop_centre) + sq(y - membrane_position(x) - y_drop_centre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - x_drop_centre) + sq(y - membrane_position(x) - y_drop_centre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    #endif
    
    /* Set fields for membrane position */
    foreach() { 
        W[] = membrane_position(x);
        Wx[] = membrane_first_derivative(x);
        Wxx[] = membrane_second_derivative(x);
        
        origc[] = c[];
    }

}

#if AMR
event refinement (i++) {
/* Adaptive grid refinement */

    // Adapts with respect to velocities and volume fraction 
    // adapt_wavelet ({u.x, u.y, f}, (double[]){1e-3, 1e-3, 1e-3},
    //     minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    refine((sq(x - x_drop_centre) + sq(y - membrane_position(x) - y_drop_centre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
            && (sq(x - x_drop_centre) + sq(y - membrane_position(x) - y_drop_centre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
            && (level < MAXLEVEL));


}
#endif


#if MOVING
event accAdjustment(i++) {
    
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


event logOutput (i++; t <= TMAX) {
    /* Find velocity norm */
    scalar un[];
    foreach()
        un[] = norm(u);

    fprintf(stderr, "t = %.5f, i = %d, v = %.8f, unMax = %g\n", t, i, \
        2 * pi * statsf(c).sum, normf(un).max);

    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, i = %d, v = %.8f\n", t, i, 2 * pi * statsf(c).sum);
    fclose(logfile);
}

event logstats (i++) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

event gfsOutput(i += 1000) {
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);
    gfs_output_no++;
}


// event output_data(i++) {
//     /* Set c to be f */
//     // scalar c = f;

//     /* Determine the curvature and heights */
//     heights(c, htest);
//     cstats s = curvature (c, kappa, sigma=1., add = false);
//     foreach() { 
//         if (c[] == 0 || c[] == 1) {
//             kappax[] = nodata;
//             kappay[] = nodata;
//         } else {
//             double kappaxVal = kappa_x(point, htest);
//             double kappayVal = kappa_y(point, htest);

//             if (fabs(kappaxVal) > 1e3) {
//                 kappax[] = nodata;
//             } else {
//                 kappax[] = fabs(kappaxVal);
//             }

//             if (fabs(kappayVal) > 1e3) {
//                 kappay[] = nodata;
//             } else {
//                 kappay[] = fabs(kappayVal);
//             }
//         }
//     }

//     /* Output gfs and curvature files */
//     char gfs_filename[80];
//     sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
//     output_gfs(file = gfs_filename);


//     /* Output the interface and curvature along it */
//     char interface_filename[80];
//     sprintf(interface_filename, "interface_%d.txt", gfs_output_no);
//     FILE *interface_file = fopen(interface_filename, "w");

//     foreach() {
//         if (c[] > 1e-6 && c[] < 1. - 1e-6) {

//             // Height function derivatives
//             double hx = (htest.y[1, 0] - htest.y[-1, 0])/2.;
//             double hxx = (htest.y[1, 0] + htest.y[-1, 0] - 2.*htest.y[])/Delta;
//             double hy = (htest.x[0, 1] - htest.x[0, -1])/2.;
//             double hyy = (htest.x[0, 1] + htest.x[0, -1] - 2.*htest.x[])/Delta;

//             // Segment info
//             coord n = interface_normal(point, c);
//             double alpha = plane_alpha(c[], n);
//             coord segment[2];
//             if (facets(n, alpha, segment) == 2) {
//                 fprintf(interface_file, \
//                     "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
//                     x + segment[0].x * Delta, y + segment[0].y * Delta,
//                     x + segment[1].x * Delta, y + segment[1].y * Delta,
//                     kappa[], kappax[], kappay[],
//                     htest.y[], hx, hxx, 
//                     htest.x[], hy, hyy,
//                     W[], Wx[], Wxx[], 
//                     x, y);
//             }
//         }
//     }
//     fclose(interface_file);

//     gfs_output_no++;

// }

event end (t = end) {
    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "Finished after %g seconds\n", end_wall_time - start_wall_time);
    fclose(logfile);

}
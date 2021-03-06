/* surfaceTensionTest.c
    Membrane at the bottom, gravity pointing down, droplet center at x=0, y>0, 
    without AMR (adaptive mesh refinement), without parallelisation
*/

#define MOVING 0 // Moving frame adjustment
#define AMR 1 // Adaptive mesh refinement
#define WALL 1 // Droplet along the wall

#include <vofi.h>
#include "test_parameters.h" // Includes all defined parameters

scalar W[], Wx[], Wxx[];

#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h"
#include "tension.h"

face vector av; // Acceleration at each face

double mag;
int gfs_output_no = 0;
double x_drop_centre;
double y_drop_centre;
double MEMBRANE_RADIUS;
double DROP_REFINED_WIDTH = 0.02;
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
FILE * fp_stats; 

/* Physical constants */
double REYNOLDS; // Reynolds number of liquid
double WEBER; // Weber number of liquid
double FROUDE; // Froude number of liquid
double RHO_R; // Density ratio
double MU_R; // Viscosity ratio

// Symmetry on left boundary (only on wall case)
#if WALL
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary
#else
u.n[left] = neumann(0.); // Neumann condition if the droplet is not at the wall
#endif

// Zero Neumann conditions at far-field boundaries
u.n[bottom] = neumann(0.);
u.n[top] = neumann(0.);
u.n[right] = neumann(0.);

vector htest[];
scalar c[];
scalar kappa[], kappax[], kappay[];
scalar origf[];
scalar avX[], avY[];

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
    return (q[1, 0] - q[-1, 0]) / (2. * Delta);
}

double y_derivative(Point point, scalar q) {
    return (q[0, 1] - q[0, -1]) / (2. * Delta);
}

double xx_derivative(Point point, scalar q) {
    return (q[1, 0] - 2 * q[] + q[-1, 0]) / (Delta * Delta);
}

double yy_derivative(Point point, scalar q) {
    return (q[0, 1] - 2 * q[] + q[0, -1]) / (Delta * Delta);
}

double xy_derivative(Point point, scalar q) {
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

    #if WALL
        x_drop_centre = 0.;
        y_drop_centre = 1.5;
        MEMBRANE_RADIUS = 1.5;
        // mag = 2.;
        mag = 0.5;
    #else
        x_drop_centre = 3.;
        y_drop_centre = 1.5;
        MEMBRANE_RADIUS = 5.;
        mag = 3.;
    #endif

    a = av;

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

    #if AMR
        init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    #else
        init_grid(1 << MAXLEVEL);
    #endif 

    size(BOX_WIDTH); // Size of the domain

    // c.refine = c.prolongation = fraction_refine;
    // c.height = htest;
    
    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    // NITERMIN = 1; // Min number of iterations (default 1)
    // NITERMAX = 300; // Max number of iterations (default 100)
    // TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

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

event init (t = 0) {

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Define the volume fraction and the curvature */
    vofi (f, MINLEVEL);
    for (int l = MINLEVEL + 1; l <= MAXLEVEL; l++) {
        refine (f[] > 0. && f[] < 1. && level < l);
        vofi (f, l);
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
        
        origf[] = f[];
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


event logOutput (t += 1e-3) {
    fprintf(stderr, "t = %.5f, i = %d, v = %.8f\n", t, i, 2 * pi * statsf(f).sum);

    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, i = %d, v = %.8f\n", t, i, 2 * pi * statsf(f).sum);
    fclose(logfile);
}

event logstats (t += 1e-3) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}


event output_data(t += 1e-2) {
    /* Set c to be f */
    scalar c = f;

    /* Determine the curvature and heights */
    heights(c, htest);
    cstats s = curvature (c, kappa, 1. / WEBER, add = false);
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
                kappax[] = fabs(kappaxVal) / WEBER;
            }

            if (fabs(kappayVal) > 1e3) {
                kappay[] = nodata;
            } else {
                kappay[] = fabs(kappayVal) / WEBER;
            }
        }
    }

    /* Output gfs and curvature files */
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);


    /* Output the interface and curvature along it */
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d.txt", gfs_output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    foreach() {
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {

            // Height function derivatives
            double hx = (htest.y[1, 0] - htest.y[-1, 0])/2.;
            double hxx = (htest.y[1, 0] + htest.y[-1, 0] - 2.*htest.y[])/Delta;
            double hy = (htest.x[0, 1] - htest.x[0, -1])/2.;
            double hyy = (htest.x[0, 1] + htest.x[0, -1] - 2.*htest.x[])/Delta;

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
                    W[], Wx[], Wxx[], 
                    x, y);
            }
        }
    }
    fclose(interface_file);

    gfs_output_no++;

}

event end (t = MAX_TIME) {
    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "Finished after %g seconds\n", end_wall_time - start_wall_time);
    fclose(logfile);

}
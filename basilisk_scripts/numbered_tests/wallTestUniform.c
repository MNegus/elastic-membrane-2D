/* test_I1.c  
    Membrane at the bottom, gravity pointing down, droplet center at x=0, y>0, 
    without AMR (adaptive mesh refinement), without parallelisation
*/

#define MOVING 1

#include <vofi.h>
#include "test_parameters.h" // Includes all defined parameters

scalar W[], Wx[], Wxx[];

#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h"
#include "tension.h"

double mag = 1.0;
int gfs_output_no = 0;
double x_drop_centre;
double y_drop_centre;
double DROP_REFINED_WIDTH = 0.05;
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
FILE * fp_stats; 

/* Physical constants */
double REYNOLDS; // Reynolds number of liquid
double WEBER; // Weber number of liquid
double FROUDE; // Froude number of liquid
double RHO_R; // Density ratio
double MU_R; // Viscosity ratio

// Symmetry on left boundary
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary

// Conditions on surface
uf.n[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);

vector htest[];
scalar c[];
scalar kappa[], kappax[], kappay[];
scalar origf[];

double membrane_position(double x) {
/* Continuous function for the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return mag * (1 - x * x / sq(MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}

double membrane_first_derivative(double x) {
/* Continuous function for the first derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return mag * (-2 * x / sq(MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}

double membrane_second_derivative(double x) {
/* Continuous function for the second derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return mag * (-2 / sq(MEMBRANE_RADIUS));
    } else {
        return 0;
    }
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

    x_drop_centre = 0.0 * BOX_WIDTH;
    y_drop_centre = 0.25 * BOX_WIDTH;

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

    init_grid(1 << MAXLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    // c.refine = c.prolongation = fraction_refine;
    // c.height = htest;

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
    
    /* Set fields for membrane position */
    foreach() { 
        W[] = membrane_position(x);
        Wx[] = membrane_first_derivative(x);
        Wxx[] = membrane_second_derivative(x);
        
        origf[] = f[];
    }

}

#if MOVING
event accAdjustment(i++) {
    face vector av = a; // Acceleration at each face
    
    // y acceleration
    foreach_face(y) {
        double v = uf.y[];
        double ut = av.x[];
        double ux = (uf.x[1, 0] - uf.x[-1, 0]) / (2. * Delta);
        double uy = (uf.x[0, 1] - uf.x[0, -1]) / (2. * Delta);
        double uxx = (uf.x[1, 0] - uf.x[] + uf.x[-1, 0]) / (Delta * Delta);
        double uyy = (uf.x[0, 1] - uf.x[] + uf.x[0, -1]) / (Delta * Delta);
        double uxy = (uf.x[1, 1] - uf.x[-1, 1] - uf.x[1, -1] + uf.x[-1, -1]) / (4. * Delta * Delta);
        double vxy = (uf.y[1, 1] - uf.y[-1, 1] - uf.y[1, -1] + uf.y[-1, -1]) / (4. * Delta * Delta);
        double vyy = (uf.y[0, 1] - uf.y[] + uf.y[0, -1]) / (Delta * Delta);

        double Wxf = interpolate(Wx, x, y);

        av.y[] += Wxf * ut + Wxf * ux * u.x[] + Wxf * uy * v \
            + (mu.y[] / rho[]) * (2. * Wxf * vxy + Wxf * Wxf * vyy \
                - (Wxf * uxx + Wxf * uyy + 2. * Wxf * Wxf * uxy \
                    + Wxf * Wxf * Wxf * uyy));
    }

    // x acceleration
    foreach_face(x) {
        double py = (p[0, 1] - p[0, -1]) / (2. * Delta);
        double uyy = (uf.x[0, 1] - uf.x[] + uf.x[0, -1]) / (Delta * Delta);
        double uxy = (uf.x[1, 1] - uf.x[-1, 1] - uf.x[1, -1] + uf.x[-1, -1]) / (4. * Delta * Delta);

        double Wxf = interpolate(Wx, x, y);
        av.x[] += (1. / rho[]) * (-Wxf * py + mu.x[] * (2. * Wxf * uxy + Wxf * Wxf * uyy));
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


event output_data(t += 1e-1) {
    /* Set c to be f */
    foreach() {
        c[] = f[];
    }

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
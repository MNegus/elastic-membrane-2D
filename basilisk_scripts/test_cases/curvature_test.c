/* curvature_test.c
    Tests the implementation of the moving frame curvature. We initialise a 
    stationary droplet with a specific membrane displacement. We do not include
    any hydrodynamics, and just check whether the curvature calculation of the
    droplet matches the stationary case. 

    Similar to the curvature.c test in the Basilisk source code, we use vofi to
    initialise the droplet, and use the curvature statistics to check the output
    matches the curvature. This time, we keep the droplet position and radius
    the same, but still vary the refinement level and perhaps the membrane 
    displacment. 

*/

#define MOVING 1

#include <vofi.h>
#include "parameters.h" // Includes all defined parameters


scalar W[], Wx[], Wxx[];
vector htest[];

#include "fractions.h"
#include "curvature.h"

double mag;
double magMax = 0.5; // Magnitude of membrane displacement
int gfs_output_no = 1;
double drop_centre;
double DROP_REFINED_WIDTH = 0.01;




double membrane_position(double x) {
/* Continuous function for the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return mag * cos(pi * x / (2 * MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}

double membrane_first_derivative(double x) {
/* Continuous function for the first derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return -mag * (pi / (2 * MEMBRANE_RADIUS)) \
            * sin(pi * x / (2 * MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}

double membrane_second_derivative(double x) {
/* Continuous function for the second derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        return -mag * pow(pi / (2 * MEMBRANE_RADIUS), 2) \
            * cos(pi * x / (2 * MEMBRANE_RADIUS));
    } else {
        return 0;
    }
}


static double droplet_phi (creal xy[2]) {
/* Level-set function for the initial position of the droplet, where xy is an 
array with xy[0] = x and xy[1] = y */
  return sq(xy[0]) + sq(xy[1] - membrane_position(xy[0]) - drop_centre) \
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

    drop_centre = INITIAL_DROP_HEIGHT + DROP_RADIUS;

    int levelmax = MAXLEVEL;

    for (mag = 0; mag <= magMax; mag += magMax) {

            init_grid(1 << MINLEVEL); // Create grid according to the minimum level
            size(BOX_WIDTH); // Size of the domain

            scalar c[], kappa[], kappax[], kappay[];
            c.refine = c.prolongation = fraction_refine;
            c.height = htest;

            /* Define the volume fraction and the curvature */
            vofi (c, MINLEVEL);
            for (int l = MINLEVEL + 1; l <= levelmax; l++) {
                refine (c[] > 0. && c[] < 1. && level < l);
                refine((sq(x) + sq(y - membrane_position(x) - drop_centre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
                    && (sq(x) + sq(y - membrane_position(x) - drop_centre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
                    && (level < l));
                vofi (c, l);
            }
            
            /* Set fields for membrane position */
            foreach() { 
                W[] = membrane_position(x);
                Wx[] = membrane_first_derivative(x);
                Wxx[] = membrane_second_derivative(x);
            }

            /* Determine the curvature and heights */
            heights(c, htest);
            cstats s = curvature (c, kappa, sigma = 2., add = false);
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
                        kappax[] = 2 * fabs(kappaxVal);
                    }

                    if (fabs(kappayVal) > 1e3) {
                        kappay[] = nodata;
                    } else {
                        kappay[] = 2 * fabs(kappayVal);
                    }
                }
            }

            /* Output gfs and curvature files */
            char gfs_filename[80];
            char curvature_filename[80];
            if (mag == 0) {
                sprintf(gfs_filename, "gfs_output_%d_flat.gfs", levelmax);
                sprintf(curvature_filename, "curvature_flat.txt");
            } else {
                sprintf(gfs_filename, "gfs_output_%d_curved.gfs", levelmax);
                sprintf(curvature_filename, "curvature_curved.txt");
            }
            output_gfs(file = gfs_filename);
            
            // Output the curvature
            FILE *curvature_file = fopen(curvature_filename, "w");
            foreach() {
                if (kappa[] != nodata) {
                    fprintf(curvature_file, "%g %g %g\n", x, y, kappa[]);
                }
            }
            fclose(curvature_file);
            
            /* Output the heights */
            char heights_x_filename[80];
            sprintf(heights_x_filename, "heights_x_mag_%g.txt", mag);
            FILE *heights_x_file = fopen(heights_x_filename, "w");

            char heights_y_filename[80];
            sprintf(heights_y_filename, "heights_y_mag_%g.txt", mag);
            FILE *heights_y_file = fopen(heights_y_filename, "w");

            foreach() {
                if (c[] > 0 && c[] < 1) {
                    if (htest.x[] != nodata) {
                        fprintf(heights_x_file, "%g %g %g\n", x, y, htest.x[]);
                    }

                    if (htest.y[] != nodata) {
                        fprintf(heights_y_file, "%g %g %g\n", x, y, htest.y[]);
                    }
                }
            }

            fclose(heights_x_file);
            fclose(heights_y_file);


            /* Output the interface and curvature along it */
            char interface_filename[80];
            sprintf(interface_filename, "interface_mag_%g.txt", mag);
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
                            "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
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


            /* Determine error in curvature */
            // int l = levelmax;
            // foreach (serial) {
            //     // fprintf (stderr, "x = %g, w = %g\n", x, membrane_position(x));
            //     if (c[] > 0. && c[] < 1.) {   
            //         double e = fabs(kappa[]/2. - 1/DROP_RADIUS)*DROP_RADIUS;
            //         n[l].volume += dv();
            //         n[l].avg += dv()*e;
            //         n[l].rms += dv()*e*e;
            //         if (e > n[l].max)
            //             n[l].max = e;
            //     }
            // }
            
            // // Output the error
            // if (n[l].volume) {
            //     n[l].avg /= n[l].volume;
            //     n[l].rms = sqrt(n[l].rms/n[l].volume);
            //     fprintf (stderr, "%g %g %g %g\n",
            //         2.*R*(1 << l),
            //         n[l].avg, n[l].rms, n[l].max);

            //     FILE *logfile = fopen("log", "a");
            //     fprintf(logfile, "%g %g %g %g %g %g %g %g\n",
            //         2.*R*(1 << l),
            //         n[l].avg, n[l].rms, n[l].max,
            //         sc[l].h/t, sc[l].f/t, sc[l].a/t, sc[l].c/t);
            //     fclose(logfile);
            // }
    }

    

    

    

}


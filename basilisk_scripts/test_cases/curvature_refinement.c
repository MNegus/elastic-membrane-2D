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

double mag = 0.5;
int gfs_output_no = 1;
double drop_centre;
double DROP_REFINED_WIDTH = 0.01;
double WEBER;




double membrane_position(double x) {
/* Continuous function for the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        // return mag * cos(pi * x / (2 * MEMBRANE_RADIUS));
        return mag * x * x;
    } else {
        return 0;
    }
}

double membrane_first_derivative(double x) {
/* Continuous function for the first derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        // return -mag * (pi / (2 * MEMBRANE_RADIUS)) \
        //     * sin(pi * x / (2 * MEMBRANE_RADIUS));
        return 2 * mag * x;
    } else {
        return 0;
    }
}

double membrane_second_derivative(double x) {
/* Continuous function for the second derivative of the membrane position */
    if (x <= MEMBRANE_RADIUS) {
        // return -mag * pow(pi / (2 * MEMBRANE_RADIUS), 2) \
        //     * cos(pi * x / (2 * MEMBRANE_RADIUS));
        return 2 * mag;
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

    WEBER = RHO_L * V * V * R / SIGMA;

    for (int refineLevel = MINLEVEL; refineLevel <= MAXLEVEL; refineLevel++) {
        fprintf(stderr, "refineLevel = %d\n", refineLevel);

        init_grid(1 << MINLEVEL); // Create grid according to the minimum level
        size(BOX_WIDTH); // Size of the domain

        scalar c[], kappa[], kappax[], kappay[];
        c.refine = c.prolongation = fraction_refine;
        c.height = htest;

        /* Define the volume fraction and the curvature */
        vofi (c, MINLEVEL);
        for (int l = MINLEVEL + 1; l <= refineLevel; l++) {
            refine (c[] > 0. && c[] < 1. && level < l);
            refine((sq(x) + sq(y - membrane_position(x) - drop_centre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
                && (sq(x) + sq(y - membrane_position(x) - drop_centre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
                && (level < refineLevel));
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
        cstats s = curvature (c, kappa, sigma = 1. / WEBER, add = false);
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
        sprintf(gfs_filename, "gfs_output_%d.gfs", refineLevel);
        output_gfs(file = gfs_filename);


        /* Output the interface and curvature along it */
        char interface_filename[80];
        sprintf(interface_filename, "interface_%d.txt", refineLevel);
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
    }

    

    

    

}


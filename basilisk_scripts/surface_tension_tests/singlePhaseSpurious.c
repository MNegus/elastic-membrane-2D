/**
Spurious currents test version...? */

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
double yCentre = 3.0; // y position of centre of drop (in lab frame)
double DROP_REFINED_WIDTH = 0.02; // Width of refined region around drop
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
double DC = 0.; // Change in domain tolerance


/* Field definitions */
face vector av; // Acceleration at each face
vector htest[]; // Height function for droplet interface
scalar c[], origc[], * interfaces = {c, origc}; // Volume fraction of droplet 
scalar kappa[], kappax[], kappay[]; // Curvature fields
scalar avX[], avY[]; // Acceleration in x and y direction

/* File names */
FILE * fp = NULL; // Output file 
FILE * fp_stats; // Stats file

/* Boundary conditions */
#if WALL    
u.n[left] = dirichlet(0.); // No flow in the x direction along boundary
#else
u.n[left] = neumann(0.);
#endif

// Zero Neumann conditions at far-field boundaries
u.n[bottom] = neumann(0.);
u.n[top] = neumann(0.);
u.n[right] = neumann(0.);


double membrane_position(double x) {
/* Continuous function for the membrane position */
    return 0.;
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
  
    DIAMETER = 2.0 * DROP_RADIUS;
    LAPLACE =  2.5 * 1200;
    MU = sqrt(DIAMETER/LAPLACE);
    TMAX = sq(DIAMETER) / MU;
    
    TOLERANCE = 1e-6;
    // stokes = true;
    c.sigma = 1;
    #if AMR
        init_grid(1 << MINLEVEL);
    #else
        init_grid(1 << MAXLEVEL);
    #endif
    size(BOX_WIDTH); // Size of the domain

    #if WALL
        xCentre = 0.0;
    #else
        xCentre = 3.0;
    #endif
    

    run();
}

/**
We allocate a field to store the previous volume fraction field (to
check for stationary solutions). */

scalar cn[];

event init (i = 0) {

  /**
  We set the constant viscosity field... */

  const face vector muc[] = {MU,MU};
  mu = muc;

  /**
  ... open a new file to store the evolution of the amplitude of
  spurious currents for the various LAPLACE, LEVEL combinations... */

  char name[80];
  sprintf (name, "La-%g-%d", LAPLACE, MAXLEVEL);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");

  /**
  ... and initialise the shape of the interface and the initial volume
  fraction field. */
  
//   fraction (c, sq(DIAMETER/2) - sq(x) - sq(y - yCentre));
    /* Define the volume fraction and the curvature */
    vofi (c, MINLEVEL);
    for (int l = MINLEVEL + 1; l <= MAXLEVEL; l++) {
        refine (c[] > 0. && c[] < 1. && level < l);
        vofi (c, l);
    }

    #if AMR
    refine((sq(x - xCentre) + sq(y - membrane_position(x) - yCentre) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - xCentre) + sq(y - membrane_position(x) - yCentre)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    #endif

  foreach()
    cn[] = c[];
  boundary ({cn});
}

event logfile (i++; t <= TMAX)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  
  double dc = change (c, cn);
  if (i > 1 && dc < DC)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g %g\n",
	   MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER), dc);
  fprintf (stderr, "%g %g %g %g\n",
	   t, MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER), dc);
    
}

event gfsOutput(i += 1000) {
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

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D spurious.gfv", "w");
  output_gfs (fp);
}
#endif

/**
We use an adaptive mesh with a constant (maximum) resolution along the
interface. */

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

/**
## Results

The maximum velocity converges toward machine zero for a wide range of
Laplace numbers on a timescale comparable to the viscous dissipation
timescale, as expected.

~~~gnuplot Evolution of the amplitude of the capillary currents $\max(|\mathbf{u}|)(D/\sigma)^{1/2}$ as a function of non-dimensional time $\tau=t\mu/D^2$ for the range of Laplace numbers indicated in the legend.
set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", \
  'La-12000-5' w l t "La=12000"
~~~

The equilibrium shape and curvature converge toward the exact shape
and curvature at close to second-order rate.

~~~gnuplot Convergence of the error on the equilibrium shape of the droplet with resolution. The diameter is given in number of grid points.
set xlabel 'D'
set ylabel 'Shape error'
set logscale x
set xtics 2
set pointsize 1
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):5 w lp t "RMS", \
            '< sort -n -k1,2 log' u (0.8*2**$1):6 w lp t "Max", \
             0.2/(x*x) t "Second order"
~~~

~~~gnuplot Convergence of the relative error on the equilibrium curvature value with resolution. The diameter is given in number of grid points.
set ylabel 'Relative curvature error'
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):($7/2.5) w lp t "Max", \
             0.6/(x*x) t "Second order"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/

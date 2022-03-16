/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact
*/

/* Computational flags */
#define CURVATUREADJUST 1 // Adjusts the curvature
#define BODYFORCEADJUST 1 // Adjusts the body force
#define IFORCEADJUST 1 // Adjusts the interfacial force
#define MEMBRANE 1 // Impose a membrane to deform the droplet
#define AMR 1 // Adaptive mesh refinement
#define WALL 0 // Droplet along the wall
#define TRANSPOSED 0 // Transposes so the membrane is along y
#define SINGLESTEP 0 // If 1, only performs one timestep 
#define REFINEMENTSTUDY 0 // Runs at multiple levels for a refinement study
#define JACOBI 1 // Jacobi relaxation

/* Dimensional fluid properties */
const double R = 1.0e-3; // Radius of droplet (metres)
const double V = 5.0; // Velocity of droplet (metres per second)
const double RHO_L = 998.0; // Density of liquid (kilograms per metre cubed)
const double RHO_G = 1.23; // Density of gas (kilograms per metre cubed)
const double MU_L = 1.0e-3; // Viscosity of liquid (Pascal seconds)
const double MU_G = 1.81e-5; // Viscosity of gas (Pascal seconds)
const double SIGMA = 7.29e-2; // Surface tension (Newtons per metre)

/* Dimensionless droplet definitions */
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double INITIAL_DROP_HEIGHT = 0.125; // Initial gap between drop and plate

/* Computational constants */
const int MINLEVEL = 5; // Minimum refinement level 
const int MAXLEVEL = 8; // Maximum refinement level
const double BOX_WIDTH = 6.0; // Width of the computational box
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-4; // Time between plate outputs
const double LOG_OUTPUT_TIMESTEP = 1e-4; // Time between log outputs
const double MAX_TIME = 1.0; // Hard maximum time 
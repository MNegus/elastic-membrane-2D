/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact
*/

/* Fluid properties */
// Water at 5 m/s, radius 1mm
const double RHO_R = 0.00120; // Density ratio
const double MU_R = 0.0183; // Viscosity ratio
const double REYNOLDS = 4990; // Reynolds number
const double WEBER = 342; // Weber number
const double FR = 50.5; // Froude number

/* Droplet definition */
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double INITIAL_DROP_HEIGHT = 0.125; // Initial gap between drop and plate

/* Membrane parameters */
const double ALPHA = 2.0;
const double BETA = 100; // Tension term
const double GAMMA = 2.0; // Bending term
const double DELTA_T = 1e-4; // Timestep for membrane solution
const int COUPLED = 1; // Set to 1 if coupled

/* Computational constants */
const int MINLEVEL = 5; // Minimum refinement level 
const int MAXLEVEL = 12; // Maximum refinement level
const double BOX_WIDTH = 6.0; // Width of the computational box
const double MEMBRANE_RADIUS = 4.0; // Width of the membrane
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-4; // Time between plate outputs
const double LOG_OUTPUT_TIMESTEP = 1e-4; // Time between log outputs
const double MAX_TIME = 0.4; // Hard maximum time 
const double MEMBRANE_START_TIME = 1e-3; // Time which to start membrane motion
const int MOVIES = 1; // Set to 1 if making movies

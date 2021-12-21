/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact
*/

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

/* Membrane parameters */
const double ALPHA = 2; // Mass ratio
const double BETA = 1; // Tension term
const double GAMMA = 2; // Bending term
const double DELTA_T = 1e-4; // Timestep for membrane solution
const int COUPLED = 1; // Set to 1 if membrane motion is coupled to the fluid
const int STATIONARY = 0; // Set to 1 if enforcing membrane to be stationary

/* Computational constants */
const int MINLEVEL = 7; // Minimum refinement level 
const int MAXLEVEL = 13; // Maximum refinement level
const double BOX_WIDTH = 24.0; // Width of the computational box
const double MEMBRANE_RADIUS = 16.0; // Width of the membrane
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-4; // Time between plate outputs
const double LOG_OUTPUT_TIMESTEP = 1e-4; // Time between log outputs
const double MAX_TIME = 0.4; // Hard maximum time 
const double MEMBRANE_START_TIME = 1e-3; // Time which to start membrane motion
const int MOVIES = 1; // Set to 1 if making movies
const int CUTOFF = 1; // Set to 1 to cutoff the pressure
const int x_min_height = 2; // Min height to cutoff in x
const int y_min_height = 3; // Min height to cutoff in y
const int FD_COARSEN_LEVEL = 0;
const int REMOVE_ENTRAPMENT = 0;
const double REMOVAL_DELAY = 0.02; // Time after pinch-off to start removal

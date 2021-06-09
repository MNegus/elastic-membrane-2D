/* Function definitions */


void initialise_membrane(double *w_previous, double *w, double *p_previous, double *p, double *p_next, int N_MEMBRANE, double DELTA_T, double L, double ALPHA, double BETA);


void multiply_matrix(double *y_arr, double *matrix_arr, double *x_arr, double scale, int ADD);


void membrane_timestep(double *w_previous, double *w, double *w_next, double *p_previous, double *p, double *p_next);
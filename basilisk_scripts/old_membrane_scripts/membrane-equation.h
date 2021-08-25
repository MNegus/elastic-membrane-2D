/* membrane-equation.h */

void initialise_membrane(double *w_previous, double *w, double *w_deriv, 
    double *p_previous, double *p, double *p_next, int M, double DELTA_X, 
    double DELTA_T, double L, double ALPHA, double BETA, double GAMMA);

void multiply_matrix(int M, double *y_arr, double *matrix_arr, double *x_arr, double scale, int ADD);

void initialise_coefficient_matrices(int M, double DELTA_X, double DELTA_T, \
    double ALPHA, double BETA, double GAMMA);


void membrane_timestep(double *w_previous, double *w, double *w_next, 
    double *w_deriv, double *p_previous, double *p, double *p_next, 
    int M, double DELTA_X, double DELTA_T) ;
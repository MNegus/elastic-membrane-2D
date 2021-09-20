function [w_next, p, w_t, d, d_t, J] = membrane_timestep_first_order(...
    xs, t, w, w_previous, p_previous, pressure_type, EPSILON, DELTA_T, ...
    DELTA_X, M, Cpressure, A_mat, B_mat)
    
    %% Solves for w_next via matrix inversion
    rhs = B_mat * w - A_mat * w_previous + Cpressure * p_previous;
    w_next = A_mat \ rhs;

    %% Solves for remaining quantities
    % Membrane derivatives
    w_t = (w_next - w_previous) / (2 * DELTA_T);
    w_tt = (w_previous - 2 * w + w_next) / (DELTA_T^2);
    
    w_x = zeros(M, 1);
    w_x(2 : M - 1) = (w(3 : M) - w(1 : M - 2)) / (2 * DELTA_X);
    w_x(M) = -w(M - 1) / (2 * DELTA_X);
    w_x_fun = @(x) interp1(xs, w_x, x, 'linear', 'extrap');
    
    % Continous functions for w and its derivatives
    w_fun = @(x) interp1(xs, w, x, 'linear', 'extrap');
    w_t_fun = @(x) interp1(xs, w_t, x, 'linear', 'extrap');
    w_tt_fun = @(x) interp1(xs, w_tt, x, 'linear', 'extrap');
    
    % Pressure, turnover points and jet thickness
    [p, d, d_t, J] = w_dependent_first_order(xs, t, w_fun, ...
        w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);
    
end
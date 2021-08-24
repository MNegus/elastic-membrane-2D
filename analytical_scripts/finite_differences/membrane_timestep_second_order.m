function [w_next, p, w_t, d, d_t, J] = membrane_timestep_second_order(xs, t, ...
    w, w_previous, w_t_previous, d_previous, d_t_previous, pressure_type, ...
    EPSILON, M, DELTA_X, DELTA_T, ...
    Cpressure, A_mat, B_mat)
    

    %% Continuous functions for w and w_x
    w_fun = @(x) interp1(xs, w, x, 'linear', 'extrap');

    w_x = zeros(M, 1);
    w_x(2 : M - 1) = (w(3 : M) - w(1 : M - 2)) / (2 * DELTA_X);
    w_x(M) = -w(M - 1) / (2 * DELTA_X);
    w_x_fun = @(x) interp1(xs, w_x, x, 'linear', 'extrap');

    %% Initial guesses for w_next
    w_next_guess = w_previous + DELTA_T * w_t_previous;
    
    
    %% Solves for w_next using fsolve
    diff = 1e4;
%     tol = 1e-6;
    tol = max(1e-14, max(abs(DELTA_T * w_t_previous)) / 100);
    maxiter = 100;
    iternum = 1;
    w_next_curr = w_next_guess;
    while ((diff > tol) && (iternum <= maxiter))
        %% Determine w_t_vals and w_tt_vals
        w_t_vals = (w_next_curr - w_previous) / (2 * DELTA_T);
        w_tt_vals = (w_previous - 2 * w + w_next_curr) / (DELTA_T^2);
        
        %% Solves for pressure
        [p_vals, ~, ~, ~] = w_dependent_quantities(xs, t, ...
            w_t_vals, w_tt_vals, w_fun, w_x_fun, d_previous, d_t_previous, pressure_type, EPSILON, DELTA_T);
        
        %% Updates guess for w_next
        rhs = B_mat * w - A_mat * w_previous + Cpressure * p_vals;
        w_next_update = A_mat \ rhs;
        
        %% Determine difference
        diff = max(abs(w_next_update - w_next_curr));
        
        %% Update w_next_curr and increments iternum
        w_next_curr = w_next_update;
        iternum = iternum + 1;
        
    end
    w_next = w_next_update;
    
    if (iternum > maxiter) 
        warning("Max iter reached in membrane timestep"); 
    end
    
    %% 
        
        
    
%     zero_fun = @(w_next_vals) full_zero_fun(w_next_vals, xs, t, w, w_previous, ...
%             w_fun, w_x_fun, A_mat, B_mat, ...
%             DELTA_T, Cpressure, pressure_type, EPSILON);
%     
%     options = optimoptions('fsolve', 'Display', 'iter');
%     w_next = fsolve(zero_fun, w_next_guess, options);
    
    %% Solves for remaining quantities
    % Membrane derivatives
    w_t = (w_next - w_previous) / (2 * DELTA_T);
    w_tt = (w_previous - 2 * w + w_next) / (DELTA_T^2);
    
    % Pressure, turnover points and jet thickness
    [p, d, d_t, J] = w_dependent_quantities(xs, t, w_t, w_tt, ...
        w_fun, w_x_fun, d_previous, d_t_previous, pressure_type, ...
        EPSILON, DELTA_T);
    
%     %% Zero function used in fsolve to find w_next
%     function res = full_zero_fun(w_next_vals, xs, t, w, w_previous, ...
%             w_fun, w_x_fun, A_mat, B_mat, ...
%             DELTA_T, Cpressure, pressure_type, EPSILON)
%         
%         %% Determine w_t_vals and w_tt_vals
%         w_t_vals = (w_next_vals - w_previous) / (2 * DELTA_T);
%         w_tt_vals = (w_previous - 2 * w + w_next_vals) / (DELTA_T^2);
%         
%         %% Solves for pressure
%         [p_vals, ~, ~, ~] = w_dependent_quantities(xs, t, ...
%             w_t_vals, w_tt_vals, w_fun, w_x_fun, pressure_type, EPSILON);
%         
%         %% Determines res, the solution to the PDE
%         res = A_mat * w_next_vals -B_mat * w + A_mat * w_previous ...
%             - Cpressure * p_vals;
%         
%     end
end
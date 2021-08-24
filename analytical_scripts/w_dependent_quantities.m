function [p, d, d_t, J] = w_dependent_quantities(xs, t, w_t, w_tt, ...
    w_fun, w_x_fun, d_previous, d_t_previous, pressure_type, EPSILON, DELTA_T)
    
    %% Determine continous functions for w_t_vals and w_tt_vals
    w_t_fun = @(x) interp1(xs, w_t, x, 'linear', 'extrap');
    w_tt_fun = @(x) interp1(xs, w_tt, x, 'linear', 'extrap');

    %% Determine d and d_t
    [d, d_t] = turnover_point_trapz(xs, t, d_previous, d_t_previous, w_fun, w_t_fun, w_x_fun, EPSILON, DELTA_T);

    % Finds d_idx such that x(d_idx) < epsilon * d but x(d_idx) >= epsilon * d
    d_idx = sum(xs < EPSILON * d);

    %% Determine m and its derivatives
    if (d_idx < 2)
        m_t_fun = @(s) zeros(size(s));
        m_tt_fun = @(s) zeros(size(s));
    else
        s_vals = xs(1 : d_idx) / EPSILON;

        m_t = cumtrapz(s_vals, w_t(1 : d_idx), 1);
        m_tt = cumtrapz(s_vals, w_tt(1 : d_idx), 1);

        m_t_fun = @(s) interp1(s_vals, m_t, s, 'linear', 'extrap');
        m_tt_fun = @(s) interp1(s_vals, m_tt, s, 'linear', 'extrap');
    end
    %% Determine time-dependent quantities
    [A, C, J] = time_dependent_quantities(d, d_t, w_t_fun, w_tt_fun, m_t_fun, EPSILON);
    
    %% Determine pressure at current timestep
%     tic;
    if pressure_type == "outer"
        p = outer_pressure(xs, m_tt_fun, w_tt_fun, d, A, EPSILON);
    elseif pressure_type == "composite"
        p = composite_pressure(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, EPSILON);
    else
        error("Invalid pressure_type");
    end
%     time = toc;
%     disp("Pressure determination time = " + time);

    % Restrict large values of p
    p(p > 1e4) = 0;

end
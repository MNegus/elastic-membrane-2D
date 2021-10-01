function [w_next, p, w_t, d, d_t, J] = membrane_timestep(xs, t, ...
    w, w_previous, p_previous, p_previous_previous, pressure_type, ...
    ALPHA, BETA, GAMMA, EPSILON, M, DELTA_X, DELTA_T, ...
    Cpressure, A_mat, B_mat)
    
        %% Determines derivatives of w at current timestep
        w_t = (w - w_previous) / DELTA_T;

        % MAKE QUICKER BY TURNING INTO A MATRIX MULTIPLICATION
        w_tt = spatial_derivatives(w, ALPHA, BETA, GAMMA, M, DELTA_X) ...
            + (1 / ALPHA) * p_previous;

        w_x = zeros(M, 1);
        w_x(2 : M - 1) = (w(3 : M) - w(1 : M - 2)) / (2 * DELTA_X);
        w_x(M) = -w(M - 1) / (2 * DELTA_X);

        %% Determines anonymous functions for w and its derivatives
        w_fun = @(x) interp1(xs, w, x, 'linear', 'extrap');
        w_t_fun = @(x) interp1(xs, w_t, x, 'linear', 'extrap');
        w_tt_fun = @(x) interp1(xs, w_tt, x, 'linear', 'extrap');
        w_x_fun = @(x) interp1(xs, w_x, x, 'linear', 'extrap');

        %% Determine d and d_t
        [d, d_t] = turnover_point(t, w_fun, w_t_fun, w_x_fun, EPSILON);

        % Finds d_idx such that x(d_idx) < epsilon * d but x(d_idx) >= epsilon * d
        d_idx = sum(xs < EPSILON * d);

        %% Determine m and its derivatives
        if (d_idx < 2)
            m_t_fun = @(s) zeros(size(s));
            m_tt_fun = @(s) zeros(size(s));
        else
            s_vals = xs(1 : d_idx) / EPSILON;

            m_t = zeros(size(s_vals));
            m_tt = zeros(size(s_vals));

            % VECTORISE THIS FOR SPEEEEED
%             tic;
%             for q = 1 : length(s_vals)
%                 m_t(q) = trapz(s_vals(1 : q), w_t(1 : q), 1);
%                 m_tt(q) = trapz(s_vals(1 : q), w_tt(1 : q), 1);
%             end
%             toc;
            
            m_t = cumtrapz(s_vals, w_t(1 : d_idx), 1);
            m_tt = cumtrapz(s_vals, w_tt(1 : d_idx), 1);

            m_t_fun = @(s) interp1(s_vals, m_t, s, 'linear', 'extrap');
            m_tt_fun = @(s) interp1(s_vals, m_tt, s, 'linear', 'extrap');
        end


        %% Determine time-dependent quantities
        [A, B, C, J] = time_dependent_quantities(d, d_t, w_t_fun, w_tt_fun, m_t_fun, EPSILON);
        
        %% Determine pressure at current timestep
        if pressure_type == "outer"
            p = outer_pressure(xs, m_tt_fun, w_tt_fun, d, A, EPSILON);
        elseif pressure_type == "composite"
            p = composite_pressure(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, EPSILON);
        else
            error("Invalid pressure_type");
        end

        % Restrict large values of p
        p(p > 1e4) = 0;
        
        %% Determine right-hand-side and solve for w_next
        rhs = B_mat * w - A_mat * w_previous ...
            + Cpressure * (p_previous_previous + 2 * p_previous + p);
        w_next = A_mat \ rhs;
end
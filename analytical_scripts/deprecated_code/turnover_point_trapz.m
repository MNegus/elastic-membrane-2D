function [d, d_t] = turnover_point_trapz(xs, t, d_previous, d_t_previous, ...
    w_fun, w_t_fun, w_x_fun, epsilon, delta_t) 

    xhats = xs / epsilon;

    if (t == 0)
        d = 0; d_t = 0;
    else
        
        %% Use iterative procedure to determine d(t)
        diff = 1e4;
        d_guess = d_previous + delta_t * d_t_previous;
        tol = max(1e-14, max(abs(delta_t * d_t_previous)) / 100);
        maxiter = 50;
        iternum = 1;
        d_curr = d_guess;
        while ((diff > tol) && (iternum <= maxiter))
            
            % Updates d
            idx = sum(xhats < d_curr);
            s_vals = xhats(1 : idx);
            integrand =  w_fun(epsilon * s_vals) ./ sqrt(d_curr^2 - s_vals.^2);
            d_update = 2 * sqrt(t - (2 / pi) * trapz(s_vals, integrand));
            
            % Determines difference
            diff = abs(d_update - d_curr);
            
            %% Update d_curr and increments iternum
            d_curr = d_update;
            iternum = iternum + 1; 
        end
        d = d_update;
    
        if (iternum > maxiter) 
            warning("Max iter in turnover_point reached, using fsolve");
            
            % Use fsolve
            options = optimoptions('fsolve', 'OptimalityTolerance', 1e-8, 'display', 'off');
            d_zero_fun = @(d) full_d_zero_fun(d, t, w_fun, epsilon);
            d = fsolve(d_zero_fun, d_guess, options);
        end

        %% Determine d'(t)
        idx = sum(xhats < d);
        s_vals = xhats(1 : idx);
        integrand_1 = epsilon * s_vals .* w_x_fun(epsilon * s_vals) ./ sqrt(d^2 - s_vals.^2);
        integrand_2 = w_t_fun(s_vals) ./ sqrt(d^2 - s_vals.^2);
        
        d_t = (1 - (2 / pi) * trapz(s_vals, integrand_2)) ...
            / (d / 2 + (2 / (pi * d)) * trapz(s_vals, integrand_1));
        
    end
    
    %% Function definition for using fsolve if needed
    function res = full_d_zero_fun(d, t, w_fun, epsilon)
        integrand = @(s) w_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        res = t - d^2 / 4 - (2 / pi) * integral(integrand, 0, d);
    end
end
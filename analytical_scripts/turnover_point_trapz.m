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
        maxiter = 10;
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
            warning("Max iter reached"); 
        end

        %% Determine d'(t)
        idx = sum(xhats < d);
        s_vals = xhats(1 : idx);
        integrand_1 = epsilon * s_vals .* w_x_fun(epsilon * s_vals) ./ sqrt(d^2 - s_vals.^2);
        integrand_2 = w_t_fun(s_vals) ./ sqrt(d^2 - s_vals.^2);
        
        d_t = (d / 2 + (2 / (pi * d)) * trapz(s_vals, integrand_1)) ...
            / (1 - (2 / pi) * trapz(s_vals, integrand_2));
        
    end
    
end
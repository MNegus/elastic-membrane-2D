 function ps = outer_pressure(xs, m_tt_fun, w_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    xhats = xs / epsilon;

    %% Finds idx such that xhat = d
    idx = sum(xhats < d);
    
    %% Determines the singular integral for 0 <= x < epsilon * d
    xhat_dependents = zeros(size(xs(1 : idx)));
    s_vals = xhats(1 : idx);
    
    for m = 1 : idx
        xhat = xhats(m);
        
        %% trapz method (USE CUMTRAPZ)
        integrand = @(s) trapz_integrand(s, xhat, m, d, m_tt_fun, w_tt_fun, epsilon);
        
        integral_value = trapz(s_vals, integrand(s_vals), 1);
        xhat_dependents(m) = xhat * m_tt_fun(xhat) - integral_value;

    end
    
    %% Returns p values
    ps(1 : idx) = (A + xhat_dependents) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    
    %% Function definition
    function vals = trapz_integrand(s, xhat, xhat_idx, d, m_tt_fun, w_tt_fun, epsilon)
        
        % Sets the bulk nodes
        vals = (2 / pi) * (sqrt(d^2 - s.^2) .* (s .* m_tt_fun(s) - xhat * m_tt_fun(xhat))) ./ (s.^2 - xhat^2);
        
        % Adjusts the node where s = xhat to avoid singularity
        if (xhat == 0)
            vals(xhat_idx) = (2 / pi) * d * w_tt_fun(0);
        else 
            vals(xhat_idx) = (2 / pi) * sqrt(d^2 - xhat^2) * (m_tt_fun(xhat) + xhat * w_tt_fun(epsilon * xhat)) / (2 * xhat);
        end
        
        
    end

end
function ps = outer_pressure_alt(xs, m_tt_fun, w_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    xhats = xs / epsilon;

    %% Finds idx such that x
    idx = sum(xhats < d);
    
    %% Determines the singular integral for 0 <= x < epsilon * d
    singular_integrals = zeros(size(xs(1 : idx)));
    
    for m = 1 : idx
        xhat = xhats(m)

        %% integral method from -d to d
%         integrand = @(s) (1 / pi) * (sqrt(d^2 - s.^2) .* (m_tt_fun(s) - m_tt_fun(xhat))) ./ (s - xhat);
%         singular_integrals(m) = xhat * m_tt_fun(xhat) - integral(integrand, -d, d);
    
        %% integral method from 0 to d
%         integrand = @(s) (2 / pi) * (sqrt(d^2 - s.^2) .* (s .* m_tt_fun(s) - xhat * m_tt_fun(xhat))) ./ (s.^2 - xhat^2);
%         singular_integrals(m) = xhat * m_tt_fun(xhat) - integral(integrand, 0, d);
        
        %% trapz method
        integrand = @(s) trapz_integrand(s, xhat, m, d, m_tt_fun, w_tt_fun, epsilon);
        s_vals = xhats(1 : idx)';
        integrand(s_vals)
        integral_value = trapz(s_vals, integrand(s_vals), 1)
%         error("integral value")
        
        singular_integrals(m) = xhat * m_tt_fun(xhat) - integral_value;

    end
    
    %% Returns p values
    ps(1 : idx) = (A + singular_integrals) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    
    
    function vals = trapz_integrand(s, xhat, xhat_idx, d, m_tt_fun, w_tt_fun, epsilon)
        vals = (2 / pi) * (sqrt(d^2 - s.^2) .* (s .* m_tt_fun(s) - xhat * m_tt_fun(xhat))) ./ (s.^2 - xhat^2);
        if (xhat == 0)
            vals(xhat_idx) = (2 / pi) * d * w_tt_fun(0);
        else 
            vals(xhat_idx) = (2 / pi) * sqrt(d^2 - xhat^2) * (m_tt_fun(xhat) + xhat * w_tt_fun(epsilon * xhat)) / (2 * xhat);
        end
        
        
    end

end
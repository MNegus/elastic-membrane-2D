function ps = outer_pressure_integral_method(xs, m_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    xhats = xs / epsilon;

    %% Finds idx such that xs(idx) < epsilon * d and xs(idx + 1) >= epsilon * d
    idx = sum(xhats < d);
    
    %% Determines the singular integral for 0 <= x < epsilon * d
    tol = 1e-10;
    singular_integrals = zeros(size(xs(1 : idx)));
    
    for m = 1 : idx
        xhat = xhats(m);
        
        %% integral/quadgk method
        integrand = @(s) sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s - xhat);
        singular_integrals(m) = (1 / pi) * (integral(integrand, -d, xhat - tol) + integral(integrand, xhat + tol, d));
        
        %% trapz method with s = xhats (kind of works)
%         integrand = @(s) 2 * s .* sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s.^2 - xhat^2);
%         
%         if (m <= 1)
%             integral_1 = 0;
%         else
%             s_vals_1 = xhats(1 : m - 1)';
%             integral_1 = trapz(s_vals_1, integrand(s_vals_1), 1);
%         end
%         
%         if (m >= idx - 1)
%             integral_2 = 0;
%         else
%             s_vals_2 = xhats(m + 1 : idx)';
%             integral_2 = trapz(s_vals_2, integrand(s_vals_2), 1);
%         end
%         
%         singular_integrals(m) = (1 / pi) * (integral_1 + integral_2);

    end
    
    %% Returns p values
    ps(1 : idx) = (1 / epsilon) * (A - singular_integrals) ./ sqrt(d^2 - xhats(1 : idx).^2);
    

end
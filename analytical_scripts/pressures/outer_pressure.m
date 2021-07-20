function ps = outer_pressure(xs, m_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    xhats = epsilon * xs;

    %% Finds idx such that xs(idx) < epsilon * d and xs(idx + 1) >= epsilon * d
    idx = sum(xhats < d);
    
    %% Determines the singular integral for 0 <= x < epsilon * d
    tol = 1e-10;
    n = 1e3;
    singular_integrals = zeros(size(xs(1 : idx)));
    
    for m = 1 : idx
        xhat = epsilon * xs(m);
        
        %% integral method
%         integrand = @(s) sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s - xhat);
%         singular_integrals(m) = (1 / pi) * (integral(integrand, -d, xhat - tol) + integral(integrand, xhat + tol, d));
        
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

        %% trapz method with svals different
%         integrand = @(s) 2 * s .* sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s.^2 - xhat^2);
%         if (m <= 1)
%             integral_1 = 0;
%         else
%             s_vals_1 = linspace(0, xhat, n + 1)';
%             s_vals_1 = s_vals_1(1 : end - 1);
%             integral_1 = trapz(s_vals_1, integrand(s_vals_1), 1);
%         end
%         
%         if (m >= idx - 1)
%             integral_2 = 0;
%         else
%             s_vals_2 = linspace(xhat, d, n + 1)';
%             s_vals_2 = s_vals_2(2 : end);
%             integral_2 = trapz(s_vals_2, integrand(s_vals_2), 1);
%         end
%         
%         singular_integrals(m) = (1 / pi) * (integral_1 + integral_2);

        %% Chebfun method
%         integrand = @(s) 2 * s .* sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s.^2 - xhat^2);
%         if ((m == 1) || (m == idx))
%             cheb_integral = 0;
%         else  
%             cheb_integrand = chebfun(integrand, [0, xhat, d], 'exps', [0, -1, 0], 'splitting', 'on');
%             cheb_integral = sum(cheb_integrand);
%         end
%       
%         singular_integrals(m) = (1 / pi) * (cheb_integral);
    end
    
    %% Returns p values
    singular_integrals
    ps(1 : idx) = (A - singular_integrals) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    

end
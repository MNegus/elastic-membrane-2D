function ps = outer_pressure(xs, w_t_fun, w_tt_fun, m_tt_fun, d, d_t, epsilon)
    ps = zeros(size(xs));

    %% Finds idx such that xs(idx) < epsilon * d and xs(idx + 1) >= epsilon * d
    idx = sum(xs < epsilon * d);
    
    %% Determines value of A(t)
    integrand_1 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
    integral_1 = (1 / pi) * integral(integrand_1, -d, d);
    
    integrand_2 = @(s) w_tt_fun(epsilon * s) .* sqrt(d^2 - s.^2);
    integral_2 = (1 / pi) * integral(integrand_2, -d, d);
    
    A = d * d_t * (1 - integral_1) - integral_2
    
    %% Determines the singular integral for 0 <= x < epsilon * d
    tol = 1e-10;
    singular_integrals = zeros(size(xs(1 : idx)));
    
    for m = 1 : idx
        xhat = epsilon * xs(m);
        integrand = @(s) sqrt(d^2 - s.^2) .* m_tt_fun(s) ./ (s - xhat);
        singular_integrals(m) = (1 / pi) * (integral(integrand, -d, xhat - tol) + integral(integrand, xhat + tol, d));
    end
    
    %% Returns p values
    ps(1 : idx) = (1 / epsilon) * (A - singular_integrals) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    

end
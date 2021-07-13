function ps = overlap_pressure(xs, w_t_fun, d, d_t, epsilon)
    ps = zeros(size(xs));
    
    %% Finds idx such that xs(idx) < epsilon * d and xs(idx + 1) >= epsilon * d
    idx = sum(xs < epsilon * d);
    
    %% Determines value of C(t)
    integrand_1 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
    integral_1 = (1 / pi) * integral(integrand_1, -d, d);
   
    C = d * d_t * (1 - integral_1)
    
    %% Returns p values
    ps(1 : idx) = (C / sqrt(2 * epsilon * d)) ./ sqrt(epsilon * d - xs(1 : idx));

end
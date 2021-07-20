function [d, d_t, A, B, C, J] = time_dependent_quantities(t, w_fun, w_t_fun, w_tt_fun, w_x_fun, m_t_fun, epsilon)

    %% Function definitions
    function res = full_d_zero_fun(d, t, w_fun, epsilon)
        integrand = @(s) w_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        res = t - d^2 / 4 - (1 / pi) * integral(integrand, -d, d);
    end

    function res = full_d_t_zero_fun(d_t, d, w_t_fun, w_x_fun, epsilon)
        integrand_1 = @(s) epsilon * s .* w_x_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        integrand_2 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        res = 1 - d_t * d / 2 - (d_t / (pi * d)) * integral(integrand_1, -d, d) ...
               - (1 / pi) * integral(integrand_2, -d, d);
    end

    %% Optimising settings
    options = optimoptions('fsolve', 'OptimalityTolerance', 1e-8);

    %% Determine d(t)
    d_zero_fun = @(d) full_d_zero_fun(d, t, w_fun, epsilon);
    d = fsolve(d_zero_fun, 2 * sqrt(t), options);
    
    %% Determine d'(t)
    d_t_zero_fun = @(d_t) full_d_t_zero_fun(d_t, d, w_t_fun, w_x_fun, epsilon);
    d_t = fsolve(d_t_zero_fun, 1 / sqrt(t), options);

    %% Determine A(t)
    A_integrand_1 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
    A_integral_1 = (1 / pi) * integral(A_integrand_1, -d, d);
    
    A_integrand_2 = @(s) w_tt_fun(epsilon * s) .* sqrt(d^2 - s.^2);
    A_integral_2 = (1 / pi) * integral(A_integrand_2, -d, d);
    
    A = d * d_t * (1 - A_integral_1) - A_integral_2;
    
    %% Determine B(t)
    B_integrand = @(s) (m_t_fun(d) - m_t_fun(s)) ./ (sqrt(d^2 - s.^2) .* (d - s));
    B = (1 / pi) * integral(B_integrand, -d, d);
    
    %% Determine C(t)
    C = A + A_integral_2;
    
    %% Determine J(t)
    J = pi * (1 - B)^2 * d / (8 * d_t^2);
    
end
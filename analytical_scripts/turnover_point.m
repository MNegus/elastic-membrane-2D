function [d, d_t] = turnover_point(t, w_fun, w_t_fun, w_x_fun, epsilon)

    %% Function definitions
    function res = full_d_zero_fun(d, t, w_fun, epsilon)
        integrand = @(s) w_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        res = t - d^2 / 4 - (2 / pi) * integral(integrand, 0, d);
    end

    function res = full_d_t_zero_fun(d_t, d, w_t_fun, w_x_fun, epsilon)
        integrand_1 = @(s) epsilon * s .* w_x_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        integrand_2 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        res = 1 - d_t * d / 2 - (2 * d_t / (pi * d)) * integral(integrand_1, 0, d) ...
               - (2 / pi) * integral(integrand_2, 0, d);
    end

    if (t == 0)
        d = 0; d_t = 0;
    else
        %% Optimising settings
        options = optimoptions('fsolve', 'OptimalityTolerance', 1e-8, ...
            'display', 'off');

        %% Determine d(t)
        d_zero_fun = @(d) full_d_zero_fun(d, t, w_fun, epsilon);
        d = fsolve(d_zero_fun, 2 * sqrt(t), options);

        %% Determine d'(t)
        d_t_zero_fun = @(d_t) full_d_t_zero_fun(d_t, d, w_t_fun, w_x_fun, epsilon);
        d_t = fsolve(d_t_zero_fun, 1 / sqrt(t), options);
    end
end
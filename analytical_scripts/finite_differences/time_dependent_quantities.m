function [A, C, J] = time_dependent_quantities(d, d_t, w_t_fun, w_tt_fun, m_t_fun, epsilon)
    
    if (d == 0)
        A = 0; C = 0; J = 0;
    else
        %% Determine A(t)
        A_integrand_1 = @(s) w_t_fun(epsilon * s) ./ sqrt(d^2 - s.^2);
        A_integral_1 = (2 / pi) * integral(A_integrand_1, 0, d);

        A_integrand_2 = @(s) w_tt_fun(epsilon * s) .* sqrt(d^2 - s.^2);
        A_integral_2 = (2 / pi) * integral(A_integrand_2, 0, d);

        A = d * d_t * (1 - A_integral_1) - A_integral_2;

        %% Determine C(t)
        C = d * d_t * (1 - A_integral_1) ;

        %% Determine J(t)
        B_integrand = @(s) (d * m_t_fun(d) - s .* m_t_fun(s)) ./ (d^2 - s.^2).^(3/2);
        B = (2 / pi) * integral(B_integrand, 0, d);
        J = pi * (1 - B)^2 * d / (8 * d_t^2);
    end
    
end
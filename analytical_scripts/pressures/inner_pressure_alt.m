function ps = inner_pressure_alt(xs, d, d_t, J, epsilon)
    
%     %% Uses xs to find etas, the pressure parameter
    zero_fun = @(etas) (xs - epsilon * d) / epsilon^3 ...
        + (J / pi) * (exp(2 * etas) + 4 * exp(etas) + 2 * etas + 1);
    etas_0 = linspace(10, -100, length(xs));
    etas = fsolve(zero_fun, etas_0);


    %% Return pressure as a function of etas
    ps = (1 / epsilon^2) * 2 * d_t^2 * exp(etas) ./ (1 + exp(etas)).^2;
end
function ps = inner_pressure(xs, d, d_t, J, epsilon)

    %% Use etas_definition and map etas back onto xs
    % Find first value in xs such that x > epsilon * d
    x_d_idx = find(xs > epsilon * d, 1, 'first');
    
    x_tildes = (xs - epsilon * d) / epsilon^3;
    
    % Negative etas part
    negative_etas = log(sqrt(-J ./ (pi * x_tildes(1 : x_d_idx - 1))));

    % Positive etas part
    positive_etas = pi * x_tildes(x_d_idx : end) / (2 * J);
    
    % Put etas arrays together
    etas = zeros(length(negative_etas) + length(positive_etas), 1);
    etas(1 : length(negative_etas)) = negative_etas;
    etas(length(negative_etas) + 1 : end) = positive_etas;
    etas = unique(etas);
    
    % Save xs as a function of etas
    xs_etas = epsilon * d + (epsilon^3 * J / pi) * (2 * etas + 1 - 4 * exp(-etas) - exp(-2 * etas));
    
    % Save ps as a function of etas
    ps_etas = (1 / epsilon^2) * 2 * d_t^2 * exp(-etas) ./ (1 + exp(-etas)).^2;
    
    % Interpolate ps_etas onto xs
    ps = interp1(xs_etas, ps_etas, xs);
    
end
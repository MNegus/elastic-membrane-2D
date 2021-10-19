function xs_adapt = adaptive_x_values(d, A, tol, EPSILON)
    
    %% Determines the theoretical pressure maximum and min
    p_max = A / sqrt(d^2 - (d - tol)^2);
    p_min = A / d;
    
    %% Determines x values adapted to the outer pressure
    p_vals = linspace(p_min, p_max, 1e3)';
    lower_xs = sqrt(d^2 - A^2 ./ p_vals.^2);
    
    %% Reflect to make upper points
    upper_xs = -lower_xs + 2 * d;
    
    %% Convert to non-outer coordinates
    xs_adapt = (1 / EPSILON) * real([lower_xs; d; upper_xs]);

end
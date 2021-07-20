function ps = overlap_pressure(xs, d, C, epsilon)
    ps = zeros(size(xs));
    
    %% Finds idx such that xs(idx) < epsilon * d and xs(idx + 1) >= epsilon * d
    idx = sum(xs < epsilon * d);
    
    %% Returns p values
    ps(1 : idx) = C ./ (sqrt(2 * epsilon * d) * sqrt(epsilon * d - xs(1 : idx)));

end
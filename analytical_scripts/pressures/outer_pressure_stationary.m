 function ps = outer_pressure_stationary(xs, d, A, epsilon)
    ps = zeros(size(xs));
    
    xhats = xs / epsilon;

    %% Finds idx such that xhat = d
    idx = sum(xhats < d);
    
    %% Returns p values
    ps(1 : idx) = A ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);

end
 function ps = outer_pressure_stationary(xs, t, epsilon)
    ps = zeros(size(xs));
    
    xhats = xs / epsilon;
    
    if (t > 0)
        %% Determine d and A
        d = 2 * sqrt(t);
        d_t = 1 / sqrt(t);
        A = d * d_t;

        %% Finds idx such that xhat = d
        idx = sum(xhats < d);

        %% Returns p values
        ps(1 : idx) = A ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
        
    end
    
    

end
function ps = composite_pressure_stationary(xs, t, epsilon)

    if (t == 0)
        ps = zeros(size(xs));
    else
        % Calculate time dependent terms
        d = 2 * sqrt(t);
        d_t = 1 / sqrt(t);
        A = d * d_t;
        C = A;
        J = pi * d / (8 * d_t^2);
        
        % Calculate outer pressure
        outer_ps = outer_pressure_stationary(xs, t, epsilon);

        % Calculate overlap pressure
        overlap_ps = overlap_pressure(xs, d, C, epsilon);

        % Calculate inner pressure
        inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

        % Calculate composite pressure
        ps = outer_ps + inner_ps - overlap_ps;
    end
    
end


function ps = composite_pressure_stationary(xs, t, d, d_t, A, C, J, epsilon)

    if (t == 0)
        ps = zeros(size(xs));
    else
        % Calculate outer pressure
        outer_ps = outer_pressure_stationary(xs, d, A, epsilon);

        % Calculate overlap pressure
        overlap_ps = overlap_pressure(xs, d, C, epsilon);

        % Calculate inner pressure
        inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

        % Calculate composite pressure
        ps = outer_ps + inner_ps - overlap_ps;
    end
    
end


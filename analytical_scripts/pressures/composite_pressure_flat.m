function ps = composite_pressure_flat(xs, t, d, d_t, J, w_t_fun, w_tt_fun, epsilon)

    if (t == 0)
        ps = zeros(size(xs));
    else
        % Calculate outer pressure
        outer_ps = outer_pressure_flat(xs, w_t_fun, w_tt_fun, d, d_t, epsilon);
        
        % Calculate overlap pressure
        C = (1 - w_t_fun(0)) * d * d_t;
        overlap_ps = overlap_pressure(xs, d, C, epsilon);

        % Calculate inner pressure
        inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

        % Calculate composite pressure
        ps = outer_ps + inner_ps - overlap_ps;
    end
    
end

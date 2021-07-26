function ps = composite_pressure(xs, t, d, d_t, A, C, J, m_tt, w_tt, epsilon)

    if (t == 0)
        ps = zeros(size(xs));
    else
        % Calculate outer pressure
        outer_ps = outer_pressure(xs, m_tt, w_tt, d, A, epsilon);

        % Calculate overlap pressure
        overlap_ps = overlap_pressure(xs, d, C, epsilon);

        % Calculate inner pressure
        inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

        % Calculate composite pressure
        ps = outer_ps + inner_ps - overlap_ps;
    end
    
end


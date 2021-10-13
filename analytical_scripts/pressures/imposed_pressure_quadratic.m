function ps = imposed_pressure_quadratic(xs, t, d, d_t, J, w_t, w_tt, b, c, epsilon, pressure_type)

    if (t == 0)
        ps = zeros(size(xs));
    else
        %% Calculate outer pressure
        idx = sum(epsilon * xs < d);
        xhats = xs(1 : idx) / epsilon;
        outer_ps = zeros(size(xs));
        outer_ps(1 : idx) = (1 / epsilon) ...
            * (-(w_tt - 0.5 * c * (d^2 + 2 * xhats.^2)) .* sqrt(d^2 - xhats.^2) ...
            + ((1 - w_t) * d * d_t + 1.5 * b * d_t * d^3) ./ sqrt(d^2 - xhats.^2));
        
        %% Either outputs outer pressure or composite
        if (pressure_type == "outer")
            ps = outer_ps;
        elseif (pressure_type == "composite")
            % Calculate overlap pressure
            C = d * d_t * (1 - w_t + 3 * b * d^2 / 2);
            overlap_ps = overlap_pressure(xs, d, C, epsilon);

            % Calculate inner pressure
            inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

            % Calculate composite pressure
            ps = outer_ps + inner_ps - overlap_ps;
        else
            error("Invalid pressure type");
        end
        
        
    end
    
end

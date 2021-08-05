function [ds, as, a_ts, a_tts, q_ts] = a_solution_t_form(ts, ...
    t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals, alpha, DELTA_T)

    % Interpolates solutions for ds, as and a_ts
    ds = interp1(t_vals_d_form, d_vals_d_form, ts);
    as = interp1(t_vals_d_form, as_d_form, ts);
    a_ts = interp1(t_vals_d_form, a_ts_d_form, ts);
    
    % Determines solution for a_tts by numerical differentiation
    a_tts = zeros(size(a_ts));
    a_tts(2 : end, :) = diff(a_ts, 1, 1) / DELTA_T;
    
    % Determines q_ts using governing equations
    q_ts = -(alpha * a_tts + kvals' .* as); 
end


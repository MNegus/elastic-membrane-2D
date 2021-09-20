function q_ts = stationary_normal_modes_solution(ts, N, L, epsilon)
    
    %% Saves d and d_t
    ds = 2 * sqrt(ts);
    d_ts = 1 ./ sqrt(ts);
    
    %% Saves q_ts
    q_ts = zeros(length(ts), N);
    for n = 1 : N
        lambda = pi * (2 * n - 1) / (2 * L);
        q_ts(:, n) = -pi * besselj(0, epsilon * lambda * ds) .* ds .* d_ts;
    end

end
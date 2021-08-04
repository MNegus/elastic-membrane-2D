function [ws, w_ts] = w_solution_normal_modes(xs, as, a_ts, L, N)
    ws = zeros(size(xs));
    w_ts = zeros(size(xs));
    
    lambda = @(n) pi * (2 * n - 1) / (2 * L);
    
    
    %% FIND AN ALTERNATIVE USING MATRIX MULTIPLICATION
    for n = 1 : N
        ws = ws + as(n) * cos(lambda(n) * xs) / sqrt(L);
        w_ts = w_ts + a_ts(n) * cos(lambda(n) * xs) / sqrt(L);
    end

end
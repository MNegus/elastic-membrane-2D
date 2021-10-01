function [ws, w_ts, ps] = w_solution_normal_modes(xs, as, a_ts, q_ts, d, L, N, EPSILON)
%     ws = zeros(size(xs));
%     w_ts = zeros(size(xs));
%     ps = zeros(size(xs));
    
%     lambda = @(n) pi * (2 * n - 1) / (2 * L);
    
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
    
    %% FIND AN ALTERNATIVE USING MATRIX MULTIPLICATION
%     for n = 1 : N
%         ws = ws + as(n) * cos(lambda(n) * xs) / sqrt(L);
%         w_ts = w_ts + a_ts(n) * cos(lambda(n) * xs) / sqrt(L);
%         ps = ps - q_ts(n) * cos(lambda(n) * xs) / sqrt(L);
%     end
    
    %% Matrix multiplication
    ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);
    w_ts = sum(a_ts .* cos(xs * lambdas), 2) / sqrt(L);
    ps = sum(-q_ts .* cos(xs * lambdas), 2) / sqrt(L);
    
    ps(xs >= EPSILON * d) = 0;
    

end
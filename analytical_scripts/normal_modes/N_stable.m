function N = N_stable(alpha, beta, gamma, L, q, delta_d)

%     lambda = @(n) pi * (2 * n - 1) / (2 * L);
% 
%     tN = @(N) 2 * pi * sqrt(alpha / (beta * lambda(N)^2 + gamma * lambda(N)^4));
%     
%     N_zero_fun = @(N) delta_d - tN(N) / q;
%     
%     N = fsolve(N_zero_fun, 2);


    if (gamma > 0)
        c = alpha * (2 * pi / (q * delta_d))^2;
        lambda2 = (-beta + sqrt(beta^2 + 4 * gamma * c)) / (2 * gamma);
        N = 0.5 + (L / pi) * sqrt(lambda2);
    else
        N = 0.5 + (2 * L / pi) * sqrt(alpha / beta);
    end
end
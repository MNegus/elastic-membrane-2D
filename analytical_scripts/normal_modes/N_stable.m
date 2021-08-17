function N = N_stable(alpha, beta, gamma, L, q, delta_d)

    lambda = @(n) pi * (2 * n - 1) / (2 * L);

    tN = @(N) 2 * pi * sqrt(alpha / (beta * lambda(N)^2 + gamma * lambda(N)^4));
    
    N_zero_fun = @(N) delta_d - tN(N) / q;
    
    N = fsolve(N_zero_fun, 2);
end
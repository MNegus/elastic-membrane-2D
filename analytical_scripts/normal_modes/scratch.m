d = 0.1;
alpha = 1;
epsilon = 1;
L = 1;
N = 4;
lambdas = pi * (2 * (1 : N) - 1) / (2 * L)
lambda = @(n) pi * (2 * n - 1) / (2 * L);

[M, S] = mass_matrix(d, alpha, epsilon, L, N)

Snn = @(n) (pi * d^2 / 2) * (besselj(0, epsilon * lambdas(n) * d)^2 + besselj(1, epsilon * lambdas(n) * d)^2);

Smn = @(m, n) (pi * d / epsilon) ...
    * (lambdas(n) * besselj(0, epsilon * lambdas(m) * d) * besselj(1, epsilon * lambdas(n) * d) ...
        - lambdas(m) * besselj(0, epsilon * lambdas(n) * d) * besselj(1, epsilon * lambdas(m) * d)) ...
        / (lambdas(n)^2  - lambdas(m)^2);

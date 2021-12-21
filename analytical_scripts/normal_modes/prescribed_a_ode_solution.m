function [t_vals, d_vals, a_vals, a_t_vals, ks] = prescribed_a_ode_solution(alpha, beta, gamma, epsilon, dd, dmax, N, L)

% Discretised dvals
d_vals = 0 : dd : dmax;

% Defines constant vectors
lambdas = pi * (2 * (1 : N)' - 1) / (2 * L);
ks = epsilon^2 * (beta * lambdas.^2 + gamma * lambdas.^4);
pns = (alpha / epsilon^2) ./ sqrt(L * lambdas);

% Solves ODE
odefun = @(d, y) full_odefun(d, y, lambdas, ks, pns, alpha, epsilon, L, N);
y0 = zeros(2 * N, 1);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4, 'Stats', 'on', 'Maxstep', dd);
[d_vals, y] = ode45(odefun, d_vals, y0, opts);
a_vals = y(:, 1 : N);
b_vals = y(:, N + 1 : 2 * N);

% Determines time values
t_vals = zeros(size(d_vals));
for k = 1 : length(t_vals)
    t_vals(k) = d_vals(k)^2 / 4 ...
        + dot(a_vals(k, :), besselj(0, epsilon * d_vals(k) * lambdas)) / sqrt(L);
end

% Determine a_t vals
a_t_vals = zeros(length(t_vals), N);
for k = 1 : length(t_vals)
    d = d_vals(k);
    Gamma_vals = besselj(0, epsilon * d * lambdas);
    q_vals = -(d^2 / 4 + dot(a_vals(k, :)', Gamma_vals) / sqrt(L)) * pns;
    a_t_vals(k, :) = (epsilon^2 / alpha) * (ks .* b_vals(k, :)' - q_vals);
end

%% ODE function definition
function dydd = full_odefun(d, y, lambdas, ks, pns, alpha, epsilon, L, N)
    d
    
    dydd = zeros(2 * N, 1);
    
    %% Load in values
    as = y(1 : N, 1);
    bs = y(N + 1 : 2 * N, 1);
    
    %% Calculate useful quantities
    
    % Useful quantities
    Gammas = besselj(0, epsilon * d * lambdas);
    Gammas_d = - epsilon * lambdas .* besselj(1, epsilon * d * lambdas);
    qs = -(d^2 / 4 + dot(as, Gammas) / sqrt(L)) * pns;
    Fs = (epsilon^2 / alpha) * (ks .* bs - qs);
    Q = (sqrt(L) * d + 2 * dot(as, Gammas_d)) / (2 * sqrt(L) - 2 * dot(Fs, Gammas));
    if Q < 0
        Q
        warning("Q is negative");
    end
    
    %% Saves derivatives
    dydd(1 : N) = Q * Fs;
    dydd(N + 1 : 2 * N) = - Q * as;
end


end
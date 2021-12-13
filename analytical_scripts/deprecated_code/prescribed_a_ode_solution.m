function [t_vals, d_vals, a_vals, a_t_vals, ks] = prescribed_a_ode_solution(alpha, beta, gamma, epsilon, dd, dmax, N, L)

% Discretised dvals
d_vals = 0 : dd : dmax;

% Solves ODE
lambdas = pi * (2 * (1 : N)' - 1) / (2 * L);
ks = beta * lambdas.^2 + gamma * lambdas.^4;
odefun = @(d, y) full_odefun(d, y, lambdas, ks, alpha, epsilon, L, N);
y0 = zeros(2 * N, 1);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4, 'Stats', 'on', 'Maxstep', dd);
[d_vals, y] = ode45(odefun, d_vals, y0, opts);
a_vals = y(:, 1 : N);
bvals = y(:, N + 1 : 2 * N);

% Determines time values
t_vals = zeros(size(d_vals));
for k = 1 : length(t_vals)
    t_vals(k) = d_vals(k)^2 / 4;
end

% Determine a_t vals
a_t_vals = zeros(length(t_vals), N);
for k = 1 : length(t_vals)
    d = d_vals(k);
    qvals = - alpha * d^2 ./ (4 * sqrt(L) * sqrt(lambdas)); 
    a_t_vals(k, :) = (1 / alpha) * (ks .* bvals(k, :)' - qvals);
end

%% ODE function definition
function dydd = full_odefun(d, y, lambdas, ks, alpha, epsilon, L, N)
    d
    
    dydd = zeros(2 * N, 1);
    
    %% Load in values
    as = y(1 : N, 1);
    bs = y(N + 1 : 2 * N, 1);
    
    %% Calculate useful quantities
    qs = - alpha * d^2 ./ (4 * sqrt(L) * sqrt(lambdas)); 
    Fs = (1 / alpha) * (ks .* bs - qs);
    
%     Gammas = besselj(0, epsilon * d * lambdas);
%     Gammas_d = - epsilon * lambdas .* besselj(1, epsilon * d * lambdas);
%     Q =(sqrt(L) * d + 2 * dot(as, Gammas_d)) / (2 * sqrt(L) - 2 * dot(Fs, Gammas));
    Q = d / 2;
    if Q < 0
        Q
        warning("Q is negative");
    end
    
    %% Saves derivatives
    dydd(1 : N) = Q * Fs;
    dydd(N + 1 : 2 * N) = - Q * as;
end


end
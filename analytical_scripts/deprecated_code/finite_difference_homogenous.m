%% finite_difference_homogenous.m
clear;
%% Parameters
L = 1;
epsilon = 1;
N = 512;
M = N - 1;
dx = L / (N - 1);
xs = (0 : dx : L - dx)';
% tmax = (L / (2 * epsilon))^2;
tmax = 0.5
t0 = 1e-6;
dt = tmax / 100;
tvals = t0 : dt : tmax;

alpha = 2 / epsilon^4;
beta = 0;
gamma = 2;

%% Initialise arrays
lambda = pi * (2 * 3 - 1) / (2 * L);
a = 0.1;

% Initial conditions for w
w0 = a * cos(lambda * xs);
w_t0 = zeros(size(xs));
w_tt0 = zeros(size(xs));
y0 = [w0; w_t0];
yp0 = [w_t0; w_tt0];

% Initialise A_mat matrix
Cbeta = beta * dx^2 / gamma;

A_upper_upper = -ones(M, 1);
A_upper_upper(3) = -2;

A_upper = (4 + Cbeta) * ones(M, 1);
A_upper(2) = 8 + 2 * Cbeta;

A_main = -(6 + 2 * Cbeta) * ones(M, 1);
A_main(2) = -(7 + 2 * Cbeta);
A_main(M) = -(5 + 2 * Cbeta);

A_lower = (4 + Cbeta) * ones(M, 1);

A_lower_lower = -ones(M, 1);

A_mat = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

%% Solves the ODE
ode_options = odeset('RelTol',1e-8, 'AbsTol', 1e-8, 'Maxstep', dt);
tol = 5e-3;
ode_fun = @(t, y, yp) full_ode_fun(t, y, yp, xs, A_mat, M, dx,  alpha,  gamma);
sol = ode15i(ode_fun, [t0 tmax], y0, yp0, ode_options);
[y, yp] = deval(sol, tvals);
    

%% Plots solution
ws = y(1 : M, :);

figure(1);
for k = 1 : length(tvals)
    
    plot(xs, ws(:, k));
    drawnow;
    pause(0.01);
end

%% ODE function definition
function res = full_ode_fun(t, y, yp, xs, A_mat, M, dx,  alpha,  gamma)
    
    res = zeros(2 * M, 1);
    
    %% Load in functions
    ws = y(1 : M, 1);
    qs = y(M + 1 : 2 * M, 1);
    w_ts = yp(1 : M, 1);
    q_ts = yp(M + 1 : 2 * M, 1);
    
    %% Governing equations
    res(1 : M) = qs - w_ts;
    res(M + 1 : 2 * M) = alpha * q_ts - (gamma / dx^4) * A_mat * ws;
    
    t
end
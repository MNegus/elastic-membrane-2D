%% finite_difference_solution.m

addpath("pressures");


%% Parameters
L = 1;
epsilon = 1;
N = 512;
M = N - 1;
dx = L / (N - 1);
xs = (0 : dx : L - dx)';
% tmax = (L / (2 * epsilon))^2;
tmax = 0.05
t0 = 1e-6;
dt = tmax / 100;
tvals = t0 : dt : tmax;

alpha = 2 / epsilon^4;
beta = 1;
gamma = 2;

%% Initialise arrays

% Initial conditions for w
w0 = zeros(size(xs));
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
% ode_options = odeset('RelTol',1e-8, 'AbsTol', 1e-8, 'Maxstep', dt);
tol = 5e-3;
ode_fun = @(t, y, yp) full_ode_fun(t, y, yp, xs, A_mat, M, L, dx, epsilon, alpha, beta, gamma);
sol = ode15i(ode_fun, tvals, y0, yp0);
[y, yp] = deval(sol, tvals);
    
%% ODE function definition
function res = full_ode_fun(t, y, yp, xs, A_mat, M, L, dx, epsilon, alpha, beta, gamma)
    
    res = zeros(2 * M, 1);
    
    %% Load in functions
    ws = y(1 : M, 1);
    qs = y(M + 1 : 2 * M, 1);
    w_ts = yp(1 : M, 1);
    q_ts = yp(M + 1 : 2 * M, 1);
    w_tts = q_ts;
    
    %% Determine d and d_t
    d = 2 * sqrt(t);
    d_t = 1 / sqrt(t);

    %% Determine time-dependent quantities
    A = d * d_t;
    C = d * d_t;
    J = pi * d / (8 * d_t^2);
    
    %% Determine composite pressures
    outer_ps = outer_pressure_stationary(xs, d, A, epsilon);
    overlap_ps = overlap_pressure(xs, d, C, epsilon);
    inner_ps = inner_pressure(xs, d, d_t, J, epsilon);
    ps = outer_ps + inner_ps - overlap_ps;
    if (t > 0)
        figure(1);
        plot(xs, ps, '-o');
        xlim([0, 1.5 * epsilon * d]);
        pause(0.1);
        drawnow;
    end
    
    %% Governing equations
    res(1 : M) = qs - w_ts;
    res(M + 1 : 2 * M) = alpha * q_ts + (gamma / dx^4) * A_mat * ws - ps;
    t
end

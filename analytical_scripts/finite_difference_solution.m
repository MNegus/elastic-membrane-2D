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
tvals = t-1 : dt : tmax;

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
    
    % Determine w_x
    w_xs = zeros(M, 1);
    w_xs(2 : M - 1) = (ws(3 : M) - ws(1 : M - 2)) / (2 * dx);
    w_xs(M) = -ws(M - 1) / (2 * dx);
    
    %% Create interpolated functions for w and its deriatives
    w_fun = @(x) interp1(xs, ws, x);
    w_t_fun = @(x) interp1(xs, w_ts, x);
    w_tt_fun = @(x) interp1(xs, w_tts, x);
    w_x_fun = @(x) interp1(xs, w_xs, x);
    
    %% Determine d and d_t
    [d, d_t] = turnover_point(t, w_fun, w_t_fun, w_x_fun, epsilon);
    
    % Finds d_idx such that x(d_idx) < epsilon * d but x(d_idx) >= epsilon * d
    d_idx = sum(xs < epsilon * d);
    
    %% Determine m and its derivatives
    if (d_idx < 2)
        m_t_fun = @(s) zeros(size(s));
        m_tt_fun = @(s) zeros(size(s));
    else
        s_vals = xs(1 : d_idx) / epsilon;
        
        m_ts = zeros(size(s_vals));
        m_tts = zeros(size(s_vals));
        
        % VECTORISE THIS FOR SPEEEEED
        for k = 1 : length(s_vals)
            m_ts(k) = trapz(s_vals(1 : k), w_ts(1 : k), 1);
            m_tts(k) = trapz(s_vals(1 : k), w_tts(1 : k), 1);
        end
        
        m_t_fun = @(s) interp1(s_vals, m_ts, s, 'linear', 'extrap');
        m_tt_fun = @(s) interp1(s_vals, m_tts, s, 'linear', 'extrap');
        
    end
    
    %% Determine time-dependent quantities
    [A, C, J] = time_dependent_quantities(d, d_t, w_t_fun, w_tt_fun, m_t_fun, epsilon);
    
    %% Determine composite pressures
    ps = composite_pressure(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, epsilon);
    if (t > 0)
        figure(1);
        plot(xs, ps);
        xlim([0, 1.5 * d]);
        pause(0.1);
        drawnow;
    end
    
    %% Governing equations
    res(1 : M) = qs - w_ts;
    res(M + 1 : 2 * M) = alpha * q_ts + (gamma / dx^4) * A_mat * ws - ps;
    t
end



%% mitchell_solution_stationary.m
% Code to use an iterative Mitchell FD method

clear;
close all;

addpath("pressures");
addpath("finite_differences");

%% Parameters
EPSILON = 1;
ALPHA = 10 / EPSILON^2; BETA = 0.0001 * EPSILON^2; GAMMA = 0.0001 * EPSILON^2; 
L = 8;
N_MEMBRANE = 2048;
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.25;
DELTA_T = 10^-4;
DELTA_X = L / (N_MEMBRANE - 1); 
T_VALS = 0 : DELTA_T : T_MAX;

% Derived parameters
if (GAMMA == 0)
    Cpressure = DELTA_X * DELTA_X / BETA;
else
    Cpressure = DELTA_X^4 / GAMMA;
end

%% Basilisk compare data
data_directory = "~/scratch/reflecting_waves/membrane_radius_4";
IMPACT_TIMESTEP = 1250;

%% Initialise arrays
xs = (0 : DELTA_X : L - DELTA_X)';
w_previous = zeros(size(xs));
w = zeros(size(xs));
w_next = zeros(size(xs));

p_previous_previous = zeros(size(xs));
p_previous = zeros(size(xs));
p = zeros(size(xs));

%% Stationary solutions to time-dependent terms
d_fun = @(t) 2 * sqrt(t);
d_t_fun = @(t) 1 / sqrt(t);
A_fun = @(t) d_fun(t) * d_t_fun(t);
C_fun = @(t) d_fun(t) * d_t_fun(t);
J_fun = @(t) pi * d_fun(t) / (8 * d_t_fun(t)^2);
    
%% Matrix definitions
[A_mat, B_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);
%% Initial conditions
t = 0
k = 0;

subplot(2, 1, 1);
plot(xs, w_previous);

subplot(2, 1, 2);
plot(xs, p_previous);
pause(0.01);


%% Solve for first w
t = t + DELTA_T
k = 1;

rhs = 0.5 * B_mat * w_previous;
w = A_mat \ rhs;

subplot(2, 1, 1);
plot(xs, w);

subplot(2, 1, 2);
plot(xs, p);
pause(0.01);


%% Loops
for k = 2 : length(T_VALS)
    
    %% Determines derivatives of w at current timestep
    w_t = (w - w_previous) / DELTA_T;
    
    % MAKE QUICKER BY TURNING INTO A MATRIX MULTIPLICATION
    w_tt = spatial_derivatives(w, ALPHA, BETA, GAMMA, M, DELTA_X) ...
        + (1 / ALPHA) * p_previous;
    
    w_x = zeros(M, 1);
    w_x(2 : M - 1) = (w(3 : M) - w(1 : M - 2)) / (2 * DELTA_X);
    w_x(M) = -w(M - 1) / (2 * DELTA_X);
    
    %% Determines anonymous functions for w and its derivatives
    w_fun = @(x) interp1(xs, w, x, 'linear', 'extrap');
    w_t_fun = @(x) interp1(xs, w_t, x, 'linear', 'extrap');
    w_tt_fun = @(x) interp1(xs, w_tt, x, 'linear', 'extrap');
    w_x_fun = @(x) interp1(xs, w_x, x, 'linear', 'extrap');
    
    %% Determine d and d_t
    [d, d_t] = turnover_point(t, w_fun, w_t_fun, w_x_fun, EPSILON);
    
    % Finds d_idx such that x(d_idx) < epsilon * d but x(d_idx) >= epsilon * d
    d_idx = sum(xs < EPSILON * d);
    
    %% Determine m and its derivatives
    if (d_idx < 2)
        m_t_fun = @(s) zeros(size(s));
        m_tt_fun = @(s) zeros(size(s));
    else
        s_vals = xs(1 : d_idx) / EPSILON;
        
        m_t = zeros(size(s_vals));
        m_tt = zeros(size(s_vals));
        
        % VECTORISE THIS FOR SPEEEEED
        for q = 1 : length(s_vals)
            m_t(q) = trapz(s_vals(1 : q), w_t(1 : q), 1);
            m_tt(q) = trapz(s_vals(1 : q), w_tt(1 : q), 1);
        end
        
        m_t_fun = @(s) interp1(s_vals, m_t, s, 'linear', 'extrap');
        m_tt_fun = @(s) interp1(s_vals, m_tt, s, 'linear', 'extrap');
    end
    
    
    %% Determine time-dependent quantities
    [A, C, J] = time_dependent_quantities(d, d_t, w_t_fun, w_tt_fun, m_t_fun, EPSILON);
    
    %% Determine pressure at current timestep
    p_stationary = composite_pressure_stationary(xs, t, d_fun(t), d_t_fun(t), A_fun(t), C_fun(t), J_fun(t), EPSILON);
    p = composite_pressure(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, EPSILON);
    
    % Restrict large values of p
    p(p > 1e4) = 0;
    
    %% Determine right-hand-side and solve for w_next
    rhs = B_mat * w - A_mat * w_previous ...
        + Cpressure * (p_previous_previous + 2 * p_previous + p);
    w_next = A_mat \ rhs;

    %% Swaps arrays
    % Swaps ws
    temp = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp;
    
    % Swaps ps
    temp = p_previous_previous;
    p_previous_previous = p_previous;
    p_previous = p;
    p = temp;

    %% Updates time
    t = T_VALS(k);
    t

    %% Plots
    % w plot
    subplot(3, 1, 1);
    plot(xs, w_next);
    
%     % Reads in Basilisk solution
%     membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", data_directory, k + IMPACT_TIMESTEP));
%     unsorted_xs = membrane_mat(:, 1);
%     unsorted_ws = membrane_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     ws = unsorted_ws(idxs);
%     hold on;
%     plot(sorted_xs, ws);
%     hold off;

    title(sprintf("t = %g, k = %d", t, k)); 
    
    % w_t plot
    subplot(3, 1, 2);
    plot(xs, w_t);
    
%     % Reads in Basilisk solution
%     membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", data_directory, k + IMPACT_TIMESTEP));
%     unsorted_xs = membrane_mat(:, 1);
%     unsorted_w_ts = membrane_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     w_ts = unsorted_w_ts(idxs);
%     hold on;
%     plot(sorted_xs, w_ts);
%     hold off;
    
    title(sprintf("t = %g, k = %d", t, k)); 
    
    
    subplot(3, 1, 3);
    plot(xs, p_previous);
    
%     % Reads in Basilisk solution
%     pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", data_directory, k + IMPACT_TIMESTEP));
%     unsorted_xs = pressure_mat(:, 1);
%     unsorted_ps = pressure_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     ps = unsorted_ps(idxs);
    hold on;
%     plot(sorted_xs, ps);
    plot(xs, p_stationary);
    hold off;
%     
%     ylim([0, 5 * p_previous(1)]);
    
    x0=400;
    y0=400;
    width=1200;
    height=800;
    
    set(gcf,'position',[x0,y0,width,height])
    
    pause(0.00001);
end

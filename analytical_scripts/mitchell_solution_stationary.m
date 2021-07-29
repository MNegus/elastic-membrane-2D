%% mitchell_solution_stationary.m
% Code to use an iterative Mitchell FD method

clear;
close all;

addpath("pressures");

%% Parameters
EPSILON = 1;
ALPHA = 2.0 / EPSILON^4; BETA = 1; GAMMA = 2; 
L = 4;
N_MEMBRANE = 1025;
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.25;
DELTA_T = 10^-4;
DELTA_X = L / (N_MEMBRANE - 1); 

% Derived parameters
Cpressure = DELTA_X^4 / GAMMA;
Calpha = ALPHA * DELTA_X^4 / (GAMMA * DELTA_T^2);
Cbeta = BETA * DELTA_X^2 / GAMMA;

%% Basilisk compare data
data_directory = "~/Desktop/decoupled_example";
IMPACT_TIMESTEP = 1250;

%% Initialise arrays
xs = (0 : DELTA_X : L - DELTA_X)';
w_previous = zeros(size(xs));
w = zeros(size(xs));
w_next = zeros(size(xs));
p_previous = zeros(size(xs));
p = zeros(size(xs));
p_next = zeros(size(xs));

%% Stationary solutions to time-dependent terms
d_fun = @(t) 2 * sqrt(t);
d_t_fun = @(t) 1 / sqrt(t);
A_fun = @(t) d_fun(t) * d_t_fun(t);
C_fun = @(t) d_fun(t) * d_t_fun(t);
J_fun = @(t) pi * d_fun(t) / (8 * d_t_fun(t)^2);
    
%% Matrix definitions
% A definition
A_upper_upper = ones(M, 1);
A_upper_upper(3) = 2;

A_upper = (-Cbeta - 4) * ones(M, 1);
A_upper(2) = -2 * Cbeta - 8;

A_main = (4 * Calpha + 2 * Cbeta + 6) * ones(M, 1);
A_main(2) = A_main(2) + 1;
A_main(M) = A_main(M) - 1;

A_lower = (-Cbeta - 4) * ones(M, 1);

A_lower_lower = ones(M, 1);

A_mat = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

% B definition
B_upper_upper = -2 * ones(M, 1);
B_upper_upper(3) = 2 * B_upper_upper(3);

B_upper = (2 * Cbeta + 8) * ones(M, 1);
B_upper(2) = 2 * B_upper(2);

B_main = (8 * Calpha - 4 * Cbeta - 12) * ones(M, 1);
B_main(2) = B_main(2) - 2;
B_main(M) = B_main(M) + 2;

B_lower = (2 * Cbeta + 8) * ones(M, 1);

B_lower_lower = -2 * ones(M, 1);

B_mat = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);

%% Initial conditions
t = 0
k = 0;

% a = 0.1;
% lambda = pi * 7 / (2 * L);
% w_previous = a * cos(lambda * xs);
w_previous = zeros(size(xs));

subplot(2, 1, 1);
plot(xs, w_previous);

subplot(2, 1, 2);
plot(xs, p_previous);
pause(0.01);


%% Solve for first w
t = t + DELTA_T
k = 1;

% Pressure at next timestep
p_next = composite_pressure_stationary(xs, t, d_fun(t), d_t_fun(t), ...
    A_fun(t), C_fun(t), J_fun(t), EPSILON);
% p_next = outer_pressure_stationary(xs, d_fun(t), A_fun(t), EPSILON);

rhs = 0.5 * B_mat * w_previous + Cpressure * p_next;
w = A_mat \ rhs;

subplot(2, 1, 1);
plot(xs, w);

subplot(2, 1, 2);
plot(xs, p);
pause(0.01);

% Swaps pressures
temp = p_previous;
p_previous = p;
p = p_next;
p_next = temp;


% Loops
while (t < T_MAX) 
    
    % Determines pressure at next timestep
    t_next = t + DELTA_T;
    p_next = composite_pressure_stationary(xs, t_next, d_fun(t_next), ...
        d_t_fun(t_next), A_fun(t_next), C_fun(t_next), J_fun(t_next), ...
        EPSILON);

%     p_next = outer_pressure_stationary(xs, d_fun(t_next), A_fun(t_next), EPSILON);
    
    % Restrict large values of p_next
    p_next(p_next > 1e4) = 0;
    
    rhs = B_mat * w - A_mat * w_previous ...
        + Cpressure * (p_previous + 2 * p + p_next);
    w_next = A_mat \ rhs;

    % Swaps ws
    temp = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp;
    
    % Swaps ps
    temp = p_previous;
    p_previous = p;
    p = p_next;
    p_next = temp;

    % Updates time
    t = t_next
    k = k + 1;

    %% Plots
    % w plot
    subplot(2, 1, 1);
    plot(xs, w_next);
    % Reads in Basilisk solution
    membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", data_directory, k + IMPACT_TIMESTEP));
    unsorted_xs = membrane_mat(:, 1);
    unsorted_ws = membrane_mat(:, 2);
    [sorted_xs, idxs] = sort(unsorted_xs);
    ws = unsorted_ws(idxs);
    hold on;
    plot(sorted_xs, ws);
    hold off;

    subplot(2, 1, 2);
    plot(xs, p_next);
    % Reads in Basilisk solution
    pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", data_directory, k + IMPACT_TIMESTEP));
    unsorted_xs = pressure_mat(:, 1);
    unsorted_ps = pressure_mat(:, 2);
    [sorted_xs, idxs] = sort(unsorted_xs);
    ps = unsorted_ps(idxs);
    hold on;
    plot(sorted_xs, ps);
    hold off;
    
    pause(0.001);



end
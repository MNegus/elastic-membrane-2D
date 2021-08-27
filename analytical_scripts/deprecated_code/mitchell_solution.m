%% mitchell_solution.m
% Code to use an iterative Mitchell FD method

clear;
close all;

addpath("pressures");

%% Parameters
EPSILON = 1;
ALPHA = 1.0 / EPSILON^4; BETA = 1; GAMMA = 1; 
L = 4;
N_MEMBRANE = 513;
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.5;
DELTA_T = 10^-4;
DELTA_X = L / (N_MEMBRANE - 1); 

% Derived parameters
Cpressure = DELTA_X^4 / GAMMA;
Calpha = ALPHA * DELTA_X^4 / (GAMMA * DELTA_T^2);
Cbeta = BETA * DELTA_X^2 / GAMMA;

%% Initialise arrays
xs = (0 : DELTA_X : L - DELTA_X)';
w_previous = zeros(size(xs));
w = zeros(size(xs));
w_next = zeros(size(xs));
p_previous = zeros(size(xs));
p = zeros(size(xs));
p_next = zeros(size(xs));
    
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
a = 0.1;
lambda = pi * 7 / (2 * L);
w_previous = a * cos(lambda * xs);

figure(1);
plot(xs, w_previous);
pause(0.01);

% Solve for first w
t = t + DELTA_T
rhs = 0.5 * B_mat * w_previous;
w = A_mat \ rhs;

plot(xs, w);
pause(0.01);

% Loops
while (t < T_MAX) 
    rhs = B_mat * w - A_mat * w_previous;
    w_next = A_mat \ rhs;

    % Swaps
    temp = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp;

    t = t + DELTA_T

    % Plots
    plot(xs, w);

    pause(0.01);

end

%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("../finite_differences");
addpath("../pressures");
addpath("../");

parent_dir = "/scratch/negus/second_order_tests";


%% Parameters
EPSILON = 1;
ALPHA = 2 / EPSILON^2; BETA = 1 * EPSILON^2; GAMMA = 2 * EPSILON^2; 
L = 4;
T_MAX = 0.05;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 1024;

%% Finite differences
fd_data_dir = sprintf("%s/finite_differences", parent_dir);
save_finite_differences_solution_second_order(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite")
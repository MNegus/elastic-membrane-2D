%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

parent_dir = "/media/michael/newarre/elastic_membrane/analytical_tests";


%% Parameters
EPSILON = 1;
ALPHA = 2 / EPSILON^2; BETA = 1 * EPSILON^2; GAMMA = 0 * EPSILON^2; 
L = 4;
T_MAX = 0.25;
DELTA_T = 1e-4;

% NM parameters
N = 128;

% FD parameters
N_MEMBRANE = 2056;

%% Normal modes
nm_data_dir = sprintf("%s/normal_modes", parent_dir);
save_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX, DELTA_T);
% save_validated_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);

%% Finite differences
% fd_data_dir = sprintf("%s/finite_differences", parent_dir);
% save_finite_differences_solution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T)
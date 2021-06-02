%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

%% Parameters
ALPHA = 0.002;
BETA = 7;
GAMMA = 0.0001;
% GAMMA = 0;
L = 4;
N_MEMBRANE = 128;
T_MAX = 0.1;
DELTA_T = 1e-4;
tvals = 0 : DELTA_T : T_MAX;
xs = linspace(0, L, N_MEMBRANE);

%% Numerically solves and compare

for t = tvals
    t
    % Exact solution
    w_exact_1 = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
    w_exact_2 = homogeneous_exact_solution(xs, t, ALPHA, BETA, 0, L);

    %  Plots full
    figure(1);
    plot(xs, w_exact_1);
    hold on
    plot(xs, w_exact_2);
    hold off;

    pause(0.01);

end



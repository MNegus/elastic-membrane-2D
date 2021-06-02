%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

clear;
close all;
%% Parameters

L = 1;
% N_MEMBRANES = [512, 1024, 2048, 4096];
% N_MEMBRANES = [512, 1024, 2048, 4096, 8192, 16384];
N_MEMBRANES = [4096, 8192, 16384];
Deltaxs =  L ./ (N_MEMBRANES - 1)
T_MAX = 0.1;
% DELTA_POWERS = linspace(-1, -5, 19)
DELTA_POWERS = -2 : -0.25 : -5;
DELTA_TS = 10.^DELTA_POWERS;
% DELTA_TS = 1e-4;

% ALPHA = 0.02; 
% BETA = ALPHA * min(Deltaxs)^2 / ((1e-2)^2)
% GAMMA = ALPHA * min(Deltaxs)^4 / ((1e-2)^2)
% 
% Dbetas = BETA * DELTA_TS.^2 ./ (ALPHA * Deltaxs'.^2)
% Dgammas = GAMMA * DELTA_TS.^2 ./ (ALPHA * Deltaxs'.^4)
ALPHA = 0.002;
BETA = 7;
GAMMA = 0.0001;


%% Numerically solves and compare

for N_MEMBRANE = N_MEMBRANES
    M = N_MEMBRANE - 1; % We ignore the end point
    xs = linspace(0, L, N_MEMBRANE);
    Deltax = L / (N_MEMBRANE - 1)
    max_errors = zeros(size(DELTA_TS));

    for k = 1 : length(DELTA_TS)
        % Reset figure
%         close(figure(1));
%         figure(1);
%         
%         close(figure(2));
%         figure(2);

        % Setting DELTA_T
        DELTA_T = DELTA_TS(k);
        Dbeta = BETA * DELTA_T^2 / (ALPHA * Deltax^2)
        Dgamma = GAMMA * DELTA_T^2 / (ALPHA * Deltax^4)

        % Initialises error term
        max_error = 0;

        % A definition
        A_upper_upper = Dgamma * ones(M, 1);
        A_upper_upper(3) = 2 * Dgamma;

        A_upper = (-Dbeta - 4 * Dgamma) * ones(M, 1);
        A_upper(2) = -2 * Dbeta - 8 * Dgamma;

        A_main = (4 + 2 * Dbeta + 6 * Dgamma) * ones(M, 1);
        A_main(2) = A_main(2) + Dgamma;
        A_main(M) = A_main(M) - Dgamma;

        A_lower = (-Dbeta - 4 * Dgamma) * ones(M, 1);

        A_lower_lower = Dgamma * ones(M, 1);

        A = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

        % B definition
        B_upper_upper = -2 * Dgamma * ones(M, 1);
        B_upper_upper(3) = 2 * B_upper_upper(3);

        B_upper = (2 * Dbeta + 8 * Dgamma) * ones(M, 1);
        B_upper(2) = 2 * B_upper(2);

        B_main = (8 - 4 * Dbeta - 12 * Dgamma) * ones(M, 1);
        B_main(2) = B_main(2) - 2 * Dgamma;
        B_main(M) = B_main(M) + 2 * Dgamma;

        B_lower = (2 * Dbeta + 8 * Dgamma) * ones(M, 1);

        B_lower_lower = -2 * Dgamma * ones(M, 1);

        B = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);

        % Initial conditions
        t = 0
        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        w_previous = w_exact(1 : end - 1);

%         plot(xs(1 : end - 1), w_previous);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
%         pause(0.01);

        % Solve for first w
        t = t + DELTA_T
        rhs = 0.5 * B * w_previous;
        w = A \ rhs;

        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);


%         plot(xs(1 : end - 1), w);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
%         pause(0.01);

%         Loops
        max_diff = 0
        while (t < T_MAX) 
            rhs = B * w - A * w_previous;
            w_next = A \ rhs;

            % Swaps
            temp = w_previous;
            w_previous = w;
            w = w_next;
            w_next = temp;

            t = t + DELTA_T

            % Exact solution
            w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);

            % Updates max errors 
            max_diff = max(abs(w - w_exact(1 : end - 1)));
            if (max_diff > max_error)
               max_error = max_diff; 
            end
            
%             % Plots both
%             figure(1);
%             plot(xs(1 : end - 1), w);
%             hold on;
%             plot(xs, w_exact);
%             hold off;
%             
%             % Plots diff
%             figure(2);
%             plot(xs(1 : end - 1), w - w_exact(1 : end - 1));
%             pause(0.01);

        end
        max_errors(k) = max_diff;
    end
    
    figure(3);
    hold on;
    plot(DELTA_TS, max_errors, '-o');
    
end

%% Change plot
figure(3);
set(gca, 'XDir','reverse');
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
legend(["N = 512", "N = 1024", "N = 2048", "N = 4096", "N = 8192"]);




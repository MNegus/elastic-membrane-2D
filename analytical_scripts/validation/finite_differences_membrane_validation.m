%% finite_differences_membrane_validation.m
% Validaes the solution for the membrane with the finite differences method

clear;
addpath("../finite_differences");
save_solutions = false;

%% Data directory
parent_dir = "/home/michael/Desktop/confirmation_validation/finite_differences_membrane_validation";

%% Loads in physical parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, ~] ...
    = parameters();

%% Physical parameters
N_imposed = 128;

tvals = 0 : DELTA_T: T_MAX;

lambda = @(n) pi * (2 * n - 1) / (2 * L);
lambdas_imposed = pi * (2 * (1 : N_imposed) - 1) / (2 * L);

% Makes array of N values
N_val = 16;
N_MEMBRANE_MAX = 8192;
N_MEMBRANES = [];
while N_val <= N_MEMBRANE_MAX
    N_MEMBRANES(end + 1) = N_val;
    N_val = N_val * 4;
end

%% Loops over parameters
if (save_solutions)
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            %% Make a directory for these parameters
            param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
            mkdir(param_dir);
            
            %% Loops over values of N_MEMBRANE
            for N_MEMBRANE = N_MEMBRANES
                DELTA_X = L / (N_MEMBRANE - 1); 
                xs = (0 : DELTA_X : L - DELTA_X)';
                
                % Makes directory for this value of N
                data_dir = sprintf("%s/N_MEMBRANE_%d", param_dir, N_MEMBRANE);
                mkdir(data_dir);
                
                % Saves numerical solution
                save_prescribed_finite_differences_solution(data_dir, ...
                    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, N_imposed)
                
            end
            
        end
    end
end
end

%% Compares each saved solution to the exact solution
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
            max_errors = [];
            ks = BETA * lambdas_imposed.^2 + GAMMA * lambdas_imposed.^4;
            
            %% Loops over values of N_MEMBRANE
            for N_MEMBRANE = N_MEMBRANES
                
                DELTA_X = L / (N_MEMBRANE - 1); 
                xs = (0 : DELTA_X : L - DELTA_X)';
                
                % Solution directory
                data_dir = sprintf("%s/N_MEMBRANE_%d", param_dir, N_MEMBRANE);
                
                %% Loops over time, saving difference between numerical and exact solutions
                max_diff = 0;
                for q = 2 : 10 : length(tvals)
                    t = tvals(q)
                    
                    % Finds exact solution for as
                    exact_as = ALPHA * (1 - cos(sqrt(ks) * t / sqrt(ALPHA))) ./ (ks .* sqrt(L * lambdas_imposed));
                    
                    % Reads in  numerical solution
                    ws_mat = matfile(sprintf("%s/w_%d.mat", data_dir, q));
                    
                    % Finds w solutions
                    exact_solution = w_solution(xs, exact_as, L, N_imposed);
                    numerical_solution = EPSILON^2 * ws_mat.w_next;
                    
                    % Finds absolute difference, and saves if largest
                    max_diff = max(max_diff, ...
                        max(abs(exact_solution - numerical_solution)));
                end
                
                % Update max_errors
                max_errors(end + 1, 1) = N_MEMBRANE;
                max_errors(end, 2) = max_diff;

            end
            
            %% Saves max_errors in directory
            save(sprintf("%s/max_errors.mat", param_dir), 'max_errors'); 
        end 
    end 
end

%% Plots errors
close(figure(1));
figure(1);
hold on;
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
            max_errors_dir = matfile(sprintf("%s/max_errors.mat", param_dir));
            max_errors = max_errors_dir.max_errors;
            plot(max_errors(:, 1), max_errors(:, 2), '-o', 'Displayname', sprintf("gamma = %g", GAMMA));
        end 
    end 
end
legend("interpreter", "latex", "fontsize", 12, "location", "southwest");
% grid on;
xlabel("$N$ (Number of $x$ grid points)", "interpreter", "latex", "fontsize", 15);
ylabel("Max norm error", "interpreter", "latex", "fontsize", 15);
set(gca, 'Ticklabelinterpreter', 'latex'); 
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
set(gca, 'fontsize', 12);
title("Finite differences, membrane validation: alpha = 1, beta = 1", "interpreter", "latex");
savefig(sprintf("%s/finite_differences_membrane_validation.fig", parent_dir));
exportgraphics(gca, sprintf("%s/finite_differences_membrane_validation.png", parent_dir), "resolution", 300);


function ws = w_solution(xs, as, L, N)
    
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
    
    % Find ws
    ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);
end

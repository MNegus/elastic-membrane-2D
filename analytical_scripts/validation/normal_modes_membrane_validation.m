%% normal_modes_membrane_validation.m
% Validates the solution for the membrane with the normal modes method

clear;
addpath("../normal_modes");
save_solutions = false;

%% Data directory
parent_dir = "/home/michael/Desktop/confirmation_validation/normal_modes_membrane_validation";

%% Loads in physical parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

tvals = 0 : DELTA_T : T_MAX;

N_imposed = 128;
lambda = @(n) pi * (2 * n - 1) / (2 * L);
lambdas_imposed = pi * (2 * (1 : N_imposed) - 1) / (2 * L);



%% Loops over parameters
if (save_solutions)
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            %% Make a directory for these parameters
            param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
            mkdir(param_dir);
            
            %% Loops over values of N, up to Nmax
            Nmax = N_stable(ALPHA, BETA, GAMMA, L, 10, DELTA_T);
            N = 2;
            while N < Nmax
                
                % Solution directory
                data_dir = sprintf("%s/N_%d", param_dir, N);
                mkdir(data_dir);
                
                % Saves numerical solution
                save_prescribed_normal_modes_solution(data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX, DELTA_T);
                
                % Doubles N
                N = N * 2;
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
            
            %% Loops over values of N, up to Nmax
            Nmax = N_stable(ALPHA, BETA, GAMMA, L, 10, DELTA_T);
            N = 2;
            while N < Nmax
                % Solution directory
                data_dir = sprintf("%s/N_%d", param_dir, N);
                
                % Opens numerical solution
                as_mat = matfile(sprintf("%s/as.mat", data_dir));
                as = as_mat.as;

                a_ts_mat = matfile(sprintf("%s/a_ts.mat", data_dir));
                a_ts = a_ts_mat.a_ts;

                q_ts_mat = matfile(sprintf("%s/q_ts.mat", data_dir));
                q_ts = q_ts_mat.q_ts;
                
                %% Loops over time, saving difference between numerical and exact solutions
                max_diff = 0;
                for q = 1 : 10 : length(tvals)
                    t = tvals(q)
                    
                    % Finds exact solution for as
                    exact_as = ALPHA * (1 - cos(sqrt(ks) * t / sqrt(ALPHA))) ./ (ks .* sqrt(L * lambdas_imposed));
                    
                    % Finds w solutions
                    exact_solution = w_solution(xs, exact_as, L, N_imposed);
                    numerical_solution = w_solution(xs, as(q, :), L, N);
                    
                    % Finds absolute difference, and saves if largest
                    max_diff = max(max_diff, ...
                        max(abs(exact_solution - numerical_solution)));
                end
                
                % Update max_errors
                max_errors(end + 1, 1) = N;
                max_errors(end, 2) = max_diff;
                
                % Doubles N
                N = N * 2
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
xline(N_imposed, 'linestyle', '--', 'Displayname', '$N_{imposed}$');
legend("interpreter", "latex", "fontsize", 12, "location", "southwest");
% grid on;
xlabel("$N$ (Number of modes)", "interpreter", "latex", "fontsize", 15);
ylabel("Max norm error", "interpreter", "latex", "fontsize", 15);
set(gca, 'Ticklabelinterpreter', 'latex'); 
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
set(gca, 'fontsize', 12);
title("Normal modes, membrane validation: alpha = 1, beta = 1", "interpreter", "latex");
savefig(sprintf("%s/normal_modes_membrane_validation.fig", parent_dir));
exportgraphics(gca, sprintf("%s/normal_modes_membrane_validation.png", parent_dir), "resolution", 300);


%% Function definitions
function ws = w_solution(xs, as, L, N)
    
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
    
    % Find ws
    ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);
end

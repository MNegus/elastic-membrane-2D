%% prescribed_finite_differences_validation.m
% Conducts validation for finite differences with prescribed pressure

clear;
close all;

addpath("../finite_differences");

parent_dir = "/media/michael/newarre/elastic_membrane/analytical_validation/finite_differences_validation";

%% Parameters
EPSILON = 1;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
DELTA_T = 1e-4;
IMPACT_TIMESTEP = 0.125 / DELTA_T
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
% N_MEMBRANE_MAX = 10924;
N_MEMBRANES = 2.^(7:14);

save_solution = false;

%% Loops over varying types and saves
figno = 0;
for varying = ["alpha", "beta", "gamma"]
    figno = figno + 1;
    
    if (save_solution)
        mkdir(sprintf("%s/%s_varying", parent_dir, varying));
    end
    
    %% Sets parameters
    if varying == "alpha"
        ALPHAS = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * ones(size(ALPHAS)) * EPSILON^2;
        Ne = 256;
    elseif varying == "beta"
        BETAS = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "gamma"
        GAMMAS = [0.5, 1, 2, 4, 8, 16, 32] * EPSILON^2;
        ALPHAS = 2 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
    no_params = length(ALPHAS);
    
    %% Loops over parameters and saves
    figure(figno);
    hold on;
    for idx = 1 : no_params
        %% Saves specific parameters
        ALPHA = ALPHAS(idx); BETA = BETAS(idx); GAMMA = GAMMAS(idx); 
        
        data_dir = sprintf("%s/%s_varying/alpha_%g-beta_%g-gamma_%g", parent_dir, varying, ALPHA, BETA, GAMMA);
        
        if (save_solution)
            mkdir(data_dir);
        end
        
        %% Finds Ne
        Nstable = N_stable(ALPHA, BETA, GAMMA, L, 10, 1e-4);
        Ne_power = floor(log2(0.75 * Nstable));
        if (varying == "alpha")
            Ne = 256;
            Ne_power = log2(Ne);
        else
            Ne = 2^Ne_power;
        end
        
        if (save_solution)
            save(sprintf("%s/Ne.mat", data_dir), 'Ne');
        end
        
        %% Determine exact solution
        lambdas = pi * (2 * (1 : Ne)' - 1) / (2 * L);
        ks = EPSILON^2 * (BETA * lambdas.^2 + GAMMA * lambdas.^4);
        ALPHAHAT = ALPHA / EPSILON^2;
        as_exact = @(t) ALPHAHAT * (1 - cos(sqrt(ks / ALPHAHAT) * t)) ./ (ks .* sqrt(L * lambdas));
        
        %% Loops over Nms and finds solution and compares to exact
        % Arrays for errors
        ws_errors = zeros(length(N_MEMBRANES), 1);
        if (save_solution)
            for N_idx = 1 : length(N_MEMBRANES)
                
                % Find N_MEMBRANE dependent variables
                N_idx
                N_MEMBRANE = N_MEMBRANES(N_idx);
                DELTA_X = L / (N_MEMBRANE - 1); 
                xs = (0 : DELTA_X : L - DELTA_X)';
                
                % Save parameter dir
                param_dir = sprintf("%s/N_%d", data_dir, N_MEMBRANE);
                mkdir(param_dir);
            
                % Find numerical solution
                freq = 100;
                save_prescribed_finite_differences_solution(param_dir, ...
                    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, Ne, freq)

                % Find exact solution with current xs
                ws_exact = @(t) w_solution(xs, as_exact(t), L, lambdas);
                
                % Loops over time and compares solutions for w
                for k = IMPACT_TIMESTEP + 1 : 100 : length(T_VALS)
                    t = T_VALS(k);

                    % Finds solution for w
                    % Reads in  numerical solution
                    ws_mat = matfile(sprintf("%s/w_%d.mat", param_dir, k - IMPACT_TIMESTEP));
                    ws = ws_mat.w_next;

                    % Compares solution to exact
                    error = max(abs(ws - ws_exact(t)));
                    ws_errors(N_idx) = max(ws_errors(N_idx), error);
                end
            end
            save(sprintf("%s/ws_errors.mat", data_dir), 'ws_errors');
        end
       
        %% Load errors back in and plot
        data_dir
        errors_mat = matfile(sprintf("%s/ws_errors.mat", data_dir));
        ws_errors = errors_mat.ws_errors;
        N_MEMBRANES
        loglog(N_MEMBRANES, ws_errors, 'Displayname', sprintf("idx = %d", idx));
        drawnow;
        
        %% OPTIONAL: Plots solutions
%         figure(1);
%         for k = IMPACT_TIMESTEP + 1 : 100 : length(T_VALS)
%             hold off;
%             t = T_VALS(k);
%             
%             for N = Nms
%                 % Load in solution
%                 ds_mat = matfile(sprintf("%s/N_%d/ds.mat", data_dir, N));
%                 ds_nm = ds_mat.ds;
%                 
%                 as_mat = matfile(sprintf("%s/N_%d/as.mat", data_dir, N));
%                 as = as_mat.as;
%         
%                 a_ts_mat = matfile(sprintf("%s/N_%d/a_ts.mat", data_dir, N));
%                 a_ts = a_ts_mat.a_ts;
%         
%                 q_ts_mat = matfile(sprintf("%s/N_%d/q_ts.mat", data_dir, N));
%                 q_ts = q_ts_mat.q_ts; 
%                 
%                 % Finds solution for w
%                 [ws, w_ts, ps] ...
%                     = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
%                         a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
%                         ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON);
%                     
%                 % Plot w
%                 plot(xs, ws, 'linewidth', 2, 'Displayname', sprintf("N = %d", N));
%                 hold on;
%             end
% 
%             plot(xs, w_exact(t), 'linewidth', 2, 'color', 'black', ...
%                 'linestyle', '--', 'Displayname', 'Exact');
%             legend;
%             drawnow;
%             title(sprintf("alpha = %g, beta = %g, gamma = %g, Ne = %d", ALPHA, BETA, GAMMA, Ne));
%             pause(0.01);
%         end
    end
    
    %% Set figure properties
    legend();
    set(gca, 'Yscale', 'log');
    set(gca, 'Xscale', 'log');
end

%% Function definitions
function ws = w_solution(xs, as, L, lambdas)
    ws = sum(cos(xs * lambdas') .* as', 2) / sqrt(L);
end


% 
% 
% %% Loops over parameters
% if (save_solutions)
% for ALPHA = ALPHAS
%     for BETA = BETAS
%         for GAMMA = GAMMAS
%             %% Make a directory for these parameters
%             param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%             mkdir(param_dir);
%             
%             %% Loops over values of N_MEMBRANE
%             for N_MEMBRANE = N_MEMBRANES
%                 DELTA_X = L / (N_MEMBRANE - 1); 
%                 xs = (0 : DELTA_X : L - DELTA_X)';
%                 
%                 % Makes directory for this value of N
%                 data_dir = sprintf("%s/N_MEMBRANE_%d", param_dir, N_MEMBRANE);
%                 mkdir(data_dir);
%                 
%                 % Saves numerical solution
%                 save_prescribed_finite_differences_solution(data_dir, ...
%                     ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, N_imposed)
%                 
%             end
%             
%         end
%     end
% end
% end
% 
% %% Compares each saved solution to the exact solution
% for ALPHA = ALPHAS
%     for BETA = BETAS
%         for GAMMA = GAMMAS
%             param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%             max_errors = [];
%             ks = BETA * lambdas_imposed.^2 + GAMMA * lambdas_imposed.^4;
%             
%             %% Loops over values of N_MEMBRANE
%             for N_MEMBRANE = N_MEMBRANES
%                 
%                 DELTA_X = L / (N_MEMBRANE - 1); 
%                 xs = (0 : DELTA_X : L - DELTA_X)';
%                 
%                 % Solution directory
%                 data_dir = sprintf("%s/N_MEMBRANE_%d", param_dir, N_MEMBRANE);
%                 
%                 %% Loops over time, saving difference between numerical and exact solutions
%                 max_diff = 0;
%                 for q = 2 : 10 : length(tvals)
%                     t = tvals(q)
%                     
%                     % Finds exact solution for as
%                     exact_as = ALPHA * (1 - cos(sqrt(ks) * t / sqrt(ALPHA))) ./ (ks .* sqrt(L * lambdas_imposed));
%                     
%                     % Reads in  numerical solution
%                     ws_mat = matfile(sprintf("%s/w_%d.mat", data_dir, q));
%                     
%                     % Finds w solutions
%                     exact_solution = w_solution(xs, exact_as, L, N_imposed);
%                     numerical_solution = EPSILON^2 * ws_mat.w_next;
%                     
%                     % Finds absolute difference, and saves if largest
%                     max_diff = max(max_diff, ...
%                         max(abs(exact_solution - numerical_solution)));
%                 end
%                 
%                 % Update max_errors
%                 max_errors(end + 1, 1) = N_MEMBRANE;
%                 max_errors(end, 2) = max_diff;
% 
%             end
%             
%             %% Saves max_errors in directory
%             save(sprintf("%s/max_errors.mat", param_dir), 'max_errors'); 
%         end 
%     end 
% end
% 
% %% Plots errors
% close(figure(1));
% figure(1);
% hold on;
% for ALPHA = ALPHAS
%     for BETA = BETAS
%         for GAMMA = GAMMAS
%             param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%             max_errors_dir = matfile(sprintf("%s/max_errors.mat", param_dir));
%             max_errors = max_errors_dir.max_errors;
%             plot(max_errors(:, 1), max_errors(:, 2), '-o', 'Displayname', sprintf("gamma = %g", GAMMA));
%         end 
%     end 
% end
% legend("interpreter", "latex", "fontsize", 12, "location", "southwest");
% % grid on;
% xlabel("$N$ (Number of $x$ grid points)", "interpreter", "latex", "fontsize", 15);
% ylabel("Max norm error", "interpreter", "latex", "fontsize", 15);
% set(gca, 'Ticklabelinterpreter', 'latex'); 
% set(gca, 'yscale', 'log');
% set(gca, 'xscale', 'log');
% set(gca, 'fontsize', 12);
% title("Finite differences, membrane validation: alpha = 1, beta = 1", "interpreter", "latex");
% savefig(sprintf("%s/finite_differences_membrane_validation.fig", parent_dir));
% exportgraphics(gca, sprintf("%s/finite_differences_membrane_validation.png", parent_dir), "resolution", 300);
% 
% 
% function ws = w_solution(xs, as, L, N)
%     
%     lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
%     
%     % Find ws
%     ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);
% end

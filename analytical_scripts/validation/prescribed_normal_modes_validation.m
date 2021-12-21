%% prescribed_normal_modes_validation.m
% Conducts validation for normal modes with prescribed pressure
clear;
close all;

addpath("../normal_modes");

parent_dir = "/media/michael/newarre/elastic_membrane/analytical_validation/normal_modes_validation";

%% Parameters
EPSILON = 1;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
DELTA_T = 1e-4;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
N_MEMBRANE = 10924;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;

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
        ALPHAS = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3] / EPSILON^2;
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
        
        %% Finds N range
        Nstable = N_stable(ALPHA, BETA, GAMMA, L, 10, 1e-4);
        Ne_power = floor(log2(0.75 * Nstable));
        if (varying == "alpha")
            Ne = 256;
            Ne_power = log2(Ne);
        else
            Ne = 2^Ne_power;
        end
        
        % Range of values of N
        Nms = 2.^(1 : Ne_power + 1);
        
        if (save_solution)
            save(sprintf("%s/Ne.mat", data_dir), 'Ne');
            save(sprintf("%s/Nms.mat", data_dir), 'Nms');
        end
        
        %% Determine exact solution
        lambdas = pi * (2 * (1 : Ne)' - 1) / (2 * L);
        ks = EPSILON^2 * (BETA * lambdas.^2 + GAMMA * lambdas.^4);
        ALPHAHAT = ALPHA / EPSILON^2;
        as_exact = @(t) ALPHAHAT * (1 - cos(sqrt(ks / ALPHAHAT) * t)) ./ (ks .* sqrt(L * lambdas));
        ws_exact = @(t) w_solution(xs, as_exact(t), L, lambdas);
        
        
        %% Loops over Nms and finds solution and compares to exact
        % Arrays for errors
        ws_errors = zeros(length(Nms), 1);
        if (save_solution)
            for N_idx = 1 : length(Nms)
                N_idx
                N = Nms(N_idx);
                param_dir = sprintf("%s/N_%d", data_dir, N);

                    mkdir(param_dir);
                    save_prescribed_normal_modes_solution(param_dir, ...
                       ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX - IMPACT_TIME, DELTA_T);


                % Loads back in solution
                ds_mat = matfile(sprintf("%s/N_%d/ds.mat", data_dir, N));
                ds_nm = ds_mat.ds;

                as_mat = matfile(sprintf("%s/N_%d/as.mat", data_dir, N));
                as = as_mat.as;

                a_ts_mat = matfile(sprintf("%s/N_%d/a_ts.mat", data_dir, N));
                a_ts = a_ts_mat.a_ts;

                q_ts_mat = matfile(sprintf("%s/N_%d/q_ts.mat", data_dir, N));
                q_ts = q_ts_mat.q_ts;

                % Loops over time and compares solutions for w
                for k = IMPACT_TIMESTEP + 1 : 100 : length(T_VALS)
                    t = T_VALS(k);

                    % Finds solution for w
                    [ws, w_ts, ps] ...
                        = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
                            a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
                            ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON);

                    % Compares solution to exact
                    error = max(abs(ws - ws_exact(t)));
                    ws_errors(N_idx) = max(ws_errors(N_idx), error);
                end
            end
            save(sprintf("%s/ws_errors.mat", data_dir), 'ws_errors');
        end
       
        %% Load errors back in and plot
        errors_mat = matfile(sprintf("%s/ws_errors.mat", data_dir));
        ws_errors = errors_mat.ws_errors;
        loglog(Nms, ws_errors);
        
        
        
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
    set(gca, 'Yscale', 'log');
    set(gca, 'Xscale', 'log');
end

%% Function definitions
function ws = w_solution(xs, as, L, lambdas)
    ws = sum(cos(xs * lambdas') .* as', 2) / sqrt(L);
end

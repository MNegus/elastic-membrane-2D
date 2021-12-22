%% prescribed_finite_differences_validation.m
% Conducts validation for finite differences with prescribed pressure

clear;
close all;

addpath("../finite_differences");
addpath("../normal_modes");
addpath("../");

parent_dir = "/media/michael/newarre/elastic_membrane/analytical_validation/finite_differences_validation";

fontsize = 28;

cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

%% Parameters
EPSILON = 1;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
N_MEMBRANE = 10924;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

DELTA_TS = [1e-2, 1e-3, 1e-4, 1e-5];
% DELTA_TS = [1e-5];

save_solution = false;

%% Loops over varying types and saves
figno = 0;
for varying = ["alpha", "beta", "gamma"]
% for varying = "alpha"
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
        
        displaynames = string(length(ALPHAS));
        for name_idx = 1 : length(ALPHAS)
           displaynames(name_idx) = "$\alpha = $" + num2str(ALPHAS(name_idx), "%.3f"); 
        end
        
        titlestr = "$\beta = 0$, $\gamma = 2$";
        
    elseif varying == "beta"
        BETAS = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
        
        displaynames = string(length(ALPHAS));
        for name_idx = 1 : length(ALPHAS)
           displaynames(name_idx) = "$\beta = $" + num2str(BETAS(name_idx)); 
        end
        
        titlestr = "$\alpha = 1$, $\gamma = 2$";
        
    elseif varying == "gamma"
        GAMMAS = [0.5, 1, 2, 4, 8, 16, 32] * EPSILON^2;
        ALPHAS = 2 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
        
        for name_idx = 1 : length(ALPHAS)
           displaynames(name_idx) = "$\gamma = $" + num2str(GAMMAS(name_idx)); 
        end
        
        titlestr = "$\alpha = 2$, $\beta = 0$";
    end
    no_params = length(ALPHAS);
    
    colors = ones(no_params, 3);
    
    color_idxs = floor(linspace(1, length(cmap), no_params));
    
    for q = 1 : no_params
        colors(q, :) = cmap(color_idxs(q), :);
    end
    
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
        elseif (varying == "gamma")
            Ne = 128;
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
        
        %% Loops over DELTA_Ts and finds solution and compares to exact
        % Arrays for errors
        ws_errors = zeros(length(DELTA_TS), 1);
        if (save_solution)
            for DELTA_idx = 1 : length(DELTA_TS)
                
                % Find DELTA_T dependent variables
                DELTA_idx
                DELTA_T = DELTA_TS(DELTA_idx);
                
                IMPACT_TIMESTEP = 0.125 / DELTA_T;
                T_VALS = 0 : DELTA_T : T_MAX - IMPACT_TIME;
                
                % Save parameter dir
                power = -log10(DELTA_T);
                param_dir = sprintf("%s/DELTA_T_1e-%d", data_dir, power);
                mkdir(param_dir);
            
                % Find numerical solution
                freq = floor(1e-2 / DELTA_T)
                save_prescribed_finite_differences_solution(param_dir, ...
                    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, Ne, freq)
                
                % Find exact solution with current xs
                ws_exact = @(t) w_solution(xs, as_exact(t), L, lambdas);
                
                max_timestep = floor((T_MAX - IMPACT_TIME) / DELTA_T)
                
                timesteps = 1 : freq : max_timestep
                
                % Loops over time and compares solutions for w
                for timestep_idx = 1 : length(timesteps)
                    timestep = timesteps(timestep_idx)
                    
                    t = T_VALS(timestep + 1);

                    % Finds solution for w
                    % Reads in  numerical solution
                    ws_mat = matfile(sprintf("%s/w_%d.mat", param_dir, timestep));
                    ws = ws_mat.w_next;

                    % Compares solution to exact
                    error = max(abs(ws - ws_exact(t)));
                    ws_errors(DELTA_idx) = max(ws_errors(DELTA_idx), error);
                end
            end
            save(sprintf("%s/ws_errors_temporal.mat", data_dir), 'ws_errors');
        end
       
        %% Load errors back in and plot
        data_dir
        errors_mat = matfile(sprintf("%s/ws_errors_temporal.mat", data_dir));
        ws_errors = errors_mat.ws_errors;
        DELTA_TS
        ws_errors
        loglog(DELTA_TS, ws_errors, '-o', 'Displayname', displaynames(idx), ...
            'color', colors(idx, :), 'markerfacecolor', colors(idx, :), ...
            'markersize', 10);
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
%     legend("interpreter", "latex", "fontsize", fontsize);
    set(gca, 'Yscale', 'log');
    set(gca, 'Xscale', 'log');
    set(gca, 'XDir', 'reverse')
    xticks([1e-5, 1e-4, 1e-3, 1e-2]);
    xlim([5e-6, 2e-2]);
    ylim([1e-7, 0.002]);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    
    xlabel("$\Delta t$", 'interpreter', 'latex');
    ylabel("$||w - w_{exact}||_\infty$", 'interpreter', 'latex');
%     ylim([2e-6, 0.015]);
    title(titlestr, "interpreter", "latex", "fontsize", fontsize);
    
    grid on;
    
    x0=400;
    y0=400;
    height=650;
    width=500;

    set(gcf,'position',[x0,y0,width,height]);
    
    name = sprintf("validation_figures/fd_temporal_validation_%s_varying", varying);
    exportgraphics(gcf, sprintf("%s.png", name), "Resolution", 300);
    savefig(gcf, sprintf("%s.fig", name));

end

%% Function definitions
function ws = w_solution(xs, as, L, lambdas)
    ws = sum(cos(xs * lambdas') .* as', 2) / sqrt(L);
end

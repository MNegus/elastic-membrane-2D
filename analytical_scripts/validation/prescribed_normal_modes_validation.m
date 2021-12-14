%% prescribed_normal_modes_validation.m
% Conducts validation for normal modes with prescribed pressure
close all;

addpath("../normal_modes");

parent_dir = "/media/michael/newarre/elastic_membrane/analytical_validation/normal_modes_validation";

%% Parameters
EPSILON = 1;
ALPHA = 2;
BETA = 1;
GAMMA = 2;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
DELTA_T = 1e-4;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
N_MEMBRANE = 10924;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;

save = false;

%% Loops over varying types and saves
for varying = ["alpha", "beta", "gamma"]
    mkdir(sprintf("%s/%s_varying", parent_dir, varying));
    
    %% Sets parameters
    if varying == "alpha"
        ALPHAS = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4] / EPSILON^2;
        BETAS = ones(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * ones(size(ALPHAS)) * EPSILON^2;
    elseif varying == "beta"
        BETAS = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
            
    elseif varying == "gamma"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 1 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
    no_params = length(ALPHAS);
    
    %% Loops over parameters and saves
    for idx = 1 : no_params
        %% Saves specific parameters
        ALPHA = ALPHAS(idx); BETA = BETAS(idx); GAMMA = GAMMAS(idx); 
        
        data_dir = sprintf("%s/%s_varying/alpha_%g-beta_%g-gamma_%g", parent_dir, varying, ALPHA, BETA, GAMMA);
        mkdir(data_dir);
        
        %% Finds N range
        Nstable = N_stable(ALPHA, BETA, GAMMA, L, 10, 1e-4);
        Ne_power = floor(log2(0.75 * Nstable));
        Ne = 2^Ne_power;
        
        % Range of values of N
        Nms = 2.^(1 : Ne_power + 1);
        
        if (save)
            save(sprintf("%s/Ne.mat", data_dir), 'Ne');
            save(sprintf("%s/Nms.mat", data_dir), 'Nms');
        end
        
        %% Loops over Nms and finds solution and compares to exact
        for N = Nms
           param_dir = sprintf("%s/N_%d", data_dir, N);
           
           if (save)
               mkdir(param_dir);
               save_prescribed_normal_modes_solution(param_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX - IMPACT_TIME, DELTA_T);
           else
               %% INSERT HERE: COMPARISON BETWEEN EXACT AND PRESCRIBED
           end
        end
    end
end

% %% Plots all solutions together
% figure(1);
% for k = IMPACT_TIMESTEP + 1 : 100 : length(T_VALS)
%     hold off;
%     for N = Nms
%         % Load in solution
%         ds_mat = matfile(sprintf("%s/N_%d/ds.mat", data_dir, N));
%         ds_nm = ds_mat.ds;
%         
%         as_mat = matfile(sprintf("%s/N_%d/as.mat", data_dir, N));
%         as = as_mat.as;
% 
%         a_ts_mat = matfile(sprintf("%s/N_%d/a_ts.mat", data_dir, N));
%         a_ts = a_ts_mat.a_ts;
% 
%         q_ts_mat = matfile(sprintf("%s/N_%d/q_ts.mat", data_dir, N));
%         q_ts = q_ts_mat.q_ts; 
%         
%         % Finds solution for w
%         [ws, w_ts, ps] ...
%             = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
%                 a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
%                 ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON);
%             
%         % Plot w
%         plot(xs, ws, 'Displayname', sprintf("N = %d", N));
%         hold on;
%     end
%     legend;
%     drawnow;
%     pause(0.01);
%     
% end
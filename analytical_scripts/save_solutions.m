%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD
close all;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

master_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";


%% Parameters
% [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
%     = parameters();

%% MANUAL PARAMETERS VARYING
EPSILON = 1;
L = 16;
T_MAX = 0.4;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

IMPACT_TIME = 0.125;

for varying = ["gamma"]
    parent_dir = sprintf("%s/%s_varying", master_dir, varying);
    
    %% Sets the parameters
    if varying == "alpha"
        ALPHAS = [1, 2, 4, 8] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "beta"
        BETAS = [0, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
            
    elseif varying == "gamma"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 1 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
%     [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
%     = parameters();
    no_params = length(ALPHAS);
%     
    parfor idx = 1 : no_params
        %% Saves specific parameters
        ALPHA = ALPHAS(idx); BETA = BETAS(idx); GAMMA = GAMMAS(idx); 
        
        data_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
        mkdir(data_dir);

        %% Finite differences
        fd_data_dir = sprintf("%s/finite_differences", data_dir);
        mkdir(fd_data_dir);

        composite_dir = sprintf("%s/composite", fd_data_dir);
        mkdir(composite_dir);
        save_finite_differences_solution(fd_data_dir, ...
            ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, ...
                "composite");
            
%         outer_dir = sprintf("%s/outer", fd_data_dir);
%         mkdir(outer_dir);
%         save_finite_differences_solution(fd_data_dir, ...
%             ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, ...
%                 "outer");
        
        %% Normal modes
%         nm_data_dir = sprintf("%s/normal_modes", data_dir);
%         mkdir(nm_data_dir);
%         N = 128;
%         save_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX, DELTA_T);
% %         save_validated_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX - IMPACT_TIME, DELTA_T);
    end
    
end



%% Normal modes
% for ALPHA = ALPHAS
%     for BETA = BETAS
%         for GAMMA = GAMMAS
%             data_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
% %             data_dir = parent_dir;
%             mkdir(data_dir);
% 
%             %% Normal modes
%             nm_data_dir = sprintf("%s/normal_modes", data_dir);
%             mkdir(nm_data_dir);
%             N = 128;
%             save_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX, DELTA_T);
% %             save_validated_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);
%         
% 
%         end
%     end
% end

            
%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD
close all;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

parent_dir = "/media/michael/newarre/elastic_membrane/scratch";

composite = true;
outer = false;
normal_modes = false;


%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters()

%% Finite differences
if ((outer) || (composite))
    for ALPHA_idx = 1 : length(ALPHAS)
        ALPHA = ALPHAS(ALPHA_idx);
        for BETA = BETAS
            for GAMMA = GAMMAS
    %         GAMMA = GAMMAS(ALPHA_idx);
                data_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
    %             data_dir = parent_dir;
                mkdir(data_dir);

                %% Finite differences
                fd_data_dir = sprintf("%s/finite_differences", data_dir);
                mkdir(fd_data_dir);

                if (composite)
                    composite_dir = sprintf("%s/composite", fd_data_dir);
                    mkdir(composite_dir);
                    tic
                    save_finite_differences_solution(fd_data_dir, ...
                        ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, ...
                        "composite")
                    toc
                end
                
                if (outer)
                    outer_dir = sprintf("%s/outer", fd_data_dir);
                    mkdir(outer_dir);
                    save_finite_differences_solution(fd_data_dir, ...
                        ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX - IMPACT_TIME, DELTA_T, ...
                        "outer")
                end

            end
        end
    end
end

%% Normal modes
if (normal_modes)
    for ALPHA = ALPHAS
        for BETA = BETAS
            for GAMMA = GAMMAS
                data_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
    %             data_dir = parent_dir;
                mkdir(data_dir);

                %% Normal modes
                nm_data_dir = sprintf("%s/normal_modes", data_dir);
                mkdir(nm_data_dir);
%                 N = floor(N_stable(ALPHA, BETA, GAMMA, L, 10, 1e-4))
%                 N = 128

%                 save_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX - IMPACT_TIME, DELTA_T);
                save_validated_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX - IMPACT_TIME, DELTA_T);


            end
        end
    end
end

            
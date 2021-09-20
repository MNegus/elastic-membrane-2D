%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

parent_dir = "/media/michael/newarre/elastic_membrane/beta_vary_test/analytical_data";


%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

%%
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            data_dir = sprintf("%s/alpha_%g_beta_%g_gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%             data_dir = parent_dir;
            mkdir(data_dir);

            %% Finite differences
            fd_data_dir = sprintf("%s/finite_differences", data_dir);
            mkdir(fd_data_dir);

            composite_dir = sprintf("%s/composite", fd_data_dir);
            mkdir(composite_dir);
            save_finite_differences_solution_first_order(fd_data_dir, ...
                ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, ...
                "composite")

        %     outer_dir = sprintf("%s/outer", fd_data_dir);
        %     mkdir(outer_dir);
        %     save_finite_differences_solution_first_order(fd_data_dir, ...
        %         ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, ...
        %         "outer")

        end
    end
end

%% Normal modes
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            data_dir = sprintf("%s/alpha_%g_beta_%g_gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%             data_dir = parent_dir;
            mkdir(data_dir);

            %% Normal modes
            nm_data_dir = sprintf("%s/normal_modes", data_dir);
            mkdir(nm_data_dir);
            % save_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, N, L, T_MAX, DELTA_T);
            save_validated_normal_modes_solution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);
        

        end
    end
end

            
% Plots the saved solutions using normal modes, FD and DNS
clear;
close all;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");


%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters()
ALPHA = ALPHAS(1);
BETA = BETAS(1);
GAMMA = GAMMAS(1);

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data";
analytical_parent_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
dns_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/dns", parent_dir, ALPHA, BETA, GAMMA);

%% Turnover matrices
fd_comp_mat = matfile(sprintf("%s/finite_differences/composite/ds.mat", analytical_parent_dir));
ds_comp = fd_comp_mat.ds;

%% Reads in DNS mass loss
mass_mat = dlmread(sprintf("%s/output.txt", dns_dir));
ts_dns = mass_mat(:, 1) - IMPACT_TIME;
mass_loss_dns = (pi / 2 - mass_mat(:, 2) / (2 * pi)) / (pi / 2);

%% Loops in time and finds the mass loss

timesteps = IMPACT_TIMESTEP + 2 : 10 : length(T_VALS);
mass_analytical = zeros(size(timesteps));

for timestep_idx = 1 : length(timesteps)
    k = timesteps(timestep_idx)
    length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t

    %% Loads in analytical solutions
    % Composites
    w_ts_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
    w_ts_comp = EPSILON^2 * w_ts_comp_mat.w_t;
            
    %% Works out fluxes
    if timestep_idx > 1
        d_idx = sum(xs < ds_comp(k - IMPACT_TIMESTEP));
        flux = trapz(xs(1 : d_idx), w_ts_comp(1 : d_idx)) * DELTA_T;
        mass_analytical(timestep_idx) = mass_analytical(timestep_idx - 1) + flux;
    end
    
end

%% 
close all;
figure(1);
hold on;
plot(ts_dns, mass_loss_dns / (2 * pi));
plot(timesteps * DELTA_T - IMPACT_TIME, mass_analytical / (pi / 2));
%% dns_pressure_validation.m
% Validates the pressure in the DNS by comparing its value at the origin
% for multiple levels

clear;
close all;

addpath("../");

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 18;

save_solutions = false;

%% Parameters
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
MAX_TIMESTEP = T_MAX / DELTA_T;
freq = 10;
timesteps = 1 : freq : MAX_TIMESTEP;
T_VALS = DELTA_T * timesteps - IMPACT_TIME;

%% Loads in directories
type = "moving";

if type == "stationary"
    parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/stationary";
else
    parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/alpha_2-beta_1-gamma_2";
end

levels = [9, 10, 11, 12, 13, 14, 15];

level_dirs = strings(length(levels), 1);

for idx = 1 : length(levels)
   level = levels(idx);
   
    if (level == 14)
        if type == "stationary"
            level_dirs(idx) = "/media/michael/newarre/elastic_membrane/stationary_membrane/membrane_outputs";
        else
            level_dirs(idx) = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/dns/membrane_outputs"
        end
    else
       level_dirs(idx) = sprintf("%s/max_level_%d/membrane_outputs", parent_dir, level);
    end
end

% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));

%% Saves pressures
if (save_solutions)
    for level_idx = 1 : length(levels)
        level = levels(level_idx);
        p0s = zeros(length(timesteps), 2);

        for t_idx = 1 : length(timesteps)
            k = timesteps(t_idx)
            t = T_VALS(t_idx);
            pressure_mat = dlmread(sprintf("%s/p_%d.txt", level_dirs(level_idx), k - 1));
            unsorted_xs = pressure_mat(:, 1);
            unsorted_ps = pressure_mat(:, 2);
            [sorted_xs, idxs] = sort(unsorted_xs);
            ps = unsorted_ps(idxs);
            p0s(t_idx, 1) = t;
            p0s(t_idx, 2) = ps(1);
        end

        save(sprintf("%s/origin_pressures_%d.mat", parent_dir, level), 'p0s');
    end
end

%% Plot solutions
figure(1);
hold on;
for level_idx = 1 : length(levels)
    level = levels(level_idx);
    
    p0s_mat = matfile(sprintf("%s/origin_pressures_%d.mat", parent_dir, level));
    p0s = p0s_mat.p0s;
    
    % Plot solution
    plot(p0s(:, 1), p0s(:, 2), 'color', cmap(color_idxs(level_idx), :), ...
        'linewidth', 2);
    drawnow; 
end

%%
grid on;
xlabel("$t$", 'interpreter', 'latex');
% Plots the saved solutions using normal modes, FD and DNS
clear;
close all;

% addpath("../finite_differences");
% addpath("../normal_modes");
% addpath("../pressures");
addpath("../");
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 16;
N_MEMBRANE = 10924;
DELTA_T = 1e-4;
T_MAX = 0.4;

EPSILON = 1;

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Data dirs
stationary_parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/stationary";
moving_parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/alpha_2-beta_1-gamma_2";

%% Levels and colors
levels = [10, 11, 12, 13, 14, 15];
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels, 3);
for q = 1 : no_levels
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Save jet areas
% jet_mat = dlmread(sprintf("%s/jet_area.txt", stationary_dir));
% stationary_jet_area = jet_mat(:, 2);
% 
% jet_mat = dlmread(sprintf("%s/jet_area.txt", moving_dir));
% moving_jet_area = jet_mat(:, 2);


%% Plot total areas
freq = 100;

figure(1);
hold on;

for level_idx = 1 : length(levels)
	level = levels(level_idx);
    
    if (level == 14)
        stationary_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane";
        moving_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/dns";
    else
        stationary_dir = sprintf("%s/max_level_%d", stationary_parent_dir, level);
        moving_dir = sprintf("%s/max_level_%d", moving_parent_dir, level);
    end
    
    area_mat = dlmread(sprintf("%s/output.txt", stationary_dir));
    ts_dns = area_mat(:, 1) - IMPACT_TIME;
    stationary_area = area_mat(:, 2) / (2 * pi);

    area_mat = dlmread(sprintf("%s/output.txt", moving_dir));
    ts_dns = area_mat(:, 1) - IMPACT_TIME;
    moving_area = area_mat(:, 2) / (2 * pi);
    
    scatter(ts_dns(1 : freq : end), stationary_area(1 : freq : end), [], colors(level_idx, :));
    scatter(ts_dns(1 : freq : end), moving_area(1 : freq : end), [], colors(level_idx, :), 'd');
    
end

yline(stationary_area(1), 'linestyle', '--', 'linewidth', 2, 'color', 'black')
grid on;
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
ylabel("Total area", 'interpreter', 'latex');
legend(["Stationary membrane", "Moving membrane", "Initial area"], 'location', 'southwest', 'interpreter', 'latex');
ylim([1.53, 1.58]);
xticks(-0.2 : 0.1 : 0.3);
xlim([-0.22, 0.3]);

x0=400;
y0=400;
height=400;
width=500;

set(gcf,'position',[x0,y0,width,height]);
% 
% exportgraphics(gcf, "validation_figures/total_area.png", "Resolution", 300);
% savefig(gcf, "validation_figures/total_area.fig");
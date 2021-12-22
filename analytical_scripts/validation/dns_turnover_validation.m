%% dns_pressure_validation.m
% Validates the pressure in the DNS by comparing its value at the origin
% for multiple levels

clear;
close all;

addpath("../");

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;


%% Parameters
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.4;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);


%% Loads in directories
for type = ["stationary", "moving"]
    close all;

    if type == "stationary"
        parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/stationary";
    else
        parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_validation/alpha_2-beta_1-gamma_2";
    end

    levels = [10, 11, 12, 13, 14, 15];

    level_dirs = strings(length(levels), 1);

    for idx = 1 : length(levels)
       level = levels(idx);

        if (level == 14)
            if type == "stationary"
                level_dirs(idx) = "/media/michael/newarre/elastic_membrane/stationary_membrane";
            else
                level_dirs(idx) = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/dns";
            end
        else
           level_dirs(idx) = sprintf("%s/max_level_%d", parent_dir, level);
        end
    end

    % Set up colors
    color_idxs = floor(linspace(1, length(cmap), length(levels)));
    no_levels = length(levels);
    colors = ones(no_levels - 1, 3);
   
    
    for q = 1 : no_levels - 1
        colors(q, :) = cmap(color_idxs(q), :);
    end

    %% Load in analytical solution
    if type == "stationary"
        ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
        ds_analytical = 2 * sqrt(ts_analytical);
    else
        analytical_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/finite_differences/composite";
        ds_mat = matfile(sprintf("%s/ds.mat", analytical_dir));
        ds_analytical = ds_mat.ds;
        ds_analytical = ds_analytical(1 : end - 1);
        ts_analytical = DELTA_T * (0 : length(ds_analytical) - 1);
    end

    %% Plot turnover points
    figure(1);
    hold on;

    dns_turnovers = zeros(NO_TIMESTEPS, length(levels));

    freq = 50;

    for level_idx = 1 : length(levels)
        level = levels(level_idx);

        turnover_mat = dlmread(sprintf("%s/turnover_points_basilisk.txt", level_dirs(level_idx)));
        ts = turnover_mat(:, 1) - IMPACT_TIME;
        size(ts)
        ds = turnover_mat(:, 2);

        dns_turnovers(:, level_idx) = ds;

        scatter(ts(1 : freq : end), ds(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
            'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
    end

    plot(ts_analytical, ds_analytical, 'linestyle', '--', 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');

    grid on;
    ylim([0, 1.2]);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
    ylabel("$d_m(t)$", 'interpreter', 'latex');
    legend('location', 'southeast', 'interpreter', 'latex');

    %% Plot L2 norm error
    axes('Position',[.32 .6 .3 .3])
    box on
    hold on;
    ds_15 = dns_turnovers(:, length(levels));


    norms = zeros(length(levels) - 1, 1);
    for level_idx = 1 : length(levels) - 1
        level = levels(level_idx);

        % Determines L2-norm difference
        diff = dns_turnovers(:, level_idx) - ds_15;
        sqrt(sum(diff.^2) / length(diff))
        norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
    end
    plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);
    
    sz = 100;
    scatter(levels(1 : end - 1), norms, sz, colors, 'filled');
    
    set(gca, 'yscale', 'log');
    xlim([min(levels) - 0.5, max(levels) - 0.5])
    ylim([10^-3, 10^(-0.5)]);
    xticks(levels(1 : end - 1));
    yticks([10^-3, 10^-2, 10^-1]);

    grid on;
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$m$", 'interpreter', 'latex');
    ylabel("$||d_m - d_{15}||_2$", 'interpreter', 'latex');

    %% Overal figure settings
    x0=400;
    y0=400;
    height=800;
    width=600;

    set(gcf,'position',[x0,y0,width,height])

    exportgraphics(gcf, sprintf("dns_turnover_%s.png", type), "Resolution", 300);
    savefig(gcf, sprintf("dns_turnover_%s.fig", type));
end
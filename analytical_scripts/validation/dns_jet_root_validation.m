close all;
clear;

%% Data loading
master_dir = "/scratch/negus/jet_root_height_validation";
MAXLEVELS = [10, 11, 12, 13];

%% Plotting parameters
markers = ['s', 'd', 'o', 'v'];
sz = 36;
freq = 50;

%% Physical parameters
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.3;

% Analytical solutions
analytical_ts = 0 : DELTA_T : T_MAX - IMPACT_TIME;
analytical_ds = 2 * sqrt(analytical_ts);
analytical_d_ts = 1 ./ sqrt(analytical_ts);
analytical_heights = pi * analytical_ds ./ ( 8 * analytical_d_ts.^2) * (1 + 4 / pi);
analytical_fluxes = pi * ones(size(analytical_ts));
analytical_energy = pi * analytical_ts;

%% Surface conditions
sub_dirs = ["no_surface", "dirichlet_surface", "neumann_surface"];
conditions = ["No applied surface condition", "No-slip surface condition", "Free slip surface condition"];

%% Loop over surface conditions and plots
figno = 1;
for condition_idx = 1 : length(sub_dirs)
    parent_dir = sprintf("%s/%s", master_dir, sub_dirs(condition_idx));

    %% Plots jet root x position validation
    close(figure(figno));
    figure(figno);
    figno = figno + 1;
    hold on;
    plot(analytical_ts, analytical_ds, 'linestyle', '--', 'linewidth', 2, 'color', 'black', 'Displayname', 'Analytical');

    for idx = 1 : length(MAXLEVELS)
        MAXLEVEL = MAXLEVELS(idx);
        filename = sprintf("%s/max_level_%d/raw_data/turnover_points_basilisk.txt", parent_dir, MAXLEVEL);
        data = readmatrix(filename);
        ts = data(:, 1) - IMPACT_TIME;
        ds = data(:, 2);
        Js = data(:, 3);
        d_ts = data(:, 4);
        fluxes = data(:, 6);
        energies = data(:, 7);


        scatter(ts(1 : freq : end), ds(1 : freq : end), ...
            sz, markers(idx), 'filled', 'Markeredgecolor', [0 0 0],  ...
            'Displayname', sprintf("Max level = %d", MAXLEVEL));

    end
    xlim([-0.01, max(ts)]);
    legend('Location', 'Northwest', 'Interpreter', 'latex', 'Fontsize', 12);
    xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 12);
    ylabel('Jet root $x$ position, $d(t)$', 'Interpreter', 'latex', 'Fontsize', 12);
    set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
    grid on;
    title(conditions(condition_idx), "Interpreter", "latex", "Fontsize", 12);
    
    
    savefig(sprintf("figures/fig_files/jet_root_x_position_%s.fig", sub_dirs(condition_idx))); 
    exportgraphics(gca, sprintf("figures/png_files/jet_root_x_position_%s.png", sub_dirs(condition_idx)));

    %% Plots jet root height validation
    close(figure(figno));
    figure(figno);
    figno = figno + 1;
    hold on;
    plot(analytical_ts, analytical_heights, 'linestyle', '--', 'linewidth', 2, 'color', 'black', 'Displayname', 'Analytical');

    for idx = 1 : length(MAXLEVELS)
        MAXLEVEL = MAXLEVELS(idx);
        filename = sprintf("%s/max_level_%d/raw_data/turnover_points_basilisk.txt", parent_dir, MAXLEVEL);
        data = readmatrix(filename);
        ts = data(:, 1) - IMPACT_TIME;
        heights = data(:, 3);

        scatter(ts(1 : freq : end), heights(1 : freq : end), ...
            sz, markers(idx), 'filled', 'Markeredgecolor', [0 0 0], ...
            'Displayname', sprintf("Max level = %d", MAXLEVEL));

    end
    xlim([-0.01, max(ts)]);
    legend('Location', 'Northwest', 'Interpreter', 'latex', 'Fontsize', 12);
    xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 12);
    ylabel('Jet root height', 'Interpreter', 'latex', 'Fontsize', 12);
    set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
    grid on;
    title(conditions(condition_idx), "Interpreter", "latex", "Fontsize", 12);
    savefig(sprintf("figures/fig_files/jet_root_height_%s.fig", sub_dirs(condition_idx))); 
    exportgraphics(gca, sprintf("figures/png_files/jet_root_height_%s.png", sub_dirs(condition_idx)));
    
    %% Plots jet root x velocity validation
    close(figure(figno));
    figure(figno);
    figno = figno + 1;
    hold on;
    plot(analytical_ts, analytical_d_ts, 'linestyle', '--', 'linewidth', 2, 'color', 'black', 'Displayname', 'Analytical');

    for idx = 1 : length(MAXLEVELS)
        MAXLEVEL = MAXLEVELS(idx);
        filename = sprintf("%s/max_level_%d/raw_data/turnover_points_basilisk.txt", parent_dir, MAXLEVEL);
        data = readmatrix(filename);
        ts = data(:, 1) - IMPACT_TIME;
        d_ts = data(:, 4);

        scatter(ts(1 : freq : end), d_ts(1 : freq : end), ...
            sz, markers(idx), 'filled', 'Markeredgecolor', [0 0 0], ...
            'Displayname', sprintf("Max level = %d", MAXLEVEL));

    end
    xlim([-0.01, max(ts)]);
    ylim([0, 10]);
    legend('Location', 'Northeast', 'Interpreter', 'latex', 'Fontsize', 12);
    xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 12);
    ylabel("Jet root $x$ velocity, $d'(t)$", 'Interpreter', 'latex', 'Fontsize', 12);
    set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
    grid on;
    title(conditions(condition_idx), "Interpreter", "latex", "Fontsize", 12);
    savefig(sprintf("figures/fig_files/jet_root_x_velocity_%s.fig", sub_dirs(condition_idx))); 
    exportgraphics(gca, sprintf("figures/png_files/jet_root_x_velocity_%s.png", sub_dirs(condition_idx)));
end




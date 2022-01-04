%% plot_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

% Options (set to 0 if don't want to plot the solution)
normal_modes = 1;
finite_differences_comp = 1;
finite_differences_outer = 0;
dns = 1;

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();
ALPHA = 2
BETA = 1
GAMMA = 2

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

fontsize = 22;

%% Timesteps to plot
% timesteps = ((IMPACT_TIME + [0.01, 0.1, 0.2]) / DELTA_T) + 1;
times = IMPACT_TIME + [0.01, 0.075, 0.2]
timesteps = ceil(times / DELTA_T) + 1

%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data";
analytical_parent_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
% dns_dir = "/media/michael/newarre/elastic_membrane/gamma_vary_test/basilisk_data/gamma_0.1";
dns_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/dns", parent_dir, ALPHA, BETA, GAMMA);

%% Loads in normal modes solutions
if (normal_modes)
    N_mat = matfile(sprintf("%s/normal_modes/N.mat", analytical_parent_dir));
    N = N_mat.N;
    
    as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
    as = as_mat.as;

    a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
    a_ts = a_ts_mat.a_ts;

    q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
    q_ts = q_ts_mat.q_ts;
end

%% Loads in turnover points
if (dns)
    % DNS turnover points
    dns_mat = dlmread(sprintf("%s/raw_data/turnover_points_basilisk.txt", dns_dir));
    ds_dns = dns_mat(:, 2);
    ts_dns = (dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;
end

if (normal_modes)
    % Normal modes turnover points
    nm_mat = matfile(sprintf("%s/normal_modes/ds.mat", analytical_parent_dir));
    ds_nm = nm_mat.ds;
end

if (finite_differences_comp)
    % FD turnover points
    fd_comp_mat = matfile(sprintf("%s/finite_differences/composite/ds.mat", analytical_parent_dir));
    ds_comp = fd_comp_mat.ds;
end

if (finite_differences_outer)
    fd_outer_mat = matfile(sprintf("%s/finite_differences/outer/ds.mat", analytical_parent_dir));
    ds_outer = fd_outer_mat.ds;
end

%% Colors of lines, depending on time
% color_mags = linspace(0, 0.5, length(timesteps));
% color_mags = [0, 0.5, 0.75];
color_mags = [0, 0, 0];
colors = ones(length(timesteps), 3);
for k = 1 : length(timesteps)
   colors(k, :) = color_mags(k) * colors(k, :); 
end

%%
% Animated lines
close all;


tiledlayout(2, 3, 'TileSpacing','Compact');


for timestep_idx = 1 : length(timesteps)
    %% Updates time
    k = timesteps(timestep_idx);
    t = T_VALS(k);
    t
    color = colors(timestep_idx, :);

    %% Loads in analytical solutions
    if (t > 0)
        if (normal_modes)
            % Normal modes
            [ws_nm, w_ts_nm, ps_nm] ...
                = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
                a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
                ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
        end
        
        if (finite_differences_comp)
            % Composites
            ws_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_comp = EPSILON^2 * ws_comp_mat.w_next;

            w_ts_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_comp = EPSILON^2 * w_ts_comp_mat.w_t;

            ps_comp_mat = matfile(sprintf("%s/finite_differences/composite/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_comp = ps_comp_mat.p;
        end
        
        if (finite_differences_outer)
            % Outers
            ws_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_outer = ws_outer_mat.w_next;

            w_ts_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_outer = w_ts_outer_mat.w_t;

            ps_outer_mat = matfile(sprintf("%s/finite_differences/outer/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_outer = ps_outer_mat.p;
        end
    end
        
    %% w plot
%     subplot(2, 3, timestep_idx);
    nexttile(timestep_idx);
    hold on

    % DNS
    if (dns)
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        plot(sorted_xs, ws, 'Color', color,'LineWidth', 2);
    end

    if (t > 0)
        if (finite_differences_comp)
            plot(xs, ws_comp, 'Color', color,'LineWidth', 2, ...
                'linestyle', '--');
        end
        
        if (normal_modes)
            plot(xs, ws_nm, 'Color', color,'LineWidth', 2, 'linestyle', ':');
        end
        
        if (finite_differences_outer)
            plot(xs, ws_outer, 'Color', color,'LineWidth', 2, ...
                'Displayname', 'FD: Outer');
        end
    end

    xlabel("$x$", "interpreter", "latex", "Fontsize", fontsize);
    ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", fontsize);
    ymax = max(ws_nm) * 1.1;
    ymin = -ymax / 5;
    ylim([ymin, ymax]);
    xlim([0, 3]);
    ax = gca;
    ax.YAxis.Exponent = 0;
%     ax.YAxis.TickLabelFormat = '%.4f';
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    grid on;
    title(sprintf("$t$ = %.2f", t), "Interpreter", "Latex", 'Fontsize', fontsize);
    
    %% w_t plot
%     subplot(2, 3, timestep_idx + 3);
    nexttile(timestep_idx + 3);
    hold on;

    % DNS
    if (dns)
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_w_ts = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        w_ts = unsorted_w_ts(idxs);
        plot(sorted_xs, w_ts, 'Color', color,'LineWidth', 2);
        
    end

    if (t > 0)
        if (finite_differences_comp)
            plot(xs, w_ts_comp, 'Color', color,'LineWidth', 2, ...
                'linestyle', '--');
        end
        
        if (normal_modes)
            plot(xs, w_ts_nm, 'Color', color,'LineWidth', 2, 'linestyle', ':');
        end
        
        if (finite_differences_outer)
            plot(xs, w_ts_outer, 'Color', color,'LineWidth', 2, ...
                'Displayname', 'FD: Outer');
        end
    end

    xlabel("$x$", "interpreter", "latex", "Fontsize", fontsize);
    ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", fontsize);
    
    ymax = max(w_ts_nm) * 1.1;
    ymin = -ymax / 2;
    ylim([ymin, ymax]);
    xlim([0, 3]);
    grid on
    
    if (timestep_idx == 1)
        yticks([-0.04, 0, 0.04, 0.08]);
    end
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize, 'fontsize', fontsize);

    %% Figure settings
    x0=400;
    y0=400;
    width=1200;
    height=700;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;


end

%% Set legend 
nexttile(2);
L = legend(["DNS ", "Analytical: Finite differences", "Analytical: Normal modes "], ...
    'Orientation','horizontal', 'Interpreter', 'latex', 'fontsize', fontsize, ...
    'Location', 'northoutside');
exportgraphics(gcf, "figures/membrane_model_comparison.png", "Resolution", 300);
savefig(gcf, "figures/membrane_model_comparison.fig");
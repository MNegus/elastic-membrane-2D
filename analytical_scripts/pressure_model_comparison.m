%% plot_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();
ALPHA = 2;
BETA = 1;
GAMMA = 2;

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
times = IMPACT_TIME + [0.01, 0.075, 0.2]
timesteps = ceil(times / DELTA_T) + 1

%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data";
analytical_parent_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
dns_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/dns", parent_dir, ALPHA, BETA, GAMMA);
stationary_dns_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane";

%% Loads in normal modes solutions
N_mat = matfile(sprintf("%s/normal_modes/N.mat", analytical_parent_dir));
N = N_mat.N;

as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
q_ts = q_ts_mat.q_ts;

ds_mat = matfile(sprintf("%s/normal_modes/ds.mat", analytical_parent_dir));
ds_nm = ds_mat.ds;

%% Colors of lines, depending on time
% color_mags = linspace(0, 0.75, length(timesteps));
% color_mags = [0, 0.5, 0.75];
color_mags = [0, 0, 0];
colors = ones(length(timesteps), 3);
for k = 1 : length(timesteps)
   colors(k, :) = color_mags(k) * colors(k, :); 
end

%% Plots lines
% Animated lines
close all;


% tiledlayout(3, 1,'TileSpacing','Compact');
figure(1);
hold on;

for timestep_idx = 1 : length(timesteps)
    %% Updates time
    k = timesteps(timestep_idx);
    t = T_VALS(k);
    t
    color = colors(timestep_idx, :);

    %% Loads in analytical solutions
    if (t > 0)
        % Normal modes
        [ws_nm, w_ts_nm, ps_nm] ...
            = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
            a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
            ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
        
        % Composite
        ps_comp_mat = matfile(sprintf("%s/finite_differences/composite/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
        ps_comp = ps_comp_mat.p;
        
        % Outer
        ps_outer_mat = matfile(sprintf("%s/finite_differences/outer/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
        ps_outer = ps_outer_mat.p;
    end
        
%     nexttile;
%     hold on;
    
    %% Moving membrane plots
    % DNS plot
    pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", dns_dir, k - 1));
    unsorted_xs = pressure_mat(:, 1);
    unsorted_ps = pressure_mat(:, 2);
    [sorted_xs, idxs] = sort(unsorted_xs);
    ps = unsorted_ps(idxs);
    plot(sorted_xs, ps, 'Color', color,'LineWidth', 2);
    
    % Analytical plot
    if (t > 0)
        plot(xs, ps_comp, 'Color', color,'LineWidth', 2, ...
            'linestyle', '--');
        
%         plot(xs, ps_outer, 'Color', color,'LineWidth', 2, ...
%             'linestyle', '-.');
        plot(xs, ps_nm, 'Color', color,'LineWidth', 2, 'linestyle', ':');
    end

    %% Stationary membrane plots
%    % DNS plot
%     pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", stationary_dns_dir, k - 1));
%     unsorted_xs = pressure_mat(:, 1);
%     unsorted_ps = pressure_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     ps = unsorted_ps(idxs);
%     plot(-sorted_xs, ps, 'Color', color,'LineWidth', 2);
    
    

end

%% Figure properties
xlabel("$x$", "interpreter", "latex", "Fontsize", fontsize);
ylabel("$p(x, -w(x, t), t)$", "interpreter", "latex", "Fontsize", fontsize);
ymax = 35;
xmax = 1;
ylim([0, ymax]);
xlim([0, xmax]);
grid on;

set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
x0=400;
y0=400;
width=1200;
height=400;
set(gcf,'position',[x0,y0,width,height])

%% Time arrow
p1 = [0.32, 0.75]
p2 = [27, 13];
x = p1 / xmax;
y = p2 / ymax;
annotation('arrow', x, y);

text(0.5, 25, '$t$', 'Interpreter', 'latex', 'fontsize', fontsize);


%% Set legend 
L = legend(["DNS ", "Analytical: Finite differences", "Analytical: Normal modes "], ...
    'Orientation','horizontal', 'Interpreter', 'latex', 'fontsize', fontsize, ...
    'Location', 'northoutside');
exportgraphics(gcf, "figures/pressure_model_comparison.png", "Resolution", 300);
savefig(gcf, "figures/pressure_model_comparison.fig");
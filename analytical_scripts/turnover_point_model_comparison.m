%% plot_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

%% Parameters
% Basilisk parameters
DELTA_T = 1e-4;
T_MAX = 0.4;
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = (0 : DELTA_T : T_MAX - IMPACT_TIME)';

%% Timesteps to plot
timesteps = IMPACT_TIMESTEP - 200 : 650 : length(T_VALS) - 500;

%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data";
stationary_dns_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane";

%% Analytical turnover points
ds_analytical = 2 * sqrt(ts_analytical);

%% Stationary DNS turnover points
dns_mat = dlmread(sprintf("%s/raw_data/turnover_points_basilisk.txt", stationary_dns_dir));
ds_dns = dns_mat(IMPACT_TIMESTEP + 1 : end, 2);
ts_dns = dns_mat(IMPACT_TIMESTEP + 1 : end, 1) - IMPACT_TIMESTEP * DELTA_T;

%% Plot lines
close all;
figure(1);
hold on;

% Plot turnover point
h(1) = plot(ts_dns, ds_dns, 'color', 'black');
h(2) = plot(ts_analytical, ds_analytical, 'linestyle', '--', 'color', 'black');
xlabel("$t$", "interpreter", "latex", "Fontsize", 18);
ylabel("$d(t)$", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);

% Plot error
yyaxis right
err_color = 0.6 * [1 1 1];
errs = (ds_analytical - ds_dns) ./ ds_dns;
h(3) = plot(ts_analytical, errs, 'color', err_color);
h(4) = yline(0.05, 'linestyle', ':');
ylabel("Error", "interpreter", "latex", "Fontsize", 18);
set(gca, 'YColor', err_color);


%% Figure properties
set(gcf, 'position', [400, 400, 1024, 400]);
legend(h(1:2), ["DNS", "Analytical"], "Interpreter", "latex", ...
    "fontsize", 15, "Location", "Northwest");


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
stationary_dns_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane";
moving_dns_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/dns";

%% Save total areas
area_mat = dlmread(sprintf("%s/output.txt", stationary_dns_dir));
ts_dns = area_mat(:, 1) - IMPACT_TIME;
stationary_area = area_mat(:, 2) / (2 * pi);

area_mat = dlmread(sprintf("%s/output.txt", moving_dns_dir));
ts_dns = area_mat(:, 1) - IMPACT_TIME;
moving_area = area_mat(:, 2) / (2 * pi);

%% Save jet areas
jet_mat = dlmread(sprintf("%s/jet_area.txt", stationary_dns_dir));
stationary_jet_area = jet_mat(:, 2);

jet_mat = dlmread(sprintf("%s/jet_area.txt", moving_dns_dir));
moving_jet_area = jet_mat(:, 2);


%% Plot total areas
freq = 100;

figure(1);
hold on;
% plot(ts_dns, stationary_area, 'linewidth', 2, 'color', cmap(1, :));
% plot(ts_dns, moving_area, 'linewidth', 2, 'color', cmap(end, :));
scatter(ts_dns(1 : freq : end), stationary_area(1 : freq : end), [], cmap(1, :));
scatter(ts_dns(1 : freq : end), moving_area(1 : freq : end), [], cmap(end, :));
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

exportgraphics(gcf, "validation_figures/total_area.png", "Resolution", 300);
savefig(gcf, "validation_figures/total_area.fig");

%% Plot total jet areas
figure(2);
hold on;
% plot(ts_dns, stationary_jet_area, 'linewidth', 2, 'color', cmap(1, :));
% plot(ts_dns, moving_jet_area, 'linewidth', 2, 'color', cmap(end, :));
scatter(ts_dns(1 : freq : end), stationary_jet_area(1 : freq : end), [], ... 
    cmap(1, :), 'd');
scatter(ts_dns(1 : freq : end), moving_jet_area(1 : freq : end), [], ...
    cmap(end, :), 'd');
grid on;
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
ylabel("Jet area", 'interpreter', 'latex');
legend(["Stationary membrane", "Moving membrane", "Initial area"], 'location', 'northwest', 'interpreter', 'latex');
xticks(-0.2 : 0.1 : 0.3);
yticks([0, 0.005, 0.01, 0.015, 0.02])
xlim([-0.22, 0.3]);
ylim([-0.001, 0.02]);
ax = gca;
ax.YAxis.Exponent = 0
ax.YAxis.TickLabelFormat = '%.3f'

x0=400;
y0=400;
height=400;
width=500;

set(gcf,'position',[x0,y0,width,height])

exportgraphics(gcf, "validation_figures/jet_area.png", "Resolution", 300);
savefig(gcf, "validation_figures/jet_area.fig");

%% Plots ratio loss
figure(3);
hold on;
total_diff = (stationary_area - moving_area)
jet_diff = (stationary_jet_area - moving_jet_area) 

total_ratio = total_diff ./ stationary_area;
jet_ratio = jet_diff ./ stationary_area;

scatter(ts_dns(1 : freq : end), total_ratio(1 : freq : end), [], ...
    'black');
scatter(ts_dns(1 : freq : end), jet_ratio(1 : freq : end), [], ...
    'black', 'd');

ylim([-0.005, 0.025]);
grid on;
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
ylabel("Area loss ratio", 'interpreter', 'latex');
legend(["Total", "Jet"], 'location', 'northwest', 'interpreter', 'latex');

x0=400;
y0=400;
height=300;
width=1000;

set(gcf,'position',[x0,y0,width,height])

exportgraphics(gcf, "validation_figures/area_ratios.png", "Resolution", 300);
savefig(gcf, "validation_figures/area_ratios.fig");
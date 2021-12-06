%% plot_stationary_pressures.m
% Plots the outer, inner and composite pressure in the stationary case for 
% examples
addpath("pressures");

xs = linspace(0, 0.4, 4096);
epsilon = 0.1;
t = 2;

% Calculate time dependent terms
d = 2 * sqrt(t);
d_t = 1 / sqrt(t);
A = d * d_t;
C = A;
J = pi * d / (8 * d_t^2);

% Calculate outer pressure
outer_ps = outer_pressure_stationary(xs, t, epsilon);

% Calculate overlap pressure
overlap_ps = overlap_pressure(xs, d, C, epsilon);

% Calculate inner pressure
inner_ps = inner_pressure(xs, d, d_t, J, epsilon);

% Calculate composite pressure
comp_ps = outer_ps + inner_ps - overlap_ps;

d_idx = sum(xs < epsilon * d);

%% Plot all
close all;
figure(1);
plot(xs(1 : d_idx), outer_ps(1 : d_idx), 'linewidth', 2, ...
    'color', 'black', 'linestyle', '--');
hold on;
plot(xs, inner_ps, 'linewidth', 2, ...
    'color', 'black', 'linestyle', ':');
plot(xs, comp_ps, 'linewidth', 2, ...
    'color', 'black');
xline(epsilon * d, 'linewidth', 1, ...
    'color', 0.55 * [1 1 1], 'linestyle', '--');

%% Figure properties
grid on;
fontsize = 22;
ylim([0, 40]);
ylabel("$p(x, 0, t)$", "interpreter", "latex", "Fontsize", fontsize);
xlabel("$x$", "interpreter", "latex", "Fontsize", fontsize);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
legend(["Outer solution", "Inner solution", "Composite solution"], ...
    'Interpreter', 'latex', 'fontsize', fontsize, ...
    'Location', 'northwest');
x0=400;
y0=400;
width=1200;
height=300;
set(gcf,'position',[x0,y0,width,height])
exportgraphics(gcf, "figures/composite_pressure_comparison.png", "Resolution", 300);
savefig(gcf, "figures/composite_pressure_comparison.fig");
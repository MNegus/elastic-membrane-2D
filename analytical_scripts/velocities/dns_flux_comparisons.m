clear;

%% Analytical functions
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);

%% Computational values
IMPACT_TIME = 0.125;
DELTA_T = 1e-4;
BOX_WIDTH = 6;
MAXLEVEL = 12;
MIN_CELL_SIZE = BOX_WIDTH / 2 ^ MAXLEVEL; 

%% Basilisk data 
data_dir = "~/scratch/multiple_slice_output/raw_data";

%% Load in fluxes text file
timestep = 3000; % Choose a timestep between about 1250 and 3000
data_mat = readmatrix(sprintf("%s/fluxes_%d.txt", data_dir, timestep));
xs = data_mat(:, 1);
fluxes = data_mat(:, 2);

%% Plots fluxes
close(figure(1));
figure(1);
scatter(xs, fluxes);
hold on;
yline(pi, 'color', 'black', 'linestyle', '--', 'linewidth', 2);

legend(["DNS result", "Analytical value"], "Interpreter", "Latex", ...
    "location", "east", "Fontsize", 12);
xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 12);
ylabel('Flux entering jet', 'Interpreter', 'latex', 'Fontsize', 12);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
grid on;
title(sprintf("Fluxes at computational time $t$ = %g", DELTA_T * timestep), "Interpreter", "Latex");

%% Plots x velocity lines
num_xs = length(xs);
close(figure(2));
figure(2);
hold on;
for k = 1 : 10 : num_xs
    x_cell = 2 * (k - 1);
    vel_data_mat = readmatrix(sprintf("%s/velocities_%d-x_cell_%d.txt", data_dir, timestep, x_cell));
    ys = vel_data_mat(:, 1);
    us = vel_data_mat(:, 3);
    vs = vel_data_mat(:, 4);
    
    x_val = xs(k) + x_cell * MIN_CELL_SIZE;
    plot(us, ys, 'Displayname', sprintf("$x$ = %.3f", x_val));
end
legend("Interpreter", "Latex", "location", "east", "Fontsize", 12);
xlabel('$u(x, y, t)$', 'Interpreter', 'latex', 'Fontsize', 12);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 12);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
grid on;
title(sprintf("$x$ velocities at $t$ = %g", DELTA_T * timestep - IMPACT_TIME), "Interpreter", "Latex");

%% Plots y velocity lines
num_xs = length(xs);
close(figure(3));
figure(3);
hold on;
for k = 1 : 10 : num_xs
    x_cell = 2 * (k - 1);
    vel_data_mat = readmatrix(sprintf("%s/velocities_%d-x_cell_%d.txt", data_dir, timestep, x_cell));
    ys = vel_data_mat(:, 1);
    us = vel_data_mat(:, 3);
    vs = vel_data_mat(:, 4);
    
    x_val = xs(k) + x_cell * MIN_CELL_SIZE;
    plot(vs, ys, 'Displayname', sprintf("$x$ = %.3f", x_val));
end
legend("Interpreter", "Latex", "location", "northeast", "Fontsize", 12);
xlabel('$v(x, y, t)$', 'Interpreter', 'latex', 'Fontsize', 12);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 12);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
grid on;
title(sprintf("$y$ velocities at $t$ = %g", DELTA_T * timestep - IMPACT_TIME), "Interpreter", "Latex");

%% Manually computes energy flux in fixed frame
num_xs = length(xs);
fluxes = zeros(size(num_xs));
close(figure(4));
figure(4);
hold on;
for k = 1 : num_xs
    x_cell = 2 * (k - 1);
    vel_data_mat = readmatrix(sprintf("%s/velocities_%d-x_cell_%d.txt", data_dir, timestep, x_cell));
    ys = vel_data_mat(:, 1);
    us = vel_data_mat(:, 3);
    vs = vel_data_mat(:, 4);
    fluxes(k) = trapz(ys, us .* (0.5 * (us.^2 + vs.^2)));
    
    x_val = xs(k) + x_cell * MIN_CELL_SIZE;
%     plot(vs, ys, 'Displayname', sprintf("$x$ = %.3f", x_val));
end
plot(xs, fluxes);
% legend("Interpreter", "Latex", "location", "northeast", "Fontsize", 12);
xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 12);
ylabel('Flux', 'Interpreter', 'latex', 'Fontsize', 12);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 12)
grid on;
title(sprintf("Flux into jet in stationary frame", DELTA_T * timestep - IMPACT_TIME), "Interpreter", "Latex");

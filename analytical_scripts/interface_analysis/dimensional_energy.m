close all;

%% Data loading
rho_l = 998;
V = 5;
R = 1e-3;

DELTA_T = 1e-4 * R / V;
IMPACT_TIME = 0.125 * R / V;

% Reads data
data = readmatrix("turnover_points_dimensional.txt");
ts = data(:, 1) - IMPACT_TIME;
ds = data(:, 2);
Js = data(:, 3);
d_ts = data(:, 4);
fluxes = data(:, 6);
energies = data(:, 7);

% Analytical solutions
analytical_ts = 0 : DELTA_T : max(ts);
analytical_ds = 2 * sqrt(analytical_ts);
analytical_d_ts = 1 ./ sqrt(analytical_ts);
analytical_Js = pi * analytical_ds ./ ( 8 * analytical_d_ts.^2) * (1 + 4 / pi);
analytical_fluxes = rho_l * V^3 * R * pi * ones(size(analytical_ts));
analytical_energy = pi * analytical_ts;

%% Plot energy flux
close(figure(6));
figure(6);
hold on;
plot(analytical_ts, analytical_fluxes, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
plot(ts, fluxes, 'linewidth', 2);
legend(["Analytical (pi)", "Basilisk"]);
xlabel("t");
ylabel("Flux");
title("Flux into jet");
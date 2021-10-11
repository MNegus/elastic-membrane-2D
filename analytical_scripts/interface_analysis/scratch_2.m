close all;

%% Data loading
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;

% Selects level
level = 11;
data = readmatrix("turnover_points.txt");
ts = DELTA_T * data(:, 1) - IMPACT_TIME;
ds = data(:, 2);
Js = data(:, 3);

% Analytical solutions
analytical_ts = 0 : DELTA_T : max(ts);
analytical_ds = 2 * sqrt(analytical_ts);
analytical_d_ts = 1 ./ sqrt(analytical_ts);
analytical_Js = pi * analytical_ds ./ ( 8 * analytical_d_ts.^2);
analytical_fluxes = pi * ones(size(analytical_ts));
analytical_energy = pi * analytical_ts;

%% Plot turnover points
close(figure(1));
figure(1);
hold on;
plot(analytical_ts, analytical_ds, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
plot(ts, ds, 'linewidth', 2);
legend(["Analytical", "Basilisk"]);
xlabel("t");
ylabel("d(t)");
title("Turnover point");

%% Plot jet thickness
close(figure(2));
figure(2);
hold on;
plot(analytical_ts, analytical_Js, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
plot(ts, Js, 'linewidth', 2);
ylim([0, 0.075]);
legend(["Analytical (J(t))", "Basilisk"]);
xlabel("t");
ylabel("J(t)");
title("Height of turnover point");


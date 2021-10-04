close all;

index = 3;

DELTA_T = 1e-4;
IMPACT_TIME = 0.125;

matlab_data = readmatrix("turnover_points.txt");
matlab_times = matlab_data(:, 1) * DELTA_T;
matlab_ds = matlab_data(:, index);

basilisk_data = readmatrix("turnover_points_basilisk.txt");
basilisk_times = basilisk_data(:, 1);
basilisk_ds = basilisk_data(:, index);

% Plot turnover points
figure(1);
hold on;
plot(matlab_times, matlab_ds);
plot(basilisk_times, basilisk_ds);
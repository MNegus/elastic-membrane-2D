addpath("../pressures");

clear;
%% Physical parameters
ALPHA = 1;
BETA = 0;
GAMMA = 12.8;
EPSILON = 1;
L = 1;
T_MAX = 0.1;
DELTA_T = 1e-4;
T_VALS = 0 : 10 * DELTA_T : T_MAX;

N_MEMBRANES = [512, 1024, 2048, 4094];

pressure_type = "outer"; % Which type of pressure solution to test

%% Data directory
data_dir = "/home/negus/Desktop/pressure_validation";

%% Imposed membrane solution
q = 0.01;
mag = 0.5;
w = @(t) mag * t^2;
w_t = @(t) 2 * mag * t;
w_tt = @(t) 2 * mag;
a = 1 - q^2 / 2;
d = @(t) 2 * sqrt(t - a * w(t));
d_t = @(t) (1 - a * w_t(t)) / sqrt(t - a * w(t));
d_tt = @(t) -(2 * a * (t - a * w(t)) * w_tt(t) + (1 - a * w_t(t))^2) / (2 * (t - a * w(t))^(3/2));
b = @(t) q^2 * (d(t) * w_t(t) - 2 * d_t(t) * w(t)) / (3 * d(t)^3);
c = @(t) q^2 * (d(t) * (d(t) * w_tt(t) - 4 * d_t(t) * w_t(t)) + w(t) * (6 * d_t(t)^2 - 2 * d(t) * d_tt(t))) / (3 * d(t)^4);
J = @(t) pi * (1 - w_t(t) + 3 * b(t) * d(t)^2 / 2)^2 * d(t) / (8 * d_t(t)^2);
full_w_fun = @(xs, t) w(t) * (1 - q^2 * xs.^2 / (EPSILON * d(t))^2);
full_w_t_fun = @(xs, t) w_t(t) - 3 * b(t) * xs.^2 / EPSILON^2;
full_w_tt_fun = @(xs, t) w_tt(t) - 3 * c(t) * xs.^2 / EPSILON^2;
full_w_x_fun = @(xs, t) w(t) * (- 2 * q^2 * xs / (EPSILON * d(t))^2);

%% Arrays for saving solutions of time-dependent quantities
ds_exact = zeros(size(T_VALS));
d_ts_exact = zeros(size(T_VALS));
Js_exact = zeros(size(T_VALS));
for k = 2 : length(T_VALS)
    t = T_VALS(k)
    ds_exact(k) = d(t);
    d_ts_exact(k) = d_t(t);
    Js_exact(k) = J(t);
end

ds_numerical = zeros(size(T_VALS));
d_ts_numerical = zeros(size(T_VALS));
Js_numerical = zeros(size(T_VALS));


%% Loops and plots both solutions
ps_diffs = zeros(size(N_MEMBRANES));
ds_diffs = zeros(size(N_MEMBRANES));
d_ts_diffs = zeros(size(N_MEMBRANES));
Js_diffs = zeros(size(N_MEMBRANES));


close(figure(1));
figure(1);

for m = 1 : length(N_MEMBRANES)
    N_MEMBRANE = N_MEMBRANES(m)
    
    DELTA_X = L / (N_MEMBRANE - 1); 
    xs = (0 : DELTA_X : L - DELTA_X)';
    
    ds_numerical = zeros(size(T_VALS));
    d_ts_numerical = zeros(size(T_VALS));
    Js_numerical = zeros(size(T_VALS));
    
    for k = 1  : length(T_VALS)   
        t = T_VALS(k)

        %% Updates functions 
        w_fun = @(xs) full_w_fun(xs, t);
        w_t_fun = @(xs) full_w_t_fun(xs, t);
        w_tt_fun = @(xs) full_w_tt_fun(xs, t);
        w_x_fun = @(xs) full_w_x_fun(xs, t);

        %% Finds numerical solution
        [ps_numerical, ds_numerical(k), d_ts_numerical(k), Js_numerical(k)] = w_dependents(xs, t, w_fun, ...
            w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);

        %% Finds exact solution
        ps_exact = imposed_pressure_quadratic(xs, t, d(t), d_t(t), J(t), w_t(t), w_tt(t), b(t), c(t), EPSILON, pressure_type);

        %% Determines difference in pressure
        idx = sum(xs < EPSILON * ds_numerical(k));
        diffs = abs(ps_exact(1 : idx - 1) - ps_numerical(1 : idx - 1));
        if (~isempty(diffs))
            ps_diffs(m) = max(ps_diffs(m), max(diffs))
        end
        
        %% Plots exact and numerical
        if (pressure_type == "outer")
            ps_stationary = outer_pressure_stationary(xs, t, EPSILON);
        else
            ps_stationary = composite_pressure_stationary(xs, t, EPSILON);
        end
        plot(xs, ps_stationary, 'linewidth', 2);
        hold on;
        plot(xs, ps_exact, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
        
        plot(xs, ps_numerical, 'linewidth', 2);
        hold off;
        legend(["Stationary", "Exact (moving)", "Numerical (moving)"]);
        xlim([0, 1]);
        ylim([0 15]);
        drawnow;
        pause(0.001);
        
        %% Plots difference
%         plot(xs, ps_exact - ps_numerical);
%         hold on;
%         xline(EPSILON * d(t));
%         hold off;
%         xlim([0, 1]);
%         drawnow;
%         pause(0.01);
    end
    
    %% Determines difference in time dependent quantitites
    ds_diffs(m) = max(abs(ds_exact - ds_numerical));
    d_ts_diffs(m) = max(abs(d_ts_exact(2 : end) - d_ts_numerical(2 : end)));
    Js_diffs(m) = max(abs(Js_exact - Js_numerical));
end

%% Plots errors
close(figure(1));
figure(1);
loglog(N_MEMBRANES, ps_diffs);
hold on
plot(N_MEMBRANES, 1e3 ./ N_MEMBRANES.^2);
title("Pressure convergence");

%% Plots ds
close(figure(2));
figure(2);
plot(T_VALS, 2 * sqrt(T_VALS));
hold on
plot(T_VALS, ds_exact);
plot(T_VALS, ds_numerical);

%%
ds_diffs
d_ts_diffs
Js_diffs

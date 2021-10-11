addpath("../pressures");

%% Data directory
data_dir = "/home/negus/Desktop/pressure_validation";

%% Physical parameters
ALPHA = 1;
BETA = 0;
GAMMA = 12.8;
EPSILON = 1;
L = 16;
T_MAX = 0.25;

% DELTA_TS = [1e-2, 1e-3, 1e-4, 1e-5];
DELTA_T = 1e-4;

N_MEMBRANES = [256, 512, 1024, 2048, 4096, 8192];
% N_MEMBRANES = 8192
T_VALS = 0 : 10 * DELTA_T : T_MAX;

pressure_type = "composite"; % Which type of pressure solution to test

%% Imposed membrane solution
a = 0.01; % Coefficient so w(t) = a * t^2
ws = a * T_VALS.^2;
w_ts = 2 * a * T_VALS;
w_tts = 2 * a * ones(size(T_VALS));
w_x_fun = @(xs) zeros(size(xs)); % Always == 0

%% Exact solutions for time dependent quantities
ds_exact = 2 * sqrt(T_VALS - ws);
d_ts_exact = (1 - w_ts) ./ sqrt(T_VALS - ws);
Js_exact = pi * (1 - w_ts).^2 .* ds_exact ./ (8 * d_ts_exact.^2);
Bs_exact = w_ts;

%% Arrays for numerical solutions
ds_numerical = zeros(size(T_VALS));
d_ts_numerical = zeros(size(T_VALS));
Js_numerical = zeros(size(T_VALS));

%% Loops
ps_diffs = zeros(size(N_MEMBRANES));
ds_diffs = zeros(size(N_MEMBRANES));
d_ts_diffs = zeros(size(N_MEMBRANES));
Js_diffs = zeros(size(N_MEMBRANES));
Bs_diffs = zeros(size(N_MEMBRANES));

close(figure(1));
figure(1);
for m = 1 : length(N_MEMBRANES)
    N_MEMBRANE = N_MEMBRANES(m)
    
    
    DELTA_X = L / (N_MEMBRANE - 1); 
    xs = (0 : DELTA_X : L - DELTA_X)';
    
    ds_numerical = zeros(size(T_VALS));
    d_ts_numerical = zeros(size(T_VALS));
    Js_numerical = zeros(size(T_VALS));
    Bs_numerical = zeros(size(T_VALS));
    
    for q = 2  : length(T_VALS)   
        t = T_VALS(q)

        %% Updates functions 
        w_fun = @(xs) ws(q) * ones(size(xs));
        w_t_fun = @(xs) w_ts(q) * ones(size(xs));
        w_tt_fun = @(xs) w_tts(q) * ones(size(xs));

        %% Finds numerical solution
        [ps_numerical, ds_numerical(q), d_ts_numerical(q), Js_numerical(q)] = w_dependents(xs, t, w_fun, ...
            w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);

        %% Finds exact solution
        if (pressure_type == "outer")
            ps_exact = outer_pressure_flat(xs, w_t_fun, w_tt_fun, ds_exact(q), d_ts_exact(q), EPSILON);
        elseif (pressure_type == "composite")
            ps_exact = composite_pressure_flat(xs, t, ds_exact(q), d_ts_exact(q), Js_exact(q), w_t_fun, w_tt_fun, EPSILON);
        end

        %% Determines difference
        max_diff = max(abs(ps_exact - ps_numerical));
        ps_diffs(m) = max(ps_diffs(m), max_diff);
        
        
        %% Plots
%         plot(xs, ps_exact);
%         hold on;
%         plot(xs, ps_numerical);
%         hold off;
%         xlim([0, 1]);
%         title(max_diff);
%         drawnow;
%         
%         pause(0.001);
    end
    
    ds_diffs(m) = max(abs(ds_exact - ds_numerical));
    d_ts_diffs(m) = max(abs(d_ts_exact(2 : end) - d_ts_numerical(2 : end)));
    Js_diffs(m) = max(abs(Js_exact - Js_numerical));
    Bs_diffs(m) = max(abs(Bs_exact - Bs_numerical));
    
end

%% Plots errors
close(figure(1));
figure(1);
loglog(N_MEMBRANES, ps_diffs);
hold on
plot(N_MEMBRANES, 1e3 ./ N_MEMBRANES.^2);
title("Pressure convergence");

%% Plots ds
ds_diffs
d_ts_diffs
Js_diffs



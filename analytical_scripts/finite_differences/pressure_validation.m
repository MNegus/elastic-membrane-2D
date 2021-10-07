addpath("../pressures");

%% Physical parameters
ALPHA = 1;
BETA = 0;
GAMMA = 12.8;
EPSILON = 1;
L = 16;
T_MAX = 0.25;
DELTA_T = 1e-4;
T_VALS = 0 : DELTA_T : T_MAX;

N_MEMBRANE = 4096;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

pressure_type = "composite"; % Which type of pressure solution to test

%% Data directory
data_dir = "/home/negus/Desktop/pressure_validation";

%% Imposed membrane solution
a = 0.01; % Coefficient so w(t) = a * t^2
full_w_fun = @(xs, t) a * t^2 * ones(size(xs));
full_w_t_fun = @(xs, t) 2 * a * t * ones(size(xs));
full_w_tt_fun = @(xs, t) 2 * a * ones(size(xs));
w_x_fun = @(xs) zeros(size(xs)); % Always == 0

%% Arrays for saving solutions of time-dependent quantities



%% Loops and plots both solutions
close(figure(1));
figure(1);



for q = 2  : 10 : length(T_VALS)   
    t = T_VALS(q)
   
    %% Updates functions 
    w_fun = @(xs) full_w_fun(xs, t);
    w_t_fun = @(xs) full_w_t_fun(xs, t);
    w_tt_fun = @(xs) full_w_tt_fun(xs, t);
    
    %% Finds numerical solution
    [ps_numerical, d, d_t, J] = w_dependents(xs, t, w_fun, ...
        w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);
    
    %% Finds exact solution
    if (pressure_type == "outer")
        ps_exact = outer_pressure_flat(xs, w_t_fun, w_tt_fun, d, d_t, EPSILON);
    elseif (pressure_type == "composite")
        ps_exact = composite_pressure_flat(xs, t, d, d_t, J, w_t_fun, w_tt_fun, EPSILON);
    end

    %% Plots exact and numerical
    plot(xs, ps_exact, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
    hold on;
    plot(xs, ps_numerical, 'linewidth', 2);
    hold off;
    legend(["Exact", "Numerical"]);
    max(abs(ps_exact - ps_numerical))
    

    xlim([0, 1]);
    ylim([0 5]);
    drawnow;
    pause(0.1);
   
    
end


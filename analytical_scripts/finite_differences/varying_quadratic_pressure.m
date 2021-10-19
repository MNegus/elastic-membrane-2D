addpath("../pressures");

clear;

%% Physical parameters
EPSILON = 1;
L = 1;
T_MAX = 0.1;
DELTA_T = 1e-4;
T_VALS = 0 : 10 * DELTA_T : T_MAX;
N_MEMBRANE = 1024;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

pressure_type = "composite"; % Which type of pressure solution to test
plotting = false;

%% Data directory
data_dir = "/home/negus/Desktop/pressure_validation";

%% Imposed membrane solution
q = 0.01;
mags = [0.01, 0.04, 0.16, 0.64, 1.28];
d_exact_vals = zeros(length(T_VALS), length(mags));
d_t_exact_vals = zeros(length(T_VALS), length(mags));
J_exact_vals = zeros(length(T_VALS), length(mags));

d_numerical_vals = zeros(length(T_VALS), length(mags));
d_t_numerical_vals = zeros(length(T_VALS), length(mags));
J_numerical_vals = zeros(length(T_VALS), length(mags));


%% Loops and plots all solutions
if plotting
    close(figure(1));
    figure(1);
end

for k = 1  : length(T_VALS)   
    t = T_VALS(k)

    %% Plots stationary
    if (pressure_type == "outer")
        ps_stationary = outer_pressure_stationary(xs, t, EPSILON);
    else
        ps_stationary = composite_pressure_stationary(xs, t, EPSILON);
    end
    
    if plotting
        plot(xs, ps_stationary, 'linestyle', '--', 'color', [0 0 0], 'linewidth', 2, 'Displayname', 'Stationary');
        hold on;
    end
    
    for mag_idx = 1 : length(mags)
        mag = mags(mag_idx);
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
        
        %% Updates functions 
        w_fun = @(xs) full_w_fun(xs, t);
        w_t_fun = @(xs) full_w_t_fun(xs, t);
        w_tt_fun = @(xs) full_w_tt_fun(xs, t);
        w_x_fun = @(xs) full_w_x_fun(xs, t);

        %% Finds numerical solution
        [ps_numerical, d_numerical_vals(k, mag_idx), d_t_numerical_vals(k, mag_idx), J_numerical_vals(k, mag_idx)] = w_dependents(xs, t, w_fun, ...
            w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);
        
        %% Finds exact solution
        ps_exact = imposed_pressure_quadratic(xs, t, d(t), d_t(t), J(t), w_t(t), w_tt(t), b(t), c(t), EPSILON, pressure_type);
        d_exact_vals(k, mag_idx) = d(t);
        d_t_exact_vals(k, mag_idx) = d_t(t);
        J_exact_vals(k, mag_idx) = J(t);
        
        %% Plots exact and numerical solution
        if plotting
            plot(xs, ps_exact, 'Displayname', sprintf("mag = %g (Exact)", mag), 'linestyle', '--', 'linewidth', 2);
            plot(xs, ps_numerical, 'Displayname', sprintf("mag = %g (Numerical)", mag));
        end
    end
    if plotting
        hold off;
        legend();
        xlim([0, 1]);
    %     ylim([0 30]);
        drawnow;
        pause(0.001);
    end


end

%% Plots d values
close(figure(2));
figure(2);
hold on;
plot(T_VALS, 2 * sqrt(T_VALS), 'Displayname', 'Stationary', 'linestyle', '--', 'color', 'black', 'linewidth', 2);

for mag_idx = 1 : length(mags)
    mag = mags(mag_idx);
    plot(T_VALS, d_exact_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Exact)", mag), 'linestyle', '--', 'linewidth', 2);
    plot(T_VALS, d_numerical_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Numerical)", mag));
end
legend();


%% Plots d_t values
close(figure(3));
figure(3);
hold on;
plot(T_VALS, 1 ./ sqrt(T_VALS), 'Displayname', 'Stationary', 'linestyle', '--', 'color', 'black', 'linewidth', 2);

for mag_idx = 1 : length(mags)
    mag = mags(mag_idx);
    plot(T_VALS, d_t_exact_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Exact)", mag), 'linestyle', '--', 'linewidth', 2);
    plot(T_VALS, d_t_numerical_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Numerical)", mag));
end
legend();

%% Plots J values
close(figure(4));
figure(4);
hold on;
% plot(T_VALS, 1 ./ sqrt(T_VALS), 'Displayname', 'Stationary', 'linestyle', '--', 'color', 'black', 'linewidth', 2);

for mag_idx = 1 : length(mags)
    mag = mags(mag_idx);
    plot(T_VALS, J_exact_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Exact)", mag), 'linestyle', '--', 'linewidth', 2);
    plot(T_VALS, J_numerical_vals(:, mag_idx), 'Displayname', sprintf("mag = %g (Numerical)", mag));
end
legend();

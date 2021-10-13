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
q = 0.5;
a = 0.01;
w = @(t) a * t^2;
w_t = @(t) 2 * a * t;
w_tt = @(t) 2 * a;
d = @(t) 2 * sqrt(t - w(t) * (1 - q^2 / 2));
d_t = @(t) (1 - w_t(t) * (1 - q^2 / 2)) / sqrt(t - w(t) * (1 - q^2 / 2));
d_tt = @(t)
b = @(t) q^2 * (d(t) * w_t(t) - 2 * d_t(t) * w(t)) / (3 * d(t)^3);
c = @(t) q^2 * (d(t) * (d(t) * w_tt(t) - 4 * d_t(t) * w_t(t)) + w(t) * (6 * w_t(t)^2 - 2 * d(t) * d_tt(t))) / (3 * d(t)^4);
full_w_fun = @(xs, t) w(t) * (1 - q^2 * xs.^2 / (EPSILON * d(t))^2);
full_w_t_fun = @(xs, t) w_t(t) - 3 * b(t) * xs.^2 / EPSILON^2;
full_w_tt_fun = @(xs, t) w_tt(t) - 3 * c(t) * xs.^2 / EPSILON^2;

% Plot solution for ws
figure(1);
for k = 2 : length(T_VALS)
   t = T_VALS(k); 
   w_vals = full_w_tt_fun(xs, t);
   plot(xs, w_vals);
   ylim([0, max(w_vals)]);
   xlim([0, 1]);
   drawnow;
   pause(0.01);
    
end

%% Arrays for saving solutions of time-dependent quantities



%% Loops and plots both solutions
% close(figure(1));
% figure(1);
% 
% 
% 
% for q = 2  : 10 : length(T_VALS)   
%     t = T_VALS(q)
%    
%     %% Updates functions 
%     w_fun = @(xs) full_w_fun(xs, t);
%     w_t_fun = @(xs) full_w_t_fun(xs, t);
%     w_tt_fun = @(xs) full_w_tt_fun(xs, t);
%     
%     %% Finds numerical solution
%     [ps_numerical, d_val, d_t_val, J_val] = w_dependents(xs, t, w_fun, ...
%         w_t_fun, w_tt_fun, w_x_fun, pressure_type, EPSILON);
%     
%     %% Finds exact solution
%     if (pressure_type == "outer")
%         ps_exact = outer_pressure_flat(xs, w_t_fun, w_tt_fun, d_val, d_t_val, EPSILON);
%     elseif (pressure_type == "composite")
%         ps_exact = composite_pressure_flat(xs, t, d_val, d_t_val, J_val, w_t_fun, w_tt_fun, EPSILON);
%     end
% 
%     %% Plots exact and numerical
%     plot(xs, ps_exact, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
%     hold on;
%     plot(xs, ps_numerical, 'linewidth', 2);
%     hold off;
%     legend(["Exact", "Numerical"]);
%     max(abs(ps_exact - ps_numerical))
%     
% 
%     xlim([0, 1]);
%     ylim([0 5]);
%     drawnow;
%     pause(0.1);
%    
%     
% end


%% Parameters (physical and ones relating to other methods)
EPSILON = 1;
ALPHA = 2 / EPSILON^2; BETA = 1 * EPSILON^2; GAMMA = 2 * EPSILON^2; 
L = 4;
N_MEMBRANE = 2056;
T_MAX = 0.25;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L)';

%% Normal modes paramters
N = 64;
D_MAX = 2 * sqrt(T_MAX);
DELTA_D = 1e-4;
D_VALS = 0 : DELTA_D : D_MAX;

%% Solves ode for avals etc
[t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form] = a_ode_solution(ALPHA, BETA, GAMMA, EPSILON, DELTA_D, D_MAX, N, L);

%% Interpolates solution onto regular t_vals
T_MAX = max(t_vals_d_form);
DELTA_T = 1e-4;
T_VALS = 0 : DELTA_T : T_MAX;
ds = interp1(t_vals_d_form, d_vals_d_form, T_VALS);
as = interp1(t_vals_d_form, as_d_form, T_VALS);
a_ts = interp1(t_vals_d_form, a_ts_d_form, T_VALS);


%% Plots solution
for k = 1 : length(T_VALS)
    t = T_VALS(k)
    [ws, w_ts] = w_solution_normal_modes(xs, as(k, :), a_ts(k, :), L, N);
    
    subplot(2, 1, 1);
    plot(xs, ws);
    hold on;
    xline(EPSILON * ds(k));
    hold off;
    
    subplot(2, 1, 2);
    plot(xs, w_ts);
    hold on;
    xline(EPSILON * ds(k));
    hold off;
    
    pause(0.01);
    
end

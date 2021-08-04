
%% Parameters (physical and ones relating to other methods)
EPSILON = 1;
ALPHA = 10 / EPSILON^2; BETA = 20 * EPSILON^2; GAMMA = 0.1 * EPSILON^2; 
L = 4;
N_MEMBRANE = 2056;
T_MAX = 0.25;
DELTA_T = 1e-4;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L)';

%% Normal modes paramters
N = 16;
D_MAX = 2 * sqrt(T_MAX);
DELTA_D = 1e-4;
D_VALS = 0 : DELTA_D : D_MAX;

%% Solves ode for avals etc
[t_vals, d_vals, as, a_ts] = a_ode_solution(ALPHA, BETA, GAMMA, EPSILON, DELTA_D, D_MAX, N, L);

%% Interpolates solution onto regular t_vals

%% Plots solution
for k = 1 : length(t_vals)
    t_vals(k)
    [ws, w_ts] = w_solution_normal_modes(xs, as(k, :), a_ts(k, :), L, N);
    
    subplot(2, 1, 1);
    plot(xs, ws);
    
    subplot(2, 1, 2);
    plot(xs, w_ts);
    
    pause(0.01);
    
   
    
end

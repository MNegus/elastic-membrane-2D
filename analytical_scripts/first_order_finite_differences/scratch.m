clear;

%% Parameters
ALPHA = 2;
BETA = 1000;
GAMMA = 2;
EPSILON = 1;
L = 4;

T_MAX = 0.25;
DELTA_T = 1e-4;
N_MEMBRANE = 20000;

%% Derived parameters
DELTA_X = L / (N_MEMBRANE - 1); 
M = N_MEMBRANE - 1;
xs = (0 : DELTA_X : L - DELTA_X)';
T_VALS = 0 : DELTA_T : T_MAX;

%% Derive matrices
[L_mat, A_mat, A0_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

%% Initialise w_previous
N_max = 4;
w_previous = zeros(size(xs));
for n = 1 : N_max
    lambda = pi * (2 * n - 1) / (2 * L);
    w_previous = w_previous + (1 / n) * cos(lambda * xs);
end

figure(1);
plot(xs, w_previous);
pause(0.001);

%% Initialise w
% w = A0_mat \ w_previous;
w = w_previous;

figure(1);
plot(xs, w);
pause(0.001);

%% Loops over time
for k = 2 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    %%
    w_next = homogeneous_membrane_timestep(w, w_previous, A_mat);
    
    %% Plots
    figure(1);
    plot(xs, w_next);
    pause(1e-9);
    drawnow;
    
    %% Swaps
    temp = w_previous;
    w_previous = w;
    w = w_next;
    w_next = temp;
    

end
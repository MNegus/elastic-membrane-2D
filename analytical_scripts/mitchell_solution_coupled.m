%% mitchell_solution_stationary.m
% Code to use an iterative Mitchell FD method

clear;
close all;

addpath("pressures");

%% Parameters
EPSILON = 1;
ALPHA = 2.0 / EPSILON^4; BETA = 1; GAMMA = 2; 
L = 4;
N_MEMBRANE = 1025;
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.25;
DELTA_T = 10^-4;
DELTA_X = L / (N_MEMBRANE - 1); 

% Derived parameters
Cpressure = DELTA_X^4 / GAMMA;
Calpha = ALPHA * DELTA_X^4 / (GAMMA * DELTA_T^2);
Cbeta = BETA * DELTA_X^2 / GAMMA;

%% Basilisk compare data
data_directory = "~/Desktop/decoupled_example";
IMPACT_TIMESTEP = 1250;

%% Initialise arrays
xs = (0 : DELTA_X : L - DELTA_X)';
w_previous = zeros(size(xs));
w = zeros(size(xs));
w_t = zeros(size(xs));
w_tt = zeros(size(xs));
w_next = zeros(size(xs));
p_previous = zeros(size(xs));
p = zeros(size(xs));
p_next = zeros(size(xs));

%% Stationary solutions to time-dependent terms
d_fun = @(t) 2 * sqrt(t);
d_t_fun = @(t) 1 / sqrt(t);
A_fun = @(t) d_fun(t) * d_t_fun(t);
C_fun = @(t) d_fun(t) * d_t_fun(t);
J_fun = @(t) pi * d_fun(t) / (8 * d_t_fun(t)^2);
    
%% Matrix definitions
% A definition
A_upper_upper = ones(M, 1);
A_upper_upper(3) = 2;

A_upper = (-Cbeta - 4) * ones(M, 1);
A_upper(2) = -2 * Cbeta - 8;

A_main = (4 * Calpha + 2 * Cbeta + 6) * ones(M, 1);
A_main(2) = A_main(2) + 1;
A_main(M) = A_main(M) - 1;

A_lower = (-Cbeta - 4) * ones(M, 1);

A_lower_lower = ones(M, 1);

A_mat = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

% B definition
B_upper_upper = -2 * ones(M, 1);
B_upper_upper(3) = 2 * B_upper_upper(3);

B_upper = (2 * Cbeta + 8) * ones(M, 1);
B_upper(2) = 2 * B_upper(2);

B_main = (8 * Calpha - 4 * Cbeta - 12) * ones(M, 1);
B_main(2) = B_main(2) - 2;
B_main(M) = B_main(M) + 2;

B_lower = (2 * Cbeta + 8) * ones(M, 1);

B_lower_lower = -2 * ones(M, 1);

B_mat = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);

%% Initial conditions
t = 0
k = 0;

w_previous = zeros(size(xs));

subplot(2, 1, 1);
plot(xs, w_previous);

subplot(2, 1, 2);
plot(xs, p_previous);
pause(0.01);


%% Solve for first w
t = t + DELTA_T
k = 1;

% Initial guess for w_next and its derivatives
w_guess = w_previous;
w_t_guess = zeros(size(w_guess));
w_tt_guess = zeros(size(w_guess));
w_x_guess = zeros(size(w_guess));
w_x_guess(2 : M - 1) = (w_guess(3 : M) - w_guess(1 : M - 2)) / (2 * DELTA_X);
w_x_guess(M) = -w_guess(M - 1) / (2 * DELTA_X);

% Stopping conditions
diff = 1e9;
tol = 1e-9;
iter_num = 1;

while (diff > tol) || (iter_num < 10)
    iter_num
    %% Create interpolant functions
    w_fun = @(x) interp1(xs, w_guess, x);
    w_t_fun = @(x) interp1(xs, w_t_guess, x);
    w_tt_fun = @(x) interp1(xs, w_tt_guess, x);
    w_x_fun = @(x) interp1(xs, w_x_guess, x);
    
    %% Update guess for d and d_t
    figure(7);
    plot(xs, w_fun(xs));
    pause(0.1);
    [d, d_t] = turnover_point(t, w_fun, w_t_fun, w_x_fun, EPSILON)
%     d = 2 * sqrt(t + DELTA_T);
%     d_t = 1 / sqrt(t + DELTA_T);
    
    
    %% Update guess for outer pressure outer pressure
    A = d * d_t;
     
    % Determine outer pressure with stationary assumption
    p_guess = outer_pressure_stationary(xs, d, A, EPSILON);
    

   %% Update guess for w
   rhs = 0.5 * B_mat * w_previous + Cpressure * p_guess;
   new_w_guess = A_mat \ rhs;
   
   %% Update guess for w_tt
   new_w_tt_guess = pde_rhs(new_w_guess, p_guess, ALPHA, BETA, GAMMA, M, DELTA_X);
   
   %% Update guess for w_t
   new_w_t_guess = w_t + (DELTA_T / 2) * (w_tt + new_w_tt_guess);
   
   %% Measures change in w
   diff = max(abs(new_w_guess - w_guess));
   
   %% Updates guesses
   w_guess = new_w_guess;
   w_t_guess = new_w_t_guess;
   w_tt_guess = new_w_tt_guess;
    
   iter_num = iter_num + 1;
end

disp("Done!")
% 
% subplot(2, 1, 1);
% plot(xs, w);
% 
% subplot(2, 1, 2);
% plot(xs, p);
% pause(0.01);
% 
% % Swaps pressures
% temp = p_previous;
% p_previous = p;
% p = p_next;
% p_next = temp;


% % Loops
% while (t < T_MAX) 
%     
%     % Determines pressure at next timestep
%     t_next = t + DELTA_T;
%     p_next = composite_pressure_stationary(xs, t_next, d_fun(t_next), ...
%         d_t_fun(t_next), A_fun(t_next), C_fun(t_next), J_fun(t_next), ...
%         EPSILON);w_tt
% 
% %     p_next = outer_pressure_stationary(xs, d_fun(t_next), A_fun(t_next), EPSILON);
%     
%     % Restrict large values of p_next
%     p_next(p_next > 1e4) = 0;
%     
%     rhs = B_mat * w - A_mat * w_previous ...
%         + Cpressure * (p_previous + 2 * p + p_next);
%     w_next = A_mat \ rhs;
% 
%     % Swaps ws
%     temp = w_previous;
%     w_previous = w;
%     w = w_next;
%     w_next = temp;
%     
%     % Swaps ps
%     temp = p_previous;
%     p_previous = p;
%     p = p_next;
%     p_next = temp;
% 
%     % Updates time
%     t = t_next
%     k = k + 1;
% 
%     %% Plots
%     % w plot
%     subplot(2, 1, 1);
%     plot(xs, w_next);
%     % Reads in Basilisk solution
%     membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", data_directory, k + IMPACT_TIMESTEP));
%     unsorted_xs = membrane_mat(:, 1);
%     unsorted_ws = membrane_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     ws = unsorted_ws(idxs);
%     hold on;
%     plot(sorted_xs, ws);
%     hold off;
% 
%     subplot(2, 1, 2);
%     plot(xs, p_next);
%     % Reads in Basilisk solution
%     pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", data_directory, k + IMPACT_TIMESTEP));
%     unsorted_xs = pressure_mat(:, 1);
%     unsorted_ps = pressure_mat(:, 2);
%     [sorted_xs, idxs] = sort(unsorted_xs);
%     ps = unsorted_ps(idxs);
%     hold on;
%     plot(sorted_xs, ps);
%     hold off;
%     
%     pause(0.001);
% 
% 
% 
% end

%% function definitions
function L = pde_rhs(ws, ps, ALPHA, BETA, GAMMA, M, DELTA_X)
    
    L = zeros(M, 1);
    
    Cbeta = BETA * DELTA_X^2 / GAMMA;
    
    %% Boundary near x = 0
    L(1) = -(6 + 2 * Cbeta) * ws(1) + 2 * (4 + Cbeta) * ws(2) - 2 * ws(3);
    L(2) = (4 + Cbeta) * ws(1) - (7 + 2 * Cbeta) * ws(2) + (4 + Cbeta) * ws(3) - ws(4);
    
    %% Bulk nodes
    L(3 : M - 2) = ...
           - ws(1 : M - 4) ...
           + (4 + Cbeta) * ws(2 : M - 3) ...
           - (6 + 2 * Cbeta) * ws(3 : M - 2) ...
           + (4 + Cbeta) * ws(4 : M - 1) ...
           - ws(5 : M);
       
   %% Boundary near x = L
   L(M - 1) = -ws(M - 3) + (4 + Cbeta) * ws(M - 2) - (6 + 2 * Cbeta) * ws(M - 1) + (4 + Cbeta) * ws(M);
   L(M) = -ws(M - 2) + (4 + Cbeta) * ws(M - 1) - (5 + 2 * Cbeta) * ws(M);
   
   %% Add ps to L and multiply by constant factor
   L = (GAMMA / (ALPHA * DELTA_X^4)) * (L + ps);

end
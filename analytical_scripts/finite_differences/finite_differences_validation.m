%% Physical parameters
ALPHA = 1;
BETA = 0;
GAMMA = 12.8;
EPSILON = 1;
L = 16;
Ne = 128;
T_MAX = 0.01;
DELTA_T = 1e-4;
tvals = 0 : delta_t : tmax;

N_MEMBRANE = 1024;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

lambda = @(n) pi * (2 * n - 1) / (2 * L);
lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
ks = beta * lambdas.^2 + gamma * lambdas.^4;

%% Data directory
data_dir = "/home/negus/Desktop/finite_differences_validation_test";

%% Saves solution
save_prescribed_finite_differences_solution(data_dir, ...
    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, Ne)

%% Plots and compares solution to exact
close(figure(1));
figure(1);

% for q = 1 : length(tvals)
for q = 2  : length(tvals)   
    t = tvals(q)
   
   
    %% Finds exact solution
    exact_as = alpha * (1 - cos(sqrt(ks) * t / sqrt(alpha))) ./ (ks .* sqrt(L * lambdas));
    exact_ws = w_solution(xs, exact_as, L, Ne);

    %% Reads in  numerical solution
    ws_mat = matfile(sprintf("%s/w_%d.mat", data_dir, q));
    numerical_ws = EPSILON^2 * ws_mat.w_next;

   %% Plots exact and numerical
   plot(xs, exact_ws, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
   hold on;
   plot(xs, numerical_ws, 'linewidth', 2);
   hold off;
   legend(["Exact", "Numerical"]);
   drawnow;
   max(abs(exact_solution - numerical_solution))
   pause(0.1);
   
    
end


function ws = w_solution(xs, as, L, N)
    
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
    
    % Find ws
    ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);
end

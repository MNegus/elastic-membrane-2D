
%% Physical parameters
alpha = 1;
beta = 0;
gamma = 12.8;
epsilon = 1;
L = 16;
N = 128;
tmax = 0.25;
delta_t = 1e-4;
tvals = 0 : delta_t : tmax;
xs = linspace(0, L, 1024);
lambda = @(n) pi * (2 * n - 1) / (2 * L);
lambdas = pi * (2 * (1 : N)' - 1) / (2 * L);
ks = beta * lambdas.^2 + gamma * lambdas.^4;

%% Data directory
data_dir = "/home/michael/Desktop/normal_modes_validation_test";



%% Saves solution
save_prescribed_normal_modes_solution(data_dir, alpha, beta, gamma, epsilon, N, L, tmax, delta_t);

%% Opens back up solution

as_mat = matfile(sprintf("%s/as.mat", data_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/a_ts.mat", data_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/q_ts.mat", data_dir));
q_ts = q_ts_mat.q_ts;

%%
size(as)

%% Plots and compares solution to exact
close(figure(1));
figure(1);

exact_as = zeros(1, N);
% for q = 1 : length(tvals)
for q = 1 : 10 : length(tvals)   
    t = tvals(q)
   
   
   %% Finds exact solution
   exact_as = alpha * (1 - cos(sqrt(ks) * t / sqrt(alpha))) ./ (ks .* sqrt(L * lambdas));
   

   %% Plots exact and numerical
   exact_solution = w_solution(xs, exact_as, L, N);
   numerical_solution = w_solution(xs, as(q, :), L, N);
   plot(xs, exact_solution, 'linewidth', 5, 'color', 0.5 * [1 1 1]);
   hold on;
   plot(xs, numerical_solution, 'linewidth', 2);
   hold off;
   legend(["Exact", "Numerical"]);
   drawnow;
   max(abs(exact_solution - numerical_solution))
   pause(0.1);
   
    
end



%% Function definitions
function ws = w_solution(xs, as, L, N)
    ws = zeros(size(xs));
    
    lambda = @(n) pi * (2 * n - 1) / (2 * L);
    
    
    %% FIND AN ALTERNATIVE USING MATRIX MULTIPLICATION
    for n = 1 : N
        ws = ws + as(n) * cos(lambda(n) * xs) / sqrt(L);
    end

end

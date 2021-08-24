%% Parameters
epsilon = 1;
alpha = 2; beta = 1; gamma = 2;
L = 4;
t_max = 0.25;
delta_t = 1e-4;
ts = 0 : delta_t : t_max;

% FD parameters
N_membrane = 2000;

%% Construct w0
N_max = 4;
delta_x = L / (N_membrane - 1); 
xs = (0 : delta_x : L - delta_x)';
w0 = zeros(size(xs));

for n = 1 : N_max
    lambda = pi * (2 * n - 1) / (2 * L);
    w0 = w0 + (1 / n) * cos(lambda * xs);
end

%%
[ws, w_ts, w_tts] = homogeneous_mol_solution(alpha, beta, gamma, epsilon, L, delta_t, t_max, ts, N_membrane, w0);

%% loops over time
figure(1);
for k = 1 : 10 :length(ts)
    k
   plot(xs, ws(:, k)); 
   drawnow;
   pause(0.01);
    
end

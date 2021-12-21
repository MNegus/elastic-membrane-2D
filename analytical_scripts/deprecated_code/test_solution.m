close all

%% Loads in parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

ALPHA = ALPHAS(1);
BETA = BETAS(1);
GAMMA = GAMMAS(1);

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Loads in normal modes solutions
parent_dir = "/media/michael/newarre/elastic_membrane/scratch";
analytical_parent_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);

N_mat = matfile(sprintf("%s/normal_modes/N.mat", analytical_parent_dir));
N = N_mat.N

as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
q_ts = q_ts_mat.q_ts;

ks = pi * (2 * (1 : N) - 1) / (2 * L);


%% Load in solution and plot them
as_analytical = zeros(1, N);
figure(1);
for k = IMPACT_TIMESTEP + 1 : 10 : length(T_VALS)
    t = T_VALS(k);
    % Plot numerical result
%     plot(1 : N, as(k - IMPACT_TIMESTEP, :));
    ws = sum(as(k - IMPACT_TIMESTEP, :) .* cos(xs * ks), 2) / sqrt(L);
    plot(xs, ws);
    hold on;
    
    % Plot analytical result
    for n = 1 : N
        as_analytical(n) = a_sol(n, t, ALPHAS, GAMMAS, L); 
    end
%     plot(1 : N, as_analytical);
    ws_analytical = sum(as_analytical .* cos(xs * ks), 2) / sqrt(L);
    plot(xs, ws_analytical);
    hold off;
    
%     xlim([5, N]);
    drawnow;
    pause(1);
end



function a = a_sol(n, t, alpha, gamma, L)
    k = pi * (2 * n - 1) / (2 * L);
    omega = sqrt(gamma / alpha) * k^2;
    mu = 2 * (alpha / gamma)^(1/4);
    
    a = (2 * pi / omega^2) * sin(omega * t + mu^2 / 4 - pi / 4) * alpha;
end
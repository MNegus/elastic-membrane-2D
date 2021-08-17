alpha = 0.1;
beta = 0;
gamma = 1e-4;
epsilon = 1;
L = 4;
q = 10;
delta_t = 1e-4;
tmax = 0.25;
% N_max = N_stable(alpha, beta, gamma, L, q, delta_d)

N_MEMBRANE = 1024;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

%%
[N, delta_d, as, a_ts, a_tts, q_ts] ...
    = validated_normal_modes_solution(alpha, beta, gamma, epsilon, N, L, tmax, delta_t, xs)

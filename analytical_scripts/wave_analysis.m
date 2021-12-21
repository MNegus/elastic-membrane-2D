[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();

%
ALPHAS
BETAS
GAMMAS
L

%% Determine omega
k = @(n) pi * (2 * n - 1) / (2 * L);
omega = @(n) sqrt(EPSILON^4 * (BETAS * k(n)^2 + GAMMAS * k(n)^4) / ALPHAS);

%% Determine phase velocity
cp = @(n) omega(n) / k(n);

%% Determine group velocity
% cg = @(n) (EPSILON^4 / ALPHAS) * (2 * BETAS * k(n) + 4 * GAMMAS * k(n)^3) / (2 * omega(n));
cg = @(n) EPSILON^2 * (BETAS + 2 * GAMMAS * k(n)^2) / sqrt(ALPHAS * (BETAS + GAMMAS * k(n)^2));

%% Determine return time
Tp = @(n) 2 * L / cp(n);
Tg = @(n) 2 * L / cg(n);

%% Output specific values
N = 512;
Tp(N)
Tg(N)
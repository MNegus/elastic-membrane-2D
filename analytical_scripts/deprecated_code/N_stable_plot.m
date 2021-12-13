%% N_stable plots
clear
delta_t = 1e-4;
delta_d = delta_t;

lambda = @(n) pi * (2 * n - 1) / (2 * L);
tfun = @(n) 2 * pi * sqrt(alpha / (beta * lambda(n)^2 + gamma * lambda(n)^4));

L = 4;
q = 10;


%%
alphas = 10.^(-6 : 1 : 6);
beta = 1;
gamma = 0;
% betas = 10.^(-6 : 1 : 2);

Ns = zeros(size(alphas));
for k = 1 : length(Ns)
   Ns(k) = N_stable(alphas(k), beta, gamma, L, q, delta_d);
end
Ns
close all
figure(1);
loglog(alphas, Ns, '-o');
xlabel("$\alpha$ / $\beta$", "interpreter", "latex", "fontsize", 15);
ylabel("$N_{max}$", "interpreter", "latex", "fontsize", 15);
title("$ \Delta t = $" + delta_t + ", $q = $" + q + ", $\gamma = 0$", "interpreter", "latex", "fontsize", 15);

exportgraphics(gca, "wave_equation_stability.png", "resolution", 300);

%%
alphas = 10.^(-10 : 1 : 8);
beta = 0;
gamma = 1;

Ns = zeros(size(alphas));
for k = 1 : length(Ns)
   Ns(k) = N_stable(alphas(k), beta, gamma, L, q, delta_d);
end
Ns
close all
figure(1);
loglog(alphas, Ns, '-o');
xlabel("$\alpha$ / $\gamma$", "interpreter", "latex", "fontsize", 15);
ylabel("$N_{max}$", "interpreter", "latex", "fontsize", 15);
title("$ \Delta t = $" + delta_t + ", $q = $" + q + ", $\beta = 0$", "interpreter", "latex", "fontsize", 15);

exportgraphics(gca, "membrane_equation_stability.png", "resolution", 300);


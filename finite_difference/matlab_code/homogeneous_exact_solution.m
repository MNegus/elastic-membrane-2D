function ws = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L)
%% homogeneous_analytical_solution
% Solves the membrane equation 
%     ALPHA * w_tt - BETA * w_xx + GAMMA * w_xxxx = 0,
% with boundary conditions
%     w_x = w_xxx at x = 0, w = w_xx = 0 at x = L,
% with an analytical solution given by
% w(x, t) = sum_{n=1}^\infinity A_n cos(l_n t) cos(lambda_n x),
% with 
% 
% lambda_n = (2n -1) pi / (2L) and 
% l_n^2 = (BETA lambda_n^2 + GAMMA lambda_n^4) / ALPHA,
% 
% and A_n given by integrating the initial condition
%
% A_n = (2/L) inon(t_0^L w(x, 0) cos(lambda_n x) dx.
%
% We consider simple initial conditions of the form 
% w(x, 0) = sum_{n=1}^{N0} A_n cos(lambda_n x),
% for a finite value of N0, so that we just know A_n exactly at the start.

%% Anonymous functions
lambda = @(n) pi * (2 * n - 1) / (2 * L);
l = @(n) sqrt(BETA * lambda(n).^2 + GAMMA * lambda(n).^4) / sqrt(ALPHA);


%% Initial condition
N0 = 3; % Take 3 terms
As = [1, 0.5, 0.25]; % Size of A_n

%% Output w
ws = zeros(size(xs));
for n = 1 : N0
   ws = ws + As(n) * cos(l(n) * t) * cos(lambda(n) * xs); 
end
ws = ws';
end
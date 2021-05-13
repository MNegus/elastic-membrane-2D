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
% A_n = (2/L) int_0^L w(x, 0) cos(lambda_n x) dx.
%
% We consider simple initial conditions of the form 
% w(x, 0) = sum_{n=1}^{N0} A_n cos(lambda_n x),
% for a finite value of N0, so that we just know A_n exactly at the start.

%% Parameters
% Physical parameters 
ALPHA = 0.1; % Mass term
BETA = 10; % Tension term
GAMMA = 1; % Bending stiffness term
L = 4; % Width of domain in x

% Computational parameters 
N_MEMBRANE = 128; % Number of grid points on the membrane
T_MAX = 1.0; % Maximum value of time
DELTA_T = 1e-4; % Timestep size

% Derived parameters
Deltax = L / (N_MEMBRANE - 1); % Grid size step
xs = linspace(0, L, N_MEMBRANE); % Array for x values

%% Anonymous functions
lambda = @(n) pi * (2 * n - 1) / (2 * L);
l = @(n) sqrt(BETA * lambda(n).^2 + GAMMA * lambda(n).^4) / sqrt(ALPHA);


%% Initial condition
N0 = 3; % Take 3 terms
As = [1, 0.5, 0.25]; % Size of A_n
w0 = zeros(size(xs)); % Initialise w0
for n = 1 : N0
   w0 = w0 + As(n) * cos(lambda(n) * xs); 
end

%% Timestepping
t = 0
w = w0;
figure(1);
plot(xs, w);

% Increments t
t = DELTA_T

while (t < T_MAX) 
    % Update w
    w = zeros(size(xs));
    for n = 1 : N0
       w = w + As(n) * cos(l(n) * t) * cos(lambda(n) * xs); 
    end
    
    % Plot w
    plot(xs, w);
    ylim([-2, 2]);
    pause(0.0001);
   
    % Increments t
    t = t + DELTA_T
end


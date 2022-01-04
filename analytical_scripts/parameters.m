function [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters()

    %% Define dimensional values (in SI units)
    EPSILON = 1;
%     rho_m = 1200; % Membrane density
%     nu = 10^-3; % Membrane thickness
%     E = 2.2 * 10^9; % Membrane Young's modulus
%     T = 17; % Membrane tension
%     rho_l = 789; % Droplet density
%     R = 1e-3; % Droplet radius
%     V = 1; % Droplet speed
    
    EPSILON = 1;
    ALPHAS = 0.5;
    BETAS = 0;
    GAMMAS = 16.7;

    %% Saran wrap
%     ALPHAS = 1.7 * 10^-2;
%     BETAS = 0;
%     GAMMAS = 1 * 0.00337;
    
    %% LDPE
%     ALPHAS = 0.9;
%     BETAS = 2.2;
%     GAMMAS = 4660;

    %% Define parameters
%     EPSILON = 1;
%     ALPHAS = rho_m * nu / (rho_l * R)
%     BETAS =  T / (rho_l * R * V)
%     GAMMAS = E * nu^3 / (3 * rho_l * R^3 * V^2)

    % alpha_varying
%     ALPHAS = [1, 2, 4, 8, 16];
%     BETAS = ones(size(ALPHAS));
%     GAMMAS = 1 * ones(size(ALPHAS));
% 
%     % beta_varying
%     BETAS = [0, 10, 20, 40, 80, 160, 320, 640, 1280];
%     ALPHAS = ones(size(BETAS));
%     GAMMAS = 2 * (ALPHAS).^3;
%     
%     % gamma_varying
%     GAMMAS = [2, 4, 8, 16, 32, 64, 128];
%     ALPHAS = 2 * ones(size(GAMMAS));
%     BETAS = zeros(size(GAMMAS));

%     ALPHAS = [1, 1.5, 2, 3, 4, 6, 8];
%     BETAS = zeros(size(ALPHAS));
%     GAMMAS = 2 * (ALPHAS).^3;

    % beta_varying
%     BETAS = [0, 10, 40, 160, 640, 2560, 10240];
%     ALPHAS = ones(size(BETAS));
%     GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3;
    
    % gamma_varying
%     GAMMAS = [2, 8, 32, 128, 512, 2048, 8192];
%     GAMMAS = 10.^[-4, -3, -2, -1, 0, 1, 2] * EPSILON^2;
%     ALPHAS = 1 * ones(size(GAMMAS));
%     BETAS = zeros(size(GAMMAS));

    L = 16;
    T_MAX = 0.4;
%     T_MAX = 0.125 + 0.01;
    DELTA_T = 1e-4;

    % FD parameters
    N_MEMBRANE = 10924;
    
    IMPACT_TIME = 0.125 / EPSILON^2;
end
function [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters()

    %% Define dimensional values (in SI units)
%     rho_m = 1200; % Membrane density
%     nu = 10^-3; % Membrane thickness
%     E = 2.2 * 10^9; % Membrane Young's modulus
%     T = 17; % Membrane tension
%     rho_l = 789; % Droplet density
%     R = 1e-3; % Droplet radius
%     V = 1; % Droplet speed
    

    %% Define parameters
    EPSILON = 1;
%     ALPHAS = rho_m * nu / (EPSILON^2 * rho_l * R)
%     BETAS =  EPSILON^2 * T / (rho_l * R * V)
%     GAMMAS = EPSILON^2 * E * nu^3 / (3 * rho_l * R^3 * V^2)

    % alpha_varying
%     ALPHAS = [1, 1.5, 2, 3, 4, 6, 8] / EPSILON^2;
%     BETAS = zeros(size(ALPHAS)) * EPSILON^2;
%     GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;

    % beta_varying
%     BETAS = [0, 10, 40, 160, 640, 2560, 10240] * EPSILON^2;
%     ALPHAS = ones(size(BETAS)) / EPSILON^2;
%     GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    
    % gamma_varying
    GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
    ALPHAS = 1 * ones(size(GAMMAS)) / EPSILON^2;
    BETAS = zeros(size(GAMMAS)) * EPSILON^2;
%  
    L = 16;
    T_MAX = 0.4;
    DELTA_T = 1e-4;

    % FD parameters
    N_MEMBRANE = 10924;
    
    IMPACT_TIME = 0.125;
end
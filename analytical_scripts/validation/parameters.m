function [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters()

    %% Define parameters
    EPSILON = 1;
    ALPHAS = 1 / EPSILON^2; 
%     ALPHAS = [0.5, 0.25, 0.125];
    BETAS =  1 * EPSILON^2; 
    GAMMAS = [0.1, 0.4, 1.6, 6.4] * EPSILON^2;
    L = 16;
    T_MAX = 0.1;
    DELTA_T = 1e-4;

    % FD parameters
    N_MEMBRANE = 1024;
    

end
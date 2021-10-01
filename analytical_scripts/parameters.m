function [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters()

    %% Define parameters
    EPSILON = 1;
    ALPHAS = 100 / EPSILON^2; 
    BETAS =  0 * EPSILON^2; 
%     GAMMAS = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4] * EPSILON;
    GAMMAS = 12.8;
    L = 16;
    T_MAX = 0.01;
    DELTA_T = 1e-4;

    % FD parameters
    N_MEMBRANE = 1024;
    

end
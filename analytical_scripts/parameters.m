function [EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters()

    %% Define parameters
    EPSILON = 1;
    ALPHAS = [1]; 
%     ALPHAS = [0.5, 0.25, 0.125];
    BETAS =  0 * EPSILON^2; 
    GAMMAS = [1.6, 3.2, 6.4] * EPSILON;
%     GAMMAS = 1;
    L = 16;
    T_MAX = 0.4;
    DELTA_T = 1e-4;

    % FD parameters
    N_MEMBRANE = 10923;
    
    IMPACT_TIME = 0.125;
end
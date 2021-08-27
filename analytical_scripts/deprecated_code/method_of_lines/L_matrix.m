function L = L_matrix(beta, gamma, M, delta_x)
    
    %% Define parameters
    Cbeta = beta / delta_x^2;
    Cgamma = gamma / delta_x^4;
    
    %% Define diagonals
    L_lower_lower = Cgamma * ones(M, 1);
    
    L_lower = -(Cbeta + 4 * Cgamma) * ones(M, 1);
    
    L_main = (2 * Cbeta + 6 * Cgamma) * ones(M, 1);
    L_main(2) = (2 * Cbeta + 7 * Cgamma);
    L_main(M) = (2 * Cbeta + 5 * Cgamma);
    
    L_upper = -(Cbeta + 4 * Cgamma) * ones(M, 1);
    L_upper(2) = -2 * (Cbeta + 4 * Cgamma);
    
    L_upper_upper = Cgamma * ones(M, 1);
    L_upper_upper(3) = 2 * Cgamma;
    
    %% Define matrix
    L = spdiags([L_lower_lower L_lower L_main L_upper L_upper_upper], ...
        -2 : 2, M, M);
end
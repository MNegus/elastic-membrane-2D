function [L, A, A0] = matrix_definitions(alpha, beta, gamma, M, delta_x, delta_t)
    
    %% Define parameters
    Cbeta = beta / delta_x^2;
    Cgamma = gamma / delta_x^4;
    
    %% Define L
    % Diagonals
    L_lower_lower = Cgamma * ones(M, 1);
    
    L_lower = -(Cbeta + 4 * Cgamma) * ones(M, 1);
    
    L_main = (2 * Cbeta + 6 * Cgamma) * ones(M, 1);
    L_main(2) = (2 * Cbeta + 7 * Cgamma);
    L_main(M) = (2 * Cbeta + 5 * Cgamma);
    
    L_upper = -(Cbeta + 4 * Cgamma) * ones(M, 1);
    L_upper(2) = -2 * (Cbeta + 4 * Cgamma);
    
    L_upper_upper = Cgamma * ones(M, 1);
    L_upper_upper(3) = 2 * Cgamma;
    
    % Sparse matrix
    L = spdiags([L_lower_lower L_lower L_main L_upper L_upper_upper], ...
        -2 : 2, M, M);
    
    %% Define A, equal to I + (delta_t^2 / alpha) * L
    % Diagonals
    A_lower_lower = (delta_t^2 / alpha) * L_lower_lower;
    A_lower = (delta_t^2 / alpha) * L_lower;
    A_main = ones(M, 1) + (delta_t^2 / alpha) * L_main;
    A_upper = (delta_t^2 / alpha) * L_upper;
    A_upper_upper = (delta_t^2 / alpha) * L_upper_upper;
    
    % Sparse matrix
    A = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], ...
        -2 : 2, M, M);
    
    %% Define A0, the matrix at t = 0
    A0 = spdiags([A_lower_lower A_lower A_main + ones(M, 1) A_upper A_upper_upper], ...
        -2 : 2, M, M);
end
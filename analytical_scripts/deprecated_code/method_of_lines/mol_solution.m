function [ws, qs, ds] = mol_solution(alpha, beta, gamma, epsilon, L, delta_t, N_membrane)

    %% Derived parameters
    M = N_membrane - 1;
    delta_x = L / (N_membrane - 1); 
    xs = (0 : delta_x : L - delta_x)';
    
    
    %% Initialise y
    
    
    
    
    
    %% Function definition
    function res = full_ode_fun(y, yp, M)
       
        %% Extract w and its derivatives
        w = y(1 : M);
        w_t = y(M + 1 : 2 * M);
        w_tt = yp(M + 1 : 2 * M);
        
        % First spacial derivative of w
        w_x = zeros(M, 1);
        w_x(2 : M - 1) = (w(3 : M) - w(1 : M - 2)) / (2 * DELTA_X);
        w_x(M) = -w(M - 1) / (2 * DELTA_X);
        
    end

end
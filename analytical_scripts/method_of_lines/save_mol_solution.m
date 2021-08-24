function save_mol_solution(data_dir, alpha, beta, gamma, epsilon, N_membrane, L, tmax, delta_t)
    
    %% Derived parameters
    ts = 0 : delta_t : tmax;
    DELTA_X = L / (N_membrane - 1); 
    M = N_membrane - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    
    

end
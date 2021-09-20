function save_finite_differences_solution_first_order(parent_dir, ...
    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, ...
    pressure_type)

    %% Saves data in sub-directory depending on pressure type
    fd_data_dir = sprintf("%s/%s", parent_dir, pressure_type);

    %% Derived parameters
    DELTA_X = L / (N_MEMBRANE - 1); 
    M = N_MEMBRANE - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : DELTA_T : T_MAX;
    
    if (GAMMA == 0)
        Cpressure = 4 * DELTA_X * DELTA_X / BETA;
    else
        Cpressure = 4 * DELTA_X^4 / GAMMA;
    end
    

    %% Matrix definitions
    [A_mat, B_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

    %% Initialise arrays
    
    % Composite solutions
    w_previous = zeros(size(xs));
    w = initialise_membrane(w_previous, A_mat, B_mat);
    w_t = zeros(size(xs));
    w_next = zeros(size(xs));
    p = zeros(size(xs));
    p_previous = zeros(size(xs));
    ds = zeros(size(T_VALS));
    d_ts = zeros(size(T_VALS));
    
    %% Save initial arrays
    save(sprintf("%s/w_%d.mat", fd_data_dir, 0), 'w_previous');
    save(sprintf("%s/w_%d.mat", fd_data_dir, 1), 'w');
    save(sprintf("%s/w_t_%d.mat", fd_data_dir, 0), 'w_t');
    save(sprintf("%s/p_%d.mat", fd_data_dir, 0), 'p');
    
    %% Loops over time
    for k = 1 : length(T_VALS) - 1
        %% Updates time
        t = T_VALS(k);
        t

        %% Composite timestep
        [w_next, p, w_t, d, d_t, J] = membrane_timestep_first_order(...
            xs, t, w, w_previous, p_previous, pressure_type, ...
            EPSILON, DELTA_T, DELTA_X, M, Cpressure, A_mat, B_mat);
        
        ds(k) = d;
        d_ts(k) = d_t;

        % Saves arrays
        save(sprintf("%s/w_%d.mat", fd_data_dir, k + 1), 'w_next');
        save(sprintf("%s/w_t_%d.mat", fd_data_dir, k), 'w_t');
        save(sprintf("%s/p_%d.mat", fd_data_dir, k), 'p');
        
        % Swaps ws
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;
        
        % Swaps pressure
        p_previous = p;

    end
    
    %% Saves ds solutions
    save(sprintf("%s/ds.mat", fd_data_dir), 'ds');
    save(sprintf("%s/d_ts.mat", fd_data_dir), 'd_ts');
    
end
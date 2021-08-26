function save_finite_differences_solution_second_order(parent_dir, ...
    alpha, beta, gamma, epsilon, N_membrane, L, tmax, delta_t, ...
    pressure_type)

    %% 
    fd_data_dir = sprintf("%s/%s", parent_dir, pressure_type);

    %% Derived parameters
    DELTA_X = L / (N_membrane - 1); 
    M = N_membrane - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : delta_t : tmax;
    
    if (gamma == 0)
        Cpressure = 4 * DELTA_X * DELTA_X / beta;
    else
        Cpressure = 4 * DELTA_X^4 / gamma;
    end
    

    %% Matrix definitions
    [A_mat, B_mat] = matrix_definitions(alpha, beta, gamma, M, DELTA_X, delta_t);

    %% Initialise arrays
    
    % Composite solutions
    w_previous = zeros(size(xs));
    w = initialise_membrane(w_previous, A_mat, B_mat);
    w_t = zeros(size(xs));
    w_t_previous = w_t;
    w_next = zeros(size(xs));
    p = zeros(size(xs));
    ds = zeros(size(T_VALS));
    d_ts = zeros(size(T_VALS));
    
    %% Save arrays
    % Composite
    save(sprintf("%s/w_%d.mat", fd_data_dir, 0), 'w_previous');
    save(sprintf("%s/w_t_%d.mat", fd_data_dir, 0), 'w_t');
    save(sprintf("%s/p_%d.mat", fd_data_dir, 0), 'p');
    
    %% Loops over time
    for k = 2 : length(T_VALS)
        %% Updates time
        t = T_VALS(k);
        t

        %% Composite timestep
        [w_next, p, w_t, d, d_t, J] ...
            = membrane_timestep_second_order(xs, t, ...
                w, w_previous, w_t_previous, ds(k - 1), d_ts(k - 1), ...
                pressure_type, ...
                epsilon, M, DELTA_X, delta_t, Cpressure, A_mat, B_mat);
        ds(k) = d;
        d_ts(k) = d_t;

        % Saves arrays
        save(sprintf("%s/w_%d.mat", fd_data_dir, k - 1), 'w_next');
        save(sprintf("%s/w_t_%d.mat", fd_data_dir, k - 1), 'w_t');
        save(sprintf("%s/p_%d.mat", fd_data_dir, k - 1), 'p');
        
        % Swaps ws
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;
        w_t_previous = w_t;

    end
    
    %% Saves ds solutions
    save(sprintf("%s/ds.mat", fd_data_dir), 'ds');
    save(sprintf("%s/d_ts.mat", fd_data_dir), 'd_ts');
    
end
function save_finite_differences_solution(data_dir, alpha, beta, gamma, epsilon, N_membrane, L, tmax, delta_t)

    %% Data dirs for each
    composite_dir = sprintf("%s/composite", data_dir);
    outer_dir = sprintf("%s/outer", data_dir);

    %% Derived parameters
    DELTA_X = L / (N_membrane - 1); 
    M = N_membrane - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : delta_t : tmax;
    
    if (gamma == 0)
        Cpressure = DELTA_X * DELTA_X / beta;
    else
        Cpressure = DELTA_X^4 / gamma;
    end
    

    %% Matrix definitions
    [A_mat, B_mat] = matrix_definitions(alpha, beta, gamma, M, DELTA_X, delta_t);

    %% Initialise arrays
    
    % Composite solutions
    w_previous_composite = zeros(size(xs));
    w_composite = initialise_membrane(w_previous_composite, A_mat, B_mat);
    w_t_composite = zeros(size(xs));
    w_next_composite = zeros(size(xs));
    p_previous_previous_composite = zeros(size(xs));
    p_previous_composite = zeros(size(xs));
    p_composite = zeros(size(xs));
    ds_composite = zeros(size(T_VALS));

    % Outer solutions
    w_previous_outer = zeros(size(xs));
    w_outer = initialise_membrane(w_previous_outer, A_mat, B_mat);
    w_t_outer = zeros(size(xs));
    w_next_outer = zeros(size(xs));
    p_previous_previous_outer = zeros(size(xs));
    p_previous_outer = zeros(size(xs));
    p_outer = zeros(size(xs));
    ds_outer = zeros(size(T_VALS));
    
    %% Save arrays
    % Composite
    save(sprintf("%s/w_%d.mat", composite_dir, 0), 'w_previous_composite');
    save(sprintf("%s/w_t_%d.mat", composite_dir, 0), 'w_t_composite');
    save(sprintf("%s/p_%d.mat", composite_dir, 0), 'p_composite');
    
    % Outer 
    save(sprintf("%s/w_%d.mat", outer_dir, 0), 'w_previous_outer');
    save(sprintf("%s/w_t_%d.mat", outer_dir, 0), 'w_t_outer');
    save(sprintf("%s/p_%d.mat", outer_dir, 0), 'p_outer');
    
    %% Loops over time
    for k = 2 : length(T_VALS)
        %% Updates time
        t = T_VALS(k);
        t

        %% Composite timestep
        [w_next_composite, p_composite, w_t_composite, d, d_t, J] = membrane_timestep(xs, t, ...
            w_composite, w_previous_composite, p_previous_composite, p_previous_previous_composite, ...
            "composite",  alpha, beta, gamma, epsilon, ...
            M, DELTA_X, delta_t, Cpressure, A_mat, B_mat);
        ds_composite(k) = d;

        % Saves arrays
        save(sprintf("%s/w_%d.mat", composite_dir, k - 1), 'w_next_composite');
        save(sprintf("%s/w_t_%d.mat", composite_dir, k - 1), 'w_t_composite');
        save(sprintf("%s/p_%d.mat", composite_dir, k - 1), 'p_composite');
        
        % Swaps ws
        temp = w_previous_composite;
        w_previous_composite = w_composite;
        w_composite = w_next_composite;
        w_next_composite = temp;

        % Swaps ps
        temp = p_previous_previous_composite;
        p_previous_previous_composite = p_previous_composite;
        p_previous_composite = p_composite;
        p_composite = temp;
        
        %% Outer timestep
        [w_next_outer, p_outer, w_t_outer, d, d_t] = membrane_timestep(xs, t, ...
            w_outer, w_previous_outer, p_previous_outer, p_previous_previous_outer, ...
            "outer",  alpha, beta, gamma, epsilon, ...
            M, DELTA_X, delta_t, Cpressure, A_mat, B_mat);
        ds_outer(k) = d;

        % Saves arrays
        save(sprintf("%s/w_%d.mat", outer_dir, k - 1), 'w_next_outer');
        save(sprintf("%s/w_t_%d.mat", outer_dir, k - 1), 'w_t_outer');
        save(sprintf("%s/p_%d.mat", outer_dir, k - 1), 'p_outer');
        
        % Swaps ws
        temp = w_previous_outer;
        w_previous_outer = w_outer;
        w_outer = w_next_outer;
        w_next_outer = temp;

        % Swaps ps
        temp = p_previous_previous_outer;
        p_previous_previous_outer = p_previous_outer;
        p_previous_outer = p_outer;
        p_outer = temp;

    end
    
    %% Saves ds solutions
    save(sprintf("%s/ds.mat", composite_dir), 'ds_composite');
    save(sprintf("%s/ds.mat", outer_dir), 'ds_outer');
    
    
end
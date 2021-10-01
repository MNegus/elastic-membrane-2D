function save_prescribed_normal_modes_solution(data_dir, alpha, beta, gamma, epsilon, N, L, tmax, delta_t)
    
    %% Derived parameters
    ts = 0 : delta_t : tmax;

    %% d dependent parameters
    d_max = 2 * sqrt(tmax);
    delta_d = delta_t;

    %% Solves the ode in d-form
    [t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals] ...
        = prescribed_a_ode_solution(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);

    %% Converts solution to t-form
    [ds, as, a_ts, a_tts, q_ts] = a_solution_t_form(ts, ...
        t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals, alpha, delta_t);
    
    %% Saves matrices in data_dir
    save(sprintf("%s/ds.mat", data_dir), 'ds');
    save(sprintf("%s/as.mat", data_dir), 'as');
    save(sprintf("%s/a_ts.mat", data_dir), 'a_ts');
    save(sprintf("%s/a_tts.mat", data_dir), 'a_tts');
    save(sprintf("%s/q_ts.mat", data_dir), 'q_ts');
end
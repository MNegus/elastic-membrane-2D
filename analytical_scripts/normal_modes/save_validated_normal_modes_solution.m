function save_validated_normal_modes_solution(data_dir, alpha, beta, gamma, epsilon, L, tmax, delta_t)

    [N, delta_d, ds, as, a_ts, a_tts, q_ts] ...
        = validated_normal_modes_solution(alpha, beta, gamma, epsilon, L, tmax, delta_t);
    
    %% Saves matrices in data_dir
    save(sprintf("%s/N.mat", data_dir), 'N');
    save(sprintf("%s/ds.mat", data_dir), 'ds');
    save(sprintf("%s/as.mat", data_dir), 'as');
    save(sprintf("%s/a_ts.mat", data_dir), 'a_ts');
    save(sprintf("%s/a_tts.mat", data_dir), 'a_tts');
    save(sprintf("%s/q_ts.mat", data_dir), 'q_ts');

end
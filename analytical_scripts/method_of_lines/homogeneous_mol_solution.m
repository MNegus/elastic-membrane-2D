function [ws, w_ts, w_tts] = homogeneous_mol_solution(alpha, beta, gamma, epsilon, L, delta_t, t_max, ts, N_membrane, w0)

    %% Derived parameters
    M = N_membrane - 1;
    delta_x = L / (N_membrane - 1); 
    xs = (0 : delta_x : L - delta_x)';
    t0 = 0;
    
    
    %% Initialise L matrix
    tic
    L = L_matrix(beta, gamma, M, delta_x);
    toc
    %% Define ODE function
    ode_fun = @(t, y, yp) full_ode_fun(t, y, yp, xs, M, L, alpha);
    
    %% Initialise y
    w_t0 = zeros(size(w0));
    w_tt0 = -L * w0 / alpha;
    y0 = [w0; w_t0];
    yp0 = [w_t0; w_tt0];
    size(y0)
    size(yp0)
    
    %% Test the output
    res = full_ode_fun(0, y0, yp0, xs, M, L, alpha)
    
    %% Solve the ode
    options = odeset('maxstep', delta_t, 'Initialstep', 1e-9);
    sol = ode15i(ode_fun, [0, t_max], y0, yp0, options);
    [y, yp] = deval(sol, ts);
    
    ws = y(1 : M, :);
    w_ts = yp(1 : M, :);
    w_tts = yp(M + 1 : 2 * M, :);
    
    
    %% Function definition
    function res = full_ode_fun(t, y, yp, xs, M, L, alpha)
        t
        res = zeros(2 * M, 1);
        
        %% Extract w and q
        w = y(1 : M, 1);
        w_t = y(M + 1 : 2 * M, 1);
        q = yp(1 : M, 1);
        q_t = yp(M + 1 : 2 * M, 1);
        
        %% Save res
        res(1 : M, 1) = q - w_t;
        res(M + 1 : 2 * M, 1) = alpha * q_t + L * w;
        plot(xs, w);
        drawnow;
        
        max(res)
    end

end
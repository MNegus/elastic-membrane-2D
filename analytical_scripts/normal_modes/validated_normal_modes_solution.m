function [N, delta_d, ds, as, a_ts, a_tts, q_ts] ...
    = validated_normal_modes_solution(alpha, beta, gamma, epsilon, L, tmax, delta_t)
    
    %% Derived parameters
    ts = 0 : delta_t : tmax;
    

    %% Parameters to be passed in
    q = 10;
    tol = 1e-4;
    N_MEMBRANE = 1024;
    DELTA_X = L / (N_MEMBRANE - 1); 
    xs = (0 : DELTA_X : L - DELTA_X)';

    %% d dependent parameters
    d_max = 2 * sqrt(tmax);

    %% Initial guess for delta_d and N
    delta_d = delta_t;

    %% Loops until found appropriate delta_d and N
    converged = 0;
    while (converged == 0)
        %% Determines N_max for current delta_d
        delta_d
        N_max = N_stable(alpha, beta, gamma, L, q, delta_d)
        
        %% Initialise N and solves ode
        N0 = min(64, floor(N_max / 4));
        N = N0;
        [t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals] ...
            = a_ode_solution(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);
        % Converts solution to t-form
        [ds, as, a_ts, a_tts, q_ts] = a_solution_t_form(ts, ...
            t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals, alpha, delta_t);
        
        
        %% Loops until we've found an N or N >= N_max
        diff = 1e6;
        while ((diff > tol) && (2 * N < N_max))
            % Doubles N
            new_N = 2 * N
            pause(2)
            % Solves ode with new N
            [new_t_vals_d_form, new_d_vals_d_form, new_as_d_form, new_a_ts_d_form, new_kvals] ...
                = a_ode_solution(alpha, beta, gamma, epsilon, delta_d, d_max, new_N, L);
            
            % New solution in t form
            [new_ds, new_as, new_a_ts, new_a_tts, new_q_ts] ...
                = a_solution_t_form(ts, new_t_vals_d_form, ...
                    new_d_vals_d_form, new_as_d_form, new_a_ts_d_form, ...
                    new_kvals, alpha, delta_t);
            
            % Compares ws solution for all time between the two
            diff = 0;
            figure(1);
            for k = 1 : length(ts)
                [ws, ~, ~] = w_solution_normal_modes(xs, as(k, :), a_ts(k, :), q_ts(k, :), ds(k), L, N, epsilon);
                [new_ws, ~, ~] = w_solution_normal_modes(xs, new_as(k, :), new_a_ts(k, :), new_q_ts(k, :), new_ds(k), L, new_N, epsilon);
                diff = max(diff, max(abs(ws - new_ws)))
                delta_d
                N
                if (diff > tol)
                    break;
                end
                plot(xs, ws);
                hold on;
                plot(xs, new_ws);
                hold off;
                title(['k = ', num2str(k), ', new_N = ', num2str(new_N), ', delta_d = ', num2str(delta_d)]);
                drawnow;
                pause(0.0000001);
            end
            
            % Update solutions
            N = new_N;
            ds = new_ds;
            as = new_as;
            a_ts = new_a_ts;
            a_tts = new_a_tts;
            q_ts = new_q_ts;
        end
       diff
        %% Checks for value of N
        if (diff < tol)
            converged = 1;
        else
%            N = N0;
           delta_d = delta_d / 10;
        end
    
        
    end

end
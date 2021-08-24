clear;

%% Parameters
ALPHA = 0.1;
BETA = 200;
GAMMA = 200;
EPSILON = 1;
L = 4;

T_MAX = 0.25;
DELTA_T = 1e-4;

N_stab = N_stable(ALPHA, BETA, GAMMA, L, 10, 1e-4)

%% N_MEMBRANE convergence
N_MEMBRANES = [16, 32, 64, 128, 256, 1024, 2056];
N_errors = zeros(size(N_MEMBRANES));

for q = 1 : length(N_MEMBRANES)
    
    N_MEMBRANE = N_MEMBRANES(q);
    
    maxdiff = 0;
    
    %% Derived parameters
    DELTA_X = L / (N_MEMBRANE - 1); 
    M = N_MEMBRANE - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : DELTA_T : T_MAX;

    %% Derive matrices
    [L_mat, A_mat, A0_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

    %% Initialise w_previous
    N_max = 8;
    as = (2.^(-1 : -1 : -N_max)).^0.1;
    as = 0.01 * as / sum(as);

    w_previous = exact_solution(xs, 0, ALPHA, BETA, GAMMA, L, as, N_max);
%     ymax = max(abs(w_previous));
%     ymin = -max(abs(w_previous));
% % 
%     figure(1);
%     plot(xs, w_previous);
%     ylim([ymin, ymax]);
%     pause(0.001);


    %% Initialise w
    % w = A0_mat \ w_previous;
    w = w_previous;

    % Loops over time
    for k = 2 : length(T_VALS)
        %% Updates time
        t = T_VALS(k);
        t

        %%
        w_next = homogeneous_membrane_timestep(w, w_previous, A_mat);
        w_exact = exact_solution(xs, t, ALPHA, BETA, GAMMA, L, as, N_max);

        maxdiff = max(maxdiff, max(abs(w - w_exact)));
        
        %% Plots
%         figure(1);
%         plot(xs, w_next);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
% 
%         ylim([ymin, ymax]);
% 
%         pause(1e-9);
%         drawnow;

        %% Swaps
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;


    end
    
    N_errors(q) = maxdiff
end

% Plots errors
figure(1);
loglog(N_MEMBRANES, N_errors, '-o');
title("N errors");


%% Timestep convergence
N_MEMBRANE = 2056;
DELTA_TS = [1e-2, 1e-3, 1e-4, 1e-5];
t_errors = zeros(size(DELTA_TS));

for q = 1 : length(DELTA_TS)
    
    DELTA_T = DELTA_TS(q);
    
    maxdiff = 0;
    
    %% Derived parameters
    DELTA_X = L / (N_MEMBRANE - 1); 
    M = N_MEMBRANE - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : DELTA_T : T_MAX;

    %% Derive matrices
    [L_mat, A_mat, A0_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

    %% Initialise w_previous
    N_max = 16;
    as = (2.^(-1 : -1 : -N_max)).^0.1;
    as = 0.01 * as / sum(as);

    w_previous = exact_solution(xs, 0, ALPHA, BETA, GAMMA, L, as, N_max);
%     ymax = max(abs(w_previous));
%     ymin = -max(abs(w_previous));
% % 
%     figure(1);
%     plot(xs, w_previous);
%     ylim([ymin, ymax]);
%     pause(0.001);


    %% Initialise w
    % w = A0_mat \ w_previous;
    w = w_previous;

    % Loops over time
    for k = 2 : length(T_VALS)
        %% Updates time
        t = T_VALS(k);
        t

        %%
        w_next = homogeneous_membrane_timestep(w, w_previous, A_mat);
        w_exact = exact_solution(xs, t, ALPHA, BETA, GAMMA, L, as, N_max);

        maxdiff = max(maxdiff, max(abs(w - w_exact)));
        
        %% Plots
%         figure(1);
%         plot(xs, w_next);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
% 
%         ylim([ymin, ymax]);
% 
%         pause(1e-9);
%         drawnow;

        %% Swaps
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;


    end
    
    t_errors(q) = maxdiff
end

%% Plots errors
figure(2);
loglog(DELTA_TS, t_errors, '-o');
title("t errors");

%% Exact solution
function ws = exact_solution(xs, t, ALPHA, BETA, GAMMA, L, as, N_max)
    ws = zeros(size(xs));
    for n = 1 : N_max
        lambda = pi * (2 * n - 1) / (2 * L);
        l = sqrt((BETA * lambda^2 + GAMMA * lambda^4) / ALPHA);
        ws = ws + as(n) * cos(lambda * xs) * cos(l * t);
    end
end

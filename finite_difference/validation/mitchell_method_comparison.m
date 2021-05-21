%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

%% Parameters
ALPHA = 0.1; BETA = 10; GAMMA = 0; 
L = 4;
N_MEMBRANE = 512; % Number of points on membrane
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.25;
DELTA_POWERS = linspace(-1, -4, 8);
DELTA_TS = 10.^DELTA_POWERS
BETA * DELTA_TS.^2 / (ALPHA * Deltax^2)
xs = linspace(0, L, N_MEMBRANE);
Deltax = L / (N_MEMBRANE - 1); 


%% Numerically solves and compares
max_errors = zeros(size(DELTA_TS));

for k = 1 : length(DELTA_TS)
    % Reset figure
    close(figure(1));
    figure(1);
    
    % Setting DELTA_T
    DELTA_T = DELTA_TS(k);
    Cbeta2 = BETA * DELTA_T^2 / (ALPHA * Deltax^2)
    
    % Initialises error term
    max_error = 0;
    
    % Matrix definitions
    A_upper = (-Cbeta2 / 4) * ones(M, 1);
    A_upper(2) = -Cbeta2 / 2;
    A_main = (1 + Cbeta2 / 2) * ones(M, 1);
    A_lower = (-Cbeta2 / 4) * ones(M, 1);
    A = spdiags([A_lower A_main A_upper], -1:1, M, M);

    B_upper = (Cbeta2 / 2) * ones(M, 1);
    B_upper(2) = Cbeta2;
    B_main = (2 - Cbeta2) * ones(M, 1);
    B_lower = (Cbeta2 / 2) * ones(M, 1);
    B = spdiags([B_lower B_main B_upper], -1:1, M, M);
    
    % Initial conditions
    t = 0
    w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
    w_previous = w_exact(1 : end - 1);

    plot(xs(1 : end - 1), w_previous);
    hold on;
    plot(xs, w_exact);
    hold off;
    pause(0.01);

    % Solve for first w
    t = t + DELTA_T
    rhs = 0.5 * B * w_previous;
    w = A \ rhs;
    
    w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
    error = max(abs(w - w_exact(1 : end - 1)));
    if error > max_error
        max_error = error;
    end

    plot(xs(1 : end - 1), w);
    hold on;
    plot(xs, w_exact);
    hold off;
    pause(0.01);

    % Loops
    while (t < T_MAX) 
        rhs = B * w - A * w_previous;
        w_next = A \ rhs;

        % Swaps
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;

        t = t + DELTA_T

        % Exact solution
        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);

        % Error saving
        error = max(abs(w - w_exact(1 : end - 1)));
        if error > max_error
            max_error = error;
        end
        
        % Plots
        plot(xs(1 : end - 1), w);
        hold on;
        plot(xs, w_exact);
        hold off;

        pause(0.01);

    end
    
    % Saves max error
    max_errors(k) = max_error
end


%% Plot max errors
loglog(DELTA_TS, max_errors, '-o');
hold on;
plot(DELTA_TS, 40 * DELTA_TS.^2);
% plot(DTS, 15 * DTS.^1.5);
% % plot(DTS, 0.5 * DTS);
set(gca, 'XDir','reverse');
%% Plots exact solution
% figure(1);
% t = 0;
% while (t <= T_MAX) 
%     t
%     ws = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
%     plot(xs, ws);
%     pause(0.01);
% 
%     t = t + DELTA_T;
% end
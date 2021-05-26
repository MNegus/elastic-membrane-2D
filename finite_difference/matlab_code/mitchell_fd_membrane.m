%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

%% Parameters
ALPHA = 0.1; BETA = 1; GAMMA = 1; 
L = 4;
N_MEMBRANES = [512, 1024, 2048, 4096];
T_MAX = 0.1;
DELTA_POWERS = linspace(-1, -5, 4);
DELTA_TS = 10.^DELTA_POWERS;

%% Numerically solves and compare
close(figure(2));
figure(2);
set(gca, 'XDir','reverse');
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');

for N_MEMBRANE = N_MEMBRANES
    M = N_MEMBRANE - 1; % We ignore the end point
    xs = linspace(0, L, N_MEMBRANE);
    Deltax = L / (N_MEMBRANE - 1); 
    max_errors = zeros(size(DELTA_TS));

    for k = 1 : length(DELTA_TS)
        % Reset figure
%         close(figure(1));
%         figure(1);

        % Setting DELTA_T
        DELTA_T = DELTA_TS(k);
        Cbeta = BETA * Deltax^2 / GAMMA;
        Calpha = ALPHA * Deltax^4 / (GAMMA * DELTA_T^2);

        % Initialises error term
        max_error = 0;

        % A definition
        A_upper_upper = ones(M, 1);
        A_upper_upper(3) = 2;

        A_upper = (-Cbeta - 4) * ones(M, 1);
        A_upper(2) = -2 * Cbeta - 8;

        A_main = (4 * Calpha + 2 * Cbeta + 6) * ones(M, 1);
        A_main(2) = A_main(2) + 1;
        A_main(M) = A_main(M) - 1;

        A_lower = (-Cbeta - 4) * ones(M, 1);

        A_lower_lower = ones(M, 1);

        A = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

        % B definition
        B_upper_upper = -2 * ones(M, 1);
        B_upper_upper(3) = 2 * B_upper_upper(3);

        B_upper = (2 * Cbeta + 8) * ones(M, 1);
        B_upper(2) = 2 * B_upper(2);

        B_main = (8 * Calpha - 4 * Cbeta - 12) * ones(M, 1);
        B_main(2) = B_main(2) - 2;
        B_main(M) = B_main(M) + 2;

        B_lower = (2 * Cbeta + 8) * ones(M, 1);

        B_lower_lower = -2 * ones(M, 1);

        B = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);

        % Initial conditions
        t = 0
        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        w_previous = w_exact(1 : end - 1)

%         plot(xs(1 : end - 1), w_previous);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
%         pause(0.01);

        % Solve for first w
        t = t + DELTA_T
        rhs = 0.5 * B * w_previous;
        w = A \ rhs;

        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        error = max(abs(w - w_exact(1 : end - 1)));
        if error > max_error
            max_error = error;
        end

%         plot(xs(1 : end - 1), w);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
%         pause(0.01);

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
    %         plot(xs(1 : end - 1), w);
    %         hold on;
    %         plot(xs, w_exact);
    %         hold off;

            pause(0.01);

        end

        % Saves max error
        max_errors(k) = max_error
    end
    figure(2);
    hold on;
    plot(DELTA_TS, max_errors, '-o');
    
end

%% Plot max errors
hold on;
plot(DELTA_TS, 40 * DELTA_TS.^2);
% plot(DTS, 15 * DTS.^1.5);
% % plot(DTS, 0.5 * DTS);
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
legend(["N = 512", "N = 1024", "N = 2048", "N = 4096", "y ~ 1 / dt^2"]);

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
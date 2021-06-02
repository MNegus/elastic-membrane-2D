%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

%% Parameters
ALPHA = 0.1; BETA = 0.001; GAMMA = 0.0001; 
L = 4;
% N_MEMBRANES = [512, 1024, 2048, 4096, 8192];
N_MEMBRANES = 128;
T_MAX = 0.2;
% DELTA_POWERS = linspace(-1, -5, 5);
% DELTA_TS = 10.^DELTA_POWERS;
DELTA_TS = 1e-4;

Deltaxs = L ./ (N_MEMBRANES - 1);
Cbetas = BETA * Deltaxs.^2 / GAMMA
Calphas = ALPHA * Deltaxs.^4 ./ (GAMMA * DELTA_TS'.^2)

%% Anonymous functions
lambda = @(n) pi * (2 * n - 1) / (2 * L);
l = @(n) sqrt(BETA * lambda(n).^2 + GAMMA * lambda(n).^4) / sqrt(ALPHA);


%% Initial condition
N0 = 3; % Take 3 terms
As = [1, 0.5, 0.25]; % Size of A_n

%% Numerically solves and compare
close(figure(2));
figure(2);
set(gca, 'XDir','reverse');
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');

for N_MEMBRANE = N_MEMBRANES
    M = N_MEMBRANE - 1; % We ignore the end point
    xs = linspace(0, L, N_MEMBRANE);
    reduced_xs = xs(1 : end - 1);
    Deltax = L / (N_MEMBRANE - 1); 
    
    L2_norms = zeros(size(DELTA_TS));

    for k = 1 : length(DELTA_TS)
        % Reset figure
%         close(figure(1));
%         figure(1);
%         
        
        % Setting DELTA_T
        DELTA_T = DELTA_TS(k);
        tvals = 0 : DELTA_T : T_MAX;
        Cbeta = BETA * Deltax^2 / GAMMA;
        Calpha = ALPHA * Deltax^4 / (GAMMA * DELTA_T^2);
        
        % Saving abs error in time
        abs_errors = zeros(size(tvals));

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

        % B definition5
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

%         
%         eigenvalues = eig(full(A));
%         figure(8);
%         plot(eigenvalues);
%         determinant = det(A)
%         prod(eigenvalues)
%         length(eigenvalues)

        % Initial conditions
        t = 0
        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        w_previous = w_exact(1 : end - 1);

        % Solve for first w
        t = t + DELTA_T
        rhs = 0.5 * B * w_previous;
        w = A \ rhs;

        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        
        % Initialise error norm
        L2_norm = sum((w - w_exact(1 : end - 1)).^2);

        % Loops
        q = 2;
        while (t <= T_MAX) 
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
            L2_norm = L2_norm + sum((w - w_exact(1 : end - 1)).^2);
            L2_norms

%             Plots full
            figure(1);
            plot(xs, w_exact);

%             title("Full");
%             legend(["Exact", "Numerical"]);
            

            % Save errors
            max(abs(w - w_exact(1 : end - 1)))
            abs_errors(q) = max(abs(w - w_exact(1 : end - 1)));
            q = q + 1;

            pause(0.01);

        end
        
        

        % Saves L2_norm
        L2_norm = sqrt(L2_norm / (length(xs) * (T_MAX / DELTA_T)));
        L2_norms(k) = L2_norm
    end
%     figure(2);
%     hold on;
%     plot(DELTA_TS, L2_norms, '-o');
    
end

%%
% close(figure(2));
% figure(2);
% plot(tvals, abs_errors);
% hold on;
% plot(tvals, 0.000009 * tvals.^2);
% legend(["Measured", "Guessed"]);
%% Plot max errors
% hold on;
% % plot(DELTA_TS, 40 * DELTA_TS.^2);
% % plot(DTS, 15 * DTS.^1.5);
% % % plot(DTS, 0.5 * DTS);
% set(gca, 'yscale','log');
% set(gca, 'xscale', 'log');
% legend(["N = 512", "N = 1024", "N = 2048", "N = 4096", "N = 8192", "y ~ 1 / dt^2"]);
% exportgraphics(gca, "timestep_validation.png");

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